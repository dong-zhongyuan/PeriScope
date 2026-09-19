"""Donor-level split generation with leakage guards.

Split Contract (docs/BENCHMARK_PROTOCOL_DRAFT.md §1): the isolation unit is
the donor; every sample/cell/slice of one donor lands in exactly one
partition. Random cell splits are structurally impossible here.
"""
from __future__ import annotations

import hashlib
import json
import random
from pathlib import Path

PARTITIONS = ["train", "calibration", "model_selection", "locked_test"]


def make_donor_split(
    donor_facts: dict[str, dict],
    *,
    ratios: dict[str, float] | None = None,
    condition_key: str = "condition",
    seed: int = 42,
    min_per_condition: dict[str, int] | None = None,
) -> dict:
    """Stratify donors by condition into partitions.

    donor_facts: {donor_id: {condition: str, ...}}. Ratios must sum to ~1.0.
    Returns {"donors": {donor_id: partition}, "meta": {...}}.
    """
    ratios = ratios or {"train": 0.6, "calibration": 0.2, "model_selection": 0.2, "locked_test": 0.0}
    assert abs(sum(ratios.values()) - 1.0) < 1e-6, f"ratios must sum to 1: {ratios}"

    by_condition: dict[str, list[str]] = {}
    for donor_id, facts in donor_facts.items():
        by_condition.setdefault(str(facts.get(condition_key, "?")), []).append(donor_id)

    rng = random.Random(seed)
    assignment: dict[str, str] = {}
    for cond, donors in sorted(by_condition.items()):
        donors = sorted(donors)
        rng.shuffle(donors)
        n = len(donors)
        # proportional base, then guarantee every partition >= 1 donor when
        # the condition has enough donors (prevents empty locked_test arms),
        # then absorb rounding drift into the last partition.
        partitions = [p for p in PARTITIONS if ratios.get(p, 0) > 0]
        cuts: dict[str, int] = {p: int(ratios[p] * n) for p in partitions}
        if n >= len(partitions):
            for p in partitions:
                if cuts[p] == 0:
                    richest = max(partitions, key=lambda q: cuts[q])
                    if cuts[richest] > 1:
                        cuts[richest] -= 1
                        cuts[p] = 1
        cuts[partitions[-1]] += n - sum(cuts.values())
        idx = 0
        for p in partitions:
            for d in donors[idx : idx + cuts[p]]:
                assignment[d] = p
            idx += cuts[p]
        if min_per_condition and cond in min_per_condition:
            # basic sanity: if a condition is starved below min, flag via exception
            counts = {p: sum(1 for d in donors if assignment[d] == p) for p in partitions}
            if any(c < min_per_condition[cond] for c in counts.values() if c > 0) and n < min_per_condition[cond]:
                raise ValueError(f"condition {cond}: only {n} donors < min {min_per_condition[cond]}")

    return {
        "donors": assignment,
        "meta": {
            "seed": seed,
            "ratios": ratios,
            "condition_key": condition_key,
            "n_donors": len(assignment),
            "condition_counts": {
                cond: len(ds) for cond, ds in sorted(by_condition.items())
            },
        },
    }


def validate_no_leakage(split: dict, cell_donor_map: dict[str, str]) -> None:
    """Every cell's donor must be assigned to exactly one partition.

    cell_donor_map: {cell_id: donor_id}. Raises on unassigned donor.
    """
    assigned = split["donors"]
    missing = {d for d in cell_donor_map.values() if d not in assigned}
    if missing:
        raise ValueError(f"leakage guard: donors present in data but absent from split: {sorted(missing)}")


def save_split(split: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(split, f, indent=2, ensure_ascii=False, sort_keys=True)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    sidecar = path.with_suffix(".sha256")
    sidecar.write_text(f"{digest}  {path.name}\n", encoding="utf-8")
    return path


def load_split(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)
