"""Evidence card creation and validation (first-class platform output).

A card is the only object allowed to carry a scientific claim; it must pass
the evidence_card contract and reference run manifests it was derived from.
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

from ..contracts.loader import validate
from ..provenance.run_manifest import git_state


def new_evidence_card(
    claim_id: str,
    *,
    claim_class: str,
    claim_text_allowed: str,
    estimand: str,
    analysis_unit: str,
    evidence_runs: list[dict],
    primary_statistic: dict,
    null_control: list[dict],
    uncertainty_status: str,
    forbidden_wording: list[str],
    open_gaps: list[str],
    target_population: str,
    validation_grade: str = "none",
    observation_status: str = "measured",
    spatial_status: str = "not_spatial",
    scenario_status: str = "observed",
    input_refs: list[str] | None = None,
    git_commit: str | None = None,
) -> dict:
    if git_commit is None:
        git_commit, _ = git_state()
    card = {
        "claim_id": claim_id,
        "claim_class": claim_class,
        "claim_text_allowed": claim_text_allowed,
        "estimand": estimand,
        "analysis_unit": analysis_unit,
        "target_population": target_population,
        "evidence_runs": evidence_runs,
        "primary_statistic": primary_statistic,
        "null_control": null_control,
        "uncertainty_status": uncertainty_status,
        "validation_grade": validation_grade,
        "forbidden_wording": forbidden_wording,
        "open_gaps": open_gaps,
        "observation_status": observation_status,
        "spatial_status": spatial_status,
        "scenario_status": scenario_status,
        "provenance": {
            "git_commit": git_commit,
            "created_at": _dt.datetime.now(_dt.timezone.utc).astimezone().isoformat(timespec="seconds"),
            "input_refs": input_refs or [],
            "claim_policy_version": "docs/CLAIM_POLICY_DRAFT.md v1",
        },
    }
    validate(card, "evidence_card")
    return card


def write_card(card: dict, out_dir: str | Path) -> Path:
    import json

    out = Path(out_dir) / f"{card['claim_id']}.json"
    if out.exists():
        raise FileExistsError(f"evidence card already exists: {out}")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        json.dump(card, f, indent=2, ensure_ascii=False, sort_keys=True)
    return out
