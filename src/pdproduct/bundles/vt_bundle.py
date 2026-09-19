"""Partial-Observation VirtualTissueBundle (pillar 3).

Four layers — Observation / Inference / Intervention-Scenario / Evidence —
each field carrying orthogonal semantic axes (observation / spatial /
scenario status). Axes never merge into one truth flag; validation enforces
this via the artifact contract.
"""
from __future__ import annotations

import json
from pathlib import Path

from ..contracts.loader import validate

AXES = ("observation_status", "spatial_status", "scenario_status")


def _axes_ok(entry: dict) -> bool:
    try:
        validate(
            {
                "artifact_id": entry.get("artifact_id", "tmp"),
                "run_id": entry.get("run_id", "tmp"),
                "path": entry.get("data_ref", "tmp"),
                "sha256": entry.get("sha256", "0" * 64),
                **{k: entry[k] for k in AXES},
                "entity_level": entry["entity_level"],
                "uncertainty_status": entry.get("uncertainty_status", "unknown"),
            },
            "artifact",
        )
        return True
    except Exception:
        return False


class VirtualTissueBundle:
    """In-memory bundle; serialize with to_json()/write()."""

    def __init__(self, bundle_id: str, disease_pack: str = "PD"):
        self.bundle_id = bundle_id
        self.disease_pack = disease_pack
        self.layers: dict[str, list[dict]] = {
            "observation": [],
            "inference": [],
            "intervention_scenario": [],
            "evidence": [],
        }

    def _add(self, layer: str, entry: dict) -> "VirtualTissueBundle":
        required = set(AXES) | {"entity_level", "modality", "data_ref"}
        missing = required - set(entry)
        if missing:
            raise ValueError(f"{layer} entry missing {sorted(missing)}")
        if layer == "observation" and entry["observation_status"] != "measured":
            raise ValueError("observation layer must be measured")
        if layer == "inference" and entry["observation_status"] != "inferred":
            raise ValueError("inference layer must be inferred")
        if layer == "intervention_scenario" and entry["scenario_status"] == "observed":
            raise ValueError("scenario layer must be model_intervention or identified_counterfactual")
        self.layers[layer].append(entry)
        return self

    def add_observation(self, **kw):
        return self._add("observation", kw)

    def add_inference(self, **kw):
        return self._add("inference", kw)

    def add_intervention_scenario(self, **kw):
        return self._add("intervention_scenario", kw)

    def add_evidence_card(self, card: dict):
        if "claim_id" not in card:
            raise ValueError("evidence entry needs claim_id")
        self.layers["evidence"].append({"claim_id": card["claim_id"], **card})
        return self

    def validate_axes(self) -> dict:
        bad = []
        for layer, entries in self.layers.items():
            if layer == "evidence":
                continue  # cards are governed by the evidence_card contract
            for e in entries:
                if not _axes_ok(e):
                    bad.append((layer, e.get("data_ref", "?")))
        return {"ok": not bad, "invalid_entries": bad}

    def summary(self) -> dict:
        return {
            "bundle_id": self.bundle_id,
            "disease_pack": self.disease_pack,
            "counts": {k: len(v) for k, v in self.layers.items()},
            "axes_valid": self.validate_axes()["ok"],
        }

    def write(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(
                {"bundle_id": self.bundle_id, "disease_pack": self.disease_pack, "layers": self.layers},
                f, indent=2, ensure_ascii=False, sort_keys=True,
            )
        return path
