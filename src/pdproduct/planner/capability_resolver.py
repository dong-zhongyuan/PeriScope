"""Capability resolver (pillar 1): structural guard against illegal bridging.

Consumes structured `io` contracts in the method registry and resolves which
methods may legally serve a task. A bridge the registry cannot type is
rejected — this is the enforcement layer behind the typed capability graph.
"""
from __future__ import annotations

from ..registry.methods import all_methods

MODALITIES = {
    "rna", "adt_protein", "plasma_protein", "spatial_rna", "spatial_protein",
    "image", "genetic", "embedding", "coordinates",
}


class BridgeRefusal(Exception):
    def __init__(self, reason: str, detail: str):
        self.reason = reason
        self.detail = detail
        super().__init__(f"BRIDGE REFUSED [{reason}]: {detail}")


def _typed_methods() -> dict[str, dict]:
    out = {}
    for mid, m in all_methods().items():
        io = m.get("io")
        if isinstance(io, dict) and io.get("input_modalities") and io.get("output_modalities"):
            out[mid] = m
    return out


def resolve(
    input_modalities: set[str],
    output_modality: str,
    *,
    require_measured_spatial_anchor: bool = False,
    entity_level: str | None = None,
) -> dict:
    """Return {legal: [...], rejected: [{method, reason}]} for a task.

    require_measured_spatial_anchor: if the output is spatial, methods whose
    io lacks a measured spatial anchor input are structurally rejected (G3b).
    """
    for m in input_modalities | {output_modality}:
        if m not in MODALITIES:
            raise BridgeRefusal("unknown_modality", f"{m!r} not in MODALITIES")
    if output_modality.startswith("spatial") and require_measured_spatial_anchor:
        if not (input_modalities & {"spatial_protein", "spatial_rna", "coordinates"}):
            raise BridgeRefusal(
                "no_measured_spatial_anchor",
                "spatial output requested without any measured spatial input (G3b)",
            )

    legal, rejected = [], []
    for mid, m in _typed_methods().items():
        io = m["io"]
        if not set(io["input_modalities"]) & input_modalities:
            rejected.append({"method": mid, "reason": "no_matching_input"})
            continue
        if output_modality not in io["output_modalities"]:
            rejected.append({"method": mid, "reason": "output_not_offered"})
            continue
        if entity_level and io.get("entity_level") and entity_level not in io["entity_level"]:
            rejected.append({"method": mid, "reason": "entity_level_mismatch"})
            continue
        legal.append(mid)
    return {"legal": sorted(legal), "rejected": rejected, "n_typed_methods": len(_typed_methods())}
