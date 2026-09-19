"""Plan IR validation: seven pre-flight gates (HANDOFF_SERVER_V1 §17.4).

A plan that fails any gate is refused BEFORE execution; refusal is structured.
Gates: 1 schema, 2 capability, 3 authorization, 4 data availability,
       5 claim legality, 6 resource, 7 privacy/license.
"""
from __future__ import annotations

import json
from pathlib import Path

from ..claim_compiler.compiler import CLAIM_CLASS_ORDER
from ..contracts.loader import repo_root, validate
from ..planner.capability_resolver import BridgeRefusal, resolve

AUTHORIZED_SCOPE_PREFIXES = ("phaseB-", "phaseA-", "user-granted")


class PlanRefusal(Exception):
    def __init__(self, gate: str, detail: str):
        self.gate = gate
        self.detail = detail
        super().__init__(f"PLAN REFUSED at gate [{gate}]: {detail}")


def validate_plan(plan: dict, *, available_datasets: set[str] | None = None) -> dict:
    """Run the seven gates; return the stamped plan or raise PlanRefusal."""
    # Gate 1: schema
    try:
        validate(plan, "plan_ir")
    except Exception as e:
        raise PlanRefusal("schema", str(e)[:200]) from e

    # Gate 2: capability — every selected method must be legal for the IO bridge
    if plan.get("input_modalities") and plan.get("output_modality"):
        need_anchor = plan["output_modality"].startswith("spatial")
        try:
            out = resolve(set(plan["input_modalities"]), plan["output_modality"], require_measured_spatial_anchor=need_anchor)
        except BridgeRefusal as e:
            raise PlanRefusal("capability", str(e)) from e
        for m in plan["selected_methods"]:
            if m not in out["legal"]:
                raise PlanRefusal("capability", f"method {m!r} not legal for {plan['input_modalities']} -> {plan['output_modality']}")

    # Gate 3: authorization
    scope = plan.get("authorization_scope", "")
    if not scope.startswith(AUTHORIZED_SCOPE_PREFIXES):
        raise PlanRefusal("authorization", f"scope {scope!r} not in authorized scope prefixes {AUTHORIZED_SCOPE_PREFIXES}")

    # Gate 4: data availability — referenced datasets must have manifests
    if available_datasets is not None:
        for ref in plan.get("observed_inputs", []):
            ds = ref.split("/")[0]
            if ds not in available_datasets:
                raise PlanRefusal("data_availability", f"dataset {ds!r} has no manifest")

    # Gate 5: claim legality — class must exist; C5 requires identified assumptions
    if plan["requested_claim_class"] not in CLAIM_CLASS_ORDER:
        raise PlanRefusal("claim_legality", "unknown claim class")
    if plan["requested_claim_class"] == "C5_identified_counterfactual":
        raise PlanRefusal("claim_legality", "C5 has no evidence source in this platform (empty set by design)")

    # Gate 6: resource
    budget = plan.get("compute_budget") or {}
    if budget.get("max_gpu_hours", 0) > 48 or budget.get("max_ram_gb", 0) > 800:
        raise PlanRefusal("resource", "budget exceeds host guardrails (48 GPU-hours / 800 GB RAM)")

    # Gate 7: privacy/license — no controlled datasets in inputs without D15 clearance
    for ref in plan.get("observed_inputs", []):
        if "ADNI" in ref or "clinical" in ref.lower():
            raise PlanRefusal("privacy_license", f"controlled asset {ref!r} requires D15 clearance")

    return {**plan, "validated": True, "gates_passed": 7}


def load_plan(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as f:
        return json.load(f)
