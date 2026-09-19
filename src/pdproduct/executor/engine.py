"""Deterministic executor: a small state machine for validated plans.

Execution is stepwise; each step maps to a registered callable. The executor
never re-runs completed steps (resume after failure) and refuses plans that
did not pass validation.
"""
from __future__ import annotations

from typing import Callable

from ..planner.plan_ir import PlanRefusal

STEP_REGISTRY: dict[str, Callable] = {}


def register_step(name: str):
    def deco(fn):
        STEP_REGISTRY[name] = fn
        return fn
    return deco


class PlanExecutor:
    STATES = ("planned", "validated", "executing", "completed", "refused", "failed")

    def __init__(self, plan: dict):
        if not plan.get("validated"):
            raise PlanRefusal("not_validated", "execute only validated plans (run validate_plan first)")
        self.plan = plan
        self.state = "validated"
        self.done: dict[str, object] = {}
        self.error: str | None = None

    def run(self) -> "PlanExecutor":
        self.state = "executing"
        remaining = {s["step_id"]: s for s in self.plan["steps"]}
        progressed = True
        while remaining and progressed:
            progressed = False
            for sid, step in list(remaining.items()):
                deps = step.get("depends_on", [])
                if not all(d in self.done for d in deps):
                    continue
                fn = STEP_REGISTRY.get(step["method"])
                if fn is None:
                    self.state, self.error = "failed", f"step method {step['method']!r} not registered"
                    return self
                try:
                    self.done[sid] = fn(plan=self.plan, upstream={d: self.done[d] for d in deps})
                except Exception as e:  # failed step keeps state for resume
                    self.state, self.error = "failed", f"step {sid}: {e}"
                    return self
                del remaining[sid]
                progressed = True
        if remaining:  # cycle or unsatisfiable deps
            self.state, self.error = "failed", f"unsatisfiable dependencies: {sorted(remaining)}"
        else:
            self.state = "completed"
        return self

    def resume(self) -> "PlanExecutor":
        """Re-enter run() keeping completed steps (failure recovery)."""
        if self.state != "failed":
            raise RuntimeError("resume only from failed state")
        self.state = "executing"
        return self.run()
