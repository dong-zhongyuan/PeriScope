"""Method capability registry: typed access to METHOD_CAPABILITY_REGISTRY.yaml.

The registry is the structural guard against illegal analysis bridging:
every method declares its inputs, outputs, forbidden uses and MVP role.
No method may execute outside what its registry entry permits.
"""
from __future__ import annotations

from functools import lru_cache
from typing import Iterable

import yaml

from ..contracts.loader import repo_root

READINESS_LEVELS = [
    "KNOWLEDGE-ONLY",
    "SOURCE-AVAILABLE",
    "CHALLENGER",
    "CONDITIONALLY-READY",
    "PRODUCTION-ELIGIBLE",
]


@lru_cache(maxsize=1)
def _load() -> dict:
    path = repo_root() / "registry" / "METHOD_CAPABILITY_REGISTRY.yaml"
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def meta() -> dict:
    return _load()["meta"]


def all_methods() -> dict[str, dict]:
    return {m["method_id"]: m for m in _load()["methods"]}


def get(method_id: str) -> dict:
    try:
        return all_methods()[method_id]
    except KeyError:
        raise KeyError(
            f"method {method_id!r} not in registry; known: {sorted(all_methods())}"
        ) from None


def query(
    category: str | None = None,
    mvp_role_contains: str | None = None,
    evidence: str | None = None,
    readiness: str | None = None,
) -> dict[str, dict]:
    out = {}
    for mid, m in all_methods().items():
        if category is not None and m.get("category") != category:
            continue
        if mvp_role_contains is not None and mvp_role_contains not in str(m.get("mvp_role", "")):
            continue
        if evidence is not None and m.get("upstream_evidence") != evidence:
            continue
        if readiness is not None and m.get("readiness") != readiness:
            continue
        out[mid] = m
    return out


def mvp_core() -> dict[str, dict]:
    """The two mandatory baselines of the 2+1 first-version scope.

    sciPENN is the conditional third method and is NOT part of the core;
    it only enters after passing its five admission gates.
    """
    core_ids = ["baselines-cell-type-mean", "baselines-ridge"]
    return {mid: get(mid) for mid in core_ids}


def readiness_summary() -> dict[str, int]:
    counts = {lvl: 0 for lvl in READINESS_LEVELS}
    for m in all_methods().values():
        lvl = m.get("readiness")
        if lvl in counts:
            counts[lvl] += 1
    return counts


def forbidden_claims(method_id: str) -> list:
    """Forbidden uses of a method; the claim compiler consults this before
    allowing any output of the method into an evidence card."""
    return list(get(method_id).get("forbidden", []))


def assert_no_forbidden_claim(method_id: str, claim_text: str) -> None:
    """Raise if a proposed claim text hits a method's forbidden-use patterns.

    Exact substring match on each forbidden entry; claim_compiler will layer
    the full seven-class policy on top of this structural check.
    """
    for pattern in forbidden_claims(method_id):
        if pattern and pattern in claim_text:
            raise ValueError(
                f"claim violates registry boundary for {method_id}: {pattern!r}"
            )
