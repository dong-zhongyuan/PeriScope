"""Claim Compiler: machine-decidable claim policy enforcement.

Implements docs/CLAIM_POLICY_DRAFT.md: the seven orthogonal claim classes,
the forbidden-wording table, evidence-dependency legality, and structured
refusal. The compiler has veto power; LLM proposals are inputs only.

Amendment 1 (ADR-0003, PI directive 2026-09-18): goal/direction wording
(therapeutic lever, treatment dose-response, biomarker, diagnostic) is
permitted WITH tier labels (model-layer / mechanism-layer); only
evidence-overclaim phrasing remains forbidden. The C1-C7 ladder and the
illegal_upgrade rule are unchanged -- the amendment opens goal language,
not evidence classes.
"""
from __future__ import annotations

import re

CLAIM_CLASS_ORDER = [
    "C1_association",
    "C2_predictive_mapping",
    "C3_spatial_colocalization",
    "C4_model_intervention_response",
    "C5_identified_counterfactual",
    "C6_genetically_proxied_effect",
    "C7_experimentally_supported_effect",
]
_CLASS_RANK = {c: i for i, c in enumerate(CLAIM_CLASS_ORDER)}

# Forbidden wording -> required replacement (CLAIM_POLICY §4, as amended by ADR-0003:
# goal/direction entries removed -- therapeutic lever / treatment dose-response /
# biomarker / diagnostic are legal with tier labels; evidence-overclaim entries kept)
FORBIDDEN_PHRASES: dict[str, str] = {
    "patient digital twin": "partial-observation virtual tissue / evidence bundle",
    "blood-to-brain causal": "cross-cohort conditional coupling candidate",
    "pharmacological synergy": "(no legal substitute: no evidence source)",
    "spatially validated": "brain-side spatial compatibility / colocalization evidence",
    "ai scientist": "evidence-aware analysis / assay-design agent",
    "self-driving lab": "assay recommendation (no prospective wet-lab loop)",
    "independent validation": "validation (list donor/study/modality independence explicitly)",
    "digital twin": "partial-observation virtual tissue",
    "完整 world model": "equilibrium simulator candidate",
    "已集成 28 个专家": "28 skill cards, checkpoint-reproduced=0",
    "28 experts": "28 skill cards",
    "causal": "(requires C5/C6/C7 evidence card; association-level wording otherwise)",
}


class ClaimRefusal(Exception):
    """Structured refusal: what is missing and which rule fired."""

    def __init__(self, rule: str, detail: str, remediation: str = ""):
        self.rule = rule
        self.detail = detail
        self.remediation = remediation
        super().__init__(f"CLAIM REFUSED [{rule}]: {detail}")


class ClaimCompiler:
    def compile(
        self,
        claim_class: str,
        claim_text: str,
        evidence_cards: list[dict],
    ) -> dict:
        """Validate and stamp a claim. Raises ClaimRefusal on any violation.

        evidence_cards: loaded card dicts (claim_id, claim_class, ...).
        """
        if claim_class not in _CLASS_RANK:
            raise ClaimRefusal(
                "unknown_class",
                f"{claim_class!r} is not one of the seven orthogonal classes",
                "use a C1..C7 class from CLAIM_POLICY §1",
            )

        text = claim_text.lower()
        hits = {p: r for p, r in FORBIDDEN_PHRASES.items() if p in text}
        if hits:
            raise ClaimRefusal(
                "forbidden_wording",
                f"claim text contains forbidden phrasing: {sorted(hits)}",
                f"use replacements: {hits}",
            )

        if not evidence_cards:
            raise ClaimRefusal(
                "no_evidence",
                "a claim cannot exist without at least one evidence card",
                "attach EC-* cards from results/cases/*/evidence_cards",
            )
        for card in evidence_cards:
            card_class = card.get("claim_class")
            if card_class not in _CLASS_RANK:
                raise ClaimRefusal("bad_card", f"card {card.get('claim_id')} has unknown class")
            # ladder rule: a claim may never exceed the strongest supporting card
            if _CLASS_RANK[claim_class] > _CLASS_RANK[card_class]:
                raise ClaimRefusal(
                    "illegal_upgrade",
                    f"claim {claim_class} exceeds evidence card {card['claim_id']} ({card_class})",
                    "claim classes never auto-upgrade; produce independent higher-class evidence first",
                )

        return {
            "claim_class": claim_class,
            "claim_text": claim_text,
            "supported_by": [c["claim_id"] for c in evidence_cards],
            "compiled": True,
        }
