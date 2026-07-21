"""Task 10 — Gate 2 conclusion gate and outward-claim honesty enforcement.

Gate 2 asks whether the evidence supports a *cross-study predictive* claim. Per
the Codex conclusion gate, "具有预测价值 / predictive value" may be written only
if **all five** conditions hold:

1. the full model beats the no-learning baseline on >=2 **independent** held-out
   studies;
2. the 95% CI of the effect does not fully cross zero;
3. known-mechanism recovery is not a training-leakage artefact;
4. structure augmentation yields an interpretable gain on the structured /
   sufficient-pLDDT subset;
5. results are not driven by a single high-homology protein cluster.

Otherwise the mandated downgrade statement applies. ``evaluate_gate2`` is a pure
mechanical function of an evidence dict; it never "rounds up". With the current
data — two studies from the *same laboratory* — condition 1 fails on
independence grounds, so the honest decision is ``GATE2_STOP``.

``verify_predictive_claims`` and ``find_unsupported_claims`` enforce that no
outward text claims predictive value under STOP and that every claim cites an
existing result table.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

GATE2_GO = "GATE2_GO"
GATE2_STOP = "GATE2_STOP"

# Codex 10.结论闸门 mandated downgrade wording (zh + en).
DOWNGRADE_STATEMENT = (
    "当前公开数据不足以证明跨研究预测能力，模型仅用于候选组织与假设生成。 "
    "Current public data are insufficient to demonstrate cross-study predictive "
    "ability; the model is used only for candidate organisation and hypothesis "
    "generation."
)

# Phrases that assert cross-study predictive value (forbidden under STOP).
_PREDICTIVE_CLAIM_PATTERNS = (
    "predictive value",
    "predicts persulfidation",
    "predict persulfidation",
    "predicts novel",
    "generalis",  # generalise / generalisation
    "generaliz",  # generalize / generalization
    "具有预测价值",
    "预测价值",
    "可预测",
)


@dataclass(frozen=True)
class GateCondition:
    name: str
    passed: bool
    detail: str


@dataclass(frozen=True)
class Gate2Decision:
    decision: str
    conditions: list[GateCondition]
    downgrade_statement: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "decision": self.decision,
            "conditions": [
                {"name": c.name, "passed": c.passed, "detail": c.detail}
                for c in self.conditions
            ],
            "downgrade_statement": (
                "" if self.decision == GATE2_GO else self.downgrade_statement
            ),
        }


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else float("nan")


def evaluate_gate2(
    evidence: dict[str, Any],
    thresholds: dict[str, Any],
) -> Gate2Decision:
    """Mechanically evaluate the five Codex conditions from an evidence dict."""
    min_studies = int(thresholds.get("min_independent_studies", 2))
    min_margin = float(thresholds.get("min_ap_margin", 0.0))
    max_p = float(thresholds.get("max_permutation_p", 0.05))
    min_gain = float(thresholds.get("min_structure_gain", 0.0))

    per_study = list(evidence.get("per_study", []))
    studies_independent = bool(evidence.get("studies_are_independent", False))
    beats = sum(
        1
        for s in per_study
        if float(s["full_ap"]) - float(s["baseline_ap"]) > min_margin
    )
    cond1 = GateCondition(
        name="independent_studies_beat_baseline",
        passed=studies_independent and beats >= min_studies,
        detail=(
            f"studies_independent={studies_independent}, "
            f"beats_baseline={beats}/{len(per_study)}, need>={min_studies}"
        ),
    )

    ci_lower = evidence.get("delta_ci_lower")
    cond2 = GateCondition(
        name="effect_ci_excludes_zero",
        passed=ci_lower is not None and float(ci_lower) > 0.0,
        detail=f"delta_ci_lower={ci_lower}",
    )

    leakage = list(evidence.get("control_leakage", []))
    units = int(evidence.get("independent_units", 0))
    cond3 = GateCondition(
        name="recovery_is_not_training_leakage",
        passed=not leakage and units >= 1,
        detail=f"control_leakage={leakage}, independent_units={units}",
    )

    gain = evidence.get("structure_gain")
    cond4 = GateCondition(
        name="structure_gain_on_structured_subset",
        passed=gain is not None and float(gain) > min_gain,
        detail=f"structure_gain={gain}, need>{min_gain}",
    )

    single_cluster = bool(evidence.get("single_cluster_driven", True))
    perm_p = evidence.get("permutation_p")
    cond5 = GateCondition(
        name="not_driven_by_single_cluster",
        passed=(
            not single_cluster
            and perm_p is not None
            and float(perm_p) <= max_p
        ),
        detail=f"single_cluster_driven={single_cluster}, permutation_p={perm_p}",
    )

    conditions = [cond1, cond2, cond3, cond4, cond5]
    decision = GATE2_GO if all(c.passed for c in conditions) else GATE2_STOP
    return Gate2Decision(
        decision=decision,
        conditions=conditions,
        downgrade_statement=DOWNGRADE_STATEMENT,
    )


def verify_predictive_claims(text: str, decision: Gate2Decision) -> list[str]:
    """Return the predictive-claim phrases present in ``text`` that are
    forbidden because the decision is STOP. Empty when the decision is GO or no
    forbidden phrase appears."""
    if decision.decision == GATE2_GO:
        return []
    lowered = text.lower()
    return [
        pattern
        for pattern in _PREDICTIVE_CLAIM_PATTERNS
        if pattern.lower() in lowered
    ]


def find_unsupported_claims(
    claims: list[dict[str, Any]],
    available_tables: set[str],
) -> list[dict[str, Any]]:
    """Return claims whose cited ``table`` is missing (or absent) from the set
    of shipped result tables."""
    return [
        claim
        for claim in claims
        if not claim.get("table") or claim["table"] not in available_tables
    ]
