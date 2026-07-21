"""RED (Task 10 conclusion gate): honesty enforcement.

Two release-level guarantees:

* The Gate-2 conclusion is a mechanical function of the five Codex conditions,
  and with only two same-laboratory studies (the current data) it must resolve
  to ``GATE2_STOP`` — condition 1 requires *independent* held-out studies.
* When the decision is STOP, no outward-facing text may claim predictive value;
  the verifier flags any such phrase. When GO, the same text is permitted.

These tests use registered evidence dicts and policy text only — no fictional
sites.

Expected RED: ``plantpersulf.evaluation.conclusion_gate`` does not exist yet.
"""

from __future__ import annotations

from plantpersulf.evaluation.conclusion_gate import (  # RED: module missing
    DOWNGRADE_STATEMENT,
    evaluate_gate2,
    verify_predictive_claims,
)

# thresholds mirror configs/gate2_v1.yaml
THRESHOLDS = {
    "min_independent_studies": 2,
    "min_ap_margin": 0.0,
    "max_permutation_p": 0.05,
    "min_structure_gain": 0.0,
}


def _two_samelab_studies_evidence() -> dict[str, object]:
    """The real situation: two studies, same lab, full model >= baseline on
    both numerically, but the studies are NOT independent."""
    return {
        "studies_are_independent": False,
        "per_study": [
            {"study": "PXD006140", "full_ap": 0.05, "baseline_ap": 0.02},
            {"study": "PXD024061", "full_ap": 0.04, "baseline_ap": 0.02},
        ],
        "delta_ci_lower": 0.01,
        "permutation_p": 0.2,
        "control_leakage": [],
        "independent_units": 2,
        "structure_gain": None,
        "single_cluster_driven": True,
    }


def test_two_samelab_studies_resolve_to_stop() -> None:
    decision = evaluate_gate2(_two_samelab_studies_evidence(), THRESHOLDS)
    assert decision.decision == "GATE2_STOP"
    # condition 1 (independent studies) must be among the failures
    failed = {c.name for c in decision.conditions if not c.passed}
    assert "independent_studies_beat_baseline" in failed


def test_stop_forbids_predictive_language() -> None:
    decision = evaluate_gate2(_two_samelab_studies_evidence(), THRESHOLDS)
    text = "Our model has predictive value for cross-study persulfidation."
    violations = verify_predictive_claims(text, decision)
    assert violations  # flagged


def test_stop_allows_downgraded_language() -> None:
    decision = evaluate_gate2(_two_samelab_studies_evidence(), THRESHOLDS)
    assert verify_predictive_claims(DOWNGRADE_STATEMENT, decision) == []


def test_go_allows_predictive_language() -> None:
    ev = {
        "studies_are_independent": True,
        "per_study": [
            {"study": "S1", "full_ap": 0.3, "baseline_ap": 0.1},
            {"study": "S2", "full_ap": 0.25, "baseline_ap": 0.1},
        ],
        "delta_ci_lower": 0.05,
        "permutation_p": 0.001,
        "control_leakage": [],
        "independent_units": 2,
        "structure_gain": 0.04,
        "single_cluster_driven": False,
    }
    decision = evaluate_gate2(ev, THRESHOLDS)
    assert decision.decision == "GATE2_GO"
    text = "The model shows predictive value on held-out studies."
    assert verify_predictive_claims(text, decision) == []
