"""RED (Phase S0, Task 6): fail-closed six-condition structure-signal decision.

The decision requires ALL six conditions to pass; any failure yields
STRUCTURE_SIGNAL_UNSTABLE. No condition can be omitted or guessed."""

from __future__ import annotations

import pytest

from plantpersulf.evaluation.structure_coverage_decision import (
    REQUIRED_CONDITIONS,
    decide_structure_signal,
)


def test_decision_requires_every_condition() -> None:
    passed: dict[str, bool] = {
        "full_exceeds_sequence_both_studies": True,
        "cluster_ci_excludes_zero_both_studies": True,
        "all_seed_directions_positive_both_studies": True,
        "full_exceeds_coverage_only_both_studies": True,
        "top_cluster_removed_gain_positive_both_studies": True,
        "all_audits_pass": True,
    }
    assert decide_structure_signal(passed).status == "STRUCTURE_SIGNAL_STABLE"
    for condition in passed:
        failed = {**passed, condition: False}
        decision = decide_structure_signal(failed)
        assert decision.status == "STRUCTURE_SIGNAL_UNSTABLE"
        assert condition in decision.failed_conditions


def test_decision_rejects_unknown_condition() -> None:
    with pytest.raises(ValueError, match="unknown condition"):
        decide_structure_signal({"fake_condition": True})


def test_decision_rejects_missing_condition() -> None:
    incomplete = {
        "full_exceeds_sequence_both_studies": True,
        "cluster_ci_excludes_zero_both_studies": True,
    }
    with pytest.raises(ValueError, match="missing"):
        decide_structure_signal(incomplete)


def test_required_conditions_are_exactly_six() -> None:
    assert len(REQUIRED_CONDITIONS) == 6
    assert "all_audits_pass" in REQUIRED_CONDITIONS
    assert "full_exceeds_sequence_both_studies" in REQUIRED_CONDITIONS
    assert "cluster_ci_excludes_zero_both_studies" in REQUIRED_CONDITIONS
    assert "all_seed_directions_positive_both_studies" in REQUIRED_CONDITIONS
    assert "full_exceeds_coverage_only_both_studies" in REQUIRED_CONDITIONS
    assert "top_cluster_removed_gain_positive_both_studies" in REQUIRED_CONDITIONS


def test_decision_gate2_status_is_always_stop() -> None:
    decision = decide_structure_signal(
        {
            "full_exceeds_sequence_both_studies": True,
            "cluster_ci_excludes_zero_both_studies": True,
            "all_seed_directions_positive_both_studies": True,
            "full_exceeds_coverage_only_both_studies": True,
            "top_cluster_removed_gain_positive_both_studies": True,
            "all_audits_pass": True,
        }
    )
    assert decision.gate2_status == "GATE2_STOP"


def test_decision_fields_are_read_only() -> None:
    decision = decide_structure_signal(
        {
            "full_exceeds_sequence_both_studies": False,
            "cluster_ci_excludes_zero_both_studies": False,
            "all_seed_directions_positive_both_studies": False,
            "full_exceeds_coverage_only_both_studies": False,
            "top_cluster_removed_gain_positive_both_studies": False,
            "all_audits_pass": False,
        }
    )
    assert decision.status == "STRUCTURE_SIGNAL_UNSTABLE"
    assert set(decision.failed_conditions) == set(REQUIRED_CONDITIONS)
