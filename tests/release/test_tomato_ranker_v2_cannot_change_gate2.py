"""Tomato development results are structurally outside Gate 2."""

from plantpersulf.evaluation.external_validation import _parse_fold_study


def test_tomato_v2_fold_names_are_not_gate2_evidence() -> None:
    assert _parse_fold_study("tomato_v2_r0f0|additive_pu") is None
