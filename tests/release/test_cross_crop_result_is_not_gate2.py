from plantpersulf.evaluation.external_validation import _parse_fold_study


def test_cross_crop_fold_names_are_not_gate2_evidence() -> None:
    assert _parse_fold_study("cross_crop_source_holdout_rice|additive_pu") is None
