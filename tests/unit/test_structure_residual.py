"""Unit tests for missing-neutral structure residual scoring."""

from plantpersulf.models.structure_residual import score_structure_residual


def test_missing_structure_gets_exact_zero_residual() -> None:
    residual = score_structure_residual(
        weights=(1.0, -1.0),
        intercept=0.5,
        structure_features=[[2.0, 1.0], [0.0, 0.0]],
        quality_mask=[True, False],
    )

    assert residual[1] == 0.0
