"""Coverage-only gain cannot admit a structure residual arm."""

from plantpersulf.evaluation.structure_admission import admit_structure_residual


def test_coverage_only_gain_blocks_structure_admission() -> None:
    decision = admit_structure_residual(
        paired_structure_deltas=(0.03, 0.02, 0.04, 0.01, 0.03),
        paired_coverage_only_deltas=(0.03, 0.02, 0.04, 0.01, 0.03),
    )

    assert decision.admitted is False
    assert decision.delta == 0
