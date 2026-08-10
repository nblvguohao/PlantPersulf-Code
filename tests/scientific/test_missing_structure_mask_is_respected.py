"""RED (Task 9): missing structure is masked out, never mean-imputed.

Scientific integrity check: a cysteine with no registered structure must not
have its structure branch silently filled with a mean/placeholder vector that
still influences the score. When the missingness mask is False, the structure
feature values for that row must have *zero* effect on the output — identical
to running with the structure branch disabled entirely.

Uses small synthetic branch features (algorithm-level property, not a
biological claim).

Expected RED: ``plantpersulf.models.structure_ranker`` does not exist yet.
"""

from __future__ import annotations

from plantpersulf.models.structure_ranker import (  # RED: module missing
    BranchFeatures,
    structure_ranker_scores,
)


def _train(
    masked_structure_fill: float = 0.0,
) -> tuple[BranchFeatures, list[str]]:
    n = 24
    sequence = [[float(i % 2), float((i + 1) % 2)] for i in range(n)]
    esm = [[float(i % 2)] * 4 for i in range(n)]
    # Half the training rows have no structure.
    structure_mask = [i % 2 == 0 for i in range(n)]
    structure = [
        [float(i % 2) * 10.0, 50.0 + (i % 2) * 40.0]
        if structure_mask[i]
        else [masked_structure_fill, masked_structure_fill]
        for i in range(n)
    ]
    labels = ["positive" if i % 2 == 0 else "unlabeled" for i in range(n)]
    return (
        BranchFeatures(
            sequence=sequence,
            esm=esm,
            structure=structure,
            structure_mask=structure_mask,
            study_ids=["S"] * n,
        ),
        labels,
    )


def _predict_row(structure_values: list[float], mask: bool) -> BranchFeatures:
    return BranchFeatures(
        sequence=[[1.0, 0.0]],
        esm=[[1.0, 1.0, 1.0, 1.0]],
        structure=[structure_values],
        structure_mask=[mask],
        study_ids=["S"],
    )


def test_masked_structure_values_do_not_change_score() -> None:
    train, train_y = _train()

    # Same row, mask=False, but wildly different structure feature values.
    row_a = _predict_row([9.0, 90.0], mask=False)
    row_b = _predict_row([-999.0, -999.0], mask=False)

    score_a = structure_ranker_scores(train, train_y, row_a, seed=0).scores[0]
    score_b = structure_ranker_scores(train, train_y, row_b, seed=0).scores[0]

    # Masked-out structure is ignored entirely, not averaged in.
    assert score_a == score_b


def test_masked_structure_values_do_not_leak_during_training() -> None:
    """Masked-out structure values must not influence the fitted model either:
    two training sets that differ *only* in the structure values of rows whose
    mask is False must yield identical predictions (no mean-imputation, no
    silent leakage of placeholder values into the structure branch)."""
    train_zeros, train_y = _train(masked_structure_fill=0.0)
    train_junk, _ = _train(masked_structure_fill=-777.0)
    predict = _predict_row([1.0, 60.0], mask=True)

    from_zeros = structure_ranker_scores(train_zeros, train_y, predict, seed=0).scores[
        0
    ]
    from_junk = structure_ranker_scores(train_junk, train_y, predict, seed=0).scores[0]

    assert from_zeros == from_junk


def _neutral_predict_row(structure_values: list[float]) -> BranchFeatures:
    # Ambiguous sequence/ESM (halfway between the two toy clusters) so the
    # head is not saturated and the structure branch's effect is observable.
    return BranchFeatures(
        sequence=[[0.5, 0.5]],
        esm=[[0.5, 0.5, 0.5, 0.5]],
        structure=[structure_values],
        structure_mask=[True],
        study_ids=["S"],
    )


def test_present_structure_values_do_change_score() -> None:
    """Sanity: when the mask is True, the structure branch is actually wired
    in (otherwise the mask test above would pass vacuously)."""
    train, train_y = _train()

    row_a = _neutral_predict_row([10.0, 90.0])
    row_b = _neutral_predict_row([0.0, 50.0])

    score_a = structure_ranker_scores(train, train_y, row_a, seed=0, epochs=40).scores[
        0
    ]
    score_b = structure_ranker_scores(train, train_y, row_b, seed=0, epochs=40).scores[
        0
    ]

    assert score_a != score_b
