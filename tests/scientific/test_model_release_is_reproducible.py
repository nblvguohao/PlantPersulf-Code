"""RED (Task 9): the ranker is bit-for-bit reproducible given a fixed seed.

A model release must be reproducible from a command with a fixed seed — no
hidden nondeterminism (thread races, unseeded dropout, unstable ordering).
Two runs with identical inputs and seed must produce identical scores and
identical uncertainty estimates.

Expected RED: ``plantpersulf.models.structure_ranker`` does not exist yet.
"""

from __future__ import annotations

from plantpersulf.models.structure_ranker import (  # RED: module missing
    BranchFeatures,
    structure_ranker_scores,
)


def _branches(n: int) -> BranchFeatures:
    sequence = [[float(i % 2), float((i + 1) % 2)] for i in range(n)]
    esm = [[float(i % 3)] * 4 for i in range(n)]
    structure = [[float(i % 2) * 10.0, 40.0 + (i % 2) * 50.0] for i in range(n)]
    structure_mask = [i % 3 != 0 for i in range(n)]
    study_ids = ["PXD000001" if i % 2 else "PXD000002" for i in range(n)]
    return BranchFeatures(
        sequence=sequence,
        esm=esm,
        structure=structure,
        structure_mask=structure_mask,
        study_ids=study_ids,
    )


def _labels(n: int) -> list[str]:
    return ["positive" if i % 2 == 0 else "unlabeled" for i in range(n)]


def test_same_seed_gives_identical_scores() -> None:
    train, train_y, predict = _branches(30), _labels(30), _branches(8)

    first = structure_ranker_scores(train, train_y, predict, seed=7)
    second = structure_ranker_scores(train, train_y, predict, seed=7)

    assert first.scores == second.scores
    assert first.uncertainty == second.uncertainty


def test_different_seed_may_differ_but_stays_finite() -> None:
    train, train_y, predict = _branches(30), _labels(30), _branches(8)

    a = structure_ranker_scores(train, train_y, predict, seed=1)
    b = structure_ranker_scores(train, train_y, predict, seed=2)

    assert len(a.scores) == len(b.scores) == 8
    for score in a.scores + b.scores:
        assert score == score  # not NaN
        assert -1e6 < score < 1e6
