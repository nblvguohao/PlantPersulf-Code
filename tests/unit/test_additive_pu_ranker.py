"""Determinism and serialization tests for the additive PU ranker."""

from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)


def test_same_seed_produces_identical_additive_scores() -> None:
    numeric_rows = [[0.0, 1.0], [1.0, 0.0], [0.2, 0.8], [0.8, 0.2]]
    labels = ["positive", "unlabeled", "positive", "unlabeled"]
    proteins = ["group-a", "group-a", "group-b", "group-b"]
    config = AdditivePuConfig(class_prior=0.5, seed=7, epochs=50)
    first = fit_additive_pu_ranker(
        numeric_rows, labels, proteins, ("f1", "f2"), config
    )
    second = fit_additive_pu_ranker(
        numeric_rows, labels, proteins, ("f1", "f2"), config
    )

    assert first.score(numeric_rows) == second.score(numeric_rows)
    restored = type(first).from_dict(first.to_dict())
    assert restored.score(numeric_rows) == first.score(numeric_rows)
