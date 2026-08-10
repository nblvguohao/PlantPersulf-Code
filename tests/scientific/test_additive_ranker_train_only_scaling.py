"""No prediction-distribution leakage into additive PU scaling."""

from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)


def test_prediction_distribution_cannot_change_fitted_scaler() -> None:
    model = fit_additive_pu_ranker(
        [[0.0], [1.0], [0.2], [0.8]],
        ["positive", "unlabeled", "positive", "unlabeled"],
        ["group-a", "group-a", "group-b", "group-b"],
        ("numeric_marker",),
        AdditivePuConfig(class_prior=0.5, seed=3, epochs=20),
    )
    frozen = (model.scaler.mean, model.scaler.std)

    model.score([[1_000_000.0], [-1_000_000.0]])

    assert (model.scaler.mean, model.scaler.std) == frozen
