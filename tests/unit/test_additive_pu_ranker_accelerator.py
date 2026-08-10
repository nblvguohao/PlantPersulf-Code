"""Device-aware execution remains explicit and scoreable."""

from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)


def test_explicit_cpu_device_and_batched_scoring_are_supported() -> None:
    features = [[0.0, 1.0], [1.0, 0.0], [0.2, 0.8], [0.8, 0.2]]
    model = fit_additive_pu_ranker(
        features,
        ["positive", "unlabeled", "positive", "unlabeled"],
        ["group-a", "group-a", "group-b", "group-b"],
        ("f1", "f2"),
        AdditivePuConfig(class_prior=0.5, seed=7, epochs=2, device="cpu"),
    )
    assert len(model.score(features, device="cpu", batch_size=2)) == len(features)
