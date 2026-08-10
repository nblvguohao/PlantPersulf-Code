"""Device-aware execution remains explicit and scoreable."""

from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)
from plantpersulf.workflows.tomato_ranker_v2 import deduplicated_additive_seed_configs


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


def test_redundant_deterministic_seed_runs_are_collapsed_per_prior() -> None:
    assert deduplicated_additive_seed_configs((0.1, 0.2), (0, 1, 2)) == (
        (0, 0.1),
        (0, 0.2),
    )


def test_additive_model_seed_does_not_change_the_fixed_training_trajectory() -> None:
    features = [[0.0, 1.0], [1.0, 0.0], [0.2, 0.8], [0.8, 0.2]]
    common = (
        features,
        ["positive", "unlabeled", "positive", "unlabeled"],
        ["group-a", "group-a", "group-b", "group-b"],
        ("f1", "f2"),
    )
    first = fit_additive_pu_ranker(
        *common, AdditivePuConfig(class_prior=0.5, seed=0, epochs=5, device="cpu")
    )
    second = fit_additive_pu_ranker(
        *common, AdditivePuConfig(class_prior=0.5, seed=4, epochs=5, device="cpu")
    )
    assert first.to_dict() == second.to_dict()
