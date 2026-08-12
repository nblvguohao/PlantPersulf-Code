"""Frozen contract for the multispecies v2 dual evaluation tracks."""

from __future__ import annotations

from pathlib import Path

import yaml


def test_multispecies_v2_config_locks_primary_species_and_two_tracks() -> None:
    """Changing this contract would alter the meaning of reported metrics."""
    path = Path("configs/experiments/multispecies_v2_global_clusters_v2.yaml")
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))

    assert cfg["version"] == 2
    assert cfg["species"]["primary"] == ["arabidopsis", "rice", "tomato"]
    assert cfg["species"]["pressure"] == ["magnaporthe"]
    assert cfg["strict_cluster_holdout"]["test_fraction"] == 0.2
    assert cfg["strict_cluster_holdout"]["development_folds"] == 5
    assert cfg["literature_random_protein"]["seeds"] == list(range(10))
    assert cfg["unlabeled_panel"]["per_positive"] == 20
    assert cfg["data_isolation"] == {
        "split_before_unlabeled_sampling": True,
        "development_fit_scope": "training_fold_only",
        "panther_role": "unlabelled_evolutionary_feature_only",
    }


def test_v4_config_locks_registered_development_positive_manifest() -> None:
    path = Path("configs/experiments/multispecies_v2_global_clusters_v4.yaml")
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))

    assert cfg["development_positive_manifest"] == {
        "path": (
            "data/processed/positives/"
            "multispecies_v2_development_positives_v4.tsv"
        ),
        "sha256": "58444a14d2c7abfcff13c8fcb6124533cb4a4f8b1f2fba38ab68ad2bbaacc296",
    }


def test_v4_config_makes_comparator_input_gaps_explicit() -> None:
    path = Path("configs/experiments/multispecies_v2_global_clusters_v4.yaml")
    cfg = yaml.safe_load(path.read_text(encoding="utf-8"))

    assert cfg["comparison_inputs"] == {
        "esm_features": None,
        "structure_features": None,
        "sul_environment_manifest": None,
        "pcysmod_scores": None,
    }
    assert cfg["comparator_policy"] == {
        "pcysmod": "quantitative_only_with_complete_frozen_scores",
        "tree_graft": "architecture_reference_only",
        "gate2_eligible": False,
    }
    assert cfg["literature_baseline_parameters"] == {
        "pu_logistic": {"holdout_fraction": 0.2},
        "random_forest": {"n_estimators": 100},
        "xgboost": {"n_estimators": 100},
    }
