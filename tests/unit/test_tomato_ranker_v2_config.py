"""Frozen configuration policy for the tomato v2 workflow."""

from pathlib import Path

import yaml

from plantpersulf.workflows.tomato_ranker_v2 import select_release_model


def test_tomato_v2_config_freezes_primary_arena_and_structure_default() -> None:
    cfg = yaml.safe_load(
        Path("configs/experiments/tomato_ranker_v2.yaml").read_text(encoding="utf-8")
    )
    assert cfg["evaluation"]["primary_arena"] == "panel"
    assert cfg["evaluation"]["split"] == "repeated_homology_cluster_cv"
    assert cfg["evaluation"]["folds"] == 5
    assert cfg["evaluation"]["repetitions"] == 5
    assert cfg["structure"]["delta_default"] == 0
    assert cfg["observation_propensity"]["enabled"] is False
    assert cfg["compute"] == {"device": "auto", "score_batch_size": 16384}
    assert cfg["models"]["fallback_on_admission_failure"] is True
    assert cfg["claim_class"] == "development_candidate_ranking_not_gate2"
    assert select_release_model(False) == "pu_logistic"
    assert select_release_model(True) == "additive_pu"
