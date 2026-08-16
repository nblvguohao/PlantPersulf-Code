"""RED/lock-in (Phase F honesty discipline): Gate 2 must never ingest
within-dataset (Split A / cluster-split) metrics as cross-study evidence.

Published cysteine-PTM predictors are commonly benchmarked with a
within-integrated-dataset split (see docs/superpowers/specs/2026-07-21-
phase-f-pu-ranker-design.md for the literature basis). ``pu_ranker_cluster_v1``
exists to give a literature-comparable number under that regime, but Gate 2
specifically requires *independent held-out studies* (Codex 10.结论闸门
condition 1). This test locks in the structural guarantee that a cluster-split
experiment's model names are never recognised as a leave-study-out fold, so
they can never silently satisfy Gate 2 evidence even if someone points
``--model-release`` at the wrong experiment by mistake.
"""

from __future__ import annotations

from scripts.validate_external import _collect_fold_metrics, _parse_fold_study


def test_leave_study_out_model_name_is_recognised() -> None:
    assert _parse_fold_study("leave_PXD006140_out|structure_ranker:full") == (
        "PXD006140"
    )


def test_cluster_split_model_name_is_not_recognised_as_a_fold() -> None:
    # This is exactly the model-name shape pu_ranker_cluster_v1 produces.
    assert _parse_fold_study("pu_ranker_cluster_v1|structure_ranker:full") is None
    assert _parse_fold_study("structure_ranker:full") is None


def test_protein_split_model_name_is_not_recognised_as_a_fold() -> None:
    # Model-name shape pu_ranker_protein_split_v1 produces (one tag per
    # repetition seed, literature-comparable regime).
    assert _parse_fold_study("protein_split_seed0|structure_ranker:full") is None
    assert _parse_fold_study("protein_split_seed9|structure_ranker:full") is None


def test_cluster_cv_model_name_is_not_recognised_as_a_fold() -> None:
    # Model-name shape pu_ranker_cluster_cv_v1 produces (development-
    # stability track, repeated cluster K-fold CV).
    assert _parse_fold_study("cluster_cv_r0f0|structure_ranker:full") is None
    assert _parse_fold_study("cluster_cv_r4f4|structure_ranker:full") is None


def test_multispecies_v2_track_tags_are_not_recognised_as_a_fold() -> None:
    # Task 9.6: neither multispecies v2 track (strict homology-cluster split
    # nor the literature-comparable random-protein split) may silently
    # satisfy Gate 2 evidence. This locks in the existing prefix guarantee
    # rather than adding a parallel isolation mechanism — both tags already
    # fail the ``leave_<study>_out`` prefix check.
    assert (
        _parse_fold_study("strict_cluster_holdout|structure_ranker:full") is None
    )
    assert (
        _parse_fold_study("literature_random_protein|structure_ranker:full")
        is None
    )


def test_cluster_split_rows_never_enter_gate2_fold_metrics() -> None:
    cluster_rows = [
        {
            "model": "pu_ranker_cluster_v1|structure_ranker:full",
            "seed": "0",
            "test_ap": "0.9",
        },
        {
            "model": "pu_ranker_cluster_v1|structure_ranker:full",
            "seed": "1",
            "test_ap": "0.88",
        },
    ]
    leave_study_out_rows = [
        {
            "model": "leave_PXD006140_out|structure_ranker:full",
            "seed": "0",
            "test_ap": "0.03",
        },
    ]

    folds = _collect_fold_metrics(
        release_rows=leave_study_out_rows,
        baseline_rows=leave_study_out_rows,
        model_tag="structure_ranker:full",
        baseline_tag="structure_ranker:full",
    )
    # Only the genuine leave-study-out row produced a fold; the (deliberately
    # inflated, 0.9-AP) cluster-split rows above never got the chance to.
    assert len(folds) == 1
    assert folds[0].holdout_study == "PXD006140"

    # And directly: feeding ONLY cluster-split rows yields no folds at all.
    cluster_only_folds = _collect_fold_metrics(
        release_rows=cluster_rows,
        baseline_rows=cluster_rows,
        model_tag="structure_ranker:full",
        baseline_tag="structure_ranker:full",
    )
    assert cluster_only_folds == []
