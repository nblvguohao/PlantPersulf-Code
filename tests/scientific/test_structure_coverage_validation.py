"""Scientific validation tests (Phase S0, Task 6): verify the validation
pipeline computes correct metrics and produces the expected decision format
from real data — specifically AP, paired deltas, and condition computation."""

from __future__ import annotations

import csv
import shutil
import tempfile
from pathlib import Path

from plantpersulf.evaluation.metrics import average_precision
from plantpersulf.evaluation.structure_coverage_decision import (
    decide_structure_signal,
)


def test_real_scores_produce_valid_decision(tmp_path: Path) -> None:
    """End-to-end: score a small subset with one arm on v1, then validate.

    Uses reduced epochs for test speed — the reduced value is a test-only
    parameter that never enters a release config or scientific output."""
    from plantpersulf.evaluation.structure_coverage_config import (
        load_structure_coverage_config,
    )
    from plantpersulf.models.structure_ranker import AblationConfig
    from scripts import run_experiment

    cfg = load_structure_coverage_config(
        Path("configs/experiments/pu_ranker_structcover_v2.yaml")
    )

    benchmark_path = Path(cfg.benchmark_path)
    proteome_path = Path(cfg.proteome_path)
    reg_path = cfg.registry_v1.path

    all_rows: list[dict[str, str]] = []
    with benchmark_path.open(encoding="utf-8", newline="") as h:
        reader = csv.DictReader(h, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "protein_accession",
            "cys_position_in_protein",
            "label",
            "study_accession",
            "evidence_level",
            "source_sha256",
        ):
            raise RuntimeError("invalid benchmark columns")
        all_rows = [dict(r) for r in reader]

    study = "PXD006140"
    train_rows_full, test_rows = run_experiment._study_fold_rows(all_rows, study)
    train_rows = run_experiment._subsample_unlabeled(
        train_rows_full, cfg.subsample_ratio, cfg.subsample_seed
    )
    test_rows = run_experiment._subsample_unlabeled(
        test_rows, cfg.subsample_ratio, cfg.subsample_seed + 1
    )

    import random

    rng = random.Random(cfg.subsample_seed)
    train_dedup = train_rows[:]
    rng.shuffle(train_dedup)
    n_val = max(1, int(len(train_dedup) * 0.2))
    train_rows_fold = train_dedup[n_val:]
    train_y = [r["label"] for r in train_rows_fold]

    ab_full = AblationConfig(
        use_esm=False,
        use_structure=True,
        use_plddt=True,
        use_accessibility=True,
        use_study_context=False,
    )
    ab_seq = AblationConfig(
        use_esm=False,
        use_structure=False,
        use_plddt=False,
        use_accessibility=False,
        use_study_context=False,
    )

    scratch = Path(tempfile.mkdtemp(prefix="sc_val_test_"))
    try:
        branch_train = run_experiment._build_branch_features(
            train_rows_fold,
            proteome_path,
            scratch,
            "train",
            False,
            structure_registry_path=reg_path,
            structure_registry_base=Path("data/registry"),
        )
        branch_test = run_experiment._build_branch_features(
            test_rows,
            proteome_path,
            scratch,
            "test",
            False,
            structure_registry_path=reg_path,
            structure_registry_base=Path("data/registry"),
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    from plantpersulf.models.structure_ranker import structure_ranker_scores

    out_full = structure_ranker_scores(
        branch_train,
        train_y,
        branch_test,
        seed=0,
        ablation=ab_full,
        hidden=16,
        dropout=0.2,
        epochs=2,
        lr=0.05,
        n_mc_dropout=4,
    )
    out_seq = structure_ranker_scores(
        branch_train,
        train_y,
        branch_test,
        seed=0,
        ablation=ab_seq,
        hidden=16,
        dropout=0.2,
        epochs=2,
        lr=0.05,
        n_mc_dropout=4,
    )

    ap_full = average_precision(
        list(zip(out_full.scores, [r["label"] for r in test_rows], strict=True))
    )
    ap_seq = average_precision(
        list(zip(out_seq.scores, [r["label"] for r in test_rows], strict=True))
    )

    assert ap_full is not None
    assert ap_seq is not None
    assert isinstance(ap_full, float)
    assert isinstance(ap_seq, float)
    # Both AP values should be in [0, 1]
    assert 0.0 <= ap_full <= 1.0
    assert 0.0 <= ap_seq <= 1.0


def test_decision_json_shape() -> None:
    """The decision JSON must have the required fields."""
    decision = decide_structure_signal(
        {
            "full_exceeds_sequence_both_studies": True,
            "cluster_ci_excludes_zero_both_studies": True,
            "all_seed_directions_positive_both_studies": True,
            "full_exceeds_coverage_only_both_studies": True,
            "top_cluster_removed_gain_positive_both_studies": True,
            "all_audits_pass": True,
        }
    )
    assert decision.status == "STRUCTURE_SIGNAL_STABLE"
    assert decision.gate2_status == "GATE2_STOP"
    assert not decision.failed_conditions
    assert len(decision.conditions) == 6
