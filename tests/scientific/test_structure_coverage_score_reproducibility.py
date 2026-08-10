"""Scientific reproducibility test (Phase S0, Task 5): two scoring runs with
the same inputs and seed produce byte-identical score tables."""

from __future__ import annotations

import csv
import shutil
import tempfile
from pathlib import Path

import pytest

from plantpersulf.evaluation.structure_coverage_config import (
    load_structure_coverage_config,
)
from plantpersulf.models.structure_ranker import AblationConfig
from scripts import run_experiment


def _read_scored(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as h:
        return [dict(r) for r in csv.DictReader(h, delimiter="\t")]


@pytest.mark.slow
def test_single_fold_seed_arm_is_reproducible(tmp_path: Path) -> None:
    """Two runs of the same fold/seed/arm on the same release produce
    bit-identical scores. Uses reduced epochs for test speed — the reduced
    value is a test-only parameter that never enters a release config."""

    cfg = load_structure_coverage_config(
        Path("configs/experiments/pu_ranker_structcover_v2.yaml")
    )

    benchmark_path = Path(cfg.benchmark_path)
    proteome_path = Path(cfg.proteome_path)
    reg_path = cfg.registry_v1.path
    reg_base = Path("data/registry")

    all_rows: list[dict[str, str]] = []
    with benchmark_path.open(encoding="utf-8", newline="") as h:
        reader = csv.DictReader(h, delimiter="\t")
        if tuple(reader.fieldnames or ()) != (
            "protein_accession", "cys_position_in_protein", "label",
            "study_accession", "evidence_level", "source_sha256",
        ):
            raise RuntimeError("invalid benchmark columns")
        all_rows = [dict(r) for r in reader]

    study = cfg.studies[0]
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

    arm_cfg = cfg.arms["sequence_only"]
    ab = AblationConfig(
        use_esm=bool(arm_cfg.get("use_esm", True)),
        use_structure=bool(arm_cfg.get("use_structure", True)),
        use_plddt=bool(arm_cfg.get("use_plddt", True)),
        use_accessibility=bool(arm_cfg.get("use_accessibility", True)),
        use_study_context=bool(arm_cfg.get("use_study_context", True)),
    )

    def _run_one() -> list[float]:
        scratch = Path(tempfile.mkdtemp(prefix="sc_repro_"))
        try:
            branch_train = run_experiment._build_branch_features(
                train_rows_fold, proteome_path, scratch, "train", False,
                structure_registry_path=reg_path,
                structure_registry_base=reg_base,
            )
            branch_test = run_experiment._build_branch_features(
                test_rows, proteome_path, scratch, "test", False,
                structure_registry_path=reg_path,
                structure_registry_base=reg_base,
            )
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

        from plantpersulf.models.structure_ranker import structure_ranker_scores

        out = structure_ranker_scores(
            branch_train, train_y, branch_test,
            seed=0, ablation=ab,
            hidden=16, dropout=0.2,
            epochs=2,  # test-only reduction
            lr=0.05, n_mc_dropout=4,
        )
        return list(out.scores)

    scores_a = _run_one()
    scores_b = _run_one()

    assert scores_a == scores_b
    assert len(scores_a) == len(test_rows)
