"""Split isolation and test-unlock policy for the v2 multispecies runner."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from plantpersulf.benchmark.multispecies_splits import (
    GlobalClusterRow,
    MultispeciesSiteRow,
    build_frozen_multispecies_split,
    global_protein_id,
    write_frozen_multispecies_split,
    write_global_cluster_table,
)
from plantpersulf.workflows.multispecies_v2 import (
    assert_test_unlocked,
    development_fold_rows,
    random_protein_train_validation_test,
    run_multispecies_experiment,
    subsample_unlabeled_by_partition,
    validate_reference_proteome_inputs,
)


def _cluster(i: int) -> GlobalClusterRow:
    species = "arabidopsis" if i % 2 == 0 else "tomato"
    accession = f"P{i:02d}"
    return GlobalClusterRow(
        species=species,
        protein_accession=accession,
        global_protein_id=global_protein_id(species, accession),
        cluster_id=f"C{i}",
        identity_threshold=0.3,
        source_proteome_sha256="b" * 64,
    )


def _rows(n: int = 10, offset: int = 0) -> list[MultispeciesSiteRow]:
    return [
        MultispeciesSiteRow(
            species="arabidopsis" if (i + offset) % 2 == 0 else "tomato",
            protein_accession=f"P{i + offset:02d}",
            cys_position=position,
            label="positive" if position == 5 else "unlabeled",
            study_accession="STUDY_A" if (i + offset) % 2 == 0 else "STUDY_T",
        )
        for i in range(n)
        for position in (5, 9)
    ]


def test_random_protein_track_keeps_validation_proteins_out_of_training() -> None:
    """Replacing protein grouping with row shuffling would leak position 9."""
    partition = random_protein_train_validation_test(_rows(), seed=3)

    train = {row.global_protein_id for row in partition["train"]}
    validation = {row.global_protein_id for row in partition["validation"]}
    test = {row.global_protein_id for row in partition["test"]}
    assert len(test) == 2
    assert len(validation) == 2
    assert len(train) == 6
    assert train.isdisjoint(validation)
    assert train.isdisjoint(test)
    assert validation.isdisjoint(test)


def test_development_fold_never_exposes_frozen_test_rows() -> None:
    """Adding test rows to a development fold would enable test tuning."""
    clusters = [_cluster(i) for i in range(10)]
    rows = _rows()
    frozen = build_frozen_multispecies_split(clusters, rows, seed=20260811)

    train, validation = development_fold_rows(rows, frozen, validation_fold=0)

    frozen_test_ids = {
        row.global_protein_id for row in frozen.rows if row.split == "test"
    }
    assert frozen_test_ids.isdisjoint({row.global_protein_id for row in train})
    assert frozen_test_ids.isdisjoint({row.global_protein_id for row in validation})


def test_test_unlock_requires_matching_config_code_and_split_hash(
    tmp_path: Path,
) -> None:
    """Changing a config after approval must invalidate the test unlock."""
    config = tmp_path / "config.yaml"
    config.write_text("version: 2\nname: v2\n", encoding="utf-8")
    config_hash = hashlib.sha256(config.read_bytes()).hexdigest()
    unlock = tmp_path / "test_unlock.json"
    unlock.write_text(
        json.dumps(
            {
                "config_sha256": config_hash,
                "code_revision": "abc123",
                "split_sha256": "c" * 64,
            }
        ),
        encoding="utf-8",
    )

    assert_test_unlocked(config, "abc123", "c" * 64, unlock)
    with pytest.raises(RuntimeError, match="test unlock mismatch"):
        assert_test_unlocked(config, "changed", "c" * 64, unlock)


def test_runner_keeps_test_locked_by_default(tmp_path: Path) -> None:
    """A default test-scoring branch would allow routine hyperparameter tuning."""
    clusters = [_cluster(i) for i in range(10)]
    rows = _rows()
    frozen = build_frozen_multispecies_split(clusters, rows, seed=20260811)
    cluster_path = tmp_path / "clusters.tsv"
    split_path = tmp_path / "split.tsv"
    reference_path = tmp_path / "reference.fasta"
    reference_path.write_text(">marker\nMARKER\n", encoding="utf-8")
    reference_sha256 = hashlib.sha256(reference_path.read_bytes()).hexdigest()
    write_global_cluster_table(cluster_path, clusters)
    write_frozen_multispecies_split(split_path, frozen)
    config = tmp_path / "config.yaml"
    config.write_text(
        "version: 2\n"
        "strict_cluster_holdout:\n"
        f"  cluster_table: {cluster_path.as_posix()}\n"
        f"  split_path: {split_path.as_posix()}\n"
        "  test_unlock_required: true\n"
        "reference_proteomes:\n"
        "  - species: marker\n"
        f"    path: {reference_path.as_posix()}\n"
        f"    sha256: {reference_sha256}\n",
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="data_isolation configuration is required"):
        run_multispecies_experiment(config, code_revision="abc123")

    config.write_text(
        config.read_text(encoding="utf-8")
        + "data_isolation:\n"
        + "  split_before_unlabeled_sampling: true\n"
        + "  development_fit_scope: training_fold_only\n"
        + "  panther_role: unlabelled_evolutionary_feature_only\n",
        encoding="utf-8",
    )
    prepared = run_multispecies_experiment(config, code_revision="abc123")

    assert prepared.test_scoring_enabled is False
    assert prepared.split.sha256 == frozen.sha256


def test_panel_sampling_happens_after_partition_without_row_reuse() -> None:
    """Sampling before the split could place one unlabeled Cys in two arenas."""
    partition = {
        "train": _rows(3),
        "validation": _rows(3, offset=3),
        "test": _rows(4, offset=6),
    }
    sampled = subsample_unlabeled_by_partition(partition, per_positive=1, seed=7)

    all_keys = [
        (row.global_protein_id, row.cys_position)
        for rows in sampled.values()
        for row in rows
        if row.label == "unlabeled"
    ]
    assert len(all_keys) == len(set(all_keys))
    assert all(
        sum(row.label == "unlabeled" for row in rows)
        <= sum(row.label == "positive" for row in rows)
        for rows in sampled.values()
    )


def test_reference_proteome_validation_rejects_mismatched_sha256(
    tmp_path: Path,
) -> None:
    """The CLI must not enumerate Cys from an unverified reference FASTA."""
    fasta = tmp_path / "reference.fasta"
    fasta.write_text(">marker\nMARKER\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="reference proteome SHA256 mismatch"):
        validate_reference_proteome_inputs(
            {
                "reference_proteomes": [
                    {"species": "marker", "path": str(fasta), "sha256": "0" * 64}
                ]
            }
        )
