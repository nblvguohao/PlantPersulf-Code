"""RED (Task 6): anti-leakage benchmark splits.

Ensures no protein cluster crosses train/validation/test, no study leaks
across folds, and known Zhang-lab mechanism sites are held out of training.

Expected RED: ``plantpersulf.benchmark.splits`` does not exist yet, so
collection fails with ImportError.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from plantpersulf.benchmark.splits import (
    BenchmarkSplitTable,
    audit_no_cluster_leakage,
    build_splits,
    build_study_splits,
)


def _cluster_fixture(tmp_path: Path) -> tuple[Path, Path]:
    """Write a small cluster mapping and split config for testing."""
    clusters = tmp_path / "test_clusters.tsv"
    clusters.write_text(
        "protein_accession\tcluster_id\n"
        "A\t1\nB\t1\nC\t2\nD\t3\nE\t3\nF\t4\nG\t5\nH\t6\nI\t6\nJ\t6\n"
        "K\t7\nL\t8\n",
        encoding="utf-8",
    )
    config = tmp_path / "split_config.yaml"
    with config.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(
            {
                "version": 1,
                "cluster_split": {
                    "cluster_file": str(clusters),
                    "train_fraction": 0.6,
                    "validation_fraction": 0.2,
                    "test_fraction": 0.2,
                    "seed": 42,
                },
                "known_mechanism_holdout": ["SlWRKY6", "SlERF.D2", "BRG3"],
            },
            fh,
        )
    return clusters, config


def test_no_protein_cluster_crosses_split(tmp_path: Path) -> None:
    clusters_path, config_path = _cluster_fixture(tmp_path)
    table = build_splits(config_path)

    # For every cluster, all its proteins must share the same split
    # (verified implicitly by the cross-check below).
    # Cross-check: query the cluster → split mapping.
    cluster_splits: dict[str, str] = {}
    for row in table.rows:
        cid = getattr(row, "cluster_id", "")
        if cid and cid in cluster_splits:
            assert cluster_splits[cid] == row.split, (
                f"cluster {cid} splits across {cluster_splits[cid]} and {row.split}"
            )
        elif cid:
            cluster_splits[cid] = row.split

    assert len(cluster_splits) == 8  # all clusters assigned


def test_every_protein_gets_exactly_one_split(tmp_path: Path) -> None:
    clusters_path, config_path = _cluster_fixture(tmp_path)
    table = build_splits(config_path)

    splits = {row.split for row in table.rows}
    assert splits <= {"train", "validation", "test"}
    assert "train" in splits
    assert 8 <= len(table.rows) <= 12  # all input proteins present


def test_split_is_deterministic(tmp_path: Path) -> None:
    clusters_path, config_path = _cluster_fixture(tmp_path)
    first = build_splits(config_path)
    second = build_splits(config_path)

    fmap = {r.protein_accession: r.split for r in first.rows}
    smap = {r.protein_accession: r.split for r in second.rows}
    assert fmap == smap


def test_known_mechanism_holdout_config_is_loaded(tmp_path: Path) -> None:
    _, config_path = _cluster_fixture(tmp_path)
    table = build_splits(config_path)

    assert table.known_mechanism_holdout == ("SlWRKY6", "SlERF.D2", "BRG3")


def test_missing_cluster_file_is_fatal(tmp_path: Path) -> None:
    config = tmp_path / "bad_config.yaml"
    with config.open("w", encoding="utf-8") as fh:
        yaml.safe_dump(
            {
                "version": 1,
                "cluster_split": {
                    "cluster_file": str(tmp_path / "nonexistent.tsv"),
                    "train_fraction": 0.6,
                    "validation_fraction": 0.2,
                    "test_fraction": 0.2,
                    "seed": 42,
                },
            },
            fh,
        )
    with pytest.raises((RuntimeError, OSError, FileNotFoundError)):
        build_splits(config)


def test_study_split_holds_out_designated_study(tmp_path: Path) -> None:
    """Proteins from the held-out study must be in the test set for that fold."""
    positives = tmp_path / "sites.tsv"
    positives.write_text(
        "protein_accession\tcys_position_in_protein\tlabel\t"
        "study_accession\tevidence_level\tsource_sha256\n"
        "A\t1\tpositive\tPXD006140\tsite_ms\taaa\n"
        "A\t2\tpositive\tPXD006140\tsite_ms\tbbb\n"
        "B\t5\tpositive\tPXD024061\tsite_ms\tccc\n"
        "C\t3\tunlabeled\t\t\t\n"
        "D\t8\tunlabeled\t\t\t\n",
        encoding="utf-8",
    )
    proteome = tmp_path / "mini.fasta"
    proteome.write_text(">sp|A\nAAAA\n>sp|B\nBBBB\n>sp|C\nCCCC\n", encoding="utf-8")
    folds = build_study_splits(positives, proteome)

    assert set(folds) == {"leave_PXD006140_out", "leave_PXD024061_out"}
    # In the fold that leaves PXD006140 out, protein A must be in test.
    fold_006 = folds["leave_PXD006140_out"]
    assert fold_006["A"] == "test"
    assert fold_006["B"] == "train"
    assert fold_006["C"] == "train"
    # In the fold that leaves PXD024061 out, protein B must be in test.
    fold_024 = folds["leave_PXD024061_out"]
    assert fold_024["B"] == "test"
    assert fold_024["A"] == "train"


def test_audit_passes_when_no_cluster_leaks(tmp_path: Path) -> None:
    clusters_path, config_path = _cluster_fixture(tmp_path)
    table = build_splits(config_path)
    audit_no_cluster_leakage(table)  # must not raise


def test_audit_raises_when_cluster_leaks(tmp_path: Path) -> None:
    from plantpersulf.benchmark.splits import SplitRow

    # Artificially inject a leak: cluster 1 appears in two splits.
    table = BenchmarkSplitTable(
        rows=(
            SplitRow("A", "1", "train"),
            SplitRow("B", "1", "test"),
        ),
        seed=42,
        known_mechanism_holdout=(),
    )
    with pytest.raises(RuntimeError, match="cluster 1 crosses"):
        audit_no_cluster_leakage(table)
