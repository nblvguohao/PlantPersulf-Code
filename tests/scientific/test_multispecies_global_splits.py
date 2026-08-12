"""Software-policy tests for the v2 global homology split.

The marker accessions in this file are non-biological test identifiers.  They
exercise split mechanics only and never enter scientific outputs.
"""

from __future__ import annotations

from pathlib import Path
from time import monotonic

import pytest

from plantpersulf.benchmark.multispecies_splits import (
    GlobalClusterRow,
    MultispeciesSiteRow,
    audit_frozen_multispecies_split,
    build_frozen_multispecies_split,
    global_protein_id,
    load_frozen_multispecies_split,
    load_global_cluster_table,
    partition_rows_by_frozen_split,
    write_frozen_multispecies_split,
    write_global_cluster_table,
)


def _cluster(
    species: str, accession: str, cluster: str
) -> GlobalClusterRow:
    return GlobalClusterRow(
        species=species,
        protein_accession=accession,
        global_protein_id=global_protein_id(species, accession),
        cluster_id=cluster,
        identity_threshold=0.3,
        source_proteome_sha256="a" * 64,
    )


def _site(species: str, accession: str, study: str) -> MultispeciesSiteRow:
    return MultispeciesSiteRow(
        species=species,
        protein_accession=accession,
        cys_position=7,
        label="positive",
        study_accession=study,
    )


def test_global_protein_key_prevents_cross_species_accession_collision() -> None:
    """Dropping species from the key would collapse these two distinct rows."""
    rows = [
        _cluster("arabidopsis", "P_SHARED", "C1"),
        _cluster("tomato", "P_SHARED", "C2"),
        _cluster("rice", "P_OTHER", "C3"),
    ]
    sites = [
        _site("arabidopsis", "P_SHARED", "ATH"),
        _site("tomato", "P_SHARED", "TOM"),
    ]

    split = build_frozen_multispecies_split(
        rows, sites, seed=20260811, n_development_folds=2
    )

    assert {
        "arabidopsis|P_SHARED",
        "tomato|P_SHARED",
    }.issubset({row.global_protein_id for row in split.rows})


def test_frozen_split_keeps_each_cluster_in_one_partition() -> None:
    """Assigning members of C1 to different partitions is homology leakage."""
    clusters = []
    sites = []
    for i in range(10):
        species = "arabidopsis" if i % 2 == 0 else "tomato"
        cluster = f"C{i}"
        for member in ("A", "B"):
            accession = f"{species}_{i}_{member}"
            clusters.append(_cluster(species, accession, cluster))
            sites.append(_site(species, accession, f"STUDY_{species}"))

    split = build_frozen_multispecies_split(clusters, sites, seed=20260811)

    audit_frozen_multispecies_split(split)
    by_cluster: dict[str, set[str]] = {}
    for row in split.rows:
        by_cluster.setdefault(row.cluster_id, set()).add(row.split)
    assert all(len(partitions) == 1 for partitions in by_cluster.values())
    assert sum(row.split == "test" for row in split.rows) == 4


def test_partition_rejects_site_without_registered_cluster() -> None:
    """A fallback-to-train branch would silently leak an unclustered protein."""
    split = build_frozen_multispecies_split(
        [
            _cluster("arabidopsis", "P1", "C1"),
            _cluster("arabidopsis", "P2", "C2"),
            _cluster("arabidopsis", "P3", "C3"),
        ],
        [],
        seed=20260811,
        n_development_folds=2,
    )
    missing = _site("arabidopsis", "NOT_REGISTERED", "ATH")

    with pytest.raises(RuntimeError, match="missing from frozen split"):
        partition_rows_by_frozen_split([missing], split)


def test_global_cluster_table_round_trip_preserves_provenance(tmp_path: Path) -> None:
    """A changed TSV header or omitted source hash would break provenance."""
    path = tmp_path / "global_clusters.tsv"
    expected = [
        _cluster("arabidopsis", "P1", "C1"),
        _cluster("tomato", "P1", "C2"),
    ]

    write_global_cluster_table(path, expected)

    assert load_global_cluster_table(path) == tuple(expected)


def test_frozen_split_writer_refuses_to_overwrite_existing_boundary(
    tmp_path: Path,
) -> None:
    """Replacing a split file would make post-hoc test tuning possible."""
    clusters = []
    sites = []
    for i in range(10):
        species = "arabidopsis" if i % 2 == 0 else "tomato"
        clusters.append(_cluster(species, f"P{i}", f"C{i}"))
        sites.append(_site(species, f"P{i}", f"S_{species}"))
    split = build_frozen_multispecies_split(clusters, sites, seed=20260811)
    path = tmp_path / "frozen.tsv"

    write_frozen_multispecies_split(path, split)

    assert load_frozen_multispecies_split(path) == split
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        write_frozen_multispecies_split(path, split)


def test_global_cluster_allocation_scales_to_realistic_cluster_counts() -> None:
    """A full global proteome must not make frozen-split construction impractical."""
    clusters = [
        _cluster("arabidopsis", f"P{i}", f"C{i}") for i in range(20_000)
    ]
    started = monotonic()
    split = build_frozen_multispecies_split(clusters, (), seed=20260811)

    assert len(split.rows) == 20_000
    assert monotonic() - started < 10.0


def test_development_folds_each_receive_a_positive_cluster() -> None:
    """An empty validation fold cannot support frozen model selection."""
    clusters = [
        _cluster("arabidopsis", f"P{i}", f"C{i}") for i in range(20)
    ]
    sites = [_site("arabidopsis", f"P{i}", "ATH") for i in range(20)]

    split = build_frozen_multispecies_split(
        clusters, sites, seed=20260811, n_development_folds=5
    )

    positive_fold_by_protein = {
        row.global_protein_id: row.development_fold
        for row in split.rows
        if row.split == "development"
    }
    folds = {
        positive_fold_by_protein[site.global_protein_id]
        for site in sites
        if site.global_protein_id in positive_fold_by_protein
    }
    assert folds == {0, 1, 2, 3, 4}


def test_development_split_requires_minimum_positive_sites_and_proteins() -> None:
    """A tiny validation fold must block training rather than look valid."""
    clusters = [
        _cluster("arabidopsis", f"P{i}", f"C{i}") for i in range(20)
    ]
    sites = [_site("arabidopsis", f"P{i}", "ATH") for i in range(20)]

    split = build_frozen_multispecies_split(
        clusters,
        sites,
        seed=20260811,
        n_development_folds=5,
        min_development_positive_sites=10,
        min_development_positive_proteins=5,
    )

    assert split.status == "DEVIATION_REQUIRES_REVIEW"
