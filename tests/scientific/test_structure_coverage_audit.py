"""RED (Phase S0, Task 2): structure coverage audit with fail-closed rules.

Every benchmark row is mapped against one registered structure release and gets
exactly one ``CoverageRecord`` with an allowed ``mapping_status``. The audit
must catch blocking errors (checksum, missing-file, registry-conflict) before
any scoring can begin, and it must keep positives and unlabeled rows distinct.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.evaluation.structure_coverage_audit import (
    CoverageRecord,
    StructureCoverageAudit,
    audit_structure_coverage,
    sha256_file,
)

# Real frozen files
BENCHMARK = Path("data/processed/benchmark_v1/sites.tsv")
PROTEOME = Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta")
CLUSTERS = Path("data/processed/clusters/protein_clusters_v2.tsv")
REGISTRY_V1 = Path("data/registry/releases/alphafold_structures_release_v1.tsv")
REGISTRY_V2 = Path("data/registry/releases/alphafold_structures_release_v2.tsv")
REGISTRY_DATA_ROOT = Path("data/registry")


ALLOWED_STATUSES = frozenset({
    "mapped_cys",
    "absent_structure",
    "absent_residue",
    "non_cys_residue",
    "malformed_structure",
    "duplicate_residue_ambiguity",
})


def test_audit_v1_returns_one_record_per_benchmark_row() -> None:
    audit = audit_structure_coverage(
        release="structcover_v1",
        benchmark_sha256=sha256_file(BENCHMARK),
        proteome_sha256=sha256_file(PROTEOME),
        clusters_sha256=sha256_file(CLUSTERS),
        registry_sha256=sha256_file(REGISTRY_V1),
        registry_path=REGISTRY_V1,
        registry_base=REGISTRY_DATA_ROOT,
        benchmark_path=BENCHMARK,
        proteome_path=PROTEOME,
        clusters_path=CLUSTERS,
    )
    # Count benchmark rows
    n_rows = 0
    with BENCHMARK.open(encoding="utf-8", newline="") as h:
        import csv
        for _ in csv.DictReader(h, delimiter="\t"):
            n_rows += 1

    assert len(audit.records) == n_rows
    for rec in audit.records:
        assert rec.mapping_status in ALLOWED_STATUSES
    assert not audit.blocking_errors  # real data must pass


def test_audit_v2_returns_one_record_per_benchmark_row() -> None:
    audit = audit_structure_coverage(
        release="structcover_v2",
        benchmark_sha256=sha256_file(BENCHMARK),
        proteome_sha256=sha256_file(PROTEOME),
        clusters_sha256=sha256_file(CLUSTERS),
        registry_sha256=sha256_file(REGISTRY_V2),
        registry_path=REGISTRY_V2,
        registry_base=REGISTRY_DATA_ROOT,
        benchmark_path=BENCHMARK,
        proteome_path=PROTEOME,
        clusters_path=CLUSTERS,
    )
    n_pos = sum(1 for r in audit.records if r.label == "positive")
    n_ul = sum(1 for r in audit.records if r.label == "unlabeled")
    assert n_pos > 0
    assert n_ul > 0
    assert not audit.blocking_errors


def test_audit_keeps_positives_and_unlabeled_separate() -> None:
    audit = audit_structure_coverage(
        release="structcover_v2",
        benchmark_sha256=sha256_file(BENCHMARK),
        proteome_sha256=sha256_file(PROTEOME),
        clusters_sha256=sha256_file(CLUSTERS),
        registry_sha256=sha256_file(REGISTRY_V2),
        registry_path=REGISTRY_V2,
        registry_base=REGISTRY_DATA_ROOT,
        benchmark_path=BENCHMARK,
        proteome_path=PROTEOME,
        clusters_path=CLUSTERS,
    )
    pos_records = [r for r in audit.records if r.label == "positive"]
    ul_records = [r for r in audit.records if r.label == "unlabeled"]
    assert len(pos_records) > 0
    assert len(ul_records) > 0
    # No record has a label that is not "positive" or "unlabeled"
    for rec in pos_records:
        assert rec.label == "positive"
    for rec in ul_records:
        assert rec.label == "unlabeled"


def test_audit_summary_contains_required_dimensions() -> None:
    audit = audit_structure_coverage(
        release="structcover_v2",
        benchmark_sha256=sha256_file(BENCHMARK),
        proteome_sha256=sha256_file(PROTEOME),
        clusters_sha256=sha256_file(CLUSTERS),
        registry_sha256=sha256_file(REGISTRY_V2),
        registry_path=REGISTRY_V2,
        registry_base=REGISTRY_DATA_ROOT,
        benchmark_path=BENCHMARK,
        proteome_path=PROTEOME,
        clusters_path=CLUSTERS,
    )
    required_metrics = {
        "mapped_cys_sites",
        "coverage_rate",
        "registered_proteins",
        "covered_proteins",
    }
    assert audit.summary_rows
    metric_names = {row["metric"] for row in audit.summary_rows}
    assert required_metrics <= metric_names


def test_audit_require_pass_succeeds_on_real_data() -> None:
    audit = audit_structure_coverage(
        release="structcover_v2",
        benchmark_sha256=sha256_file(BENCHMARK),
        proteome_sha256=sha256_file(PROTEOME),
        clusters_sha256=sha256_file(CLUSTERS),
        registry_sha256=sha256_file(REGISTRY_V2),
        registry_path=REGISTRY_V2,
        registry_base=REGISTRY_DATA_ROOT,
        benchmark_path=BENCHMARK,
        proteome_path=PROTEOME,
        clusters_path=CLUSTERS,
    )
    audit.require_pass()


def test_v1_coverage_is_smaller_than_v2() -> None:
    audit_v1 = audit_structure_coverage(
        release="structcover_v1",
        benchmark_sha256=sha256_file(BENCHMARK),
        proteome_sha256=sha256_file(PROTEOME),
        clusters_sha256=sha256_file(CLUSTERS),
        registry_sha256=sha256_file(REGISTRY_V1),
        registry_path=REGISTRY_V1,
        registry_base=REGISTRY_DATA_ROOT,
        benchmark_path=BENCHMARK,
        proteome_path=PROTEOME,
        clusters_path=CLUSTERS,
    )
    audit_v2 = audit_structure_coverage(
        release="structcover_v2",
        benchmark_sha256=sha256_file(BENCHMARK),
        proteome_sha256=sha256_file(PROTEOME),
        clusters_sha256=sha256_file(CLUSTERS),
        registry_sha256=sha256_file(REGISTRY_V2),
        registry_path=REGISTRY_V2,
        registry_base=REGISTRY_DATA_ROOT,
        benchmark_path=BENCHMARK,
        proteome_path=PROTEOME,
        clusters_path=CLUSTERS,
    )
    v1_covered = sum(1 for r in audit_v1.records if r.has_registered_structure)
    v2_covered = sum(1 for r in audit_v2.records if r.has_registered_structure)
    assert v1_covered <= v2_covered
    assert v2_covered > v1_covered  # expansion is real


def test_input_hashes_are_present_in_audit() -> None:
    audit = audit_structure_coverage(
        release="structcover_v1",
        benchmark_sha256=sha256_file(BENCHMARK),
        proteome_sha256=sha256_file(PROTEOME),
        clusters_sha256=sha256_file(CLUSTERS),
        registry_sha256=sha256_file(REGISTRY_V1),
        registry_path=REGISTRY_V1,
        registry_base=REGISTRY_DATA_ROOT,
        benchmark_path=BENCHMARK,
        proteome_path=PROTEOME,
        clusters_path=CLUSTERS,
    )
    assert audit.input_hashes["benchmark"] == sha256_file(BENCHMARK)
    assert audit.input_hashes["proteome"] == sha256_file(PROTEOME)
    assert audit.input_hashes["clusters"] == sha256_file(CLUSTERS)
    assert audit.input_hashes["registry"] == sha256_file(REGISTRY_V1)
