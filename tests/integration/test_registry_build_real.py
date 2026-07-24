import csv
import hashlib
from pathlib import Path

import pytest


@pytest.mark.network
def test_registry_builder_registers_all_official_metadata(tmp_path: Path) -> None:
    from plantpersulf.provenance.registry import (
        audit_registry,
        fetch_registered_metadata,
    )

    summary = fetch_registered_metadata(
        config_path=Path("configs/data_sources.yaml"),
        registry_dir=tmp_path,
    )

    assert summary.dataset_count == 12
    assert summary.sample_count > 0
    for file_name in (
        "datasets.tsv",
        "files.tsv",
        "samples.tsv",
        "publications.tsv",
    ):
        assert (tmp_path / file_name).is_file()

    with (tmp_path / "datasets.tsv").open(encoding="utf-8", newline="") as handle:
        datasets = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["accession"] for row in datasets} == {
        "PXD006140",
        "PXD024061",
        "PXD035795",
        "PXD039999",
        "PXD051570",
        "PXD063170",
        "PXD038309",
        "PXD072089",
        "GSE163745",
        "GSE142713",
        "GSE142712",
        "GSE267238",
    }
    geo_datasets = [row for row in datasets if row["repository"] == "GEO"]
    assert all(row["bioproject_accession"].startswith("PRJNA") for row in geo_datasets)
    assert all(row["sra_study_accession"].startswith("SRP") for row in geo_datasets)

    with (tmp_path / "files.tsv").open(encoding="utf-8", newline="") as handle:
        files = list(csv.DictReader(handle, delimiter="\t"))
    cache_rows = [row for row in files if row["record_type"] == "metadata_cache"]
    assert len(cache_rows) == 12
    sra_cache_rows = [
        row for row in files if row["record_type"] == "sra_metadata_cache"
    ]
    assert len(sra_cache_rows) == 4
    for row in cache_rows + sra_cache_rows:
        cache_path = tmp_path / row["path"]
        assert cache_path.is_file()
        assert hashlib.sha256(cache_path.read_bytes()).hexdigest() == row["sha256"]

    with (tmp_path / "samples.tsv").open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    assert samples
    assert all(row["biosample_accession"].startswith("SAMN") for row in samples)
    assert all(row["sra_experiment_accession"].startswith("SRX") for row in samples)
    assert all(row["sra_run_accessions"].startswith("SRR") for row in samples)

    assert audit_registry(
        registry_dir=tmp_path,
        config_path=Path("configs/data_sources.yaml"),
    ) == summary
