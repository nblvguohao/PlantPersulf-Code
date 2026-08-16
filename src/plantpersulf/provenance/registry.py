"""Build tabular provenance registries from official metadata caches."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.download.geo import GeoClient, GeoSeries
from plantpersulf.download.iprox import IproxClient, IproxDataset
from plantpersulf.download.pride import PrideClient, PrideProject
from plantpersulf.download.sra import SraClient, SraStudy
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.schema import DatasetRecord


@dataclass(frozen=True)
class DataSource:
    repository: str
    accession: str
    scientific_role: str


@dataclass(frozen=True)
class RegistrySummary:
    dataset_count: int
    file_count: int
    sample_count: int
    publication_count: int


Row = dict[str, str]

DATASET_FIELDS = (
    "accession",
    "repository",
    "scientific_role",
    "title",
    "organism",
    "bioproject_accession",
    "sra_study_accession",
    "publication_date",
    "source_url",
    "metadata_retrieved_at",
    "metadata_cache_path",
    "metadata_sha256",
    "license_or_usage",
)
FILE_FIELDS = (
    "dataset_accession",
    "repository",
    "record_type",
    "file_name",
    "file_category",
    "source_url",
    "size_bytes",
    "remote_checksum",
    "remote_checksum_algorithm",
    "path",
    "sha256",
    "status",
    "retrieved_at",
)
SAMPLE_FIELDS = (
    "dataset_accession",
    "sample_accession",
    "title",
    "organism",
    "source_name",
    "characteristics_json",
    "biosample_accession",
    "sra_experiment_accession",
    "sra_run_accessions",
)
PUBLICATION_FIELDS = (
    "dataset_accession",
    "repository",
    "identifier_type",
    "identifier",
    "doi",
    "citation",
)


def load_data_sources(config_path: Path) -> tuple[DataSource, ...]:
    try:
        loaded: object = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid data source config: {config_path}") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError("data source config must be a mapping")
    config = cast(dict[str, Any], loaded)
    sources: list[DataSource] = []
    for key, repository, prefix in (
        ("pride", "PRIDE", "PXD"),
        ("iprox", "iProX", "PXD"),
        ("geo", "GEO", "GSE"),
    ):
        entries = config.get(key)
        if not isinstance(entries, list) or not entries:
            raise RuntimeError(f"data source config requires non-empty {key} list")
        for raw_entry in entries:
            if not isinstance(raw_entry, dict):
                raise RuntimeError(f"data source entry in {key} must be a mapping")
            entry = cast(dict[str, Any], raw_entry)
            accession = str(entry.get("accession", "")).strip().upper()
            role = str(entry.get("role", "")).strip()
            if not accession.startswith(prefix) or not role:
                raise RuntimeError(f"invalid {repository} data source entry: {entry}")
            sources.append(DataSource(repository, accession, role))
    identities = {(source.repository, source.accession) for source in sources}
    if len(identities) != len(sources):
        raise RuntimeError("data source config contains duplicate accessions")
    return tuple(sources)


def _write_tsv(path: Path, fieldnames: tuple[str, ...], rows: list[Row]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(f"{path.suffix}.tmp")
    with temporary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary_path.replace(path)


def _checksum_algorithm(checksum: str) -> str:
    if re.fullmatch(r"[0-9a-fA-F]{40}", checksum):
        return "SHA1"
    if re.fullmatch(r"[0-9a-fA-F]{32}", checksum):
        return "MD5"
    return "repository_reported" if checksum else ""


def _preferred_url(urls: tuple[str, ...]) -> str:
    for prefix in ("https://", "ftp://", "http://"):
        for url in urls:
            if url.startswith(prefix):
                return url
    return urls[0] if urls else ""


def _cache_file_row(
    accession: str,
    repository: str,
    source_url: str,
    cache_path: Path,
    cache_sha256: str,
    retrieved_at: str,
    registry_dir: Path,
) -> Row:
    return {
        "dataset_accession": accession,
        "repository": repository,
        "record_type": "metadata_cache",
        "file_name": cache_path.name,
        "file_category": "official_metadata",
        "source_url": source_url,
        "size_bytes": str(cache_path.stat().st_size),
        "remote_checksum": "",
        "remote_checksum_algorithm": "",
        "path": cache_path.relative_to(registry_dir).as_posix(),
        "sha256": cache_sha256,
        "status": "cached",
        "retrieved_at": retrieved_at,
    }


def _append_pride_rows(
    source: DataSource,
    project: PrideProject,
    registry_dir: Path,
    datasets: list[Row],
    files: list[Row],
    publications: list[Row],
) -> None:
    cache_path = project.cache_path.relative_to(registry_dir).as_posix()
    datasets.append(
        {
            "accession": project.accession,
            "repository": "PRIDE",
            "scientific_role": source.scientific_role,
            "title": project.title,
            "organism": "; ".join(project.organisms),
            "bioproject_accession": "",
            "sra_study_accession": "",
            "publication_date": project.publication_date,
            "source_url": project.source_url,
            "metadata_retrieved_at": project.retrieved_at,
            "metadata_cache_path": cache_path,
            "metadata_sha256": project.cache_sha256,
            "license_or_usage": project.license_name,
        }
    )
    files.append(
        _cache_file_row(
            project.accession,
            "PRIDE",
            project.source_url,
            project.cache_path,
            project.cache_sha256,
            project.retrieved_at,
            registry_dir,
        )
    )
    for project_file in project.files:
        files.append(
            {
                "dataset_accession": project.accession,
                "repository": "PRIDE",
                "record_type": "source_file",
                "file_name": project_file.file_name,
                "file_category": project_file.category,
                "source_url": _preferred_url(project_file.public_urls),
                "size_bytes": str(project_file.size_bytes),
                "remote_checksum": project_file.checksum,
                "remote_checksum_algorithm": _checksum_algorithm(project_file.checksum),
                "path": "",
                "sha256": "",
                "status": "remote_only",
                "retrieved_at": project.retrieved_at,
            }
        )
    for reference in project.references:
        if reference.pubmed_id:
            publications.append(
                {
                    "dataset_accession": project.accession,
                    "repository": "PRIDE",
                    "identifier_type": "PubMed",
                    "identifier": reference.pubmed_id,
                    "doi": reference.doi,
                    "citation": reference.citation,
                }
            )


def _append_geo_rows(
    source: DataSource,
    series: GeoSeries,
    sra_study: SraStudy,
    registry_dir: Path,
    datasets: list[Row],
    files: list[Row],
    samples: list[Row],
    publications: list[Row],
) -> None:
    if series.bioproject_accession != sra_study.bioproject_accession:
        raise RuntimeError(f"GEO/SRA BioProject mismatch for {series.accession}")
    if (
        series.sra_study_accession
        and series.sra_study_accession != sra_study.sra_study_accession
    ):
        raise RuntimeError(f"GEO/SRA study mismatch for {series.accession}")
    samples_by_experiment = {
        sample.sra_experiment_accession: sample for sample in series.samples
    }
    runs_by_experiment: dict[str, list[str]] = {}
    for run in sra_study.runs:
        sample = samples_by_experiment.get(run.experiment_accession)
        if sample is None:
            raise RuntimeError(
                f"SRA run {run.run_accession} is not a GEO sample in {series.accession}"
            )
        if (
            run.biosample_accession != sample.biosample_accession
            or run.sample_name != sample.accession
        ):
            raise RuntimeError(
                f"GEO/SRA sample identity mismatch for {sample.accession}"
            )
        runs_by_experiment.setdefault(run.experiment_accession, []).append(
            run.run_accession
        )
    if set(runs_by_experiment) != set(samples_by_experiment):
        raise RuntimeError(f"GEO samples lack SRA runs for {series.accession}")

    cache_path = series.cache_path.relative_to(registry_dir).as_posix()
    datasets.append(
        {
            "accession": series.accession,
            "repository": "GEO",
            "scientific_role": source.scientific_role,
            "title": series.title,
            "organism": series.organism,
            "bioproject_accession": series.bioproject_accession,
            "sra_study_accession": sra_study.sra_study_accession,
            "publication_date": series.publication_date,
            "source_url": series.source_url,
            "metadata_retrieved_at": series.retrieved_at,
            "metadata_cache_path": cache_path,
            "metadata_sha256": series.cache_sha256,
            "license_or_usage": "NCBI GEO data disclaimer",
        }
    )
    files.append(
        _cache_file_row(
            series.accession,
            "GEO",
            series.source_url,
            series.cache_path,
            series.cache_sha256,
            series.retrieved_at,
            registry_dir,
        )
    )
    files.append(
        {
            **_cache_file_row(
                series.accession,
                "GEO",
                sra_study.source_url,
                sra_study.cache_path,
                sra_study.cache_sha256,
                sra_study.retrieved_at,
                registry_dir,
            ),
            "record_type": "sra_metadata_cache",
            "file_category": "official_sra_metadata",
        }
    )
    for source_url in series.supplementary_files:
        if source_url and source_url != "NONE":
            files.append(
                {
                    "dataset_accession": series.accession,
                    "repository": "GEO",
                    "record_type": "source_file",
                    "file_name": source_url.rsplit("/", 1)[-1],
                    "file_category": "supplementary",
                    "source_url": source_url,
                    "size_bytes": "",
                    "remote_checksum": "",
                    "remote_checksum_algorithm": "",
                    "path": "",
                    "sha256": "",
                    "status": "remote_only",
                    "retrieved_at": series.retrieved_at,
                }
            )
    for sample in series.samples:
        samples.append(
            {
                "dataset_accession": series.accession,
                "sample_accession": sample.accession,
                "title": sample.title,
                "organism": sample.organism,
                "source_name": sample.source_name,
                "characteristics_json": json.dumps(
                    sample.characteristics,
                    ensure_ascii=False,
                    separators=(",", ":"),
                ),
                "biosample_accession": sample.biosample_accession,
                "sra_experiment_accession": sample.sra_experiment_accession,
                "sra_run_accessions": ";".join(
                    sorted(runs_by_experiment[sample.sra_experiment_accession])
                ),
            }
        )
    for pubmed_id in series.pubmed_ids:
        publications.append(
            {
                "dataset_accession": series.accession,
                "repository": "GEO",
                "identifier_type": "PubMed",
                "identifier": pubmed_id,
                "doi": "",
                "citation": "",
            }
        )


def _append_iprox_rows(
    source: DataSource,
    dataset: IproxDataset,
    registry_dir: Path,
    datasets: list[Row],
    files: list[Row],
    publications: list[Row],
) -> None:
    cache_path = dataset.cache_path.relative_to(registry_dir).as_posix()
    datasets.append(
        {
            "accession": dataset.accession,
            "repository": "iProX",
            "scientific_role": source.scientific_role,
            "title": dataset.title,
            "organism": "; ".join(dataset.organisms),
            "bioproject_accession": "",
            "sra_study_accession": "",
            "publication_date": dataset.publication_date,
            "source_url": dataset.source_url,
            "metadata_retrieved_at": dataset.retrieved_at,
            "metadata_cache_path": cache_path,
            "metadata_sha256": dataset.cache_sha256,
            "license_or_usage": "ProteomeXchange public data policies",
        }
    )
    files.append(
        _cache_file_row(
            dataset.accession,
            "iProX",
            dataset.source_url,
            dataset.cache_path,
            dataset.cache_sha256,
            dataset.retrieved_at,
            registry_dir,
        )
    )
    for dataset_file in dataset.files:
        files.append(
            {
                "dataset_accession": dataset.accession,
                "repository": "iProX",
                "record_type": "source_file",
                "file_name": dataset_file.file_name,
                "file_category": dataset_file.category,
                "source_url": dataset_file.source_url,
                "size_bytes": "",
                "remote_checksum": "",
                "remote_checksum_algorithm": "",
                "path": "",
                "sha256": "",
                "status": "remote_only",
                "retrieved_at": dataset.retrieved_at,
            }
        )
    for reference in dataset.references:
        if reference.pubmed_id:
            publications.append(
                {
                    "dataset_accession": dataset.accession,
                    "repository": "iProX",
                    "identifier_type": "PubMed",
                    "identifier": reference.pubmed_id,
                    "doi": "",
                    "citation": reference.citation,
                }
            )


def fetch_registered_metadata(
    config_path: Path = Path("configs/data_sources.yaml"),
    registry_dir: Path = Path("data/registry"),
) -> RegistrySummary:
    """Fetch configured official metadata and rebuild the four registry tables."""
    sources = load_data_sources(config_path)
    pride_client = PrideClient(cache_dir=registry_dir / "cache/pride")
    iprox_client = IproxClient(cache_dir=registry_dir / "cache/iprox")
    geo_client = GeoClient(cache_dir=registry_dir / "cache/geo")
    sra_client = SraClient(cache_dir=registry_dir / "cache/sra")
    datasets: list[Row] = []
    files: list[Row] = []
    samples: list[Row] = []
    publications: list[Row] = []
    for source in sources:
        if source.repository == "PRIDE":
            _append_pride_rows(
                source,
                pride_client.get_project(source.accession),
                registry_dir,
                datasets,
                files,
                publications,
            )
        elif source.repository == "iProX":
            _append_iprox_rows(
                source,
                iprox_client.get_dataset(source.accession),
                registry_dir,
                datasets,
                files,
                publications,
            )
        else:
            series = geo_client.get_series(source.accession)
            sra_query = series.sra_study_accession or series.bioproject_accession
            _append_geo_rows(
                source,
                series,
                sra_client.get_runs(sra_query),
                registry_dir,
                datasets,
                files,
                samples,
                publications,
            )

    datasets.sort(key=lambda row: (row["repository"], row["accession"]))
    files.sort(
        key=lambda row: (
            row["repository"],
            row["dataset_accession"],
            row["record_type"],
            row["file_name"],
        )
    )
    samples.sort(key=lambda row: (row["dataset_accession"], row["sample_accession"]))
    publications.sort(
        key=lambda row: (
            row["dataset_accession"],
            row["identifier_type"],
            row["identifier"],
        )
    )
    _write_tsv(registry_dir / "datasets.tsv", DATASET_FIELDS, datasets)
    _write_tsv(registry_dir / "files.tsv", FILE_FIELDS, files)
    _write_tsv(registry_dir / "samples.tsv", SAMPLE_FIELDS, samples)
    _write_tsv(
        registry_dir / "publications.tsv",
        PUBLICATION_FIELDS,
        publications,
    )
    return RegistrySummary(
        dataset_count=len(datasets),
        file_count=len(files),
        sample_count=len(samples),
        publication_count=len(publications),
    )


def _read_tsv(path: Path, expected_fields: tuple[str, ...]) -> list[Row]:
    if not path.is_file():
        raise RuntimeError(f"registry table missing: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != expected_fields:
            raise RuntimeError(f"registry table has invalid columns: {path}")
        rows: list[Row] = []
        for raw_row in reader:
            if None in raw_row or any(value is None for value in raw_row.values()):
                raise RuntimeError(f"registry table has malformed row: {path}")
            rows.append({key: str(value) for key, value in raw_row.items()})
    return rows


def audit_registry(
    registry_dir: Path = Path("data/registry"),
    config_path: Path = Path("configs/data_sources.yaml"),
) -> RegistrySummary:
    """Fail unless registry tables fully trace configured official metadata."""
    sources = load_data_sources(config_path)
    datasets = _read_tsv(registry_dir / "datasets.tsv", DATASET_FIELDS)
    files = _read_tsv(registry_dir / "files.tsv", FILE_FIELDS)
    samples = _read_tsv(registry_dir / "samples.tsv", SAMPLE_FIELDS)
    publications = _read_tsv(
        registry_dir / "publications.tsv",
        PUBLICATION_FIELDS,
    )
    expected_identities = {(source.repository, source.accession) for source in sources}
    observed_identities = {(row["repository"], row["accession"]) for row in datasets}
    if observed_identities != expected_identities or len(datasets) != len(sources):
        raise RuntimeError("dataset registry does not exactly match data source config")

    dataset_accessions = {row["accession"] for row in datasets}
    for row in datasets:
        DatasetRecord(
            accession=row["accession"],
            repository=row["repository"],
            source_url=row["source_url"],
            scientific_role=row["scientific_role"],
            metadata_retrieved_at=row["metadata_retrieved_at"],
            metadata_sha256=row["metadata_sha256"],
        )
        for required_field in (
            "title",
            "organism",
            "publication_date",
            "metadata_cache_path",
            "license_or_usage",
        ):
            if not row[required_field].strip():
                raise RuntimeError(f"dataset {row['accession']} lacks {required_field}")
        if row["repository"] == "GEO" and (
            not row["bioproject_accession"].startswith("PRJNA")
            or not row["sra_study_accession"].startswith("SRP")
        ):
            raise RuntimeError(
                f"GEO dataset {row['accession']} lacks BioProject or SRA study"
            )
        matching_cache_rows = [
            file_row
            for file_row in files
            if file_row["dataset_accession"] == row["accession"]
            and file_row["repository"] == row["repository"]
            and file_row["record_type"] == "metadata_cache"
        ]
        if len(matching_cache_rows) != 1:
            raise RuntimeError(
                f"dataset {row['accession']} requires exactly one metadata cache"
            )
        cache_row = matching_cache_rows[0]
        if (
            cache_row["path"] != row["metadata_cache_path"]
            or cache_row["sha256"] != row["metadata_sha256"]
            or cache_row["status"] != "cached"
        ):
            raise RuntimeError(
                f"dataset cache registration mismatch: {row['accession']}"
            )
        assert_registered_input(
            registry_dir / row["metadata_cache_path"],
            registry_path=registry_dir / "files.tsv",
        )
        matching_sra_cache_rows = [
            file_row
            for file_row in files
            if file_row["dataset_accession"] == row["accession"]
            and file_row["repository"] == row["repository"]
            and file_row["record_type"] == "sra_metadata_cache"
        ]
        expected_sra_cache_count = 1 if row["repository"] == "GEO" else 0
        if len(matching_sra_cache_rows) != expected_sra_cache_count:
            raise RuntimeError(
                f"dataset {row['accession']} has invalid SRA cache count"
            )
        for sra_cache_row in matching_sra_cache_rows:
            assert_registered_input(
                registry_dir / sra_cache_row["path"],
                registry_path=registry_dir / "files.tsv",
            )

    for row in files:
        identity = (row["repository"], row["dataset_accession"])
        if identity not in expected_identities:
            raise RuntimeError(f"file row references unknown dataset: {identity}")
        if row["record_type"] in {"metadata_cache", "sra_metadata_cache"}:
            if not row["path"] or not row["sha256"] or row["status"] != "cached":
                raise RuntimeError("metadata cache row lacks local path or SHA256")
        elif row["record_type"] == "source_file":
            if not row["source_url"] or row["status"] != "remote_only":
                raise RuntimeError(
                    "remote source file row lacks official URL or status"
                )
            if row["path"] or row["sha256"]:
                raise RuntimeError("remote-only source file cannot claim local SHA256")
        else:
            raise RuntimeError(f"unsupported file record type: {row['record_type']}")

    for row in samples:
        if row["dataset_accession"] not in dataset_accessions:
            raise RuntimeError("sample row references unknown dataset")
        required_sample_fields = (
            "sample_accession",
            "organism",
            "biosample_accession",
            "sra_experiment_accession",
            "sra_run_accessions",
        )
        if any(not row[field] for field in required_sample_fields):
            raise RuntimeError("sample row lacks official accession or organism")
    for row in publications:
        if row["dataset_accession"] not in dataset_accessions:
            raise RuntimeError("publication row references unknown dataset")
        if not row["identifier_type"] or not row["identifier"]:
            raise RuntimeError("publication row lacks identifier")

    return RegistrySummary(
        dataset_count=len(datasets),
        file_count=len(files),
        sample_count=len(samples),
        publication_count=len(publications),
    )
