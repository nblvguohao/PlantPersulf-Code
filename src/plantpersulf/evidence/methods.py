"""Fail-closed acquisition and audit of experimental method sources."""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.download.base import DownloadRequest, download_verified_file
from plantpersulf.provenance.hashing import hash_file

METHOD_FIELDS = (
    "study_accession",
    "identifier_type",
    "identifier",
    "repository",
    "repository_record_id",
    "official_url",
    "local_path",
    "retrieved_at",
    "size_bytes",
    "sha256",
    "method_scope",
    "data_use_status",
    "downloader_version",
)
DOWNLOADER_VERSION = "plantpersulf/0.0.0"


@dataclass(frozen=True)
class MethodSpec:
    study_accession: str
    identifier_type: str
    identifier: str
    repository: str
    repository_record_id: str
    official_url: str
    destination: Path
    expected_size_bytes: int
    method_scope: str
    data_use_status: str


@dataclass(frozen=True)
class MethodSource:
    study_accession: str
    identifier_type: str
    identifier: str
    repository: str
    repository_record_id: str
    official_url: str
    local_path: Path
    retrieved_at: str
    size_bytes: int
    sha256: str
    method_scope: str
    data_use_status: str
    downloader_version: str


def _required_string(mapping: dict[str, Any], field: str) -> str:
    value = mapping.get(field)
    if not isinstance(value, str) or not value.strip():
        raise RuntimeError(f"method source has invalid {field}")
    return value


def _load_specs(path: Path) -> tuple[MethodSpec, ...]:
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid evidence method config: {path}") from exc
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise RuntimeError("evidence method config requires version 1")
    raw_sources = loaded.get("sources")
    if not isinstance(raw_sources, list) or not raw_sources:
        raise RuntimeError("evidence method config requires sources")
    specs: list[MethodSpec] = []
    for raw_source in raw_sources:
        if not isinstance(raw_source, dict):
            raise RuntimeError("evidence method source must be a mapping")
        source = cast(dict[str, Any], raw_source)
        size = source.get("expected_size_bytes")
        if not isinstance(size, int) or size <= 0:
            raise RuntimeError("method source has invalid expected_size_bytes")
        specs.append(
            MethodSpec(
                study_accession=_required_string(source, "study_accession"),
                identifier_type=_required_string(source, "identifier_type"),
                identifier=_required_string(source, "identifier"),
                repository=_required_string(source, "repository"),
                repository_record_id=_required_string(source, "repository_record_id"),
                official_url=_required_string(source, "official_url"),
                destination=Path(_required_string(source, "destination")),
                expected_size_bytes=size,
                method_scope=_required_string(source, "method_scope"),
                data_use_status=_required_string(source, "data_use_status"),
            )
        )
    if len({spec.study_accession for spec in specs}) != len(specs):
        raise RuntimeError("duplicate evidence method study")
    return tuple(specs)


def _read_registry(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != METHOD_FIELDS:
            raise RuntimeError("evidence method registry has invalid columns")
        rows = [dict(row) for row in reader]
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise RuntimeError("evidence method registry has malformed row")
    return rows


def _write_registry(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tsv.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=METHOD_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _source_from_row(row: dict[str, str], registry_path: Path) -> MethodSource:
    local_path = (registry_path.parent / row["local_path"]).resolve()
    try:
        size = int(row["size_bytes"])
    except ValueError as exc:
        raise RuntimeError("evidence method registry has invalid size") from exc
    return MethodSource(
        study_accession=row["study_accession"],
        identifier_type=row["identifier_type"],
        identifier=row["identifier"],
        repository=row["repository"],
        repository_record_id=row["repository_record_id"],
        official_url=row["official_url"],
        local_path=local_path,
        retrieved_at=row["retrieved_at"],
        size_bytes=size,
        sha256=row["sha256"],
        method_scope=row["method_scope"],
        data_use_status=row["data_use_status"],
        downloader_version=row["downloader_version"],
    )


def _identity(spec: MethodSpec) -> dict[str, str]:
    return {
        "study_accession": spec.study_accession,
        "identifier_type": spec.identifier_type,
        "identifier": spec.identifier,
        "repository": spec.repository,
        "repository_record_id": spec.repository_record_id,
        "official_url": spec.official_url,
        "size_bytes": str(spec.expected_size_bytes),
        "method_scope": spec.method_scope,
        "data_use_status": spec.data_use_status,
        "downloader_version": DOWNLOADER_VERSION,
    }


def audit_method_sources(
    config_path: Path = Path("configs/evidence_method_sources_v1.yaml"),
    registry_path: Path = Path("data/registry/evidence_methods.tsv"),
) -> tuple[MethodSource, ...]:
    """Verify every reviewed method source against its local SHA256."""
    specs = _load_specs(config_path)
    rows = _read_registry(registry_path)
    if {row["study_accession"] for row in rows} != {
        spec.study_accession for spec in specs
    }:
        raise RuntimeError("evidence method registry does not match config")
    sources: list[MethodSource] = []
    for spec in specs:
        matches = [
            row for row in rows if row["study_accession"] == spec.study_accession
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"method source does not resolve exactly once: {spec.study_accession}"
            )
        row = matches[0]
        if any(row[field] != value for field, value in _identity(spec).items()):
            raise RuntimeError(
                f"method source provenance mismatch: {spec.study_accession}"
            )
        source = _source_from_row(row, registry_path)
        if (
            not source.retrieved_at
            or re.fullmatch(r"[0-9a-f]{64}", source.sha256) is None
            or not source.local_path.is_file()
            or source.local_path.stat().st_size != source.size_bytes
            or hash_file(source.local_path, "sha256") != source.sha256
        ):
            raise RuntimeError(
                f"method source failed local audit: {spec.study_accession}"
            )
        sources.append(source)
    return tuple(sources)


def acquire_method_source(
    accession: str,
    config_path: Path = Path("configs/evidence_method_sources_v1.yaml"),
    registry_path: Path = Path("data/registry/evidence_methods.tsv"),
    raw_root: Path = Path("data/raw"),
) -> MethodSource:
    """Download one exact allow-listed method source and register its SHA256."""
    normalized = accession.strip().upper()
    matches = [
        spec for spec in _load_specs(config_path) if spec.study_accession == normalized
    ]
    if len(matches) != 1:
        raise RuntimeError(f"evidence method accession is not approved: {normalized}")
    spec = matches[0]
    rows = _read_registry(registry_path)
    existing = [row for row in rows if row["study_accession"] == normalized]
    if existing:
        if len(existing) != 1:
            raise RuntimeError(f"duplicate method source: {normalized}")
        return next(
            source
            for source in audit_method_sources(config_path, registry_path)
            if source.study_accession == normalized
        )
    destination = raw_root / spec.destination
    result = download_verified_file(
        DownloadRequest(
            source_url=spec.official_url,
            destination=destination,
            expected_size=spec.expected_size_bytes,
            timeout_seconds=120.0,
        )
    )
    row = {
        **_identity(spec),
        "local_path": Path(
            os.path.relpath(result.destination, registry_path.parent)
        ).as_posix(),
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "sha256": result.sha256,
    }
    rows.append(row)
    rows.sort(key=lambda item: item["study_accession"])
    _write_registry(registry_path, rows)
    return next(
        source
        for source in audit_method_sources(config_path, registry_path)
        if source.study_accession == normalized
    )
