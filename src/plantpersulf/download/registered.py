"""Resolve explicitly approved downloads against the official file registry."""

from __future__ import annotations

import csv
import os
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import cast

import yaml

from plantpersulf.download.base import (
    DownloadRequest,
    canonical_download_url,
    download_verified_file,
)
from plantpersulf.provenance.hashing import hash_file

DOWNLOAD_FIELDS = (
    "dataset_accession",
    "repository",
    "file_name",
    "file_class",
    "source_url",
    "download_url",
    "retrieved_at",
    "size_bytes",
    "remote_checksum_algorithm",
    "remote_checksum",
    "path",
    "sha256",
    "license_or_usage",
    "downloader_version",
)


@dataclass(frozen=True)
class RegisteredSelection:
    dataset_accession: str
    repository: str
    record_type: str
    file_name: str
    file_class: str
    source_url: str
    size_bytes: int
    remote_checksum: str
    remote_checksum_algorithm: str
    path: str
    sha256: str
    status: str


@dataclass(frozen=True)
class DownloadSummary:
    selected_count: int
    downloaded_count: int
    cached_count: int


def _load_approved(path: Path) -> dict[str, dict[str, list[str]]]:
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid download selection config: {path}") from exc
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise RuntimeError("download selection config requires version 1")
    approved = loaded.get("approved")
    if not isinstance(approved, dict):
        raise RuntimeError("download selection config requires approved mapping")
    return cast(dict[str, dict[str, list[str]]], approved)


def resolve_registered_selection(
    accession: str,
    file_classes: tuple[str, ...],
    selection_path: Path,
    files_registry_path: Path,
) -> tuple[RegisteredSelection, ...]:
    """Resolve exact approved names; never infer a class from file metadata."""
    normalized_accession = accession.strip().upper()
    approved = _load_approved(selection_path)
    accession_config = approved.get(normalized_accession)
    if not isinstance(accession_config, dict):
        raise RuntimeError(f"accession is not approved: {normalized_accession}")
    with files_registry_path.open(encoding="utf-8", newline="") as handle:
        registry_rows = list(csv.DictReader(handle, delimiter="\t"))
    selected: list[RegisteredSelection] = []
    for file_class in file_classes:
        names = accession_config.get(file_class)
        if not isinstance(names, list) or not names:
            raise RuntimeError(
                f"file class is not approved for {normalized_accession}: {file_class}"
            )
        for file_name in names:
            matches = [
                row
                for row in registry_rows
                if row["dataset_accession"] == normalized_accession
                and row["file_name"] == file_name
            ]
            if len(matches) != 1:
                raise RuntimeError(
                    f"approved file does not resolve exactly once: {file_name}"
                )
            row = matches[0]
            selected.append(
                RegisteredSelection(
                    dataset_accession=normalized_accession,
                    repository=row["repository"],
                    record_type=row["record_type"],
                    file_name=file_name,
                    file_class=file_class,
                    source_url=row["source_url"],
                    size_bytes=int(row["size_bytes"]),
                    remote_checksum=row["remote_checksum"],
                    remote_checksum_algorithm=row["remote_checksum_algorithm"],
                    path=row["path"],
                    sha256=row["sha256"],
                    status=row["status"],
                )
            )
    return tuple(selected)


def _dataset_license(path: Path, accession: str) -> str:
    with path.open(encoding="utf-8", newline="") as handle:
        matches = [
            row
            for row in csv.DictReader(handle, delimiter="\t")
            if row["accession"] == accession
        ]
    if len(matches) != 1 or not matches[0]["license_or_usage"].strip():
        raise RuntimeError(f"dataset license is not registered: {accession}")
    return matches[0]["license_or_usage"]


def _write_download_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(".tsv.tmp")
    with temporary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=DOWNLOAD_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary_path.replace(path)


def download_registered_files(
    accession: str,
    file_classes: tuple[str, ...],
    selection_path: Path,
    files_registry_path: Path,
    datasets_registry_path: Path,
    downloads_registry_path: Path,
    raw_dir: Path,
) -> DownloadSummary:
    """Download exact approved source rows and write local SHA256 provenance."""
    selected = resolve_registered_selection(
        accession,
        file_classes,
        selection_path,
        files_registry_path,
    )
    license_or_usage = _dataset_license(datasets_registry_path, accession)
    rows: list[dict[str, str]] = []
    if downloads_registry_path.exists():
        with downloads_registry_path.open(encoding="utf-8", newline="") as handle:
            rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    downloaded_count = 0
    cached_count = 0
    for item in selected:
        if item.record_type == "metadata_cache":
            metadata_path = files_registry_path.parent / item.path
            if (
                item.status != "cached"
                or not metadata_path.is_file()
                or metadata_path.stat().st_size != item.size_bytes
                or hash_file(metadata_path, "sha256") != item.sha256
            ):
                raise RuntimeError(
                    f"registered metadata cache failed audit: {item.file_name}"
                )
            cached_count += 1
            continue
        if item.record_type != "source_file":
            raise RuntimeError(
                f"selected record is not downloadable: {item.file_name}"
            )
        existing_matches = [
            row
            for row in rows
            if row["dataset_accession"] == item.dataset_accession
            and row["file_name"] == item.file_name
        ]
        if existing_matches:
            if len(existing_matches) != 1:
                raise RuntimeError(
                    f"duplicate download manifest rows: {item.file_name}"
                )
            existing = existing_matches[0]
            expected_identity = {
                "repository": item.repository,
                "file_class": item.file_class,
                "source_url": item.source_url,
                "download_url": canonical_download_url(item.source_url),
                "size_bytes": str(item.size_bytes),
                "remote_checksum_algorithm": item.remote_checksum_algorithm,
                "remote_checksum": item.remote_checksum,
                "license_or_usage": license_or_usage,
                "downloader_version": "plantpersulf/0.0.0",
            }
            if any(
                existing[field] != value
                for field, value in expected_identity.items()
            ):
                raise RuntimeError(
                    f"download provenance mismatch: {item.file_name}"
                )
            existing_path = downloads_registry_path.parent / existing["path"]
            if (
                not existing_path.is_file()
                or existing_path.stat().st_size != item.size_bytes
                or hash_file(existing_path, "sha256") != existing["sha256"]
                or (
                    bool(item.remote_checksum)
                    and hash_file(
                        existing_path,
                        item.remote_checksum_algorithm,
                    )
                    != item.remote_checksum
                )
            ):
                raise RuntimeError(
                    f"existing download failed audit: {item.file_name}"
                )
            cached_count += 1
            continue
        destination = raw_dir / item.dataset_accession / item.file_name
        result = download_verified_file(
            DownloadRequest(
                source_url=item.source_url,
                destination=destination,
                expected_size=item.size_bytes,
                expected_checksum=item.remote_checksum,
                expected_checksum_algorithm=item.remote_checksum_algorithm,
            )
        )
        rows = [
            row
            for row in rows
            if not (
                row["dataset_accession"] == item.dataset_accession
                and row["file_name"] == item.file_name
            )
        ]
        rows.append(
            {
                "dataset_accession": item.dataset_accession,
                "repository": item.repository,
                "file_name": item.file_name,
                "file_class": item.file_class,
                "source_url": item.source_url,
                "download_url": canonical_download_url(item.source_url),
                "retrieved_at": datetime.now(timezone.utc).isoformat(),
                "size_bytes": str(result.size_bytes),
                "remote_checksum_algorithm": item.remote_checksum_algorithm,
                "remote_checksum": item.remote_checksum,
                "path": Path(
                    os.path.relpath(
                        result.destination,
                        downloads_registry_path.parent,
                    )
                ).as_posix(),
                "sha256": result.sha256,
                "license_or_usage": license_or_usage,
                "downloader_version": "plantpersulf/0.0.0",
            }
        )
        rows.sort(key=lambda row: (row["dataset_accession"], row["file_name"]))
        _write_download_manifest(downloads_registry_path, rows)
        downloaded_count += 1
    rows.sort(key=lambda row: (row["dataset_accession"], row["file_name"]))
    _write_download_manifest(downloads_registry_path, rows)
    return DownloadSummary(len(selected), downloaded_count, cached_count)


def audit_downloaded_files(
    accession: str,
    selection_path: Path,
    files_registry_path: Path,
    datasets_registry_path: Path,
    downloads_registry_path: Path,
) -> DownloadSummary:
    """Fail unless every approved local file matches both provenance registries."""
    normalized_accession = accession.strip().upper()
    approved = _load_approved(selection_path)
    accession_config = approved.get(normalized_accession)
    if not isinstance(accession_config, dict):
        raise RuntimeError(f"accession is not approved: {normalized_accession}")
    selected = resolve_registered_selection(
        normalized_accession,
        tuple(accession_config),
        selection_path,
        files_registry_path,
    )
    license_or_usage = _dataset_license(
        datasets_registry_path,
        normalized_accession,
    )
    if not downloads_registry_path.is_file():
        raise RuntimeError(f"download registry missing: {downloads_registry_path}")
    with downloads_registry_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != DOWNLOAD_FIELDS:
            raise RuntimeError("download registry has invalid columns")
        manifest_rows = [dict(row) for row in reader]
    selected_source_names = {
        item.file_name for item in selected if item.record_type == "source_file"
    }
    observed_source_names = {
        row["file_name"]
        for row in manifest_rows
        if row["dataset_accession"] == normalized_accession
    }
    if observed_source_names != selected_source_names:
        raise RuntimeError("download registry does not match approved selection")
    downloaded_count = 0
    cached_count = 0
    for item in selected:
        if item.record_type == "metadata_cache":
            metadata_path = files_registry_path.parent / item.path
            if (
                item.status != "cached"
                or not metadata_path.is_file()
                or metadata_path.stat().st_size != item.size_bytes
                or hash_file(metadata_path, "sha256") != item.sha256
            ):
                raise RuntimeError(
                    f"registered metadata cache failed audit: {item.file_name}"
                )
            cached_count += 1
            continue
        matches = [
            row
            for row in manifest_rows
            if row["dataset_accession"] == normalized_accession
            and row["file_name"] == item.file_name
        ]
        if len(matches) != 1:
            raise RuntimeError(f"download manifest row missing: {item.file_name}")
        row = matches[0]
        expected_identity = {
            "repository": item.repository,
            "file_class": item.file_class,
            "source_url": item.source_url,
            "download_url": canonical_download_url(item.source_url),
            "size_bytes": str(item.size_bytes),
            "remote_checksum_algorithm": item.remote_checksum_algorithm,
            "remote_checksum": item.remote_checksum,
            "license_or_usage": license_or_usage,
            "downloader_version": "plantpersulf/0.0.0",
        }
        if any(row[field] != value for field, value in expected_identity.items()):
            raise RuntimeError(f"download provenance mismatch: {item.file_name}")
        local_path = downloads_registry_path.parent / row["path"]
        if not local_path.is_file() or local_path.stat().st_size != item.size_bytes:
            raise RuntimeError(f"downloaded file size mismatch: {item.file_name}")
        if hash_file(local_path, "sha256") != row["sha256"]:
            raise RuntimeError(f"downloaded file SHA256 mismatch: {item.file_name}")
        if item.remote_checksum and hash_file(
            local_path, item.remote_checksum_algorithm
        ) != item.remote_checksum:
            raise RuntimeError(
                f"downloaded file remote checksum mismatch: {item.file_name}"
            )
        downloaded_count += 1
    return DownloadSummary(len(selected), downloaded_count, cached_count)
