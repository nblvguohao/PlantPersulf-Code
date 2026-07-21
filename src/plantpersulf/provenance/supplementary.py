"""Fail-closed registry of persulfidation-site supplement/search-output sources.

Some Gate-1 site-level evidence lives in files that are not the primary PRIDE
result files: a paper supplement (PXD006140 Dataset S3) and deposited MaxQuant
site tables extracted from a large archive (PXD024061). This module registers
each such source with full provenance and audits it against its local SHA256, so
the site parsers only ever consume hash-verified, allow-listed inputs.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.provenance.hashing import hash_file

SUPPLEMENTARY_FIELDS = (
    "study_accession",
    "member",
    "publication_doi",
    "source_kind",
    "container",
    "container_url",
    "local_path",
    "retrieved_at",
    "size_bytes",
    "sha256",
    "data_level",
    "scientific_use",
    "downloader_version",
)
DOWNLOADER_VERSION = "plantpersulf/0.0.0"


@dataclass(frozen=True)
class SupplementarySpec:
    study_accession: str
    member: str
    publication_doi: str
    source_kind: str
    container: str
    container_url: str
    destination: Path
    expected_size_bytes: int
    data_level: str
    scientific_use: str


@dataclass(frozen=True)
class SupplementarySource:
    study_accession: str
    member: str
    publication_doi: str
    source_kind: str
    container: str
    container_url: str
    local_path: Path
    retrieved_at: str
    size_bytes: int
    sha256: str
    data_level: str
    scientific_use: str
    downloader_version: str


def _required_string(mapping: dict[str, Any], field: str) -> str:
    value = mapping.get(field)
    if not isinstance(value, str) or not value.strip():
        raise RuntimeError(f"supplementary source has invalid {field}")
    return value


def _load_specs(path: Path) -> tuple[SupplementarySpec, ...]:
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid supplementary source config: {path}") from exc
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise RuntimeError("supplementary source config requires version 1")
    raw_sources = loaded.get("sources")
    if not isinstance(raw_sources, list) or not raw_sources:
        raise RuntimeError("supplementary source config requires sources")
    specs: list[SupplementarySpec] = []
    for raw_source in raw_sources:
        if not isinstance(raw_source, dict):
            raise RuntimeError("supplementary source must be a mapping")
        source = cast(dict[str, Any], raw_source)
        size = source.get("expected_size_bytes")
        if not isinstance(size, int) or size <= 0:
            raise RuntimeError("supplementary source has invalid expected_size_bytes")
        specs.append(
            SupplementarySpec(
                study_accession=_required_string(source, "study_accession"),
                member=_required_string(source, "member"),
                publication_doi=_required_string(source, "publication_doi"),
                source_kind=_required_string(source, "source_kind"),
                container=_required_string(source, "container"),
                container_url=_required_string(source, "container_url"),
                destination=Path(_required_string(source, "destination")),
                expected_size_bytes=size,
                data_level=_required_string(source, "data_level"),
                scientific_use=_required_string(source, "scientific_use"),
            )
        )
    keys = [(spec.study_accession, spec.member) for spec in specs]
    if len(set(keys)) != len(keys):
        raise RuntimeError("duplicate supplementary source key")
    return tuple(specs)


def _read_registry(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != SUPPLEMENTARY_FIELDS:
            raise RuntimeError("supplementary registry has invalid columns")
        rows = [dict(row) for row in reader]
    if any(
        None in row or any(value is None for value in row.values()) for row in rows
    ):
        raise RuntimeError("supplementary registry has malformed row")
    return rows


def _identity(spec: SupplementarySpec) -> dict[str, str]:
    return {
        "study_accession": spec.study_accession,
        "member": spec.member,
        "publication_doi": spec.publication_doi,
        "source_kind": spec.source_kind,
        "container": spec.container,
        "container_url": spec.container_url,
        "size_bytes": str(spec.expected_size_bytes),
        "data_level": spec.data_level,
        "scientific_use": spec.scientific_use,
        "downloader_version": DOWNLOADER_VERSION,
    }


def audit_supplementary_sources(
    config_path: Path = Path("configs/supplementary_sources_v1.yaml"),
    registry_path: Path = Path("data/registry/supplementary_sources.tsv"),
) -> tuple[SupplementarySource, ...]:
    """Verify every registered supplement source against its local SHA256."""
    specs = _load_specs(config_path)
    rows = _read_registry(registry_path)
    registry_keys = {(row["study_accession"], row["member"]) for row in rows}
    if registry_keys != {(spec.study_accession, spec.member) for spec in specs}:
        raise RuntimeError("supplementary registry does not match config")

    sources: list[SupplementarySource] = []
    for spec in specs:
        matches = [
            row
            for row in rows
            if row["study_accession"] == spec.study_accession
            and row["member"] == spec.member
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"supplementary source does not resolve once: {spec.member}"
            )
        row = matches[0]
        if any(row[field] != value for field, value in _identity(spec).items()):
            raise RuntimeError(
                f"supplementary source provenance mismatch: {spec.member}"
            )
        local_path = (registry_path.parent / row["local_path"]).resolve()
        try:
            size = int(row["size_bytes"])
        except ValueError as exc:
            raise RuntimeError("supplementary source has invalid size") from exc
        if (
            not row["retrieved_at"]
            or re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is None
            or not local_path.is_file()
            or local_path.stat().st_size != size
            or hash_file(local_path, "sha256") != row["sha256"]
        ):
            raise RuntimeError(
                f"supplementary source failed local audit: {spec.member}"
            )
        sources.append(
            SupplementarySource(
                study_accession=row["study_accession"],
                member=row["member"],
                publication_doi=row["publication_doi"],
                source_kind=row["source_kind"],
                container=row["container"],
                container_url=row["container_url"],
                local_path=local_path,
                retrieved_at=row["retrieved_at"],
                size_bytes=size,
                sha256=row["sha256"],
                data_level=row["data_level"],
                scientific_use=row["scientific_use"],
                downloader_version=row["downloader_version"],
            )
        )
    return tuple(sources)
