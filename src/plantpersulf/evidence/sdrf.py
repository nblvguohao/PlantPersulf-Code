"""Positional SDRF audit without silent filename normalization."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.provenance.hashing import hash_file

EXPECTED_SDRF_HEADER = (
    "source name",
    "characteristics[organism]",
    "characteristics[organism part]",
    "characteristics[cell type]",
    "characteristics[ancestry category]",
    "characteristics[age]",
    "characteristics[sex]",
    "characteristics[disease]",
    "characteristics[biological replicate]",
    "characteristics[individual]",
    "assay name",
    "technology type",
    "comment[data file]",
    "comment[file uri]",
    "comment[label]",
    "comment[instrument]",
    "comment[cleavage agent details]",
    "comment[precursor mass tolerance]",
    "comment[fragment mass tolerance]",
    "comment[modification parameters]",
    "comment[modification parameters]",
    "comment[modification parameters]",
    "comment[modification parameters]",
    "comment[modification parameters]",
    "comment[fraction identifier]",
    "Factor Value[Treatment]",
)


@dataclass(frozen=True)
class SdrfRecord:
    row_number: int
    source_name: str
    assay_name: str
    replicate: str
    treatment: str
    label: str
    instrument: str
    data_file_raw: str
    modifications: tuple[str, ...]
    raw_fields: tuple[str, ...]


@dataclass(frozen=True)
class SdrfFileMapping:
    row_number: int
    source_name: str
    assay_name: str
    data_file: str
    mapping_status: str
    matched_registry_file: str
    source_sha256: str
    source_locator: str


@dataclass(frozen=True)
class SdrfMappingConflict:
    row_number: int
    source_name: str
    assay_name: str
    data_file: str
    conflict_type: str
    detail: str
    source_sha256: str
    source_locator: str


@dataclass(frozen=True)
class SdrfAudit:
    header: tuple[str, ...]
    rows: tuple[SdrfRecord, ...]
    file_mappings: tuple[SdrfFileMapping, ...]
    conflicts: tuple[SdrfMappingConflict, ...]
    duplicate_header_positions: dict[str, tuple[int, ...]]
    source_sha256: str


def _registry_names(files_registry_path: Path) -> set[str]:
    with files_registry_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {
            "dataset_accession",
            "file_name",
        }.issubset(reader.fieldnames):
            raise RuntimeError("files registry has invalid columns")
        return {
            row["file_name"]
            for row in reader
            if row["dataset_accession"] == "PXD035795"
        }


def _duplicate_positions(header: tuple[str, ...]) -> dict[str, tuple[int, ...]]:
    positions: dict[str, list[int]] = {}
    for index, name in enumerate(header):
        positions.setdefault(name, []).append(index)
    return {
        name: tuple(indices)
        for name, indices in positions.items()
        if len(indices) > 1
    }


def parse_sdrf(
    path: Path,
    source_sha256: str,
    files_registry_path: Path,
) -> SdrfAudit:
    """Parse the registered PXD035795 SDRF as positional string rows."""
    if hash_file(path, "sha256") != source_sha256:
        raise RuntimeError(f"SDRF source SHA256 mismatch: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        parsed = [tuple(row) for row in csv.reader(handle, delimiter="\t")]
    if not parsed or parsed[0] != EXPECTED_SDRF_HEADER:
        raise RuntimeError("PXD035795 SDRF header mismatch")
    raw_rows = parsed[1:]
    if len(raw_rows) != 6 or any(
        len(row) != len(EXPECTED_SDRF_HEADER) for row in raw_rows
    ):
        raise RuntimeError("PXD035795 SDRF row shape mismatch")

    registry_names = _registry_names(files_registry_path)
    records: list[SdrfRecord] = []
    mappings: list[SdrfFileMapping] = []
    conflicts: list[SdrfMappingConflict] = []
    for row_number, row in enumerate(raw_rows, start=2):
        record = SdrfRecord(
            row_number=row_number,
            source_name=row[0],
            assay_name=row[10],
            replicate=row[8],
            treatment=row[25],
            label=row[14],
            instrument=row[15],
            data_file_raw=row[12],
            modifications=tuple(row[19:24]),
            raw_fields=row,
        )
        records.append(record)
        components = tuple(part.strip() for part in record.data_file_raw.split(","))
        if len(components) != 2 or any(not component for component in components):
            raise RuntimeError(f"SDRF data-file list is invalid at row {row_number}")
        for component in components:
            matched = component if component in registry_names else ""
            status = "exact" if matched else "conflict"
            locator = f"SDRF.txt:row={row_number}:data_file={component}"
            mappings.append(
                SdrfFileMapping(
                    row_number=row_number,
                    source_name=record.source_name,
                    assay_name=record.assay_name,
                    data_file=component,
                    mapping_status=status,
                    matched_registry_file=matched,
                    source_sha256=source_sha256,
                    source_locator=locator,
                )
            )
            if not matched:
                conflicts.append(
                    SdrfMappingConflict(
                        row_number=row_number,
                        source_name=record.source_name,
                        assay_name=record.assay_name,
                        data_file=component,
                        conflict_type="unmatched_exact_filename",
                        detail="no exact PXD035795 file_name exists in files.tsv",
                        source_sha256=source_sha256,
                        source_locator=locator,
                    )
                )
    return SdrfAudit(
        header=parsed[0],
        rows=tuple(records),
        file_mappings=tuple(mappings),
        conflicts=tuple(conflicts),
        duplicate_header_positions=_duplicate_positions(parsed[0]),
        source_sha256=source_sha256,
    )
