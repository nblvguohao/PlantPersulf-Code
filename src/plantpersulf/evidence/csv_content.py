"""Exact-schema readers for registered PXD035795 result CSV files."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.provenance.hashing import hash_file

PEPTIDE_HEADER = (
    "Peptide",
    "Quality",
    "Significance",
    "Avg. ppm",
    "Avg. Area",
    "NPC_rep1",
    "NPC_rep2",
    "NPC_rep3",
    "APC_rep1",
    "APC_rep2",
    "APC_rep3",
    "Sample Profile (Ratio)",
    "0h",
    "3days",
    "Group Profile (Ratio)",
    "Max Ratio",
    "#Vector",
    "Accession",
    "PTM",
)
PROTEIN_HEADER = (
    "Protein Group",
    "Protein ID",
    "Accession",
    "Significance",
    "Coverage (%)",
    "#Peptides",
    "#Unique",
    "PTM",
    "NPC_rep1 Area",
    "NPC_rep2 Area",
    "NPC_rep3 Area",
    "APC_rep1 Area",
    "APC_rep2 Area",
    "APC_rep3 Area",
    "NPC_rep1 Area (top-3 peptides)",
    "NPC_rep2 Area (top-3 peptides)",
    "NPC_rep3 Area (top-3 peptides)",
    "APC_rep1 Area (top-3 peptides)",
    "APC_rep2 Area (top-3 peptides)",
    "APC_rep3 Area (top-3 peptides)",
    "Sample Profile (Ratio)",
    "0h Area",
    "3days Area",
    "0h Area (top-3 peptides)",
    "3days Area (top-3 peptides)",
    "Group Profile (Ratio)",
    "Avg. Mass",
    "Description",
)


@dataclass(frozen=True)
class CsvContentRecord:
    row_number: int
    raw_fields: tuple[str, ...]
    declared_modifications: tuple[str, ...]
    evidence_class: str
    site_localization_status: str
    source_sha256: str
    source_locator: str


@dataclass(frozen=True)
class CsvAudit:
    header: tuple[str, ...]
    rows: tuple[CsvContentRecord, ...]
    record_type: str
    source_sha256: str
    encoding: str
    delimiter: str


def _declared_modifications(value: str) -> tuple[str, ...]:
    return tuple(part.strip() for part in value.split(";") if part.strip())


def _parse_csv(
    path: Path,
    source_sha256: str,
    expected_header: tuple[str, ...],
    expected_rows: int,
    record_type: str,
    ptm_index: int,
    site_localization_status: str,
) -> CsvAudit:
    if hash_file(path, "sha256") != source_sha256:
        raise RuntimeError(f"CSV source SHA256 mismatch: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        parsed = [tuple(row) for row in csv.reader(handle)]
    if not parsed or parsed[0] != expected_header:
        raise RuntimeError(f"{record_type} header mismatch")
    raw_rows = parsed[1:]
    if len(raw_rows) != expected_rows or any(
        len(row) != len(expected_header) for row in raw_rows
    ):
        raise RuntimeError(f"{record_type} row shape mismatch")
    records = tuple(
        CsvContentRecord(
            row_number=row_number,
            raw_fields=row,
            declared_modifications=_declared_modifications(row[ptm_index]),
            evidence_class="identification_only",
            site_localization_status=site_localization_status,
            source_sha256=source_sha256,
            source_locator=f"{path.name}:row={row_number}",
        )
        for row_number, row in enumerate(raw_rows, start=2)
    )
    return CsvAudit(
        header=parsed[0],
        rows=records,
        record_type=record_type,
        source_sha256=source_sha256,
        encoding="utf-8-sig",
        delimiter=",",
    )


def parse_peptide_csv(path: Path, source_sha256: str) -> CsvAudit:
    """Preserve every PXD035795 peptide CSV field without site inference."""
    return _parse_csv(
        path=path,
        source_sha256=source_sha256,
        expected_header=PEPTIDE_HEADER,
        expected_rows=9326,
        record_type="peptide_csv",
        ptm_index=18,
        site_localization_status="not_available",
    )


def parse_protein_csv(path: Path, source_sha256: str) -> CsvAudit:
    """Preserve every PXD035795 protein CSV field without site promotion."""
    return _parse_csv(
        path=path,
        source_sha256=source_sha256,
        expected_header=PROTEIN_HEADER,
        expected_rows=1461,
        record_type="protein_csv",
        ptm_index=7,
        site_localization_status="not_applicable",
    )
