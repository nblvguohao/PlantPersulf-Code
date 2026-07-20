"""Immutable records for provenance-preserving proteomics parsing."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

PARSER_VERSION = "omssa_csv_v1"
EVIDENCE_LEVEL = "psm_coordinate_only"

PSM_FIELDS = (
    "study_accession",
    "sample_id",
    "source_file",
    "spectrum_id",
    "peptide_sequence",
    "modified_sequence",
    "search_engine_accession_raw",
    "protein_accession_raw",
    "protein_accession_canonical",
    "start_one_based",
    "stop_one_based",
    "defline_raw",
    "modification_name_raw",
    "evidence_level",
    "quant_value",
    "quant_unit",
    "parser_version",
    "source_sha256",
)
SITE_FIELDS = (
    "study_accession",
    "sample_id",
    "source_file",
    "spectrum_id",
    "peptide_sequence",
    "modified_sequence",
    "protein_accession_raw",
    "protein_accession_canonical",
    "cys_position_in_peptide",
    "cys_position_in_protein",
    "modification_name_raw",
    "evidence_level",
    "quant_value",
    "quant_unit",
    "parser_version",
    "source_sha256",
)
ISSUE_FIELDS = (
    "source_file",
    "spectrum_id",
    "protein_accession_raw",
    "modified_sequence",
    "reason",
    "detail",
)
OUTPUT_FILENAMES = (
    "psms.tsv",
    "sites.tsv",
    "conflicts.tsv",
    "excluded.tsv",
    "missing_metadata.tsv",
    "manifest.json",
)


@dataclass(frozen=True)
class ProteomicsSource:
    study_accession: str
    source_file: Path
    registry_file: Path
    source_sha256: str
    parser_version: str = PARSER_VERSION
    evidence_scope: str = "site_coordinates"
    sample_id: str = ""


@dataclass(frozen=True)
class PeptideSpectrumMatch:
    study_accession: str
    sample_id: str
    source_file: str
    spectrum_id: str
    peptide_sequence: str
    modified_sequence: str
    search_engine_accession_raw: str
    protein_accession_raw: str
    protein_accession_canonical: str
    start_one_based: int
    stop_one_based: int
    defline_raw: str
    modification_name_raw: str
    evidence_level: str
    quant_value: str
    quant_unit: str
    parser_version: str
    source_sha256: str


@dataclass(frozen=True)
class ParseIssue:
    source_file: str
    spectrum_id: str
    protein_accession_raw: str
    modified_sequence: str
    reason: str
    detail: str


@dataclass(frozen=True)
class ParsedProteomics:
    psms: tuple[PeptideSpectrumMatch, ...]
    conflicts: tuple[ParseIssue, ...]
    excluded: tuple[ParseIssue, ...]
    missing: tuple[ParseIssue, ...]


@dataclass(frozen=True)
class ReferenceSequence:
    accession: str
    sequence: str
    sequence_version: str
    source_file: Path
    source_sha256: str


@dataclass(frozen=True)
class SiteEvidence:
    study_accession: str
    sample_id: str
    source_file: str
    spectrum_id: str
    peptide_sequence: str
    modified_sequence: str
    protein_accession_raw: str
    protein_accession_canonical: str
    cys_position_in_peptide: int
    cys_position_in_protein: int
    modification_name_raw: str
    evidence_level: str
    quant_value: str
    quant_unit: str
    parser_version: str
    source_sha256: str
    reference_source_sha256: str


@dataclass(frozen=True)
class SiteNormalization:
    sites: tuple[SiteEvidence, ...]
    conflicts: tuple[ParseIssue, ...]


@dataclass(frozen=True)
class ParseSummary:
    psm_count: int
    site_count: int
    conflict_count: int
    excluded_count: int
    missing_metadata_count: int
    missing_sequence_count: int
