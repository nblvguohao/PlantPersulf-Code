"""Residue-resolved persulfidation sites from PXD006140 Dataset S3.

Dataset S3 (`ident_peptides`) records the authors' own custom persulfidation PTMs
on cysteine — `MOD:99998` (Sulfide), `MOD:99997`/`MOD:99996` (CN-Biotin-Sulfide) —
each with a peptide-internal position (`MOD:99998 C26`) and protein coordinates
(`inferredCoords`, e.g. `O03042{259-285}`). This parser emits a `site_ms` record
only for a cysteine that (1) carries an author-designated persulfidation PTM,
(2) is not from a decoy row, (3) maps to a single protein whose registered
sequence reproduces the reported peptide slice, and (4) is a cysteine at the
computed protein coordinate. Everything else — MSBT/`MOD:00110`/unmodified
cysteines, decoys — yields no site; isoform multi-mapped spectra become conflicts.
"""

from __future__ import annotations

import csv
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.proteomics.metadata import (
    ParseIssue,
    ReferenceSequence,
    SiteEvidence,
    SiteNormalization,
)

PARSER_VERSION = "dataset_s3_persulfidation_v1"
STUDY_ACCESSION = "PXD006140"
EVIDENCE_LEVEL = "site_ms"
PERSULFIDATION_MODS: dict[str, str] = {
    "99998": "Sulfide",
    "99997": "CN-Biotin-Sulfide",
    "99996": "CN-Biotin-Na-Sulfide",
}
DECOY_FLAGS = frozenset({"1", "true", "yes"})
REQUIRED_COLUMNS = (
    "scan_id",
    "peptide",
    "mods",
    "is_decoy",
    "inferredIDs",
    "inferredCoords",
    "FDR_peptide",
)

_MOD = re.compile(r"MOD:(\d+)\s+([A-Z])(\d+)")
_COORD = re.compile(r"([^{;]+)\{(\d+)-(\d+)\}")


def _parse_persulfidation_cys(mods: str) -> list[tuple[str, int]]:
    """Return (mod_code, peptide_position) for persulfidation PTMs on cysteine."""
    hits: list[tuple[str, int]] = []
    for code, residue, position in _MOD.findall(mods):
        if code in PERSULFIDATION_MODS and residue == "C":
            hits.append((code, int(position)))
    return hits


def _parse_coordinates(inferred_coords: str) -> list[tuple[str, int, int]]:
    return [
        (protein, int(start), int(stop))
        for protein, start, stop in _COORD.findall(inferred_coords)
    ]


def _issue(row: Mapping[str, str], reason: str, detail: str) -> ParseIssue:
    return ParseIssue(
        source_file="",
        spectrum_id=row["scan_id"],
        protein_accession_raw=row["inferredIDs"],
        modified_sequence=row["mods"],
        reason=reason,
        detail=detail,
    )


def parse_dataset_s3_sites(
    source_tsv: Path,
    references: Mapping[str, ReferenceSequence],
    source_sha256: str,
    study_accession: str = STUDY_ACCESSION,
) -> SiteNormalization:
    """Emit author-designated persulfidation sites verified against sequences."""
    with source_tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REQUIRED_COLUMNS:
            raise RuntimeError(f"Dataset S3 subset has invalid columns: {source_tsv}")
        rows = [dict(row) for row in reader]
    if any(
        None in row or any(value is None for value in row.values()) for row in rows
    ):
        raise RuntimeError(f"Dataset S3 subset has malformed row: {source_tsv}")

    sites: list[SiteEvidence] = []
    conflicts: list[ParseIssue] = []
    for row in rows:
        if row["is_decoy"].strip().lower() in DECOY_FLAGS:
            continue
        cys_mods = _parse_persulfidation_cys(row["mods"])
        if not cys_mods:
            continue
        coordinates = _parse_coordinates(row["inferredCoords"])
        if len(coordinates) != 1:
            conflicts.append(
                _issue(
                    row,
                    "isoform_multi_mapping",
                    "persulfidation spectrum maps to multiple proteins",
                )
            )
            continue
        protein, start, stop = coordinates[0]
        peptide = row["peptide"]
        reference = references.get(protein)
        if reference is None:
            conflicts.append(
                _issue(
                    row,
                    (
                        "isoform_sequence_unavailable"
                        if "-" in protein
                        else "sequence_unavailable"
                    ),
                    "exact registered protein sequence is unavailable",
                )
            )
            continue
        if stop - start + 1 != len(peptide) or (
            reference.sequence[start - 1 : stop] != peptide
        ):
            conflicts.append(
                _issue(
                    row,
                    "peptide_sequence_conflict",
                    "reported peptide does not match the registered sequence slice",
                )
            )
            continue
        for code, peptide_position in cys_mods:
            protein_position = start + peptide_position - 1
            if (
                peptide_position < 1
                or peptide_position > len(peptide)
                or peptide[peptide_position - 1] != "C"
                or reference.sequence[protein_position - 1] != "C"
            ):
                conflicts.append(
                    _issue(
                        row,
                        "cysteine_coordinate_conflict",
                        "persulfidation position is not a verified cysteine",
                    )
                )
                continue
            sites.append(
                SiteEvidence(
                    study_accession=study_accession,
                    sample_id="",
                    source_file=source_tsv.as_posix(),
                    spectrum_id=row["scan_id"],
                    peptide_sequence=peptide,
                    modified_sequence=row["mods"],
                    protein_accession_raw=protein,
                    protein_accession_canonical=(
                        "" if "-" in protein else protein
                    ),
                    cys_position_in_peptide=peptide_position,
                    cys_position_in_protein=protein_position,
                    modification_name_raw=PERSULFIDATION_MODS[code],
                    evidence_level=EVIDENCE_LEVEL,
                    quant_value="",
                    quant_unit="",
                    parser_version=PARSER_VERSION,
                    source_sha256=source_sha256,
                    reference_source_sha256=reference.source_sha256,
                )
            )
    return SiteNormalization(tuple(sites), tuple(conflicts))


MAXQUANT_PARSER_VERSION = "maxquant_sites_persulfidation_v1"
MAXQUANT_STUDY_ACCESSION = "PXD024061"
MAXQUANT_CLASS_I_THRESHOLD = 0.75
MAXQUANT_LOCALIZATION_UNIT = "maxquant_localization_prob"
MAXQUANT_REQUIRED_COLUMNS = (
    "Proteins",
    "Leading proteins",
    "Positions within proteins",
    "Position",
    "Amino acid",
    "Localization prob",
    "Reverse",
    "Potential contaminant",
    "Sequence window",
)


@dataclass(frozen=True)
class MaxquantPersulfidationSites:
    sites: tuple[SiteEvidence, ...]
    low_confidence: tuple[SiteEvidence, ...]
    conflicts: tuple[ParseIssue, ...]
    excluded: tuple[ParseIssue, ...]


def _maxquant_issue(row: Mapping[str, str], reason: str, detail: str) -> ParseIssue:
    return ParseIssue(
        source_file="",
        spectrum_id="",
        protein_accession_raw=row["Leading proteins"],
        modified_sequence=row["Sequence window"],
        reason=reason,
        detail=detail,
    )


def _sequence_window_matches(
    window: str,
    position: int,
    reference: ReferenceSequence,
) -> bool:
    """Verify a MaxQuant sequence window against the registered sequence.

    The window is centred on the modified residue; out-of-range flanks are
    padded with ``_`` and are skipped. Every non-padded residue must equal the
    reference residue at its computed protein coordinate.
    """
    if not window:
        return False
    center = len(window) // 2
    if window[center] != "C":
        return False
    for index, residue in enumerate(window):
        if residue == "_":
            continue
        protein_position = position + index - center
        if (
            protein_position < 1
            or protein_position > len(reference.sequence)
            or reference.sequence[protein_position - 1] != residue
        ):
            return False
    return True


def parse_maxquant_persulfidation_sites(
    source_tsv: Path,
    references: Mapping[str, ReferenceSequence],
    source_sha256: str,
    modification_name: str,
    study_accession: str = MAXQUANT_STUDY_ACCESSION,
) -> MaxquantPersulfidationSites:
    """Extract persulfidation sites from one MaxQuant PTM-site table.

    Class-I sites (localization prob >= threshold) are emitted as ``site_ms``;
    lower-probability sites are retained in a flagged tier; razor multi-protein
    rows are conflicts; decoy/contaminant rows are excluded.
    """
    with source_tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        if not set(MAXQUANT_REQUIRED_COLUMNS).issubset(fieldnames):
            raise RuntimeError(
                f"MaxQuant site table lacks required columns: {source_tsv}"
            )
        rows = [
            {column: row[column] for column in MAXQUANT_REQUIRED_COLUMNS}
            for row in reader
        ]
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise RuntimeError(f"MaxQuant site table has malformed row: {source_tsv}")

    sites: list[SiteEvidence] = []
    low_confidence: list[SiteEvidence] = []
    conflicts: list[ParseIssue] = []
    excluded: list[ParseIssue] = []
    for row in rows:
        if row["Amino acid"] != "C":
            continue
        if row["Reverse"] == "+":
            excluded.append(_maxquant_issue(row, "reverse_decoy", "decoy hit"))
            continue
        if row["Potential contaminant"] == "+":
            excluded.append(
                _maxquant_issue(row, "potential_contaminant", "contaminant hit")
            )
            continue
        if ";" in row["Proteins"]:
            conflicts.append(
                _maxquant_issue(
                    row,
                    "razor_multi_mapping",
                    "site maps to multiple proteins",
                )
            )
            continue
        protein = row["Proteins"]
        reference = references.get(protein)
        if reference is None:
            conflicts.append(
                _maxquant_issue(
                    row,
                    (
                        "isoform_sequence_unavailable"
                        if "-" in protein
                        else "sequence_unavailable"
                    ),
                    "exact registered protein sequence is unavailable",
                )
            )
            continue
        try:
            position = int(row["Position"])
            localization = float(row["Localization prob"])
        except ValueError:
            conflicts.append(
                _maxquant_issue(
                    row,
                    "invalid_site_field",
                    "Position and Localization prob must be numeric",
                )
            )
            continue
        window = row["Sequence window"]
        if (
            position < 1
            or position > len(reference.sequence)
            or reference.sequence[position - 1] != "C"
            or not _sequence_window_matches(window, position, reference)
        ):
            conflicts.append(
                _maxquant_issue(
                    row,
                    "cysteine_coordinate_conflict",
                    "site position/window does not verify against the sequence",
                )
            )
            continue
        site = SiteEvidence(
            study_accession=study_accession,
            sample_id="",
            source_file=source_tsv.as_posix(),
            spectrum_id="",
            peptide_sequence=window,
            modified_sequence=window,
            protein_accession_raw=protein,
            protein_accession_canonical="" if "-" in protein else protein,
            cys_position_in_peptide=len(window) // 2 + 1,
            cys_position_in_protein=position,
            modification_name_raw=modification_name,
            evidence_level=EVIDENCE_LEVEL,
            quant_value=row["Localization prob"],
            quant_unit=MAXQUANT_LOCALIZATION_UNIT,
            parser_version=MAXQUANT_PARSER_VERSION,
            source_sha256=source_sha256,
            reference_source_sha256=reference.source_sha256,
        )
        if localization >= MAXQUANT_CLASS_I_THRESHOLD:
            sites.append(site)
        else:
            low_confidence.append(site)
    return MaxquantPersulfidationSites(
        tuple(sites),
        tuple(low_confidence),
        tuple(conflicts),
        tuple(excluded),
    )
