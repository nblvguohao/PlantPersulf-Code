"""Conservative cysteine coordinates verified against registered FASTA bytes."""

from __future__ import annotations

import csv
import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any, cast

from plantpersulf.proteomics.metadata import (
    ISSUE_FIELDS,
    OUTPUT_FILENAMES,
    PARSER_VERSION,
    PSM_FIELDS,
    SITE_FIELDS,
    ParseIssue,
    ParseSummary,
    PeptideSpectrumMatch,
    ReferenceSequence,
    SiteEvidence,
    SiteNormalization,
)
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file

_FASTA_HEADER = re.compile(
    r"^>(?:sp|tr)\|([^|]+)\|.*(?:^|\s)SV=([0-9]+)(?:\s|$)"
)
_SEQUENCE = re.compile(r"^[A-Z]+$")


def load_reference_sequence(
    path: Path,
    registry_path: Path,
) -> ReferenceSequence:
    """Load exactly one registered UniProt FASTA record."""
    assert_registered_input(path, registry_path)
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2 or any(line.startswith(">") for line in lines[1:]):
        raise RuntimeError(f"reference FASTA requires exactly one record: {path}")
    match = _FASTA_HEADER.match(lines[0])
    if match is None:
        raise RuntimeError(f"reference FASTA has invalid UniProt header: {path}")
    sequence = "".join(line for line in lines[1:] if line)
    if not sequence or _SEQUENCE.fullmatch(sequence) is None:
        raise RuntimeError(f"reference FASTA has invalid sequence: {path}")
    return ReferenceSequence(
        accession=match.group(1),
        sequence=sequence,
        sequence_version=match.group(2),
        source_file=path,
        source_sha256=hash_file(path, "sha256"),
    )


def _conflict(psm: PeptideSpectrumMatch, reason: str, detail: str) -> ParseIssue:
    return ParseIssue(
        source_file=psm.source_file,
        spectrum_id=psm.spectrum_id,
        protein_accession_raw=psm.protein_accession_raw,
        modified_sequence=psm.modified_sequence,
        reason=reason,
        detail=detail,
    )


def normalize_sites(
    psm: PeptideSpectrumMatch,
    references: Mapping[str, ReferenceSequence],
    evidence_scope: str,
) -> SiteNormalization:
    """Emit sites only after the complete reported peptide slice verifies."""
    if evidence_scope == "protein_level_only":
        return SiteNormalization((), ())
    if evidence_scope != "site_coordinates":
        raise RuntimeError(f"unsupported evidence scope: {evidence_scope}")

    reference = references.get(psm.protein_accession_raw)
    if reference is None:
        reason = (
            "isoform_sequence_unavailable"
            if "-" in psm.protein_accession_raw
            else "sequence_unavailable"
        )
        return SiteNormalization(
            (),
            (
                _conflict(
                    psm,
                    reason,
                    "exact registered protein sequence is unavailable",
                ),
            ),
        )

    expected_length = psm.stop_one_based - psm.start_one_based + 1
    if expected_length != len(psm.peptide_sequence):
        return SiteNormalization(
            (),
            (
                _conflict(
                    psm,
                    "peptide_length_conflict",
                    "inclusive OMSSA coordinates do not equal peptide length",
                ),
            ),
        )
    observed = reference.sequence[psm.start_one_based - 1 : psm.stop_one_based]
    if observed != psm.peptide_sequence:
        return SiteNormalization(
            (),
            (
                _conflict(
                    psm,
                    "peptide_sequence_conflict",
                    "reported peptide does not match the registered sequence slice",
                ),
            ),
        )

    sites = tuple(
        SiteEvidence(
            study_accession=psm.study_accession,
            sample_id=psm.sample_id,
            source_file=psm.source_file,
            spectrum_id=psm.spectrum_id,
            peptide_sequence=psm.peptide_sequence,
            modified_sequence=psm.modified_sequence,
            protein_accession_raw=psm.protein_accession_raw,
            protein_accession_canonical=psm.protein_accession_canonical,
            cys_position_in_peptide=offset + 1,
            cys_position_in_protein=psm.start_one_based + offset,
            modification_name_raw=psm.modification_name_raw,
            evidence_level=psm.evidence_level,
            quant_value=psm.quant_value,
            quant_unit=psm.quant_unit,
            parser_version=psm.parser_version,
            source_sha256=psm.source_sha256,
            reference_source_sha256=reference.source_sha256,
        )
        for offset, residue in enumerate(psm.peptide_sequence)
        if residue == "C"
    )
    return SiteNormalization(sites, ())


def _read_tsv(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise RuntimeError(f"site output has invalid columns: {path}")
        rows = [dict(row) for row in reader]
    if any(
        None in row or any(value is None for value in row.values())
        for row in rows
    ):
        raise RuntimeError(f"site output has malformed row: {path}")
    return cast(list[dict[str, str]], rows)


def _manifest_mapping(value: object, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise RuntimeError(f"site manifest requires mapping: {context}")
    return cast(dict[str, Any], value)


def _manifest_items(value: object, context: str) -> list[Any]:
    if not isinstance(value, list):
        raise RuntimeError(f"site manifest requires list: {context}")
    return value


def audit_site_output(
    accession: str,
    output_root: Path = Path("data/interim"),
    registry_dir: Path = Path("data/registry"),
) -> ParseSummary:
    """Audit deterministic parser output and recheck every emitted site."""
    del registry_dir  # Registries are fixed in the manifest for exact replay.
    normalized = accession.strip().upper()
    output_directory = output_root / normalized / "proteomics_parser_v1"
    observed_files = {
        path.name for path in output_directory.iterdir() if path.is_file()
    }
    if observed_files != set(OUTPUT_FILENAMES):
        raise RuntimeError("site output file set differs from parser schema")
    manifest_path = output_directory / "manifest.json"
    try:
        parsed: object = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"invalid site manifest: {manifest_path}") from exc
    manifest = _manifest_mapping(parsed, "root")
    if (
        manifest.get("schema_version") != 1
        or manifest.get("study_accession") != normalized
        or manifest.get("parser_version") != PARSER_VERSION
    ):
        raise RuntimeError("site manifest identity mismatch")
    if (
        manifest.get("biological_values_modified") is not False
        or manifest.get("labels_created") is not False
    ):
        raise RuntimeError("site manifest violates scientific integrity policy")

    source_hashes: set[str] = set()
    for raw_source in _manifest_items(manifest.get("sources"), "sources"):
        source = _manifest_mapping(raw_source, "sources[]")
        path = Path(str(source.get("path", "")))
        registry_path = Path(str(source.get("registry_path", "")))
        expected_sha256 = str(source.get("sha256", ""))
        assert_registered_input(path, registry_path)
        if hash_file(path, "sha256") != expected_sha256:
            raise RuntimeError(f"source SHA256 mismatch: {path}")
        source_hashes.add(expected_sha256)

    references: dict[str, ReferenceSequence] = {}
    for raw_reference in _manifest_items(
        manifest.get("references"),
        "references",
    ):
        entry = _manifest_mapping(raw_reference, "references[]")
        path = Path(str(entry.get("path", "")))
        registry_path = Path(str(entry.get("registry_path", "")))
        reference = load_reference_sequence(path, registry_path)
        if (
            reference.accession != entry.get("accession")
            or reference.sequence_version != entry.get("sequence_version")
            or reference.source_sha256 != entry.get("sha256")
        ):
            raise RuntimeError(f"reference sequence identity mismatch: {path}")
        references[reference.accession] = reference

    output_entries = _manifest_items(manifest.get("outputs"), "outputs")
    declared_output_names: set[str] = set()
    for raw_output in output_entries:
        entry = _manifest_mapping(raw_output, "outputs[]")
        name = str(entry.get("file", ""))
        declared_output_names.add(name)
        path = output_directory / name
        if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
            raise RuntimeError(f"output SHA256 mismatch: {path}")
    if declared_output_names != set(OUTPUT_FILENAMES) - {"manifest.json"}:
        raise RuntimeError("site manifest output list mismatch")

    psms = _read_tsv(output_directory / "psms.tsv", PSM_FIELDS)
    sites = _read_tsv(output_directory / "sites.tsv", SITE_FIELDS)
    conflicts = _read_tsv(output_directory / "conflicts.tsv", ISSUE_FIELDS)
    excluded = _read_tsv(output_directory / "excluded.tsv", ISSUE_FIELDS)
    missing = _read_tsv(
        output_directory / "missing_metadata.tsv",
        ISSUE_FIELDS,
    )
    summary = ParseSummary(
        psm_count=len(psms),
        site_count=len(sites),
        conflict_count=len(conflicts),
        excluded_count=len(excluded),
        missing_metadata_count=len(missing),
        missing_sequence_count=sum(
            row["reason"]
            in {"sequence_unavailable", "isoform_sequence_unavailable"}
            for row in conflicts
        ),
    )
    if _manifest_mapping(manifest.get("counts"), "counts") != {
        key: value for key, value in summary.__dict__.items()
    }:
        raise RuntimeError("site manifest counts mismatch")

    for row in sites:
        if (
            not row["source_sha256"]
            or row["source_sha256"] not in source_hashes
            or not row["parser_version"]
        ):
            raise RuntimeError("site lacks registered parser provenance")
        if row["evidence_level"] != "psm_coordinate_only" or (
            "negative" in row["evidence_level"].lower()
        ):
            raise RuntimeError("site output contains forbidden evidence level")
        site_reference = references.get(row["protein_accession_raw"])
        if site_reference is None:
            raise RuntimeError("emitted site lacks exact registered sequence")
        try:
            peptide_position = int(row["cys_position_in_peptide"])
            protein_position = int(row["cys_position_in_protein"])
        except ValueError as exc:
            raise RuntimeError("site output has invalid coordinate") from exc
        peptide = row["peptide_sequence"]
        if (
            peptide_position < 1
            or peptide_position > len(peptide)
            or peptide[peptide_position - 1] != "C"
        ):
            raise RuntimeError("site peptide cysteine coordinate mismatch")
        start_one_based = protein_position - peptide_position + 1
        observed = site_reference.sequence[
            start_one_based - 1 : start_one_based - 1 + len(peptide)
        ]
        if (
            observed != peptide
            or site_reference.sequence[protein_position - 1] != "C"
        ):
            raise RuntimeError("site protein coordinate mismatch")
    return summary
