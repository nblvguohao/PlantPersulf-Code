"""Streaming parser for registered OMSSA CSV search results."""

from __future__ import annotations

import csv
import json
import re
import shutil
import tempfile
from dataclasses import asdict
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.proteomics.metadata import (
    EVIDENCE_LEVEL,
    ISSUE_FIELDS,
    OUTPUT_FILENAMES,
    PARSER_VERSION,
    PSM_FIELDS,
    SITE_FIELDS,
    ParsedProteomics,
    ParseIssue,
    ParseSummary,
    PeptideSpectrumMatch,
    ProteomicsSource,
    SiteEvidence,
)
from plantpersulf.proteomics.site_normalizer import (
    load_reference_sequence,
    normalize_sites,
)
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file

OMSSA_COLUMNS = (
    "Spectrum number",
    "Filename/id",
    "Peptide",
    "E-value",
    "Mass",
    "gi",
    "Accession",
    "Start",
    "Stop",
    "Defline",
    "Mods",
    "Charge",
    "Theo Mass",
    "P-value",
    "NIST score",
)

_TARGET_DEFLINE = re.compile(r"^(?:sp|tr)\|([^|]+)\|")
_ANY_DEFLINE = re.compile(r"^(?:rev_)?(?:sp|tr)\|([^|]+)\|")


def _accession_from_defline(defline: str) -> str:
    match = _ANY_DEFLINE.match(defline)
    return match.group(1) if match is not None else ""


def _issue(
    source: ProteomicsSource,
    spectrum_id: str,
    modified_sequence: str,
    defline: str,
    reason: str,
    detail: str,
) -> ParseIssue:
    return ParseIssue(
        source_file=source.source_file.as_posix(),
        spectrum_id=spectrum_id,
        protein_accession_raw=_accession_from_defline(defline),
        modified_sequence=modified_sequence,
        reason=reason,
        detail=detail,
    )


def _exclusion_reason(defline: str) -> str:
    lowered = defline.lower()
    if lowered.startswith(("rev_sp|", "rev_tr|")) or " reversed " in lowered:
        return "decoy"
    if lowered.startswith(("con_", "cont_", "contaminant")):
        return "contaminant"
    return ""


def parse_omssa(source: ProteomicsSource) -> ParsedProteomics:
    """Parse one registered OMSSA export without deduplicating mappings."""
    assert_registered_input(source.source_file, source.registry_file)
    actual_sha256 = hash_file(source.source_file, "sha256")
    if actual_sha256 != source.source_sha256:
        raise RuntimeError(
            f"OMSSA source SHA256 differs from approved source: {source.source_file}"
        )

    psms: list[PeptideSpectrumMatch] = []
    conflicts: list[ParseIssue] = []
    excluded: list[ParseIssue] = []
    missing: list[ParseIssue] = []
    with source.source_file.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        if tuple(reader.fieldnames or ()) != OMSSA_COLUMNS:
            raise RuntimeError("OMSSA source has invalid columns")
        if not source.sample_id:
            missing.append(
                _issue(
                    source,
                    "",
                    "",
                    "",
                    "missing_sample_metadata",
                    "no registered sample maps to this source file",
                )
            )
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                conflicts.append(
                    _issue(
                        source,
                        str(row.get("Spectrum number", "")),
                        str(row.get("Peptide", "")),
                        str(row.get("Defline", "")),
                        "invalid_row",
                        "OMSSA row does not match the registered header",
                    )
                )
                continue
            spectrum_id = row["Spectrum number"]
            modified_sequence = row["Peptide"]
            defline = row["Defline"]
            exclusion_reason = _exclusion_reason(defline)
            if exclusion_reason:
                excluded.append(
                    _issue(
                        source,
                        spectrum_id,
                        modified_sequence,
                        defline,
                        exclusion_reason,
                        "record excluded before biological output",
                    )
                )
                continue
            accession_match = _TARGET_DEFLINE.match(defline)
            if accession_match is None:
                conflicts.append(
                    _issue(
                        source,
                        spectrum_id,
                        modified_sequence,
                        defline,
                        "unparseable_accession",
                        "defline lacks an exact sp or tr UniProt accession",
                    )
                )
                continue
            protein_accession = accession_match.group(1)
            try:
                start = int(row["Start"])
                stop = int(row["Stop"])
            except ValueError:
                conflicts.append(
                    _issue(
                        source,
                        spectrum_id,
                        modified_sequence,
                        defline,
                        "invalid_coordinate",
                        "Start and Stop must be integers",
                    )
                )
                continue
            psms.append(
                PeptideSpectrumMatch(
                    study_accession=source.study_accession,
                    sample_id=source.sample_id,
                    source_file=source.source_file.as_posix(),
                    spectrum_id=spectrum_id,
                    peptide_sequence=modified_sequence.upper(),
                    modified_sequence=modified_sequence,
                    search_engine_accession_raw=row["Accession"],
                    protein_accession_raw=protein_accession,
                    protein_accession_canonical=(
                        "" if "-" in protein_accession else protein_accession
                    ),
                    start_one_based=start,
                    stop_one_based=stop,
                    defline_raw=defline,
                    modification_name_raw=row["Mods"],
                    evidence_level=EVIDENCE_LEVEL,
                    quant_value="",
                    quant_unit="",
                    parser_version=source.parser_version,
                    source_sha256=actual_sha256,
                )
            )
    return ParsedProteomics(
        psms=tuple(psms),
        conflicts=tuple(conflicts),
        excluded=tuple(excluded),
        missing=tuple(missing),
    )


def _write_tsv(
    path: Path,
    fieldnames: tuple[str, ...],
    records: tuple[PeptideSpectrumMatch | SiteEvidence | ParseIssue, ...],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hash_file(path, "sha256")
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _summary(
    psms: tuple[PeptideSpectrumMatch, ...],
    sites: tuple[SiteEvidence, ...],
    conflicts: tuple[ParseIssue, ...],
    excluded: tuple[ParseIssue, ...],
    missing: tuple[ParseIssue, ...],
) -> ParseSummary:
    return ParseSummary(
        psm_count=len(psms),
        site_count=len(sites),
        conflict_count=len(conflicts),
        excluded_count=len(excluded),
        missing_metadata_count=len(missing),
        missing_sequence_count=sum(
            issue.reason
            in {"sequence_unavailable", "isoform_sequence_unavailable"}
            for issue in conflicts
        ),
    )


def publish_parsed_proteomics(
    sources: tuple[ProteomicsSource, ...],
    reference_paths: tuple[Path, ...],
    reference_registry: Path,
    output_directory: Path,
) -> ParseSummary:
    """Parse registered sources and atomically publish deterministic TSV files."""
    if not sources:
        raise RuntimeError("proteomics publication requires at least one source")
    accessions = {source.study_accession for source in sources}
    if len(accessions) != 1:
        raise RuntimeError("proteomics sources require one study accession")
    references = tuple(
        load_reference_sequence(path, reference_registry)
        for path in sorted(reference_paths)
    )
    reference_map = {reference.accession: reference for reference in references}
    if len(reference_map) != len(references):
        raise RuntimeError("duplicate reference sequence accession")

    psms: list[PeptideSpectrumMatch] = []
    sites: list[SiteEvidence] = []
    conflicts: list[ParseIssue] = []
    excluded: list[ParseIssue] = []
    missing: list[ParseIssue] = []
    for source in sources:
        parsed = parse_omssa(source)
        psms.extend(parsed.psms)
        conflicts.extend(parsed.conflicts)
        excluded.extend(parsed.excluded)
        missing.extend(parsed.missing)
        for psm in parsed.psms:
            if "C" not in psm.peptide_sequence:
                continue
            normalized = normalize_sites(
                psm,
                reference_map,
                source.evidence_scope,
            )
            sites.extend(normalized.sites)
            conflicts.extend(normalized.conflicts)

    psm_records = tuple(psms)
    site_records = tuple(sites)
    conflict_records = tuple(conflicts)
    excluded_records = tuple(excluded)
    missing_records = tuple(missing)
    summary = _summary(
        psm_records,
        site_records,
        conflict_records,
        excluded_records,
        missing_records,
    )

    output_directory.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = Path(
        tempfile.mkdtemp(
            prefix=f".{output_directory.name}.",
            dir=output_directory.parent,
        )
    )
    temporary_directory: Path | None = temporary_path
    try:
        _write_tsv(temporary_path / "psms.tsv", PSM_FIELDS, psm_records)
        _write_tsv(temporary_path / "sites.tsv", SITE_FIELDS, site_records)
        _write_tsv(
            temporary_path / "conflicts.tsv",
            ISSUE_FIELDS,
            conflict_records,
        )
        _write_tsv(
            temporary_path / "excluded.tsv",
            ISSUE_FIELDS,
            excluded_records,
        )
        _write_tsv(
            temporary_path / "missing_metadata.tsv",
            ISSUE_FIELDS,
            missing_records,
        )
        output_hashes = [
            {
                "file": name,
                "sha256": hash_file(temporary_path / name, "sha256"),
            }
            for name in OUTPUT_FILENAMES
            if name != "manifest.json"
        ]
        manifest: dict[str, object] = {
            "schema_version": 1,
            "study_accession": next(iter(accessions)),
            "parser_version": PARSER_VERSION,
            "sources": [
                {
                    "path": source.source_file.as_posix(),
                    "registry_path": source.registry_file.as_posix(),
                    "sha256": source.source_sha256,
                    "evidence_scope": source.evidence_scope,
                }
                for source in sources
            ],
            "references": [
                {
                    "accession": reference.accession,
                    "path": reference.source_file.as_posix(),
                    "registry_path": reference_registry.as_posix(),
                    "sequence_version": reference.sequence_version,
                    "sha256": reference.source_sha256,
                }
                for reference in references
            ],
            "outputs": output_hashes,
            "counts": asdict(summary),
            "biological_values_modified": False,
            "labels_created": False,
        }
        serialized = json.dumps(
            manifest,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
        (temporary_path / "manifest.json").write_bytes(
            f"{serialized}\n".encode()
        )
        if output_directory.exists():
            if _tree_hashes(output_directory) != _tree_hashes(temporary_path):
                raise RuntimeError(
                    f"existing parser output differs: {output_directory}"
                )
            return summary
        temporary_path.replace(output_directory)
        temporary_directory = None
        return summary
    finally:
        if temporary_directory is not None:
            shutil.rmtree(temporary_directory)


def _load_mapping(path: Path) -> dict[str, Any]:
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid proteomics selection config: {path}") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError("proteomics selection config must be a mapping")
    return cast(dict[str, Any], loaded)


def _resolve_sources(
    accession: str,
    registry_dir: Path,
    selection_path: Path,
) -> tuple[ProteomicsSource, ...]:
    config = _load_mapping(selection_path)
    if config.get("version") != 1:
        raise RuntimeError("proteomics selection config requires version 1")
    approved = config.get("approved")
    if not isinstance(approved, dict):
        raise RuntimeError("proteomics selection config requires approved mapping")
    accession_config = approved.get(accession)
    if not isinstance(accession_config, dict):
        raise RuntimeError(f"proteomics accession is not approved: {accession}")
    names = accession_config.get("results")
    if not isinstance(names, list) or not names or any(
        not isinstance(name, str) or not name for name in names
    ):
        raise RuntimeError(f"proteomics results are not approved: {accession}")

    downloads_path = registry_dir / "downloads.tsv"
    with downloads_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sources: list[ProteomicsSource] = []
    for name in names:
        matches = [
            row
            for row in rows
            if row["dataset_accession"] == accession
            and row["file_name"] == name
        ]
        if len(matches) != 1:
            raise RuntimeError(f"proteomics source does not resolve once: {name}")
        row = matches[0]
        sources.append(
            ProteomicsSource(
                study_accession=accession,
                source_file=downloads_path.parent / row["path"],
                registry_file=downloads_path,
                source_sha256=row["sha256"],
            )
        )
    return tuple(sources)


def _resolve_references(
    accession: str,
    registry_dir: Path,
) -> tuple[Path, tuple[Path, ...]]:
    registry_path = registry_dir / "reference_sequences.tsv"
    with registry_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "study_accession",
            "protein_accession",
            "repository",
            "source_url",
            "retrieved_at",
            "sequence_version",
            "path",
            "size_bytes",
            "sha256",
            "scientific_use",
        }
        if set(reader.fieldnames or ()) != required:
            raise RuntimeError("reference sequence registry has invalid columns")
        matches = [
            row
            for row in reader
            if row["study_accession"] == accession
            and row["scientific_use"] == "coordinate_validation_only"
        ]
    if not matches:
        raise RuntimeError(f"reference sequences are not registered: {accession}")
    matches.sort(key=lambda row: row["protein_accession"])
    return registry_path, tuple(
        registry_path.parent / row["path"] for row in matches
    )


def parse_proteomics_accession(
    accession: str,
    output_root: Path = Path("data/interim"),
    registry_dir: Path = Path("data/registry"),
    selection_path: Path = Path("configs/download_selection.yaml"),
) -> ParseSummary:
    """Resolve and publish one accession using only approved registries."""
    normalized = accession.strip().upper()
    sources = _resolve_sources(normalized, registry_dir, selection_path)
    reference_registry, reference_paths = _resolve_references(
        normalized,
        registry_dir,
    )
    return publish_parsed_proteomics(
        sources,
        reference_paths,
        reference_registry,
        output_root / normalized / "proteomics_parser_v1",
    )
