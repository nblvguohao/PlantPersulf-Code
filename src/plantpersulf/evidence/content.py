"""Deterministic content audit of registered real evidence files."""

from __future__ import annotations

import json
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.download.registered import DOWNLOAD_FIELDS
from plantpersulf.evidence.csv_content import (
    CsvAudit,
    parse_peptide_csv,
    parse_protein_csv,
)
from plantpersulf.evidence.methods import audit_method_sources
from plantpersulf.evidence.mzidentml import MzidAudit, parse_mzidentml
from plantpersulf.evidence.preflight import (
    FILE_DECISION_FIELDS,
    GAP_FIELDS,
    OUTPUT_NAMES,
    QUEUE_FIELDS,
    STUDY_FIELDS,
    _read_tsv,
    _tree_hashes,
    _write_tsv,
    audit_metadata_preflight,
    build_metadata_preflight,
)
from plantpersulf.evidence.sdrf import SdrfAudit, parse_sdrf
from plantpersulf.provenance.hashing import hash_file

CONTENT_AUDIT_VERSION = 1
PARSER_VERSION = "plantpersulf/0.0.0"
CONTENT_EVIDENCE_FIELDS = (
    "dataset_accession",
    "record_id",
    "record_type",
    "source_file",
    "source_sha256",
    "source_locator",
    "evidence_class",
    "method",
    "method_source_sha256",
    "protein_accession",
    "peptide_sequence",
    "modification_raw",
    "raw_values_json",
    "sample_mapping_status",
    "quantification_status",
    "site_localization_status",
)
CANDIDATE_FIELDS = (
    "dataset_accession",
    "source_file",
    "source_sha256",
    "spectrum_result_id",
    "spectrum_id",
    "spectrum_identification_item_id",
    "rank",
    "pass_threshold",
    "peptide_id",
    "peptide_sequence",
    "peptide_evidence_id",
    "db_sequence_id",
    "protein_accession",
    "modification_index",
    "location",
    "residue",
    "protein_position",
    "monoisotopic_mass_delta",
    "cv_accession",
    "cv_name",
    "cv_value",
    "evidence_class",
    "method_mapping_status",
    "conflict_status",
    "source_locator",
)
SAMPLE_MAPPING_FIELDS = (
    "dataset_accession",
    "source_file",
    "source_sha256",
    "row_number",
    "source_name",
    "assay_name",
    "data_file",
    "mapping_status",
    "matched_registry_file",
    "source_locator",
)
MAPPING_CONFLICT_FIELDS = (
    "dataset_accession",
    "source_file",
    "source_sha256",
    "record_id",
    "conflict_type",
    "detail",
    "source_locator",
)
CONTENT_SCHEMA_FIELDS = (
    "dataset_accession",
    "source_file",
    "source_sha256",
    "format",
    "encoding",
    "delimiter_or_namespace",
    "header_json",
    "record_count",
    "structural_counts_json",
    "parser_version",
    "parse_status",
    "unsupported_reason",
)
CONTENT_OUTPUT_NAMES = OUTPUT_NAMES + (
    "candidate_modifications.tsv",
    "sample_file_mapping.tsv",
    "mapping_conflicts.tsv",
    "content_schema.tsv",
)

Row = dict[str, str]


@dataclass(frozen=True)
class ContentAuditSummary:
    study_count: int
    evidence_record_count: int
    candidate_modification_count: int
    sample_file_mapping_count: int
    mapping_conflict_count: int
    content_schema_count: int
    site_ms_count: int
    metadata_gap_count: int


def _load_empty_mappings(path: Path) -> None:
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid evidence method mappings: {path}") from exc
    if loaded != {"version": 1, "mappings": []}:
        raise RuntimeError("content audit v1 requires an empty reviewed mapping list")


def _source_rows(registry_dir: Path) -> dict[str, Row]:
    downloads_path = registry_dir / "downloads.tsv"
    rows = _read_tsv(downloads_path, DOWNLOAD_FIELDS)
    selected = {
        row["file_name"]: row
        for row in rows
        if row["dataset_accession"] == "PXD035795"
        and row["file_name"]
        in {"SDRF.txt", "peptide.csv", "proteins.csv", "peptides_1_1_0.mzid.gz"}
    }
    if set(selected) != {
        "SDRF.txt",
        "peptide.csv",
        "proteins.csv",
        "peptides_1_1_0.mzid.gz",
    }:
        raise RuntimeError("PXD035795 content sources are not fully registered")
    return selected


def _source_path(registry_dir: Path, row: Row) -> Path:
    path = (registry_dir / row["path"]).resolve()
    if not path.is_file() or hash_file(path, "sha256") != row["sha256"]:
        raise RuntimeError(f"content source failed SHA256 audit: {path}")
    return path


def _compact_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=isinstance(value, dict),
    )


def _csv_evidence(
    audit: CsvAudit,
    source_file: str,
    method_sha256: str,
    evidence_class: str,
    protein_index: int,
    peptide_index: int | None,
    modification_index: int,
) -> list[Row]:
    rows: list[Row] = []
    for record in audit.rows:
        values = dict(zip(audit.header, record.raw_fields, strict=True))
        rows.append(
            {
                "dataset_accession": "PXD035795",
                "record_id": (
                    f"PXD035795:{source_file}:{record.row_number}"
                ),
                "record_type": audit.record_type,
                "source_file": source_file,
                "source_sha256": record.source_sha256,
                "source_locator": record.source_locator,
                "evidence_class": evidence_class,
                "method": "dimedone_switch_lc_ms_ms",
                "method_source_sha256": method_sha256,
                "protein_accession": record.raw_fields[protein_index],
                "peptide_sequence": (
                    record.raw_fields[peptide_index]
                    if peptide_index is not None
                    else ""
                ),
                "modification_raw": record.raw_fields[modification_index],
                "raw_values_json": _compact_json(values),
                "sample_mapping_status": "partial_exact_raw_mapping",
                "quantification_status": (
                    "unresolved_no_registered_column_normalization"
                ),
                "site_localization_status": record.site_localization_status,
            }
        )
    return rows


def _candidate_rows(audit: MzidAudit) -> list[Row]:
    rows: list[Row] = []
    for candidate in audit.candidates:
        row = {
            field: str(getattr(candidate, field))
            for field in CANDIDATE_FIELDS
            if field not in {"dataset_accession", "source_file", "source_sha256"}
        }
        row.update(
            {
                "dataset_accession": "PXD035795",
                "source_file": "peptides_1_1_0.mzid.gz",
                "source_sha256": candidate.source_sha256,
            }
        )
        rows.append(row)
    return rows


def _sample_mapping_rows(audit: SdrfAudit) -> list[Row]:
    return [
        {
            "dataset_accession": "PXD035795",
            "source_file": "SDRF.txt",
            "source_sha256": row.source_sha256,
            "row_number": str(row.row_number),
            "source_name": row.source_name,
            "assay_name": row.assay_name,
            "data_file": row.data_file,
            "mapping_status": row.mapping_status,
            "matched_registry_file": row.matched_registry_file,
            "source_locator": row.source_locator,
        }
        for row in audit.file_mappings
    ]


def _conflict_rows(sdrf: SdrfAudit, mzid: MzidAudit) -> list[Row]:
    rows = [
        {
            "dataset_accession": "PXD035795",
            "source_file": "SDRF.txt",
            "source_sha256": conflict.source_sha256,
            "record_id": f"SDRF.txt:row={conflict.row_number}",
            "conflict_type": conflict.conflict_type,
            "detail": conflict.detail,
            "source_locator": conflict.source_locator,
        }
        for conflict in sdrf.conflicts
    ]
    rows.extend(
        {
            "dataset_accession": "PXD035795",
            "source_file": "peptides_1_1_0.mzid.gz",
            "source_sha256": conflict.source_sha256,
            "record_id": conflict.spectrum_identification_item_id,
            "conflict_type": conflict.conflict_type,
            "detail": conflict.detail,
            "source_locator": conflict.source_locator,
        }
        for conflict in mzid.conflicts
    )
    rows.sort(key=lambda row: (row["source_file"], row["source_locator"]))
    return rows


def _schema_rows(
    sdrf: SdrfAudit,
    peptide: CsvAudit,
    protein: CsvAudit,
    mzid: MzidAudit,
) -> list[Row]:
    return [
        {
            "dataset_accession": "PXD035795",
            "source_file": "SDRF.txt",
            "source_sha256": sdrf.source_sha256,
            "format": "SDRF",
            "encoding": "utf-8-sig",
            "delimiter_or_namespace": "tab",
            "header_json": _compact_json(list(sdrf.header)),
            "record_count": str(len(sdrf.rows)),
            "structural_counts_json": _compact_json(
                {
                    "duplicate_header_positions": {
                        name: list(indices)
                        for name, indices in sdrf.duplicate_header_positions.items()
                    }
                }
            ),
            "parser_version": PARSER_VERSION,
            "parse_status": "parsed",
            "unsupported_reason": "",
        },
        {
            "dataset_accession": "PXD035795",
            "source_file": "peptide.csv",
            "source_sha256": peptide.source_sha256,
            "format": "CSV",
            "encoding": peptide.encoding,
            "delimiter_or_namespace": peptide.delimiter,
            "header_json": _compact_json(list(peptide.header)),
            "record_count": str(len(peptide.rows)),
            "structural_counts_json": "{}",
            "parser_version": PARSER_VERSION,
            "parse_status": "parsed",
            "unsupported_reason": "",
        },
        {
            "dataset_accession": "PXD035795",
            "source_file": "proteins.csv",
            "source_sha256": protein.source_sha256,
            "format": "CSV",
            "encoding": protein.encoding,
            "delimiter_or_namespace": protein.delimiter,
            "header_json": _compact_json(list(protein.header)),
            "record_count": str(len(protein.rows)),
            "structural_counts_json": "{}",
            "parser_version": PARSER_VERSION,
            "parse_status": "parsed",
            "unsupported_reason": "",
        },
        {
            "dataset_accession": "PXD035795",
            "source_file": "peptides_1_1_0.mzid.gz",
            "source_sha256": mzid.source_sha256,
            "format": "mzIdentML",
            "encoding": "utf-8-gzip",
            "delimiter_or_namespace": mzid.namespace,
            "header_json": "[]",
            "record_count": str(
                mzid.structural_counts["SpectrumIdentificationItem"]
            ),
            "structural_counts_json": _compact_json(mzid.structural_counts),
            "parser_version": PARSER_VERSION,
            "parse_status": "parsed",
            "unsupported_reason": "",
        },
    ]


def _updated_studies(path: Path) -> list[Row]:
    rows = _read_tsv(path, STUDY_FIELDS)
    statuses = {
        "PXD006140": "retained_task4_coordinate_only",
        "PXD024061": "unsupported_stage_a_checksum_only",
        "PXD035795": "complete_cycle_2",
        "PXD039999": "unsupported_stage_a_checksum_only",
    }
    for row in rows:
        accession = row["dataset_accession"]
        row["content_audit_status"] = statuses[accession]
        if accession == "PXD035795":
            row["evidence_class"] = "protein_level_only"
    return rows


def _updated_decisions(path: Path) -> list[Row]:
    rows = _read_tsv(path, FILE_DECISION_FIELDS)
    classes = {
        "peptide.csv": "identification_only",
        "proteins.csv": "protein_level_only",
    }
    for row in rows:
        if row["dataset_accession"] == "PXD035795":
            row["evidence_class"] = classes.get(row["file_name"], "unresolved")
    return rows


def _updated_gaps(path: Path) -> list[Row]:
    rows = [
        row
        for row in _read_tsv(path, GAP_FIELDS)
        if row["gap_type"] != "content_audit_pending"
        and not (
            row["dataset_accession"] == "PXD035795"
            and row["gap_type"] == "missing_sample_mapping"
        )
    ]
    rows.extend(
        [
            {
                "dataset_accession": "PXD035795",
                "file_name": "SDRF.txt",
                "gap_type": "partial_sample_file_mapping",
                "detail": (
                    "six RAW names map exactly and six encoded MGF names conflict"
                ),
            },
            {
                "dataset_accession": "PXD035795",
                "file_name": "peptide.csv;proteins.csv",
                "gap_type": "quantification_mapping_unresolved",
                "detail": (
                    "normalization and exact CSV-to-assay mapping are not registered"
                ),
            },
            {
                "dataset_accession": "PXD035795",
                "file_name": "peptides_1_1_0.mzid.gz",
                "gap_type": "site_method_mapping_absent",
                "detail": "reviewed method mapping list is empty in version 1",
            },
        ]
    )
    study_order = {
        accession: index
        for index, accession in enumerate(
            ("PXD006140", "PXD024061", "PXD035795", "PXD039999")
        )
    }
    rows.sort(
        key=lambda row: (
            study_order[row["dataset_accession"]],
            row["gap_type"],
            row["file_name"],
        )
    )
    return rows


def _summary(
    studies: list[Row],
    evidence: list[Row],
    candidates: list[Row],
    mappings: list[Row],
    conflicts: list[Row],
    schemas: list[Row],
    gaps: list[Row],
) -> ContentAuditSummary:
    return ContentAuditSummary(
        study_count=len(studies),
        evidence_record_count=len(evidence),
        candidate_modification_count=len(candidates),
        sample_file_mapping_count=len(mappings),
        mapping_conflict_count=len(conflicts),
        content_schema_count=len(schemas),
        site_ms_count=sum(
            row["evidence_class"] == "site_ms" for row in evidence + candidates
        ),
        metadata_gap_count=len(gaps),
    )


def _publish_staged(staged: Path, output: Path) -> None:
    if not output.exists():
        staged.replace(output)
        return
    try:
        audit_content_output(output)
    except RuntimeError as content_error:
        audit_metadata_preflight(output)
        backup = output.with_name(f".{output.name}.cycle1-backup")
        if backup.exists():
            raise RuntimeError(
                f"content audit backup already exists: {backup}"
            ) from content_error
        output.replace(backup)
        try:
            staged.replace(output)
        except Exception:
            backup.replace(output)
            raise
        shutil.rmtree(backup)
        return
    if _tree_hashes(output) != _tree_hashes(staged):
        raise RuntimeError(f"existing evidence content audit differs: {output}")


def build_content_audit(
    policy_path: Path = Path("configs/evidence_preflight_v1.yaml"),
    selection_path: Path = Path("configs/download_selection.yaml"),
    method_config_path: Path = Path("configs/evidence_method_sources_v1.yaml"),
    method_mapping_path: Path = Path("configs/evidence_method_mappings_v1.yaml"),
    registry_dir: Path = Path("data/registry"),
    output_directory: Path = Path("data/interim/evidence_preflight_v1"),
) -> ContentAuditSummary:
    """Build a label-free content inventory from registered real files."""
    _load_empty_mappings(method_mapping_path)
    method_sources = audit_method_sources(
        method_config_path,
        registry_dir / "evidence_methods.tsv",
    )
    method_matches = [
        source
        for source in method_sources
        if source.study_accession == "PXD035795"
    ]
    if len(method_matches) != 1:
        raise RuntimeError("PXD035795 method source is not registered exactly once")
    method = method_matches[0]

    source_rows = _source_rows(registry_dir)
    sdrf_row = source_rows["SDRF.txt"]
    peptide_row = source_rows["peptide.csv"]
    protein_row = source_rows["proteins.csv"]
    mzid_row = source_rows["peptides_1_1_0.mzid.gz"]
    sdrf = parse_sdrf(
        _source_path(registry_dir, sdrf_row),
        sdrf_row["sha256"],
        registry_dir / "files.tsv",
    )
    peptide = parse_peptide_csv(
        _source_path(registry_dir, peptide_row), peptide_row["sha256"]
    )
    protein = parse_protein_csv(
        _source_path(registry_dir, protein_row), protein_row["sha256"]
    )
    mzid = parse_mzidentml(
        _source_path(registry_dir, mzid_row), mzid_row["sha256"]
    )

    evidence = _csv_evidence(
        peptide,
        "peptide.csv",
        method.sha256,
        "identification_only",
        protein_index=17,
        peptide_index=0,
        modification_index=18,
    )
    evidence.extend(
        _csv_evidence(
            protein,
            "proteins.csv",
            method.sha256,
            "protein_level_only",
            protein_index=2,
            peptide_index=None,
            modification_index=7,
        )
    )
    candidates = _candidate_rows(mzid)
    mappings = _sample_mapping_rows(sdrf)
    conflicts = _conflict_rows(sdrf, mzid)
    schemas = _schema_rows(sdrf, peptide, protein, mzid)

    output_directory.parent.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(
        tempfile.mkdtemp(
            prefix=f".{output_directory.name}.content.",
            dir=output_directory.parent,
        )
    )
    staged = temporary_root / output_directory.name
    try:
        build_metadata_preflight(
            policy_path=policy_path,
            selection_path=selection_path,
            registry_dir=registry_dir,
            output_directory=staged,
        )
        studies = _updated_studies(staged / "study_inventory.tsv")
        decisions = _updated_decisions(staged / "file_decisions.tsv")
        gaps = _updated_gaps(staged / "metadata_gaps.tsv")
        queue = _read_tsv(staged / "large_file_queue.tsv", QUEUE_FIELDS)
        summary = _summary(
            studies, evidence, candidates, mappings, conflicts, schemas, gaps
        )
        _write_tsv(staged / "study_inventory.tsv", STUDY_FIELDS, studies)
        _write_tsv(
            staged / "file_decisions.tsv", FILE_DECISION_FIELDS, decisions
        )
        _write_tsv(
            staged / "evidence_records.tsv", CONTENT_EVIDENCE_FIELDS, evidence
        )
        _write_tsv(staged / "metadata_gaps.tsv", GAP_FIELDS, gaps)
        _write_tsv(staged / "large_file_queue.tsv", QUEUE_FIELDS, queue)
        _write_tsv(
            staged / "candidate_modifications.tsv", CANDIDATE_FIELDS, candidates
        )
        _write_tsv(
            staged / "sample_file_mapping.tsv", SAMPLE_MAPPING_FIELDS, mappings
        )
        _write_tsv(
            staged / "mapping_conflicts.tsv", MAPPING_CONFLICT_FIELDS, conflicts
        )
        _write_tsv(staged / "content_schema.tsv", CONTENT_SCHEMA_FIELDS, schemas)

        manifest_path = staged / "manifest.json"
        parsed: object = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not isinstance(parsed, dict):
            raise RuntimeError("cycle-1 manifest is not a mapping")
        manifest = cast(dict[str, Any], parsed)
        additional_inputs = (
            method_config_path,
            method_mapping_path,
            registry_dir / "evidence_methods.tsv",
        )
        manifest["schema_version"] = 2
        manifest["content_audit_version"] = CONTENT_AUDIT_VERSION
        manifest["inputs"] = cast(list[object], manifest["inputs"]) + [
            {"path": path.as_posix(), "sha256": hash_file(path, "sha256")}
            for path in additional_inputs
        ]
        manifest["registered_method_sources"] = [
            {
                "path": method.local_path.as_posix(),
                "sha256": method.sha256,
            }
        ]
        manifest["counts"] = asdict(summary)
        manifest["content_audit_complete"] = True
        manifest["labels_created"] = False
        manifest["nondetection_labeled_negative"] = False
        manifest["biological_values_modified"] = False
        manifest["outputs"] = [
            {"file": name, "sha256": hash_file(staged / name, "sha256")}
            for name in CONTENT_OUTPUT_NAMES
            if name != "manifest.json"
        ]
        manifest_path.write_bytes(
            (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
        )
        _publish_staged(staged, output_directory)
        return summary
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def _manifest_mapping(value: object, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise RuntimeError(f"content manifest requires mapping: {context}")
    return cast(dict[str, Any], value)


def _manifest_items(value: object, context: str) -> list[Any]:
    if not isinstance(value, list):
        raise RuntimeError(f"content manifest requires list: {context}")
    return value


def audit_content_output(
    output_directory: Path = Path("data/interim/evidence_preflight_v1"),
) -> ContentAuditSummary:
    """Rehash all Cycle 2 inputs, real sources, method sources, and outputs."""
    observed = {path.name for path in output_directory.iterdir() if path.is_file()}
    if observed != set(CONTENT_OUTPUT_NAMES):
        raise RuntimeError("evidence content output file set mismatch")
    manifest_path = output_directory / "manifest.json"
    try:
        parsed: object = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"invalid evidence content manifest: {manifest_path}"
        ) from exc
    manifest = _manifest_mapping(parsed, "root")
    if (
        manifest.get("schema_version") != 2
        or manifest.get("preflight_version") != "v1"
        or manifest.get("content_audit_version") != CONTENT_AUDIT_VERSION
    ):
        raise RuntimeError("evidence content manifest identity mismatch")
    if (
        manifest.get("content_audit_complete") is not True
        or manifest.get("labels_created") is not False
        or manifest.get("nondetection_labeled_negative") is not False
        or manifest.get("biological_values_modified") is not False
    ):
        raise RuntimeError("evidence content integrity flags are invalid")
    for collection in ("inputs", "registered_sources", "registered_method_sources"):
        for raw_entry in _manifest_items(manifest.get(collection), collection):
            entry = _manifest_mapping(raw_entry, f"{collection}[]")
            path = Path(str(entry.get("path", "")))
            if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
                raise RuntimeError(f"evidence content source SHA256 mismatch: {path}")
    declared: set[str] = set()
    for raw_entry in _manifest_items(manifest.get("outputs"), "outputs"):
        entry = _manifest_mapping(raw_entry, "outputs[]")
        name = str(entry.get("file", ""))
        declared.add(name)
        path = output_directory / name
        if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
            raise RuntimeError(f"evidence content output SHA256 mismatch: {path}")
    if declared != set(CONTENT_OUTPUT_NAMES) - {"manifest.json"}:
        raise RuntimeError("evidence content manifest output list mismatch")

    studies = _read_tsv(output_directory / "study_inventory.tsv", STUDY_FIELDS)
    _read_tsv(output_directory / "file_decisions.tsv", FILE_DECISION_FIELDS)
    evidence = _read_tsv(
        output_directory / "evidence_records.tsv", CONTENT_EVIDENCE_FIELDS
    )
    gaps = _read_tsv(output_directory / "metadata_gaps.tsv", GAP_FIELDS)
    _read_tsv(output_directory / "large_file_queue.tsv", QUEUE_FIELDS)
    candidates = _read_tsv(
        output_directory / "candidate_modifications.tsv", CANDIDATE_FIELDS
    )
    mappings = _read_tsv(
        output_directory / "sample_file_mapping.tsv", SAMPLE_MAPPING_FIELDS
    )
    conflicts = _read_tsv(
        output_directory / "mapping_conflicts.tsv", MAPPING_CONFLICT_FIELDS
    )
    schemas = _read_tsv(
        output_directory / "content_schema.tsv", CONTENT_SCHEMA_FIELDS
    )
    summary = _summary(
        studies, evidence, candidates, mappings, conflicts, schemas, gaps
    )
    if summary.site_ms_count != 0:
        raise RuntimeError("content audit v1 cannot contain site_ms evidence")
    if any(row["evidence_class"] != "unresolved" for row in candidates):
        raise RuntimeError("content audit candidates must remain unresolved")
    if len(schemas) != 4 or any(row["parse_status"] != "parsed" for row in schemas):
        raise RuntimeError("content audit schema coverage is incomplete")
    if _manifest_mapping(manifest.get("counts"), "counts") != asdict(summary):
        raise RuntimeError("evidence content counts mismatch")
    return summary
