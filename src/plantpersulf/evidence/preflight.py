"""Deterministic metadata inventory before biological evidence classification."""

from __future__ import annotations

import csv
import json
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.download.registered import (
    DOWNLOAD_FIELDS,
    audit_downloaded_files,
)
from plantpersulf.evidence.policy import load_preflight_policy
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file

PREFLIGHT_VERSION = "v1"
OUTPUT_NAMES = (
    "study_inventory.tsv",
    "file_decisions.tsv",
    "evidence_records.tsv",
    "metadata_gaps.tsv",
    "large_file_queue.tsv",
    "manifest.json",
)
STUDY_FIELDS = (
    "dataset_accession",
    "scientific_role",
    "publication_date",
    "official_file_count",
    "selected_source_count",
    "downloaded_source_count",
    "sample_mapping_count",
    "evidence_class",
    "content_audit_status",
)
FILE_DECISION_FIELDS = (
    "dataset_accession",
    "file_name",
    "record_type",
    "file_category",
    "source_url",
    "registry_size_bytes",
    "remote_checksum_algorithm",
    "remote_checksum",
    "decision",
    "reason",
    "evidence_class",
    "local_path",
    "local_size_bytes",
    "local_sha256",
)
EVIDENCE_FIELDS = (
    "dataset_accession",
    "source_file",
    "source_sha256",
    "source_locator",
    "evidence_class",
    "method",
    "sample_mapping_status",
    "quantification_status",
    "site_localization_status",
)
GAP_FIELDS = (
    "dataset_accession",
    "file_name",
    "gap_type",
    "detail",
)
QUEUE_FIELDS = (
    "dataset_accession",
    "file_name",
    "file_category",
    "registry_size_bytes",
    "reason",
    "required_authorization",
)


@dataclass(frozen=True)
class PreflightSummary:
    study_count: int
    file_count: int
    registered_metadata_count: int
    downloaded_count: int
    downloaded_existing_scope_count: int
    deferred_large_count: int
    deferred_count: int
    evidence_record_count: int
    metadata_gap_count: int


Row = dict[str, str]


def _read_tsv(path: Path, expected_fields: tuple[str, ...] | None = None) -> list[Row]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if (
            expected_fields is not None
            and tuple(reader.fieldnames or ()) != expected_fields
        ):
            raise RuntimeError(f"preflight input has invalid columns: {path}")
        rows = [dict(row) for row in reader]
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise RuntimeError(f"preflight input has malformed row: {path}")
    return cast(list[Row], rows)


def _write_tsv(
    path: Path,
    fields: tuple[str, ...],
    rows: list[Row],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _selection_mapping(path: Path) -> dict[str, dict[str, tuple[str, ...]]]:
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid preflight selection: {path}") from exc
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise RuntimeError("preflight selection requires version 1")
    approved = loaded.get("approved")
    if not isinstance(approved, dict):
        raise RuntimeError("preflight selection requires approved mapping")
    result: dict[str, dict[str, tuple[str, ...]]] = {}
    for raw_accession, raw_classes in approved.items():
        if not isinstance(raw_accession, str) or not isinstance(raw_classes, dict):
            raise RuntimeError("preflight selection has invalid study mapping")
        result[raw_accession] = {}
        for raw_class, raw_names in raw_classes.items():
            if (
                not isinstance(raw_class, str)
                or not isinstance(raw_names, list)
                or any(not isinstance(name, str) or not name for name in raw_names)
            ):
                raise RuntimeError("preflight selection has invalid file class")
            result[raw_accession][raw_class] = tuple(raw_names)
    return result


def _selection_names(
    mapping: dict[str, dict[str, tuple[str, ...]]],
    accession: str,
) -> dict[str, str]:
    classes = mapping.get(accession)
    if classes is None:
        raise RuntimeError(f"preflight study lacks selection: {accession}")
    names: dict[str, str] = {}
    for file_class, class_names in classes.items():
        for name in class_names:
            if name in names:
                raise RuntimeError(f"preflight selection duplicates file: {name}")
            names[name] = file_class
    return names


def _unique_by(
    rows: list[Row],
    fields: tuple[str, ...],
    context: str,
) -> dict[tuple[str, ...], Row]:
    result: dict[tuple[str, ...], Row] = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        if key in result:
            raise RuntimeError(f"duplicate {context}: {key}")
        result[key] = row
    return result


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hash_file(path, "sha256")
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _summary(decisions: list[Row], gaps: list[Row]) -> PreflightSummary:
    counts = {
        decision: sum(row["decision"] == decision for row in decisions)
        for decision in {
            "registered_metadata",
            "downloaded",
            "downloaded_existing_scope",
            "deferred_large_not_stage_a",
            "deferred_not_stage_a",
        }
    }
    return PreflightSummary(
        study_count=len({row["dataset_accession"] for row in decisions}),
        file_count=len(decisions),
        registered_metadata_count=counts["registered_metadata"],
        downloaded_count=counts["downloaded"],
        downloaded_existing_scope_count=counts["downloaded_existing_scope"],
        deferred_large_count=counts["deferred_large_not_stage_a"],
        deferred_count=counts["deferred_not_stage_a"],
        evidence_record_count=0,
        metadata_gap_count=len(gaps),
    )


def _decision(
    file_row: Row,
    selected_class: str,
    download_row: Row | None,
    threshold: int,
    downloads_path: Path,
) -> Row:
    accession = file_row["dataset_accession"]
    if file_row["record_type"] == "metadata_cache":
        decision = "registered_metadata"
        reason = "official metadata cache is registered"
    elif selected_class and download_row is not None:
        decision = (
            "downloaded_existing_scope" if accession == "PXD006140" else "downloaded"
        )
        reason = "exact approved source is downloaded and checksummed"
    elif int(file_row["size_bytes"]) >= threshold:
        decision = "deferred_large_not_stage_a"
        reason = "file exceeds the reviewed Stage A threshold"
    else:
        decision = "deferred_not_stage_a"
        reason = "file is not in the reviewed Stage A allow-list"

    local_path = ""
    local_size = ""
    local_sha256 = ""
    if download_row is not None:
        path = downloads_path.parent / download_row["path"]
        assert_registered_input(path, downloads_path)
        actual_sha256 = hash_file(path, "sha256")
        if actual_sha256 != download_row["sha256"] or path.stat().st_size != int(
            download_row["size_bytes"]
        ):
            raise RuntimeError(f"preflight downloaded source mismatch: {path}")
        local_path = path.as_posix()
        local_size = download_row["size_bytes"]
        local_sha256 = actual_sha256
    return {
        "dataset_accession": accession,
        "file_name": file_row["file_name"],
        "record_type": file_row["record_type"],
        "file_category": file_row["file_category"],
        "source_url": file_row["source_url"],
        "registry_size_bytes": file_row["size_bytes"],
        "remote_checksum_algorithm": file_row["remote_checksum_algorithm"],
        "remote_checksum": file_row["remote_checksum"],
        "decision": decision,
        "reason": reason,
        "evidence_class": "unresolved",
        "local_path": local_path,
        "local_size_bytes": local_size,
        "local_sha256": local_sha256,
    }


def build_metadata_preflight(
    policy_path: Path = Path("configs/evidence_preflight_v1.yaml"),
    selection_path: Path = Path("configs/download_selection.yaml"),
    registry_dir: Path = Path("data/registry"),
    output_directory: Path = Path("data/interim/evidence_preflight_v1"),
) -> PreflightSummary:
    """Publish a metadata-only inventory from registered real inputs."""
    policy = load_preflight_policy(policy_path)
    selection = _selection_mapping(selection_path)
    files_path = registry_dir / "files.tsv"
    datasets_path = registry_dir / "datasets.tsv"
    downloads_path = registry_dir / "downloads.tsv"
    samples_path = registry_dir / "samples.tsv"
    files = _read_tsv(files_path)
    datasets = _read_tsv(datasets_path)
    downloads = _read_tsv(downloads_path, DOWNLOAD_FIELDS)
    samples = _read_tsv(samples_path)
    download_map = _unique_by(
        downloads,
        ("dataset_accession", "file_name"),
        "download row",
    )
    dataset_map = _unique_by(datasets, ("accession",), "dataset row")
    study_order = {accession: index for index, accession in enumerate(policy.studies)}

    decisions: list[Row] = []
    studies: list[Row] = []
    gaps: list[Row] = []
    selected_classes: dict[tuple[str, str], str] = {}
    for accession in policy.studies:
        names = _selection_names(selection, accession)
        for name, file_class in names.items():
            selected_classes[(accession, name)] = file_class
        audit_downloaded_files(
            accession=accession,
            selection_path=selection_path,
            files_registry_path=files_path,
            datasets_registry_path=datasets_path,
            downloads_registry_path=downloads_path,
        )

    target_files = [row for row in files if row["dataset_accession"] in study_order]
    target_files.sort(
        key=lambda row: (study_order[row["dataset_accession"]], row["file_name"])
    )
    for file_row in target_files:
        key = (file_row["dataset_accession"], file_row["file_name"])
        decisions.append(
            _decision(
                file_row,
                selected_classes.get(key, ""),
                download_map.get(key),
                policy.large_file_threshold_bytes,
                downloads_path,
            )
        )

    for accession in policy.studies:
        dataset = dataset_map.get((accession,))
        if dataset is None:
            raise RuntimeError(f"preflight dataset is not registered: {accession}")
        accession_files = [
            row for row in decisions if row["dataset_accession"] == accession
        ]
        accession_samples = [
            row for row in samples if row["dataset_accession"] == accession
        ]
        selected_count = sum(
            bool(selected_classes.get((accession, row["file_name"])))
            and row["record_type"] == "source_file"
            for row in accession_files
        )
        downloaded_count = sum(
            row["decision"] in {"downloaded", "downloaded_existing_scope"}
            for row in accession_files
        )
        studies.append(
            {
                "dataset_accession": accession,
                "scientific_role": dataset["scientific_role"],
                "publication_date": dataset["publication_date"],
                "official_file_count": str(len(accession_files)),
                "selected_source_count": str(selected_count),
                "downloaded_source_count": str(downloaded_count),
                "sample_mapping_count": str(len(accession_samples)),
                "evidence_class": "unresolved",
                "content_audit_status": "pending_cycle_2",
            }
        )
        if not accession_samples:
            gaps.append(
                {
                    "dataset_accession": accession,
                    "file_name": "",
                    "gap_type": "missing_sample_mapping",
                    "detail": "no official sample mapping is registered",
                }
            )

    for (accession, name), file_class in sorted(selected_classes.items()):
        if accession != "PXD006140" and file_class in {"design", "results"}:
            gaps.append(
                {
                    "dataset_accession": accession,
                    "file_name": name,
                    "gap_type": "content_audit_pending",
                    "detail": "registered real file awaits cycle-2 schema audit",
                }
            )
    gaps.append(
        {
            "dataset_accession": "PXD006140",
            "file_name": "",
            "gap_type": "site_evidence_policy_pending",
            "detail": "Task 4 coordinates are not positive evidence labels",
        }
    )
    for accession in ("PXD024061", "PXD039999"):
        gaps.append(
            {
                "dataset_accession": accession,
                "file_name": "",
                "gap_type": "processed_result_not_selected_stage_a",
                "detail": "no small processed result is available in Stage A",
            }
        )
    gaps.sort(
        key=lambda row: (
            study_order[row["dataset_accession"]],
            row["gap_type"],
            row["file_name"],
        )
    )
    queue = [
        {
            "dataset_accession": row["dataset_accession"],
            "file_name": row["file_name"],
            "file_category": row["file_category"],
            "registry_size_bytes": row["registry_size_bytes"],
            "reason": row["reason"],
            "required_authorization": "separate reviewed large-file acquisition",
        }
        for row in decisions
        if row["decision"] == "deferred_large_not_stage_a"
    ]
    summary = _summary(decisions, gaps)

    output_directory.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = Path(
        tempfile.mkdtemp(
            prefix=f".{output_directory.name}.",
            dir=output_directory.parent,
        )
    )
    temporary_directory: Path | None = temporary_path
    try:
        _write_tsv(temporary_path / "study_inventory.tsv", STUDY_FIELDS, studies)
        _write_tsv(
            temporary_path / "file_decisions.tsv",
            FILE_DECISION_FIELDS,
            decisions,
        )
        _write_tsv(temporary_path / "evidence_records.tsv", EVIDENCE_FIELDS, [])
        _write_tsv(temporary_path / "metadata_gaps.tsv", GAP_FIELDS, gaps)
        _write_tsv(temporary_path / "large_file_queue.tsv", QUEUE_FIELDS, queue)
        inputs = [
            policy_path,
            selection_path,
            datasets_path,
            files_path,
            downloads_path,
            samples_path,
        ]
        outputs = [name for name in OUTPUT_NAMES if name != "manifest.json"]
        manifest: dict[str, object] = {
            "schema_version": 1,
            "preflight_version": PREFLIGHT_VERSION,
            "inputs": [
                {"path": path.as_posix(), "sha256": hash_file(path, "sha256")}
                for path in inputs
            ],
            "registered_sources": [
                {
                    "path": row["local_path"],
                    "sha256": row["local_sha256"],
                }
                for row in decisions
                if row["local_path"]
            ],
            "outputs": [
                {
                    "file": name,
                    "sha256": hash_file(temporary_path / name, "sha256"),
                }
                for name in outputs
            ],
            "counts": asdict(summary),
            "content_audit_complete": False,
            "labels_created": False,
            "nondetection_labeled_negative": False,
            "biological_values_modified": False,
        }
        serialized = json.dumps(manifest, indent=2, sort_keys=True)
        (temporary_path / "manifest.json").write_bytes(f"{serialized}\n".encode())
        if output_directory.exists():
            if _tree_hashes(output_directory) != _tree_hashes(temporary_path):
                raise RuntimeError(
                    f"existing evidence preflight differs: {output_directory}"
                )
            return summary
        temporary_path.replace(output_directory)
        temporary_directory = None
        return summary
    finally:
        if temporary_directory is not None:
            shutil.rmtree(temporary_directory)


def _mapping(value: object, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise RuntimeError(f"preflight manifest requires mapping: {context}")
    return cast(dict[str, Any], value)


def _items(value: object, context: str) -> list[Any]:
    if not isinstance(value, list):
        raise RuntimeError(f"preflight manifest requires list: {context}")
    return value


def audit_metadata_preflight(
    output_directory: Path = Path("data/interim/evidence_preflight_v1"),
) -> PreflightSummary:
    """Rehash deterministic preflight outputs and registered source bytes."""
    observed = {path.name for path in output_directory.iterdir() if path.is_file()}
    if observed != set(OUTPUT_NAMES):
        raise RuntimeError("evidence preflight output file set mismatch")
    manifest_path = output_directory / "manifest.json"
    try:
        parsed: object = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"invalid evidence preflight manifest: {manifest_path}"
        ) from exc
    manifest = _mapping(parsed, "root")
    if (
        manifest.get("schema_version") != 1
        or manifest.get("preflight_version") != PREFLIGHT_VERSION
    ):
        raise RuntimeError("evidence preflight manifest identity mismatch")
    if (
        manifest.get("content_audit_complete") is not False
        or manifest.get("labels_created") is not False
        or manifest.get("nondetection_labeled_negative") is not False
        or manifest.get("biological_values_modified") is not False
    ):
        raise RuntimeError("evidence preflight integrity flags are invalid")
    for raw_input in _items(manifest.get("inputs"), "inputs"):
        entry = _mapping(raw_input, "inputs[]")
        path = Path(str(entry.get("path", "")))
        if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
            raise RuntimeError(f"preflight input SHA256 mismatch: {path}")
    for raw_source in _items(manifest.get("registered_sources"), "registered_sources"):
        entry = _mapping(raw_source, "registered_sources[]")
        path = Path(str(entry.get("path", "")))
        if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
            raise RuntimeError(f"preflight source SHA256 mismatch: {path}")
    declared: set[str] = set()
    for raw_output in _items(manifest.get("outputs"), "outputs"):
        entry = _mapping(raw_output, "outputs[]")
        name = str(entry.get("file", ""))
        declared.add(name)
        path = output_directory / name
        if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
            raise RuntimeError(f"output SHA256 mismatch: {path}")
    if declared != set(OUTPUT_NAMES) - {"manifest.json"}:
        raise RuntimeError("evidence preflight manifest output list mismatch")
    studies = _read_tsv(output_directory / "study_inventory.tsv", STUDY_FIELDS)
    decisions = _read_tsv(output_directory / "file_decisions.tsv", FILE_DECISION_FIELDS)
    evidence = _read_tsv(output_directory / "evidence_records.tsv", EVIDENCE_FIELDS)
    gaps = _read_tsv(output_directory / "metadata_gaps.tsv", GAP_FIELDS)
    queue = _read_tsv(output_directory / "large_file_queue.tsv", QUEUE_FIELDS)
    if evidence:
        raise RuntimeError("cycle-1 preflight cannot contain evidence records")
    if len(queue) != sum(
        row["decision"] == "deferred_large_not_stage_a" for row in decisions
    ):
        raise RuntimeError("evidence preflight large-file queue mismatch")
    summary = _summary(decisions, gaps)
    if len(studies) != summary.study_count:
        raise RuntimeError("evidence preflight study count mismatch")
    if _mapping(manifest.get("counts"), "counts") != asdict(summary):
        raise RuntimeError("evidence preflight counts mismatch")
    return summary
