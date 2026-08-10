"""Deterministic Task 5 readiness audit without benchmark labels."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.evidence.content import (
    CANDIDATE_FIELDS,
    CONTENT_EVIDENCE_FIELDS,
    audit_content_output,
)
from plantpersulf.evidence.preflight import (
    STUDY_FIELDS,
    _read_tsv,
    _tree_hashes,
    _write_tsv,
)
from plantpersulf.proteomics.metadata import SITE_FIELDS
from plantpersulf.proteomics.site_normalizer import audit_site_output
from plantpersulf.provenance.hashing import hash_file

READINESS_VERSION = 1
OUTPUT_NAMES = (
    "study_readiness.tsv",
    "blockers.tsv",
    "manifest.json",
)
STUDY_READINESS_FIELDS = (
    "study_accession",
    "eligible_site_count",
    "coordinate_only_count",
    "non_site_evidence_count",
    "unresolved_candidate_count",
    "content_audit_status",
    "readiness_status",
)
BLOCKER_FIELDS = (
    "scope",
    "study_accession",
    "blocker_code",
    "observed_value",
    "required_value",
    "source_locator",
)
STUDIES = (
    "PXD006140",
    "PXD024061",
    "PXD035795",
    "PXD039999",
)
ELIGIBLE_LEVELS = (
    "site_ms",
    "site_mutagenesis",
    "site_biochemical",
)
COORDINATE_ONLY_LEVELS = ("psm_coordinate_only",)
PROHIBITED_LABELS = ("positive", "unlabeled", "negative")

Row = dict[str, str]


@dataclass(frozen=True)
class ReadinessPolicy:
    studies: tuple[str, ...]
    eligible_levels: tuple[str, ...]
    coordinate_only_levels: tuple[str, ...]
    required_study_count: int


@dataclass(frozen=True)
class BenchmarkReadinessSummary:
    decision: str
    study_count: int
    eligible_study_count: int
    eligible_site_count: int
    coordinate_only_count: int
    non_site_evidence_count: int
    unresolved_candidate_count: int
    required_study_count: int
    blocker_count: int


def _load_policy(path: Path) -> ReadinessPolicy:
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid benchmark readiness policy: {path}") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError("benchmark readiness policy requires a mapping")
    expected = {
        "version": 1,
        "studies": list(STUDIES),
        "eligible_site_evidence_levels": list(ELIGIBLE_LEVELS),
        "coordinate_only_evidence_levels": list(COORDINATE_ONLY_LEVELS),
        "minimum_distinct_site_evidence_studies": 2,
        "prohibited_record_labels": list(PROHIBITED_LABELS),
    }
    if loaded != expected:
        raise RuntimeError("benchmark readiness policy differs from frozen v1")
    return ReadinessPolicy(
        studies=STUDIES,
        eligible_levels=ELIGIBLE_LEVELS,
        coordinate_only_levels=COORDINATE_ONLY_LEVELS,
        required_study_count=2,
    )


def _relative_to_output(path: Path, output_directory: Path) -> str:
    return Path(os.path.relpath(path, output_directory)).as_posix()


def _input_entry(
    name: str,
    path: Path,
    output_directory: Path,
    resolution: str,
) -> dict[str, str]:
    stored_path = (
        path.as_posix()
        if resolution == "repository"
        else _relative_to_output(path, output_directory)
    )
    return {
        "name": name,
        "path": stored_path,
        "resolution": resolution,
        "sha256": hash_file(path, "sha256"),
    }


def _readiness_status(row: Row) -> str:
    if int(row["eligible_site_count"]):
        return "eligible_site_evidence"
    if int(row["coordinate_only_count"]):
        return "coordinate_only_not_site_evidence"
    if int(row["non_site_evidence_count"]) or int(row["unresolved_candidate_count"]):
        return "non_site_evidence_only"
    return row["content_audit_status"]


def _study_rows(
    policy: ReadinessPolicy,
    parser_rows: list[Row],
    content_studies: list[Row],
    evidence_rows: list[Row],
    candidate_rows: list[Row],
) -> list[Row]:
    content_by_study = {row["dataset_accession"]: row for row in content_studies}
    if tuple(content_by_study) != policy.studies:
        raise RuntimeError("content audit studies differ from readiness policy")

    rows: list[Row] = []
    for accession in policy.studies:
        parser_for_study = [
            row for row in parser_rows if row["study_accession"] == accession
        ]
        evidence_for_study = [
            row for row in evidence_rows if row["dataset_accession"] == accession
        ]
        candidates_for_study = [
            row for row in candidate_rows if row["dataset_accession"] == accession
        ]
        eligible_count = sum(
            row["evidence_level"] in policy.eligible_levels for row in parser_for_study
        ) + sum(
            row["evidence_class"] in policy.eligible_levels
            for row in evidence_for_study + candidates_for_study
        )
        row = {
            "study_accession": accession,
            "eligible_site_count": str(eligible_count),
            "coordinate_only_count": str(
                sum(
                    item["evidence_level"] in policy.coordinate_only_levels
                    for item in parser_for_study
                )
            ),
            "non_site_evidence_count": str(
                sum(
                    item["evidence_class"] not in policy.eligible_levels
                    for item in evidence_for_study
                )
            ),
            "unresolved_candidate_count": str(
                sum(
                    item["evidence_class"] not in policy.eligible_levels
                    for item in candidates_for_study
                )
            ),
            "content_audit_status": content_by_study[accession]["content_audit_status"],
            "readiness_status": "",
        }
        row["readiness_status"] = _readiness_status(row)
        rows.append(row)
    return rows


def _blockers(
    studies: list[Row],
    policy: ReadinessPolicy,
) -> list[Row]:
    rows: list[Row] = []
    for study in studies:
        if study["readiness_status"] == "eligible_site_evidence":
            continue
        accession = study["study_accession"]
        source = (
            "../PXD006140/proteomics_parser_v1/sites.tsv"
            if accession == "PXD006140"
            else "../evidence_preflight_v1/study_inventory.tsv"
        )
        rows.append(
            {
                "scope": "study",
                "study_accession": accession,
                "blocker_code": study["readiness_status"],
                "observed_value": study["eligible_site_count"],
                "required_value": "site-specific evidence",
                "source_locator": source,
            }
        )
    eligible_studies = sum(int(row["eligible_site_count"]) > 0 for row in studies)
    if eligible_studies == 0:
        rows.append(
            {
                "scope": "global",
                "study_accession": "",
                "blocker_code": "no_eligible_site_evidence",
                "observed_value": "0",
                "required_value": "at least one eligible site",
                "source_locator": "configs/benchmark_readiness_v1.yaml",
            }
        )
    if eligible_studies < policy.required_study_count:
        rows.append(
            {
                "scope": "global",
                "study_accession": "",
                "blocker_code": ("insufficient_distinct_site_evidence_studies"),
                "observed_value": str(eligible_studies),
                "required_value": str(policy.required_study_count),
                "source_locator": "configs/benchmark_readiness_v1.yaml",
            }
        )
    return sorted(
        rows,
        key=lambda row: (
            row["scope"],
            row["study_accession"],
            row["blocker_code"],
        ),
    )


def _summary(
    studies: list[Row],
    blockers: list[Row],
    required_study_count: int,
) -> BenchmarkReadinessSummary:
    eligible_studies = sum(int(row["eligible_site_count"]) > 0 for row in studies)
    return BenchmarkReadinessSummary(
        decision=("GO" if eligible_studies >= required_study_count else "STOP"),
        study_count=len(studies),
        eligible_study_count=eligible_studies,
        eligible_site_count=sum(int(row["eligible_site_count"]) for row in studies),
        coordinate_only_count=sum(int(row["coordinate_only_count"]) for row in studies),
        non_site_evidence_count=sum(
            int(row["non_site_evidence_count"]) for row in studies
        ),
        unresolved_candidate_count=sum(
            int(row["unresolved_candidate_count"]) for row in studies
        ),
        required_study_count=required_study_count,
        blocker_count=len(blockers),
    )


def _publish(staged: Path, output: Path) -> None:
    if not output.exists():
        staged.replace(output)
        return
    audit_benchmark_readiness(output)
    if _tree_hashes(output) != _tree_hashes(staged):
        raise RuntimeError(f"existing benchmark readiness differs: {output}")


def build_benchmark_readiness(
    policy_path: Path,
    parser_output_root: Path,
    content_output_directory: Path,
    output_directory: Path,
    registry_dir: Path,
) -> BenchmarkReadinessSummary:
    """Build a readiness decision from audited real evidence outputs."""
    policy = _load_policy(policy_path)
    audit_site_output("PXD006140", parser_output_root, registry_dir)
    audit_content_output(content_output_directory)

    parser_directory = parser_output_root / "PXD006140" / "proteomics_parser_v1"
    parser_rows = _read_tsv(parser_directory / "sites.tsv", SITE_FIELDS)
    content_studies = _read_tsv(
        content_output_directory / "study_inventory.tsv", STUDY_FIELDS
    )
    evidence_rows = _read_tsv(
        content_output_directory / "evidence_records.tsv",
        CONTENT_EVIDENCE_FIELDS,
    )
    candidate_rows = _read_tsv(
        content_output_directory / "candidate_modifications.tsv",
        CANDIDATE_FIELDS,
    )
    studies = _study_rows(
        policy,
        parser_rows,
        content_studies,
        evidence_rows,
        candidate_rows,
    )
    blockers = _blockers(studies, policy)
    summary = _summary(studies, blockers, policy.required_study_count)

    output_directory.parent.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(
        tempfile.mkdtemp(
            prefix=f".{output_directory.name}-",
            dir=output_directory.parent,
        )
    )
    staged = temporary_root / output_directory.name
    staged.mkdir()
    try:
        _write_tsv(
            staged / "study_readiness.tsv",
            STUDY_READINESS_FIELDS,
            studies,
        )
        _write_tsv(staged / "blockers.tsv", BLOCKER_FIELDS, blockers)
        inputs = [
            _input_entry(
                "readiness_policy",
                policy_path,
                output_directory,
                "repository",
            ),
            _input_entry(
                "proteomics_parser_manifest",
                parser_directory / "manifest.json",
                output_directory,
                "output_relative",
            ),
            _input_entry(
                "evidence_content_manifest",
                content_output_directory / "manifest.json",
                output_directory,
                "output_relative",
            ),
        ]
        manifest: dict[str, object] = {
            "schema_version": 1,
            "readiness_version": READINESS_VERSION,
            "decision": summary.decision,
            "benchmark_created": False,
            "labels_created": False,
            "nondetection_labeled_negative": False,
            "biological_values_modified": False,
            "counts": asdict(summary),
            "inputs": inputs,
            "outputs": [
                {"file": name, "sha256": hash_file(staged / name, "sha256")}
                for name in OUTPUT_NAMES
                if name != "manifest.json"
            ],
        }
        (staged / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        _publish(staged, output_directory)
        return summary
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def _mapping(value: object, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise RuntimeError(f"readiness manifest requires mapping: {context}")
    return cast(dict[str, Any], value)


def _items(value: object, context: str) -> list[Any]:
    if not isinstance(value, list):
        raise RuntimeError(f"readiness manifest requires list: {context}")
    return value


def audit_benchmark_readiness(
    output_directory: Path,
) -> BenchmarkReadinessSummary:
    """Rehash the readiness inputs and verify a label-free decision."""
    observed = {path.name for path in output_directory.iterdir() if path.is_file()}
    if observed != set(OUTPUT_NAMES):
        raise RuntimeError("benchmark readiness output file set mismatch")
    manifest_path = output_directory / "manifest.json"
    try:
        parsed: object = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"invalid readiness manifest: {manifest_path}") from exc
    manifest = _mapping(parsed, "root")
    if (
        manifest.get("schema_version") != 1
        or manifest.get("readiness_version") != READINESS_VERSION
    ):
        raise RuntimeError("benchmark readiness manifest identity mismatch")
    if (
        manifest.get("benchmark_created") is not False
        or manifest.get("labels_created") is not False
        or manifest.get("nondetection_labeled_negative") is not False
        or manifest.get("biological_values_modified") is not False
    ):
        raise RuntimeError("benchmark readiness integrity flags are invalid")

    input_paths: dict[str, Path] = {}
    for raw_entry in _items(manifest.get("inputs"), "inputs"):
        entry = _mapping(raw_entry, "inputs[]")
        resolution = str(entry.get("resolution", ""))
        stored_path = Path(str(entry.get("path", "")))
        path = (
            stored_path
            if resolution == "repository"
            else output_directory / stored_path
        )
        if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
            raise RuntimeError(f"readiness input SHA256 mismatch: {path}")
        input_paths[str(entry.get("name", ""))] = path
    if set(input_paths) != {
        "readiness_policy",
        "proteomics_parser_manifest",
        "evidence_content_manifest",
    }:
        raise RuntimeError("benchmark readiness input list mismatch")
    policy = _load_policy(input_paths["readiness_policy"])
    parser_directory = input_paths["proteomics_parser_manifest"].parent
    audit_site_output("PXD006140", parser_directory.parent.parent)
    audit_content_output(input_paths["evidence_content_manifest"].parent)

    declared: set[str] = set()
    for raw_entry in _items(manifest.get("outputs"), "outputs"):
        entry = _mapping(raw_entry, "outputs[]")
        name = str(entry.get("file", ""))
        declared.add(name)
        path = output_directory / name
        if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
            raise RuntimeError(f"readiness output SHA256 mismatch: {path}")
    if declared != set(OUTPUT_NAMES) - {"manifest.json"}:
        raise RuntimeError("benchmark readiness output list mismatch")

    studies = _read_tsv(
        output_directory / "study_readiness.tsv", STUDY_READINESS_FIELDS
    )
    blockers = _read_tsv(output_directory / "blockers.tsv", BLOCKER_FIELDS)
    if tuple(row["study_accession"] for row in studies) != policy.studies:
        raise RuntimeError("benchmark readiness study order mismatch")
    if any(
        "label" in field.lower() for field in STUDY_READINESS_FIELDS + BLOCKER_FIELDS
    ):
        raise RuntimeError("benchmark readiness output contains a label field")
    summary = _summary(studies, blockers, policy.required_study_count)
    if _mapping(manifest.get("counts"), "counts") != asdict(summary):
        raise RuntimeError("benchmark readiness counts mismatch")
    if manifest.get("decision") != summary.decision:
        raise RuntimeError("benchmark readiness decision mismatch")
    if summary.decision == "STOP" and not blockers:
        raise RuntimeError("benchmark STOP requires explicit blockers")
    return summary
