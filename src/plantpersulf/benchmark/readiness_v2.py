"""Benchmark readiness gate v2 — reads published persulfidation site evidence.

V2 replaced the frozen V1 gate (which saw only psm_coordinate_only/unresolved
evidence and correctly returned STOP) after Phase A's evidence-acquisition
cycles. It consumes the deterministic persulfidation site outputs published by
``persulfidation_publish.py``, counts eligible ``site_ms`` records per study,
and issues GO only when at least two distinct studies hold such evidence.

The V1 gate and its policy remain frozen and unaltered.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import cast

import yaml

from plantpersulf.proteomics.persulfidation_publish import (
    PERSULF_SITE_FIELDS,
    audit_persulfidation_sites,
)
from plantpersulf.provenance.hashing import hash_file

READINESS_VERSION = 2
OUTPUT_NAMES = (
    "study_readiness.tsv",
    "blockers.tsv",
    "manifest.json",
)
STUDY_READINESS_FIELDS = (
    "study_accession",
    "eligible_site_count",
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
ELIGIBLE_LEVELS = ("site_ms",)
PROHIBITED_LABELS = ("positive", "unlabeled", "negative")

Row = dict[str, str]


@dataclass(frozen=True)
class ReadinessV2Policy:
    version: int
    studies: tuple[str, ...]
    eligible_levels: tuple[str, ...]
    required_study_count: int


@dataclass(frozen=True)
class BenchmarkReadinessV2Summary:
    decision: str
    study_count: int
    eligible_study_count: int
    eligible_site_count: int
    required_study_count: int
    blocker_count: int


def _load_policy(path: Path) -> ReadinessV2Policy:
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid benchmark readiness v2 policy: {path}") from exc
    if not isinstance(loaded, dict):
        raise RuntimeError("benchmark readiness v2 policy requires a mapping")
    expected = {
        "version": 2,
        "studies": ["PXD006140", "PXD024061"],
        "eligible_site_evidence_levels": ["site_ms"],
        "minimum_distinct_site_evidence_studies": 2,
        "prohibited_record_labels": ["positive", "unlabeled", "negative"],
    }
    if loaded != expected:
        raise RuntimeError("benchmark readiness v2 policy differs from frozen v2")
    return ReadinessV2Policy(
        version=2,
        studies=("PXD006140", "PXD024061"),
        eligible_levels=ELIGIBLE_LEVELS,
        required_study_count=2,
    )


def _read_tsv(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    import csv

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise RuntimeError(f"readiness source has invalid columns: {path}")
        rows = [dict(row) for row in reader]
    if any(
        None in row or any(value is None for value in row.values()) for row in rows
    ):
        raise RuntimeError(f"readiness source has malformed row: {path}")
    return cast(list[dict[str, str]], rows)


def _write_tsv(
    path: Path,
    fieldnames: tuple[str, ...],
    records: list[Row],
) -> None:
    import csv

    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(records)


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hash_file(path, "sha256")
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _relative_to_output(path: Path, output_directory: Path) -> str:
    try:
        return Path(os.path.relpath(path, output_directory)).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _blockers(studies: list[Row], required_count: int) -> list[Row]:
    rows: list[Row] = []
    for study in studies:
        if study["readiness_status"] == "eligible_site_evidence":
            continue
        rows.append(
            {
                "scope": "study",
                "study_accession": study["study_accession"],
                "blocker_code": study["readiness_status"],
                "observed_value": study["eligible_site_count"],
                "required_value": "site-specific evidence",
                "source_locator": (
                    f"../{study['study_accession']}/"
                    f"persulfidation_sites_v1/sites.tsv"
                ),
            }
        )
    eligible = sum(int(row["eligible_site_count"]) > 0 for row in studies)
    if eligible == 0:
        rows.append(
            {
                "scope": "global",
                "study_accession": "",
                "blocker_code": "no_eligible_site_evidence",
                "observed_value": "0",
                "required_value": "at least one eligible site",
                "source_locator": "configs/benchmark_readiness_v2.yaml",
            }
        )
    if eligible < required_count:
        rows.append(
            {
                "scope": "global",
                "study_accession": "",
                "blocker_code": (
                    "insufficient_distinct_site_evidence_studies"
                ),
                "observed_value": str(eligible),
                "required_value": str(required_count),
                "source_locator": "configs/benchmark_readiness_v2.yaml",
            }
        )
    return sorted(
        rows,
        key=lambda row: (row["scope"], row["study_accession"], row["blocker_code"]),
    )


def _summary(
    studies: list[Row],
    blockers: list[Row],
    required: int,
) -> BenchmarkReadinessV2Summary:
    eligible = sum(int(row["eligible_site_count"]) > 0 for row in studies)
    return BenchmarkReadinessV2Summary(
        decision="GO" if eligible >= required else "STOP",
        study_count=len(studies),
        eligible_study_count=eligible,
        eligible_site_count=sum(int(row["eligible_site_count"]) for row in studies),
        required_study_count=required,
        blocker_count=len(blockers),
    )


def _input_entry(
    name: str,
    path: Path,
    output_directory: Path,
) -> dict[str, str]:
    stored_path = _relative_to_output(path, output_directory)
    return {
        "name": name,
        "path": stored_path,
        "sha256": hash_file(path, "sha256"),
    }


def build_benchmark_readiness_v2(
    policy_path: Path,
    site_output_root: Path,
    output_directory: Path,
) -> BenchmarkReadinessV2Summary:
    """Build the v2 readiness decision from published persulfidation site outputs."""
    policy = _load_policy(policy_path)

    studies: list[Row] = []
    for accession in policy.studies:
        site_dir = site_output_root / accession / "persulfidation_sites_v1"
        try:
            audit_persulfidation_sites(accession, site_dir)
        except (RuntimeError, OSError):
            studies.append(
                {
                    "study_accession": accession,
                    "eligible_site_count": "0",
                    "readiness_status": "site_audit_failed",
                }
            )
            continue
        sites = _read_tsv(site_dir / "sites.tsv", PERSULF_SITE_FIELDS)
        eligible_count = sum(
            row["evidence_level"] in policy.eligible_levels for row in sites
        )
        studies.append(
            {
                "study_accession": accession,
                "eligible_site_count": str(eligible_count),
                "readiness_status": (
                    "eligible_site_evidence"
                    if eligible_count > 0
                    else "no_site_evidence"
                ),
            }
        )

    blockers = _blockers(studies, policy.required_study_count)
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
                "readiness_policy_v2",
                policy_path,
                output_directory,
            ),
        ]
        for accession in policy.studies:
            site_manifest = (
                site_output_root
                / accession
                / "persulfidation_sites_v1"
                / "manifest.json"
            )
            if site_manifest.is_file():
                inputs.append(
                    _input_entry(
                        f"site_manifest_{accession}",
                        site_manifest,
                        output_directory,
                    )
                )
        manifest: dict[str, object] = {
            "schema_version": 2,
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
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        if output_directory.exists():
            if _tree_hashes(output_directory) != _tree_hashes(staged):
                raise RuntimeError(
                    f"existing benchmark readiness v2 differs: {output_directory}"
                )
            return summary
        staged.replace(output_directory)
        return summary
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)
