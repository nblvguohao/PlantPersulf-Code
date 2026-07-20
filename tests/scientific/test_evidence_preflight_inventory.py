import csv
import hashlib
import json
from pathlib import Path

import pytest

POLICY_PATH = Path("configs/evidence_preflight_v1.yaml")
SELECTION_PATH = Path("configs/download_selection.yaml")
REGISTRY_DIR = Path("data/registry")
OUTPUT_NAMES = {
    "study_inventory.tsv",
    "file_decisions.tsv",
    "evidence_records.tsv",
    "metadata_gaps.tsv",
    "large_file_queue.tsv",
    "manifest.json",
}


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _build(output: Path):  # type: ignore[no-untyped-def]
    from plantpersulf.evidence.preflight import build_metadata_preflight

    return build_metadata_preflight(
        policy_path=POLICY_PATH,
        selection_path=SELECTION_PATH,
        registry_dir=REGISTRY_DIR,
        output_directory=output,
    )


def test_real_registry_preflight_is_deterministic_and_unlabeled(
    tmp_path: Path,
) -> None:
    first_output = tmp_path / "first/evidence_preflight_v1"
    second_output = tmp_path / "second/evidence_preflight_v1"

    first = _build(first_output)
    second = _build(second_output)

    assert first == second
    assert first.study_count == 4
    assert first.file_count == 77
    assert first.registered_metadata_count == 4
    assert first.downloaded_count == 7
    assert first.downloaded_existing_scope_count == 2
    assert first.deferred_large_count == 56
    assert first.deferred_count == 8
    assert first.evidence_record_count == 0
    assert {path.name for path in first_output.iterdir()} == OUTPUT_NAMES
    assert _tree_hashes(first_output) == _tree_hashes(second_output)

    decisions = _rows(first_output / "file_decisions.tsv")
    assert len(decisions) == 77
    assert all(row["evidence_class"] == "unresolved" for row in decisions)
    assert [
        row["decision"]
        for row in decisions
        if row["dataset_accession"] == "PXD024061"
        and row["file_name"] == "txt_persulfproject.zip"
    ] == ["deferred_large_not_stage_a"]
    assert [
        row["registry_size_bytes"]
        for row in decisions
        if row["dataset_accession"] == "PXD024061"
        and row["file_name"] == "txt_persulfproject.zip"
    ] == ["2417963957"]
    assert all(
        row["decision"].startswith("deferred")
        for row in decisions
        if row["file_category"] in {"RAW", "PEAK"}
        and row["file_name"] != "checksum.txt"
    )
    assert [
        row["decision"]
        for row in decisions
        if row["dataset_accession"] == "PXD039999"
        and row["file_name"] == "checksum.txt"
    ] == ["downloaded"]
    mzid = next(
        row
        for row in decisions
        if row["file_name"] == "peptides_1_1_0.mzid.gz"
    )
    assert mzid["registry_size_bytes"] == "1967088"
    assert mzid["local_size_bytes"] == "184107"
    assert mzid["local_sha256"] == (
        "62105dfde42d9675bc5dc7c83969d971eec57c9e02983659cd2df96ec4d1b5fe"
    )

    evidence_rows = _rows(first_output / "evidence_records.tsv")
    assert evidence_rows == []
    assert "label" not in _rows(first_output / "study_inventory.tsv")[0]
    manifest = json.loads(
        (first_output / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["content_audit_complete"] is False
    assert manifest["labels_created"] is False
    assert manifest["nondetection_labeled_negative"] is False
    assert manifest["biological_values_modified"] is False


def test_preflight_reports_real_metadata_gaps(tmp_path: Path) -> None:
    output = tmp_path / "evidence_preflight_v1"
    _build(output)

    gaps = _rows(output / "metadata_gaps.tsv")
    assert any(
        row["dataset_accession"] == "PXD006140"
        and row["gap_type"] == "missing_sample_mapping"
        for row in gaps
    )
    assert {
        row["file_name"]
        for row in gaps
        if row["gap_type"] == "content_audit_pending"
    } == {
        "SDRF.txt",
        "peptide.csv",
        "peptides_1_1_0.mzid.gz",
        "proteins.csv",
    }


def test_preflight_audit_detects_declared_hash_corruption(tmp_path: Path) -> None:
    from plantpersulf.evidence.preflight import audit_metadata_preflight

    output = tmp_path / "evidence_preflight_v1"
    clean = _build(output)
    assert audit_metadata_preflight(output) == clean

    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outputs"][0]["sha256"] = "0" * 64
    manifest_path.write_bytes(
        (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    )

    with pytest.raises(RuntimeError, match="output SHA256 mismatch"):
        audit_metadata_preflight(output)


def test_cli_exposes_evidence_preflight_commands() -> None:
    from plantpersulf.cli import build_parser

    build = build_parser().parse_args(
        ["build-evidence-preflight", "--version", "v1"]
    )
    audit = build_parser().parse_args(
        ["audit-evidence-preflight", "--version", "v1"]
    )

    assert build.command == "build-evidence-preflight"
    assert audit.command == "audit-evidence-preflight"
