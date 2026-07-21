import csv
import hashlib
import json
from pathlib import Path

OUTPUT_NAMES = {
    "study_inventory.tsv",
    "file_decisions.tsv",
    "evidence_records.tsv",
    "metadata_gaps.tsv",
    "large_file_queue.tsv",
    "candidate_modifications.tsv",
    "sample_file_mapping.tsv",
    "mapping_conflicts.tsv",
    "content_schema.tsv",
    "manifest.json",
}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(
            path.read_bytes()
        ).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _build(output: Path):  # type: ignore[no-untyped-def]
    from plantpersulf.evidence.content import build_content_audit

    return build_content_audit(output_directory=output)


def test_real_content_audit_is_deterministic_and_label_free(
    tmp_path: Path,
) -> None:
    first_output = tmp_path / "first/evidence_preflight_v1"
    second_output = tmp_path / "second/evidence_preflight_v1"

    first = _build(first_output)
    second = _build(second_output)

    assert first == second
    assert first.study_count == 4
    assert first.evidence_record_count == 10_787
    assert first.candidate_modification_count == 25
    assert first.sample_file_mapping_count == 12
    assert first.mapping_conflict_count == 6
    assert first.content_schema_count == 4
    assert first.site_ms_count == 0
    assert {path.name for path in first_output.iterdir()} == OUTPUT_NAMES
    assert _tree_hashes(first_output) == _tree_hashes(second_output)

    manifest = json.loads(
        (first_output / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["content_audit_complete"] is True
    assert manifest["labels_created"] is False
    assert manifest["nondetection_labeled_negative"] is False
    assert manifest["biological_values_modified"] is False
    assert manifest["counts"]["site_ms_count"] == 0


def test_content_records_preserve_evidence_levels_without_site_promotion(
    tmp_path: Path,
) -> None:
    output = tmp_path / "evidence_preflight_v1"
    _build(output)

    evidence = _rows(output / "evidence_records.tsv")
    assert len(evidence) == 10_787
    assert sum(row["record_type"] == "peptide_csv" for row in evidence) == 9326
    assert sum(row["record_type"] == "protein_csv" for row in evidence) == 1461
    assert {
        row["evidence_class"]
        for row in evidence
        if row["record_type"] == "peptide_csv"
    } == {"identification_only"}
    assert {
        row["evidence_class"]
        for row in evidence
        if row["record_type"] == "protein_csv"
    } == {"protein_level_only"}
    assert not any(row["evidence_class"] == "site_ms" for row in evidence)
    assert "label" not in evidence[0]
    assert all(len(row["source_sha256"]) == 64 for row in evidence)
    assert all(len(row["method_source_sha256"]) == 64 for row in evidence)
    assert json.loads(evidence[0]["raw_values_json"])["Avg. Area"] == (
        "1.2568E9"
    )


def test_candidates_and_file_mapping_conflicts_remain_explicit(
    tmp_path: Path,
) -> None:
    output = tmp_path / "evidence_preflight_v1"
    _build(output)

    candidates = _rows(output / "candidate_modifications.tsv")
    mappings = _rows(output / "sample_file_mapping.tsv")
    conflicts = _rows(output / "mapping_conflicts.tsv")
    schema = _rows(output / "content_schema.tsv")

    assert len(candidates) == 25
    assert {row["evidence_class"] for row in candidates} == {"unresolved"}
    assert {row["method_mapping_status"] for row in candidates} == {
        "absent_v1"
    }
    assert len(mappings) == 12
    assert sum(row["mapping_status"] == "exact" for row in mappings) == 6
    assert sum(row["mapping_status"] == "conflict" for row in mappings) == 6
    assert len(conflicts) == 6
    assert {row["conflict_type"] for row in conflicts} == {
        "unmatched_exact_filename"
    }
    assert len(schema) == 4
    assert {row["parse_status"] for row in schema} == {"parsed"}


def test_content_audit_replaces_pending_gaps_with_observed_limits(
    tmp_path: Path,
) -> None:
    output = tmp_path / "evidence_preflight_v1"
    _build(output)

    gaps = _rows(output / "metadata_gaps.tsv")
    assert not any(row["gap_type"] == "content_audit_pending" for row in gaps)
    pxd035795_gaps = {
        row["gap_type"]
        for row in gaps
        if row["dataset_accession"] == "PXD035795"
    }
    assert pxd035795_gaps == {
        "partial_sample_file_mapping",
        "quantification_mapping_unresolved",
        "site_method_mapping_absent",
    }


def test_content_audit_cli_and_hash_audit(tmp_path: Path) -> None:
    from plantpersulf.cli import build_parser
    from plantpersulf.evidence.content import audit_content_output

    output = tmp_path / "evidence_preflight_v1"
    summary = _build(output)
    assert audit_content_output(output) == summary

    build = build_parser().parse_args(
        ["build-evidence-content-audit", "--version", "v1"]
    )
    audit = build_parser().parse_args(
        ["audit-evidence-content", "--version", "v1"]
    )
    assert build.command == "build-evidence-content-audit"
    assert audit.command == "audit-evidence-content"
