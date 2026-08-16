import csv
import hashlib
import json
from pathlib import Path

import pytest

POLICY = Path("configs/benchmark_readiness_v1.yaml")
REGISTRY = Path("data/registry")


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _build_real_prerequisites(tmp_path: Path) -> tuple[Path, Path]:
    from plantpersulf.evidence.content import build_content_audit
    from plantpersulf.proteomics.peptide_parser import (
        parse_proteomics_accession,
    )

    parser_root = tmp_path / "interim"
    content_directory = parser_root / "evidence_preflight_v1"
    parse_proteomics_accession(
        accession="PXD006140",
        output_root=parser_root,
        registry_dir=REGISTRY,
    )
    build_content_audit(
        registry_dir=REGISTRY,
        output_directory=content_directory,
    )
    return parser_root, content_directory


def _build_readiness(tmp_path: Path, name: str = "readiness") -> Path:
    from plantpersulf.benchmark.readiness import build_benchmark_readiness

    parser_root, content_directory = _build_real_prerequisites(tmp_path)
    output = tmp_path / name
    build_benchmark_readiness(
        policy_path=POLICY,
        parser_output_root=parser_root,
        content_output_directory=content_directory,
        output_directory=output,
        registry_dir=REGISTRY,
    )
    return output


def test_current_real_evidence_forces_benchmark_stop(tmp_path: Path) -> None:
    from plantpersulf.benchmark.readiness import build_benchmark_readiness

    parser_root, content_directory = _build_real_prerequisites(tmp_path)
    summary = build_benchmark_readiness(
        policy_path=POLICY,
        parser_output_root=parser_root,
        content_output_directory=content_directory,
        output_directory=tmp_path / "readiness",
        registry_dir=REGISTRY,
    )

    assert summary.decision == "STOP"
    assert summary.eligible_site_count == 0
    assert summary.eligible_study_count == 0
    assert summary.coordinate_only_count == 3
    assert summary.non_site_evidence_count == 10_787
    assert summary.unresolved_candidate_count == 25
    assert summary.required_study_count == 2


def test_readiness_output_is_deterministic_and_has_no_labels(
    tmp_path: Path,
) -> None:
    first = _build_readiness(tmp_path / "first")
    second = _build_readiness(tmp_path / "second")

    assert {path.name for path in first.iterdir()} == {
        "study_readiness.tsv",
        "blockers.tsv",
        "manifest.json",
    }
    assert _tree_hashes(first) == _tree_hashes(second)
    for tsv_path in first.glob("*.tsv"):
        header = tsv_path.read_text(encoding="utf-8").splitlines()[0]
        assert "label" not in header.lower()
    manifest = json.loads((first / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["decision"] == "STOP"
    assert manifest["benchmark_created"] is False
    assert manifest["labels_created"] is False
    assert manifest["nondetection_labeled_negative"] is False
    assert manifest["biological_values_modified"] is False
    assert not (tmp_path / "data/processed/benchmark_v1").exists()


def test_readiness_preserves_study_specific_blockers(tmp_path: Path) -> None:
    output = _build_readiness(tmp_path)
    studies = _rows(output / "study_readiness.tsv")
    blockers = _rows(output / "blockers.tsv")

    assert [row["study_accession"] for row in studies] == [
        "PXD006140",
        "PXD024061",
        "PXD035795",
        "PXD039999",
    ]
    assert {row["readiness_status"] for row in studies} == {
        "coordinate_only_not_site_evidence",
        "non_site_evidence_only",
        "unsupported_stage_a_checksum_only",
    }
    assert any(
        row["blocker_code"] == "insufficient_distinct_site_evidence_studies"
        and row["observed_value"] == "0"
        and row["required_value"] == "2"
        for row in blockers
    )


def test_readiness_audit_rejects_output_tampering(tmp_path: Path) -> None:
    from plantpersulf.benchmark.readiness import audit_benchmark_readiness

    output = _build_readiness(tmp_path)
    audit_benchmark_readiness(output)
    with (output / "blockers.tsv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    rows[-1][-1] = "changed"
    with (output / "blockers.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)

    with pytest.raises(RuntimeError, match="output SHA256 mismatch"):
        audit_benchmark_readiness(output)


def test_cli_exposes_benchmark_readiness_without_building_benchmark() -> None:
    from plantpersulf.cli import build_parser

    build = build_parser().parse_args(["build-benchmark-readiness", "--version", "v1"])
    audit = build_parser().parse_args(["audit-benchmark-readiness", "--version", "v1"])

    assert build.command == "build-benchmark-readiness"
    assert audit.command == "audit-benchmark-readiness"
