"""RED (Phase A3): v2 benchmark readiness gate — official GO from public data.

The frozen v1 gate (``configs/benchmark_readiness_v1.yaml`` + ``readiness.py``)
was designed to fail-closed on the pre-persulfidation evidence state (0 eligible
sites, STOP). The v2 gate reads the two study outputs published by
``persulfidation_publish.py`` — PXD006140 (Dataset S3, 320 class-I sites) and
PXD024061 (MaxQuant, 73 class-I sites) — and must output GO.

Expected RED: ``build_benchmark_readiness_v2`` / ``audit_benchmark_readiness_v2``
do not exist yet, so collection fails with ImportError (or the policy file is
missing).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.benchmark.readiness_v2 import (  # RED: module missing
    build_benchmark_readiness_v2,
)


def _publish_sites(tmp_path: Path) -> tuple[Path, Path]:
    """Publish both study outputs into tmp_path so the test is self-contained."""
    from plantpersulf.proteomics.persulfidation_publish import (
        publish_persulfidation_sites,
    )

    publish_persulfidation_sites("PXD006140", output_root=tmp_path)
    publish_persulfidation_sites("PXD024061", output_root=tmp_path)
    site_root = tmp_path
    output = tmp_path / "readiness_v2"
    return site_root, output


def test_real_evidence_outputs_yield_go(tmp_path: Path) -> None:
    site_root, output_dir = _publish_sites(tmp_path)
    summary = build_benchmark_readiness_v2(
        policy_path=Path("configs/benchmark_readiness_v2.yaml"),
        site_output_root=site_root,
        output_directory=output_dir,
    )

    assert summary.decision == "GO"
    assert summary.eligible_study_count == 2
    assert summary.eligible_site_count >= 320 + 73
    assert summary.study_count == 2  # only the two confirmed studies


def test_go_manifest_has_no_labels_and_is_fail_safe(tmp_path: Path) -> None:
    site_root, output_dir = _publish_sites(tmp_path)
    build_benchmark_readiness_v2(
        policy_path=Path("configs/benchmark_readiness_v2.yaml"),
        site_output_root=site_root,
        output_directory=output_dir,
    )

    import json

    manifest = json.loads(
        (output_dir / "manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["decision"] == "GO"
    assert manifest["benchmark_created"] is False  # benchmark not yet built
    assert manifest["labels_created"] is False
    assert manifest["nondetection_labeled_negative"] is False
    assert manifest["biological_values_modified"] is False
    assert {path.name for path in output_dir.iterdir()} == {
        "study_readiness.tsv",
        "blockers.tsv",
        "manifest.json",
    }


def test_insufficient_studies_still_goes_stop(tmp_path: Path) -> None:
    """When fewer than the minimum studies are available, the gate must fail-closed."""
    site_root, output_dir = _publish_sites(tmp_path)
    # Tamper: remove PXD024061 output to simulate only one study.
    import shutil

    shutil.rmtree(site_root / "PXD024061")
    summary = build_benchmark_readiness_v2(
        policy_path=Path("configs/benchmark_readiness_v2.yaml"),
        site_output_root=site_root,
        output_directory=output_dir / "single_study",
    )

    assert summary.decision == "STOP"
    assert summary.eligible_study_count == 1
    assert summary.blocker_count > 0


def test_go_output_is_deterministic(tmp_path: Path) -> None:
    site_root, output_dir = _publish_sites(tmp_path)
    first = build_benchmark_readiness_v2(
        policy_path=Path("configs/benchmark_readiness_v2.yaml"),
        site_output_root=site_root,
        output_directory=output_dir,
    )
    # Re-publish into the same directory must produce byte-identical output.
    second = build_benchmark_readiness_v2(
        policy_path=Path("configs/benchmark_readiness_v2.yaml"),
        site_output_root=site_root,
        output_directory=output_dir,
    )
    assert first == second


def test_missing_policy_rejects_without_defaulting(tmp_path: Path) -> None:
    site_root, _ = _publish_sites(tmp_path)
    with pytest.raises((RuntimeError, OSError)):
        build_benchmark_readiness_v2(
            policy_path=Path("configs/benchmark_readiness_v2_nonexistent.yaml"),
            site_output_root=site_root,
            output_directory=tmp_path / "no_policy",
        )
