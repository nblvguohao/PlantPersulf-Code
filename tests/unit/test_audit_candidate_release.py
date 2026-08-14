"""Gate 0: audit of the candidate release package.

Runs the real audit against the release package: artifact existence and
pinned hashes, bundle load + deterministic scoring, predictive-claims
discipline ([]), JSON parsing, and the co-signature gate.

The co-signature gate is checked from both sides. Against the live package
it must PASS with ``--expect-signed``, because Gate 0 was co-signed and
frozen on 2026-08-13 and must stay that way; a future edit that reverts a
protocol document to draft has to break this suite. The failure branch is
checked against a synthetic draft package instead of the live one, so the
test asserts the gate's behaviour rather than the repository's state on the
day it was written.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

_RELEASE_DIR = (
    Path(__file__).resolve().parents[2]
    / "results"
    / "candidates"
    / "multispecies_v2_candidate_release_v1"
)
_AUDIT_SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts" / "audit_candidate_release.py"
)

pytestmark = pytest.mark.skipif(
    not (_RELEASE_DIR / "manifest.json").is_file() or not _AUDIT_SCRIPT.is_file(),
    reason="release package or audit script not present",
)


def _run_audit(release_dir: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(_AUDIT_SCRIPT), "--release-dir", str(release_dir), *extra],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def test_audit_passes_on_the_frozen_release() -> None:
    result = _run_audit(_RELEASE_DIR)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "AUDIT PASS" in result.stdout


def test_audit_expect_signed_passes_because_gate0_is_co_signed() -> None:
    result = _run_audit(_RELEASE_DIR, "--expect-signed")
    assert result.returncode == 0, result.stdout + result.stderr
    assert "AUDIT PASS" in result.stdout


def test_audit_expect_signed_rejects_a_draft_package(tmp_path: Path) -> None:
    draft_dir = tmp_path / "draft_release"
    draft_dir.mkdir()
    (draft_dir / "sap_protocol.json").write_text(
        json.dumps({"version": "0.1"}), encoding="utf-8"
    )
    (draft_dir / "manifest.json").write_text(
        json.dumps(
            {
                "produced_artifacts": [
                    {
                        "artifact": "sap_protocol.json",
                        "status": "draft_pending_co_signature",
                        "sha256": "draft - hash recorded at final freeze",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )

    result = _run_audit(draft_dir, "--expect-signed")
    assert result.returncode == 1
    assert "DRAFT (unsigned): sap_protocol.json" in result.stdout
    assert "still draft" in result.stdout


def test_audit_checks_predictive_claims() -> None:
    result = _run_audit(_RELEASE_DIR)
    assert "predictive claims: clean ([])" in result.stdout
