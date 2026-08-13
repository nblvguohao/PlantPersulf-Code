"""Gate 0: pre-freeze audit of the candidate release package.

Runs the real audit against the release package: artifact existence and
pinned hashes, bundle load + deterministic scoring, predictive-claims
discipline ([]), JSON parsing, and the co-signature gate (which must
FAIL with --expect-signed while protocol documents are unsigned).
"""

from __future__ import annotations

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
_AUDIT_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "audit_candidate_release.py"

pytestmark = pytest.mark.skipif(
    not (_RELEASE_DIR / "manifest.json").is_file()
    or not _AUDIT_SCRIPT.is_file(),
    reason="release package or audit script not present",
)


def _run_audit(*extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(_AUDIT_SCRIPT), "--release-dir", str(_RELEASE_DIR), *extra],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


def test_audit_passes_before_signing() -> None:
    result = _run_audit()
    assert result.returncode == 0, result.stdout + result.stderr
    assert "AUDIT PASS" in result.stdout


def test_audit_fails_when_expects_signed_but_draft() -> None:
    result = _run_audit("--expect-signed")
    assert result.returncode == 1
    assert "still draft" in result.stdout


def test_audit_checks_predictive_claims() -> None:
    result = _run_audit()
    assert "predictive claims: clean ([])" in result.stdout
