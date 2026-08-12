"""Fail-closed runtime contracts for multispecies v2 task execution."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from plantpersulf.workflows import multispecies_v2 as runtime


def test_checkpoint_resumes_only_with_identical_fingerprint(tmp_path: Path) -> None:
    """Changing a frozen input/config/code identity must reject resume."""
    fingerprint = runtime.TaskFingerprint(
        track="strict_cluster_holdout",
        fold=0,
        seed=0,
        model="pu_logistic",
        input_sha256={"development_sites": "a" * 64},
        config_sha256="b" * 64,
        code_revision="commit-marker",
    )
    checkpoint = tmp_path / "fold0.checkpoint.json"

    runtime.write_task_checkpoint(checkpoint, fingerprint, {"epoch": 3})

    assert runtime.load_resumable_checkpoint(checkpoint, fingerprint) == {"epoch": 3}
    with pytest.raises(RuntimeError, match="checkpoint fingerprint mismatch"):
        runtime.load_resumable_checkpoint(
            checkpoint,
            replace(fingerprint, config_sha256="c" * 64),
        )


def test_checkpoint_rejects_tampered_state(tmp_path: Path) -> None:
    """A changed training state must not be mistaken for resumable output."""
    fingerprint = runtime.TaskFingerprint(
        track="strict_cluster_holdout",
        fold=1,
        seed=0,
        model="pu_logistic",
        input_sha256={"development_sites": "a" * 64},
        config_sha256="b" * 64,
        code_revision="commit-marker",
    )
    checkpoint = tmp_path / "fold1.checkpoint.json"
    runtime.write_task_checkpoint(checkpoint, fingerprint, {"epoch": 3})
    payload = checkpoint.read_text(encoding="utf-8").replace('"epoch":3', '"epoch":4')
    checkpoint.write_text(payload, encoding="utf-8")

    with pytest.raises(RuntimeError, match="checkpoint state hash mismatch"):
        runtime.load_resumable_checkpoint(checkpoint, fingerprint)
