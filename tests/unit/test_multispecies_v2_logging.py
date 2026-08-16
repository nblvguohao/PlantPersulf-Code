"""Reproducible training-log and manifest contracts for v2."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from plantpersulf.workflows.multispecies_v2 import (
    TrainingEvent,
    append_training_event,
    write_run_manifest,
)


def test_training_event_writer_emits_required_jsonl_fields(tmp_path: Path) -> None:
    """Dropping fold or validation AP would make a training run unauditable."""
    path = tmp_path / "training.jsonl"
    append_training_event(
        path,
        TrainingEvent(
            timestamp="2026-08-11T00:00:00Z",
            track="strict_cluster_holdout",
            fold=2,
            seed=3,
            model="pu_logistic",
            epoch=1,
            loss=0.4,
            val_ap=0.2,
            lr=0.01,
            device="cuda:0",
            gpu_memory=1024,
            wall_seconds=5.0,
        ),
    )

    assert json.loads(path.read_text(encoding="utf-8")) == {
        "device": "cuda:0",
        "epoch": 1,
        "fold": 2,
        "gpu_memory": 1024,
        "loss": 0.4,
        "lr": 0.01,
        "model": "pu_logistic",
        "seed": 3,
        "timestamp": "2026-08-11T00:00:00Z",
        "track": "strict_cluster_holdout",
        "val_ap": 0.2,
        "wall_seconds": 5.0,
    }


def test_manifest_rejects_missing_input_hash(tmp_path: Path) -> None:
    """A result without every input hash cannot be used for a scientific claim."""
    with pytest.raises(ValueError, match="input_sha256"):
        write_run_manifest(
            tmp_path / "manifest.json",
                config_sha256="a" * 64,
                split_sha256="b" * 64,
                code_revision="abc123",
                code_sha256="c" * 64,
                input_sha256={},
            command=["python", "script.py"],
            device="cpu",
        )


def test_training_event_rejects_non_finite_metric(tmp_path: Path) -> None:
    """NaN metrics would make a JSONL training audit non-reproducible."""
    with pytest.raises(ValueError, match="finite"):
        append_training_event(
            tmp_path / "training.jsonl",
            TrainingEvent(
                timestamp="2026-08-11T00:00:00Z",
                track="strict_cluster_holdout",
                fold=0,
                seed=0,
                model="pu_logistic",
                epoch=0,
                loss=float("nan"),
                val_ap=None,
                lr=0.01,
                device="cpu",
                gpu_memory=None,
                wall_seconds=0.0,
            ),
        )
