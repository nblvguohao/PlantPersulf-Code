"""Release audits for the multispecies v2 runtime manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from plantpersulf.workflows.multispecies_v2 import (
    FROZEN_TEST_MODELS,
    FROZEN_TEST_SEED,
    FROZEN_TEST_TRACK,
    TaskFingerprint,
    TrainingEvent,
    append_training_event,
    audit_v2_run_manifest,
    collect_runtime_environment,
    write_run_manifest,
    write_task_checkpoint,
)


def _config(tmp_path: Path, *, folds: int = 1) -> Path:
    path = tmp_path / "runtime.yaml"
    path.write_text(
        "strict_cluster_holdout:\n"
        f"  development_folds: {folds}\n"
        "literature_random_protein:\n"
        "  seeds: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]\n"
        "models: [pu_logistic, random_forest, xgboost, "
        "esm_linear_head, structure_ranker]\n"
        "comparison_inputs:\n  sul_environment_manifest: environment.json\n",
        encoding="utf-8",
    )
    return path


def _task_fingerprint(config_sha256: str = "a" * 64) -> TaskFingerprint:
    return TaskFingerprint(
        track="strict_cluster_holdout",
        fold=0,
        seed=0,
        model="pu_logistic",
        input_sha256={"development_sites": "c" * 64},
        config_sha256=config_sha256,
        code_revision="commit-marker",
        code_sha256="d" * 64,
    )


def _training_log(path: Path) -> Path:
    append_training_event(
        path,
        TrainingEvent(
            timestamp="2026-08-12T00:00:00+00:00",
            track="strict_cluster_holdout",
            fold=0,
            seed=0,
            model="pu_logistic",
            epoch=0,
            loss=None,
            val_ap=0.5,
            lr=None,
            device="cpu",
            gpu_memory=None,
            wall_seconds=1.0,
        ),
    )
    return path


def test_release_manifest_rejects_tampered_checkpoint_hash(tmp_path: Path) -> None:
    """A changed checkpoint makes a scientific execution manifest invalid."""
    checkpoint = tmp_path / "fold0.checkpoint.json"
    checkpoint.write_text('{"policy":"marker"}\n', encoding="utf-8")
    log = _training_log(tmp_path / "training.jsonl")
    manifest_path = tmp_path / "manifest.json"
    config = _config(tmp_path)
    config_sha = hashlib.sha256(config.read_bytes()).hexdigest()
    write_run_manifest(
        manifest_path,
        config_sha256=config_sha,
        config_path=config,
        split_sha256="b" * 64,
        code_revision="commit-marker",
        code_sha256="d" * 64,
        input_sha256={"development_sites": "c" * 64},
        command=["python", "script.py"],
        device="cpu",
        environment=collect_runtime_environment(),
        gpu_map=(),
        artifacts={"training_log": log, "checkpoint": checkpoint},
        checkpoint_paths={"fold0": checkpoint},
        task_roster=(_task_fingerprint(config_sha),),
    )
    checkpoint.write_text('{"policy":"changed"}\n', encoding="utf-8")

    with pytest.raises(RuntimeError, match="artifact SHA256 mismatch"):
        audit_v2_run_manifest(manifest_path)


def test_release_manifest_rejects_missing_required_field(tmp_path: Path) -> None:
    """A release manifest cannot silently omit its environment provenance."""
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 3,
                "config_sha256": "a" * 64,
                "split_sha256": "b" * 64,
                "code_revision": "commit-marker",
                "code_sha256": "d" * 64,
                "input_sha256": {"development_sites": "c" * 64},
                "command": ["python", "script.py"],
                "device": "cpu",
                "gpu_map": [],
                "dirty": False,
                "dirty_paths": [],
                "artifacts": {},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="required fields"):
        audit_v2_run_manifest(manifest)


def test_manifest_refuses_to_replace_nonidentical_existing_output(
    tmp_path: Path,
) -> None:
    """A second invocation cannot silently replace a different release record."""
    path = tmp_path / "manifest.json"
    path.write_text('{"policy":"old"}\n', encoding="utf-8")

    with pytest.raises(RuntimeError, match="non-identical manifest"):
        write_run_manifest(
            path,
            config_sha256="a" * 64,
            split_sha256="b" * 64,
            code_revision="commit-marker",
            code_sha256="d" * 64,
            input_sha256={"development_sites": "c" * 64},
            command=["python", "script.py"],
            device="cpu",
            environment=collect_runtime_environment(),
        )


def test_release_manifest_rejects_log_without_every_roster_task(tmp_path: Path) -> None:
    """Release audit requires one complete JSONL event for every declared task."""
    checkpoint = tmp_path / "fold0.checkpoint.json"
    checkpoint.write_text('{"policy":"marker"}\n', encoding="utf-8")
    log = _training_log(tmp_path / "training.jsonl")
    config = _config(tmp_path, folds=2)
    config_sha = hashlib.sha256(config.read_bytes()).hexdigest()
    roster = _task_fingerprint(config_sha)
    missing = TaskFingerprint(
        track="strict_cluster_holdout",
        fold=1,
        seed=1,
        model="pu_logistic",
        input_sha256=roster.input_sha256,
        config_sha256=roster.config_sha256,
        code_revision=roster.code_revision,
        code_sha256=roster.code_sha256,
    )
    manifest_path = tmp_path / "manifest.json"
    write_run_manifest(
        manifest_path,
        config_sha256=config_sha,
        config_path=config,
        split_sha256="b" * 64,
        code_revision="commit-marker",
        code_sha256="d" * 64,
        input_sha256={"development_sites": "c" * 64},
        command=["python", "script.py"],
        device="cpu",
        environment=collect_runtime_environment(),
        artifacts={"training_log": log, "checkpoint": checkpoint},
        task_roster=(roster, missing),
    )

    with pytest.raises(RuntimeError, match="training log task coverage"):
        audit_v2_run_manifest(manifest_path)


def test_release_manifest_rejects_semantically_invalid_training_event(
    tmp_path: Path,
) -> None:
    """Hash-consistent JSONL still fails when event values violate the schema."""
    checkpoint = tmp_path / "fold0.checkpoint.json"
    checkpoint.write_text('{"policy":"marker"}\n', encoding="utf-8")
    log = _training_log(tmp_path / "training.jsonl")
    event = json.loads(log.read_text(encoding="utf-8"))
    event["wall_seconds"] = -1.0
    log.write_text(json.dumps(event) + "\n", encoding="utf-8")
    manifest_path = tmp_path / "manifest.json"
    config = _config(tmp_path)
    config_sha = hashlib.sha256(config.read_bytes()).hexdigest()
    write_run_manifest(
        manifest_path,
        config_sha256=config_sha,
        config_path=config,
        split_sha256="b" * 64,
        code_revision="commit-marker",
        code_sha256="d" * 64,
        input_sha256={"development_sites": "c" * 64},
        command=["python", "script.py"],
        device="cpu",
        environment=collect_runtime_environment(),
        artifacts={"training_log": log, "checkpoint": checkpoint},
        task_roster=(_task_fingerprint(config_sha),),
    )

    with pytest.raises(RuntimeError, match="training event fields are invalid"):
        audit_v2_run_manifest(manifest_path)


def test_release_manifest_requires_checkpoint_for_every_roster_task(
    tmp_path: Path,
) -> None:
    """A log plus a self-declared roster cannot substitute for task state."""
    log = _training_log(tmp_path / "training.jsonl")
    manifest_path = tmp_path / "manifest.json"
    config = _config(tmp_path)
    config_sha = hashlib.sha256(config.read_bytes()).hexdigest()
    write_run_manifest(
        manifest_path,
        config_sha256=config_sha,
        config_path=config,
        split_sha256="b" * 64,
        code_revision="commit-marker",
        code_sha256="d" * 64,
        input_sha256={"development_sites": "c" * 64},
        command=["python", "script.py"],
        device="cpu",
        environment=collect_runtime_environment(),
        artifacts={"training_log": log},
        task_roster=(_task_fingerprint(config_sha),),
    )

    with pytest.raises(RuntimeError, match="checkpoint coverage"):
        audit_v2_run_manifest(manifest_path)


def test_release_manifest_derives_strict_roster_from_configuration(
    tmp_path: Path,
) -> None:
    """Omitting a configured fold must fail even if roster/log/checkpoint agree."""
    config = _config(tmp_path, folds=2)
    config_sha = hashlib.sha256(config.read_bytes()).hexdigest()
    roster = _task_fingerprint(config_sha)
    checkpoint = tmp_path / "fold0.json"
    from plantpersulf.workflows.multispecies_v2 import write_task_checkpoint

    write_task_checkpoint(checkpoint, roster, {"validation_ap": 0.5})
    log = _training_log(tmp_path / "training.jsonl")
    manifest_path = tmp_path / "manifest.json"
    write_run_manifest(
        manifest_path,
        config_sha256=config_sha,
        config_path=config,
        split_sha256="b" * 64,
        code_revision="commit-marker",
        code_sha256="d" * 64,
        input_sha256={"development_sites": "c" * 64},
        command=["python", "script.py", "--prepare-development"],
        device="cpu",
        environment=collect_runtime_environment(),
        artifacts={"training_log": log},
        checkpoint_paths={"fold0": checkpoint},
        task_roster=(roster,),
    )

    with pytest.raises(RuntimeError, match="configured task roster"):
        audit_v2_run_manifest(manifest_path)


def _frozen_test_fingerprint(config_sha256: str, model: str) -> TaskFingerprint:
    return TaskFingerprint(
        track=FROZEN_TEST_TRACK,
        fold=0,
        seed=FROZEN_TEST_SEED,
        model=model,
        input_sha256={"development_sites": "c" * 64},
        config_sha256=config_sha256,
        code_revision="commit-marker",
        code_sha256="d" * 64,
    )


def _write_frozen_test_manifest(
    tmp_path: Path, roster: list[TaskFingerprint]
) -> Path:
    config = _config(tmp_path)
    config_sha = hashlib.sha256(config.read_bytes()).hexdigest()
    log = tmp_path / "training.jsonl"
    artifacts: dict[str, Path] = {"training_log": log}
    checkpoint_paths: dict[str, Path] = {}
    for fingerprint in roster:
        append_training_event(
            log,
            TrainingEvent(
                timestamp="2026-08-13T00:00:00+00:00",
                track=FROZEN_TEST_TRACK,
                fold=0,
                seed=FROZEN_TEST_SEED,
                model=fingerprint.model,
                epoch=0,
                loss=None,
                val_ap=None,
                lr=None,
                device="cpu",
                gpu_memory=None,
                wall_seconds=1.0,
            ),
        )
        checkpoint = tmp_path / f"{fingerprint.model}.checkpoint.json"
        write_task_checkpoint(checkpoint, fingerprint, {"arm": fingerprint.model})
        checkpoint_paths[fingerprint.model] = checkpoint
    manifest_path = tmp_path / "manifest.json"
    write_run_manifest(
        manifest_path,
        config_sha256=config_sha,
        config_path=config,
        split_sha256="b" * 64,
        code_revision="commit-marker",
        code_sha256="d" * 64,
        input_sha256={"development_sites": "c" * 64},
        command=["python", "script.py", "--score-test"],
        device="cpu",
        environment=collect_runtime_environment(),
        artifacts=artifacts,
        checkpoint_paths=checkpoint_paths,
        task_roster=tuple(roster),
    )
    return manifest_path


def test_frozen_test_roster_passes_with_both_frozen_arms(tmp_path: Path) -> None:
    """The one-shot frozen-test track is audited against its pinned roster."""
    config = _config(tmp_path)
    config_sha = hashlib.sha256(config.read_bytes()).hexdigest()
    roster = [
        _frozen_test_fingerprint(config_sha, model) for model in FROZEN_TEST_MODELS
    ]
    manifest_path = _write_frozen_test_manifest(tmp_path, roster)

    audit_v2_run_manifest(manifest_path)


def test_frozen_test_roster_rejects_missing_baseline_arm(tmp_path: Path) -> None:
    """Dropping the baseline arm from the one-shot roster must fail the audit."""
    config = _config(tmp_path)
    config_sha = hashlib.sha256(config.read_bytes()).hexdigest()
    roster = [_frozen_test_fingerprint(config_sha, "structure_ranker")]
    manifest_path = _write_frozen_test_manifest(tmp_path, roster)

    with pytest.raises(RuntimeError, match="configured task roster"):
        audit_v2_run_manifest(manifest_path)
