"""Fail-closed runtime contracts for multispecies v2 task execution."""

from __future__ import annotations

import time
from dataclasses import replace
from pathlib import Path
from threading import get_ident

import pytest

from plantpersulf.workflows import multispecies_v2 as runtime


def test_task_fingerprint_hashes_physical_inputs_and_dirty_code(tmp_path: Path) -> None:
    """Resume must reject either changed input bytes or changed dirty code bytes."""
    input_path = tmp_path / "development.tsv"
    dirty_path = tmp_path / "runtime.py"
    input_path.write_text("registered-input-marker\n", encoding="utf-8")
    dirty_path.write_text("dirty-code-marker\n", encoding="utf-8")

    fingerprint = runtime.build_task_fingerprint(
        track="strict_cluster_holdout",
        fold=0,
        seed=0,
        model="pu_logistic",
        input_paths={"development_rows": input_path},
        config_path=input_path,
        code_revision="commit-marker",
        dirty_paths=(dirty_path,),
    )

    assert fingerprint.input_sha256["development_rows"] == runtime._sha256_file(
        input_path
    )
    assert fingerprint.config_sha256 == runtime._sha256_file(input_path)
    assert fingerprint.code_sha256 == runtime._sha256_file(dirty_path)


def test_runtime_environment_records_dependencies_and_gpu_probe() -> None:
    """Release provenance needs executable, dependency, and GPU identities."""
    environment = runtime.collect_runtime_environment()

    assert environment["python_version"]
    assert environment["platform"]
    assert "PyYAML" in environment["dependencies"]
    assert "gpu_probe" in environment


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
        code_sha256="c" * 64,
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
        code_sha256="c" * 64,
    )
    checkpoint = tmp_path / "fold1.checkpoint.json"
    runtime.write_task_checkpoint(checkpoint, fingerprint, {"epoch": 3})
    payload = checkpoint.read_text(encoding="utf-8").replace('"epoch":3', '"epoch":4')
    checkpoint.write_text(payload, encoding="utf-8")

    with pytest.raises(RuntimeError, match="checkpoint state hash mismatch"):
        runtime.load_resumable_checkpoint(checkpoint, fingerprint)


def test_task_device_assignment_round_robins_configured_gpus() -> None:
    """Independent fold tasks must be mapped deterministically, not by DDP."""
    assert runtime.assign_task_device(0, ("cuda:0", "cuda:1"), "cpu") == "cuda:0"
    assert runtime.assign_task_device(3, ("cuda:0", "cuda:1"), "cpu") == "cuda:1"
    assert runtime.assign_task_device(2, (), "cpu") == "cpu"


def test_oom_deviation_only_allows_smaller_batch() -> None:
    """OOM recovery must not alter a model or increase the evaluation batch."""
    assert runtime.record_oom_batch_deviation(128, 64) == {
        "original_batch_size": 128,
        "replacement_batch_size": 64,
    }
    with pytest.raises(ValueError, match="smaller"):
        runtime.record_oom_batch_deviation(128, 128)


def test_parallel_task_group_uses_multiple_workers_when_configured() -> None:
    """Independent fold jobs must really dispatch concurrently, not just map devices."""
    observed: set[int] = set()

    def runner(task: str, device: str) -> str:
        observed.add(get_ident())
        time.sleep(0.05)
        return f"{task}:{device}"

    results = runtime.run_task_group(
        ("fold0", "fold1"),
        max_parallel_tasks=2,
        gpu_map=("cuda:0", "cuda:1"),
        default_device="cpu",
        runner=runner,
    )

    assert results == ["fold0:cuda:0", "fold1:cuda:1"]
    assert len(observed) == 2


def test_oom_retry_reduces_only_batch_size_and_records_deviation() -> None:
    """A real OOM path retries the same operation with a smaller score batch."""
    attempted: list[int] = []

    def score(batch_size: int) -> int:
        attempted.append(batch_size)
        if batch_size == 8:
            raise RuntimeError("CUDA out of memory")
        return batch_size

    result, effective_batch, deviations = runtime.run_with_oom_batch_retry(
        score, initial_batch_size=8
    )

    assert result == 4
    assert effective_batch == 4
    assert attempted == [8, 4]
    assert deviations == [{"original_batch_size": 8, "replacement_batch_size": 4}]
