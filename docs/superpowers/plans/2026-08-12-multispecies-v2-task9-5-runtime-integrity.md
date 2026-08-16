# Task 9.5 Runtime Integrity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add hash-bound task logging, atomic checkpointing, safe resume, GPU assignment and release-grade execution manifests to multispecies v2.

**Architecture:** Extend `plantpersulf.workflows.multispecies_v2` with small, pure runtime-contract types and file writers. The v2 CLI creates one contract per `(track, fold, seed, model)` and persists only development execution state. Manifest validation is a separate release boundary that re-hashes every declared file.

**Tech Stack:** Python 3.10, stdlib `dataclasses/json/hashlib/tempfile/subprocess`, PyYAML, pytest, Ruff, mypy.

## Global Constraints

- Task scope is Task 9.5 only; do not start Task 9.6 statistics or change frozen splits.
- All scientific inputs remain registered and SHA256 verified; tests use only temporary policy markers.
- A checkpoint may resume only when input, config and code SHA256 are exactly equal.
- OOM recovery may only lower `score_batch_size`, with both sizes recorded.
- Each JSONL event has exactly: `timestamp, track, fold, seed, model, epoch, loss, val_ap, lr, device, gpu_memory, wall_seconds`.
- Every write is atomic; existing non-identical outputs fail closed.

---

### Task 1: Task fingerprints and checkpoint resume gate

**Files:**
- Modify: `src/plantpersulf/workflows/multispecies_v2.py`
- Test: `tests/unit/test_multispecies_v2_runtime.py`

**Interfaces:**
- Produces `TaskFingerprint(track: str, fold: int, seed: int, model: str, input_sha256: dict[str, str], config_sha256: str, code_revision: str)`.
- Produces `write_task_checkpoint(path: Path, fingerprint: TaskFingerprint, state: dict[str, object]) -> dict[str, object]`.
- Produces `load_resumable_checkpoint(path: Path, fingerprint: TaskFingerprint) -> dict[str, object]`.

- [ ] **Step 1: Write the failing atomic checkpoint test**

```python
def test_checkpoint_resumes_only_with_identical_fingerprint(tmp_path: Path) -> None:
    fingerprint = TaskFingerprint("strict_cluster_holdout", 0, 0, "pu_logistic", {"sites": "a" * 64}, "b" * 64, "commit")
    path = tmp_path / "fold0.checkpoint.json"
    write_task_checkpoint(path, fingerprint, {"epoch": 3})
    assert load_resumable_checkpoint(path, fingerprint) == {"epoch": 3}
    mismatched = replace(fingerprint, config_sha256="c" * 64)
    with pytest.raises(RuntimeError, match="checkpoint fingerprint mismatch"):
        load_resumable_checkpoint(path, mismatched)
```

- [ ] **Step 2: Verify RED**

Run: `pytest tests/unit/test_multispecies_v2_runtime.py::test_checkpoint_resumes_only_with_identical_fingerprint -q`

Expected: FAIL because `TaskFingerprint` and checkpoint functions are absent.

- [ ] **Step 3: Implement the minimal checkpoint contract**

```python
@dataclass(frozen=True)
class TaskFingerprint:
    track: str
    fold: int
    seed: int
    model: str
    input_sha256: dict[str, str]
    config_sha256: str
    code_revision: str

def write_task_checkpoint(path: Path, fingerprint: TaskFingerprint, state: dict[str, object]) -> dict[str, object]:
    payload = {"fingerprint": asdict(fingerprint), "state": state}
    return _atomic_json_write_new_or_identical(path, payload)
```

Validate SHA256 fields before writing; serialize sorted JSON; re-read and validate the same payload in `load_resumable_checkpoint`.

- [ ] **Step 4: Verify GREEN**

Run: `pytest tests/unit/test_multispecies_v2_runtime.py::test_checkpoint_resumes_only_with_identical_fingerprint -q`

Expected: PASS.

- [ ] **Step 5: Add tamper and incomplete-payload tests**

Write tests that remove `code_revision` or alter persisted `state_sha256`, then require `RuntimeError`. Run the test file.

- [ ] **Step 6: Commit**

```powershell
git add src/plantpersulf/workflows/multispecies_v2.py tests/unit/test_multispecies_v2_runtime.py
git commit -m "feat: add hash-bound v2 task checkpoints"
```

### Task 2: Event schema, deterministic GPU allocation, and OOM deviation

**Files:**
- Modify: `src/plantpersulf/workflows/multispecies_v2.py`
- Modify: `tests/unit/test_multispecies_v2_logging.py`
- Test: `tests/unit/test_multispecies_v2_runtime.py`

**Interfaces:**
- Produces `assign_task_device(task_index: int, gpu_map: tuple[str, ...], default_device: str) -> str`.
- Produces `record_oom_batch_deviation(original_batch_size: int, replacement_batch_size: int) -> dict[str, int]`.
- `append_training_event` rejects incomplete, invalid or non-finite values.

- [ ] **Step 1: Write the failing deterministic device test**

```python
def test_task_device_assignment_round_robins_configured_gpus() -> None:
    assert assign_task_device(0, ("cuda:0", "cuda:1"), "cpu") == "cuda:0"
    assert assign_task_device(3, ("cuda:0", "cuda:1"), "cpu") == "cuda:1"
```

- [ ] **Step 2: Verify RED**

Run: `pytest tests/unit/test_multispecies_v2_runtime.py::test_task_device_assignment_round_robins_configured_gpus -q`

Expected: FAIL because the device allocator is absent.

- [ ] **Step 3: Implement minimum task-level allocation**

```python
def assign_task_device(task_index: int, gpu_map: tuple[str, ...], default_device: str) -> str:
    if task_index < 0 or not default_device:
        raise ValueError("task index and default device are required")
    return gpu_map[task_index % len(gpu_map)] if gpu_map else default_device
```

Reject empty GPU names and non-`cuda:<integer>` entries. Do not create a DDP process group.

- [ ] **Step 4: Verify GREEN**

Run: `pytest tests/unit/test_multispecies_v2_runtime.py::test_task_device_assignment_round_robins_configured_gpus -q`

Expected: PASS.

- [ ] **Step 5: Write and verify OOM/event tests**

```python
def test_oom_deviation_only_allows_smaller_batch() -> None:
    assert record_oom_batch_deviation(128, 64) == {"original_batch_size": 128, "replacement_batch_size": 64}
    with pytest.raises(ValueError, match="smaller"):
        record_oom_batch_deviation(128, 128)
```

Extend the event test so `float("nan")` in `loss` raises `ValueError`. Run both tests RED, implement, then rerun GREEN.

- [ ] **Step 6: Commit**

```powershell
git add src/plantpersulf/workflows/multispecies_v2.py tests/unit/test_multispecies_v2_logging.py tests/unit/test_multispecies_v2_runtime.py
git commit -m "feat: add v2 runtime event and device contracts"
```

### Task 3: Release-grade manifest and verification

**Files:**
- Modify: `src/plantpersulf/workflows/multispecies_v2.py`
- Modify: `tests/unit/test_multispecies_v2_logging.py`
- Create: `tests/release/test_multispecies_v2_runtime_manifest.py`

**Interfaces:**
- Extend `write_run_manifest(path: Path, *, config_sha256: str, split_sha256: str, code_revision: str, input_sha256: dict[str, str], command: list[str], device: str, environment: dict[str, str], gpu_map: tuple[str, ...], artifacts: dict[str, Path], checkpoint_paths: dict[str, Path]) -> dict[str, object]`.
- Produce `audit_v2_run_manifest(path: Path) -> dict[str, object]`.

- [ ] **Step 1: Write the failing release test**

```python
def test_release_manifest_rejects_tampered_checkpoint_hash(tmp_path: Path) -> None:
    checkpoint = tmp_path / "fold0.checkpoint.json"
    checkpoint.write_text('{"policy":"marker"}\n', encoding="utf-8")
    manifest_path = tmp_path / "manifest.json"
    write_run_manifest(
        manifest_path,
        config_sha256="a" * 64,
        split_sha256="b" * 64,
        code_revision="commit",
        input_sha256={"sites": "c" * 64},
        command=["python", "script.py"],
        device="cpu",
        environment={"python": "3.10"},
        gpu_map=(),
        artifacts={"checkpoint": checkpoint},
        checkpoint_paths={"fold0": checkpoint},
    )
    checkpoint.write_text('{"policy":"changed"}\n', encoding="utf-8")
    with pytest.raises(RuntimeError, match="artifact SHA256 mismatch"):
        audit_v2_run_manifest(manifest_path)
```

- [ ] **Step 2: Verify RED**

Run: `pytest tests/release/test_multispecies_v2_runtime_manifest.py::test_release_manifest_rejects_tampered_checkpoint_hash -q`

Expected: FAIL because the audit function or extended schema is absent.

- [ ] **Step 3: Implement a stable manifest schema**

Extend the existing atomic writer to record: command, environment, dependency versions, `device`, `gpu_map`, `code_revision`, `dirty`, `dirty_paths`, sorted input hashes, and every artifact path/SHA256. `audit_v2_run_manifest` must parse exact schema, rehash each path, and reject missing required fields.

- [ ] **Step 4: Verify GREEN**

Run: `pytest tests/release/test_multispecies_v2_runtime_manifest.py -q`

Expected: PASS for complete manifest and tamper/missing-field failures.

- [ ] **Step 5: Commit**

```powershell
git add src/plantpersulf/workflows/multispecies_v2.py tests/unit/test_multispecies_v2_logging.py tests/release/test_multispecies_v2_runtime_manifest.py
git commit -m "feat: audit multispecies v2 runtime manifests"
```

### Task 4: Integrate runtime contract into the v2 development entrypoint

**Files:**
- Modify: `scripts/train_multispecies_v2.py`
- Modify: `configs/experiments/multispecies_v2_global_clusters_v5.yaml`
- Modify: `tests/scientific/test_train_multispecies_v2_cli.py`

**Interfaces:**
- `--prepare-development` writes events/checkpoints/manifest under the configured output directory.
- `--resume` is accepted only for exact task fingerprints.

- [ ] **Step 1: Write the failing CLI behavior test**

```python
def test_prepare_development_writes_hash_bound_runtime_artifacts(tmp_path: Path, monkeypatch) -> None:
    # Use existing real-code fold preparation with policy-marker filesystem inputs.
    cli.main(["--config", str(config), "--prepare-development"])
    assert (output_dir / "training.jsonl").is_file()
    assert (output_dir / "manifest.json").is_file()
```

The fixture must assert manifest audit succeeds and that no test labels are supplied to the development training function.

- [ ] **Step 2: Verify RED**

Run: `pytest tests/scientific/test_train_multispecies_v2_cli.py::test_prepare_development_writes_hash_bound_runtime_artifacts -q`

Expected: FAIL because the current CLI does not persist runtime artifacts.

- [ ] **Step 3: Implement minimal integration**

Read `compute.device`, `compute.score_batch_size`, `compute.max_parallel_tasks`, and `compute.gpu_map` from v5. For each of five folds, create a fingerprint before training; load a matching checkpoint only when `--resume` is set; append a complete event; atomically checkpoint fold result; write and audit the manifest after all folds. The CLI must never invoke the test path.

- [ ] **Step 4: Verify GREEN**

Run: `pytest tests/scientific/test_train_multispecies_v2_cli.py -q`

Expected: PASS.

- [ ] **Step 5: Run Task 9.5 completion gate**

```powershell
pytest tests/unit -q
pytest tests/scientific -q
pytest tests/release -q
python -m plantpersulf.cli audit-registry
python -m plantpersulf.cli audit-files
python -m plantpersulf.cli audit-leakage --split-version multispecies_strict_v4
ruff check .
mypy src/plantpersulf
```

- [ ] **Step 6: Commit**

```powershell
git add scripts/train_multispecies_v2.py configs/experiments/multispecies_v2_global_clusters_v5.yaml tests/scientific/test_train_multispecies_v2_cli.py
git commit -m "feat: integrate v2 runtime integrity controls"
```

## Self-review

The plan covers task-level GPU scheduling, JSONL logging, atomic hash-bound checkpoints, exact resume, OOM-only batch-size deviations, manifest provenance and release rejection. It does not modify scientific labels, frozen splits, evaluation claims or Task 9.6 statistics. Each production behavior has a named RED test and a direct GREEN command.
