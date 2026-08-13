"""Leakage-safe split helpers for the multispecies v2 experiment runner."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import math
import platform
import random
import re
import subprocess
import sys
import tempfile
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from functools import cache
from pathlib import Path
from typing import TypeVar

import yaml

from plantpersulf.benchmark.multispecies_splits import (
    DEVELOPMENT_SPLIT,
    TEST_SPLIT,
    FrozenMultispeciesSplit,
    MultispeciesSiteRow,
    load_frozen_multispecies_split,
    load_global_cluster_table,
)
from plantpersulf.models.traditional import pu_logistic_regression_scores
from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow
from plantpersulf.provenance.audit import assert_registered_input


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def _sha256_path(path: Path) -> str:
    """Hash a physical file or a directory tree with stable relative framing."""
    resolved = path.resolve()
    if resolved.is_file():
        return _sha256_file(resolved)
    if not resolved.is_dir():
        raise FileNotFoundError(resolved)
    digest = hashlib.sha256()
    for child in sorted(item for item in resolved.rglob("*") if item.is_file()):
        digest.update(child.relative_to(resolved).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(_sha256_file(child)))
    return digest.hexdigest()


@dataclass(frozen=True)
class TrainingEvent:
    timestamp: str
    track: str
    fold: int
    seed: int
    model: str
    epoch: int
    loss: float | None
    val_ap: float | None
    lr: float | None
    device: str
    gpu_memory: int | None
    wall_seconds: float


@dataclass(frozen=True)
class TaskFingerprint:
    """Immutable identity required before a task checkpoint may be resumed."""

    track: str
    fold: int
    seed: int
    model: str
    input_sha256: dict[str, str]
    config_sha256: str
    code_revision: str
    code_sha256: str


def _canonical_json_bytes(payload: object) -> bytes:
    try:
        serialized = json.dumps(payload, sort_keys=True, separators=(",", ":"))
        return f"{serialized}\n".encode()
    except (TypeError, ValueError) as exc:
        raise ValueError("runtime payload must be JSON serializable") from exc


def _validate_task_fingerprint(fingerprint: TaskFingerprint) -> None:
    if (
        not fingerprint.track
        or fingerprint.fold < 0
        or fingerprint.seed < 0
        or not fingerprint.model
        or not fingerprint.code_revision
        or not _valid_sha256(fingerprint.code_sha256)
        or not _valid_sha256(fingerprint.config_sha256)
        or not fingerprint.input_sha256
        or any(
            not name or not _valid_sha256(value)
            for name, value in fingerprint.input_sha256.items()
        )
    ):
        raise ValueError("task fingerprint is invalid")


def _hash_existing_paths(paths: tuple[Path, ...]) -> str:
    """Hash one dirty file directly or many files with stable path framing."""
    if not paths:
        return hashlib.sha256(b"").hexdigest()
    if len(paths) == 1:
        if not paths[0].is_file():
            raise RuntimeError(f"dirty code path is not a file: {paths[0]}")
        return _sha256_file(paths[0])
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: item.as_posix()):
        if not path.is_file():
            raise RuntimeError(f"dirty code path is not a file: {path}")
        digest.update(path.as_posix().encode())
        digest.update(b"\0")
        digest.update(bytes.fromhex(_sha256_file(path)))
    return digest.hexdigest()


def build_task_fingerprint(
    *,
    track: str,
    fold: int,
    seed: int,
    model: str,
    input_paths: dict[str, Path],
    config_path: Path,
    code_revision: str,
    dirty_paths: tuple[Path, ...],
) -> TaskFingerprint:
    """Bind a runtime task to the exact bytes consumed by that process."""
    if not input_paths or any(
        not name or not path.exists() for name, path in input_paths.items()
    ):
        raise ValueError("task inputs must name existing files")
    if not config_path.is_file():
        raise ValueError("task config must be an existing file")
    fingerprint = TaskFingerprint(
        track=track,
        fold=fold,
        seed=seed,
        model=model,
        input_sha256={
            name: _sha256_path(path) for name, path in sorted(input_paths.items())
        },
        config_sha256=_sha256_file(config_path),
        code_revision=code_revision,
        code_sha256=_hash_existing_paths(dirty_paths),
    )
    _validate_task_fingerprint(fingerprint)
    return fingerprint


def write_task_checkpoint(
    path: Path, fingerprint: TaskFingerprint, state: dict[str, object]
) -> dict[str, object]:
    """Atomically persist state only when the target is new or identical."""
    _validate_task_fingerprint(fingerprint)
    state_sha256 = hashlib.sha256(_canonical_json_bytes(state)).hexdigest()
    payload: dict[str, object] = {
        "fingerprint": asdict(fingerprint),
        "state": state,
        "state_sha256": state_sha256,
    }
    encoded = _canonical_json_bytes(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != encoded:
            raise RuntimeError(
                f"refusing to overwrite non-identical checkpoint: {path}"
            )
        return state
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(encoded)
    temporary.replace(path)
    return state


def load_resumable_checkpoint(
    path: Path, fingerprint: TaskFingerprint
) -> dict[str, object]:
    """Return state only after exact fingerprint and state-hash verification."""
    _validate_task_fingerprint(fingerprint)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"checkpoint is invalid: {path}") from exc
    if not isinstance(payload, dict) or set(payload) != {
        "fingerprint",
        "state",
        "state_sha256",
    }:
        raise RuntimeError("checkpoint payload is incomplete")
    stored_fingerprint = payload["fingerprint"]
    state = payload["state"]
    state_sha256 = payload["state_sha256"]
    if stored_fingerprint != asdict(fingerprint):
        raise RuntimeError("checkpoint fingerprint mismatch")
    if (
        not isinstance(state, dict)
        or not isinstance(state_sha256, str)
        or hashlib.sha256(_canonical_json_bytes(state)).hexdigest() != state_sha256
    ):
        raise RuntimeError("checkpoint state hash mismatch")
    return state


def assign_task_device(
    task_index: int, gpu_map: tuple[str, ...], default_device: str
) -> str:
    """Map independent tasks round-robin onto configured GPUs without DDP."""
    if task_index < 0 or not default_device:
        raise ValueError("task index and default device are required")
    if any(re.fullmatch(r"cuda:\d+", device) is None for device in gpu_map):
        raise ValueError("gpu_map entries must be cuda:<integer>")
    return gpu_map[task_index % len(gpu_map)] if gpu_map else default_device


def record_oom_batch_deviation(
    original_batch_size: int, replacement_batch_size: int
) -> dict[str, int]:
    """Record the sole permitted OOM deviation: a strictly smaller batch."""
    if (
        original_batch_size <= 0
        or replacement_batch_size <= 0
        or replacement_batch_size >= original_batch_size
    ):
        raise ValueError("OOM replacement batch size must be positive and smaller")
    return {
        "original_batch_size": original_batch_size,
        "replacement_batch_size": replacement_batch_size,
    }


TaskValue = TypeVar("TaskValue")
TaskResult = TypeVar("TaskResult")


def run_task_group(
    tasks: tuple[TaskValue, ...],
    *,
    max_parallel_tasks: int,
    gpu_map: tuple[str, ...],
    default_device: str,
    runner: Callable[[TaskValue, str], TaskResult],
    on_complete: Callable[[int, TaskResult], None] | None = None,
) -> list[TaskResult]:
    """Run independent tasks concurrently and preserve their declared order."""
    if max_parallel_tasks < 1:
        raise ValueError("max_parallel_tasks must be positive")
    if gpu_map and max_parallel_tasks > len(gpu_map):
        raise ValueError(
            "max_parallel_tasks must not exceed the number of configured GPUs"
        )
    assignments = [
        (task, assign_task_device(index, gpu_map, default_device))
        for index, task in enumerate(tasks)
    ]
    with ThreadPoolExecutor(max_workers=max_parallel_tasks) as executor:
        futures = {
            executor.submit(runner, task, device): index
            for index, (task, device) in enumerate(assignments)
        }
        results: dict[int, TaskResult] = {}
        for future in as_completed(futures):
            index = futures[future]
            result = future.result()
            results[index] = result
            if on_complete is not None:
                on_complete(index, result)
    return [results[index] for index in range(len(tasks))]


def run_with_oom_batch_retry(
    operation: Callable[[int], TaskResult], *, initial_batch_size: int
) -> tuple[TaskResult, int, list[dict[str, int]]]:
    """Retry an unchanged scoring operation after CUDA OOM by halving its batch."""
    if initial_batch_size <= 0:
        raise ValueError("positive initial batch size is required")
    batch_size = initial_batch_size
    deviations: list[dict[str, int]] = []
    while True:
        try:
            return operation(batch_size), batch_size, deviations
        except RuntimeError as exc:
            if "out of memory" not in str(exc).lower() or batch_size == 1:
                raise
            replacement = max(1, batch_size // 2)
            deviations.append(record_oom_batch_deviation(batch_size, replacement))
            batch_size = replacement


def collect_runtime_environment(
    python_executable: Path | None = None,
) -> dict[str, object]:
    """Capture exact Python packages and current GPU inventory for provenance."""
    executable = python_executable or Path(sys.executable)
    if (
        python_executable is None
        or executable.resolve() == Path(sys.executable).resolve()
    ):
        dependencies: dict[str, str] = {}
        for distribution in importlib.metadata.distributions():
            try:
                name = distribution.metadata["Name"]
            except KeyError:
                continue
            if name:
                dependencies[name] = distribution.version
        python_version = platform.python_version()
        platform_name = platform.platform()
    else:
        probe_code = (
            "import importlib.metadata,json,platform;"
            "print(json.dumps({'python_version':platform.python_version(),"
            "'platform':platform.platform(),'dependencies':"
            "{d.metadata['Name']:d.version for d in importlib.metadata.distributions() "
            "if d.metadata['Name']}}))"
        )
        try:
            isolated = subprocess.run(
                [str(executable), "-c", probe_code],
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            )
            isolated_payload = json.loads(isolated.stdout)
        except (OSError, subprocess.SubprocessError, json.JSONDecodeError) as exc:
            raise RuntimeError("isolated runtime environment probe failed") from exc
        dependencies = isolated_payload.get("dependencies")
        python_version = isolated_payload.get("python_version")
        platform_name = isolated_payload.get("platform")
        if (
            not isinstance(dependencies, dict)
            or not isinstance(python_version, str)
            or not isinstance(platform_name, str)
        ):
            raise RuntimeError("isolated runtime environment probe is invalid")
    try:
        probe = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index,name,driver_version,memory.total",
                "--format=csv,noheader",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=10,
        )
        gpu_probe: object = [line.strip() for line in probe.stdout.splitlines() if line]
    except (FileNotFoundError, subprocess.SubprocessError):
        gpu_probe = []
    return {
        "python_executable": str(executable.resolve()),
        "python_version": python_version,
        "platform": platform_name,
        "dependencies": dict(sorted(dependencies.items())),
        "gpu_probe": gpu_probe,
    }


@dataclass(frozen=True)
class MultispeciesExperimentPreparation:
    config_path: Path
    config_sha256: str
    split: FrozenMultispeciesSplit
    cluster_count: int
    test_scoring_enabled: bool


@dataclass(frozen=True)
class V2DevelopmentPipeline:
    """Preprocessing state fitted exclusively from one development fit fold."""

    imputation_values: tuple[float, ...]
    selected_feature_indices: tuple[int, ...]
    feature_mean: tuple[float, ...]
    feature_std: tuple[float, ...]
    pu_prior: float

    def transform(
        self,
        rows: Iterable[MultispeciesV2SiteRow],
        features: dict[tuple[str, int], tuple[float | None, ...]],
    ) -> list[list[float]]:
        transformed: list[list[float]] = []
        for row in rows:
            raw = features[(row.global_protein_id, row.cys_position)]
            filled = [
                value if value is not None else self.imputation_values[index]
                for index, value in enumerate(raw)
            ]
            selected = [filled[index] for index in self.selected_feature_indices]
            transformed.append(
                [
                    (value - mean) / std
                    for value, mean, std in zip(
                        selected, self.feature_mean, self.feature_std, strict=True
                    )
                ]
            )
        return transformed


@dataclass(frozen=True)
class V2DevelopmentFoldResult:
    """Development-only model-selection result; frozen test is absent by type."""

    pipeline: V2DevelopmentPipeline
    validation_ap: float


def fit_v2_development_pipeline(
    fit_rows: Iterable[MultispeciesV2SiteRow],
    features: dict[tuple[str, int], tuple[float | None, ...]],
) -> V2DevelopmentPipeline:
    """Fit imputation, feature selection, scaling, and PU prior on fit rows."""
    rows = tuple(fit_rows)
    if not rows:
        raise ValueError("fit_rows must not be empty")
    raw_rows = [features[(row.global_protein_id, row.cys_position)] for row in rows]
    width = len(raw_rows[0])
    if not width or any(len(values) != width for values in raw_rows):
        raise ValueError("feature vectors must have one non-zero shared width")
    imputation = tuple(
        sum(value for value in column if value is not None)
        / max(sum(value is not None for value in column), 1)
        for column in zip(*raw_rows, strict=True)
    )
    filled = [
        [
            value if value is not None else imputation[index]
            for index, value in enumerate(raw)
        ]
        for raw in raw_rows
    ]
    selected = tuple(
        index
        for index in range(width)
        if any(row[index] != filled[0][index] for row in filled[1:])
    )
    if not selected:
        raise RuntimeError("fit rows contain no selectable feature")
    selected_rows = [[row[index] for index in selected] for row in filled]
    mean = tuple(
        sum(row[index] for row in selected_rows) / len(selected_rows)
        for index in range(len(selected))
    )
    std = tuple(
        (
            calculated
            if (
                calculated := (
                    sum((row[index] - mean[index]) ** 2 for row in selected_rows)
                    / len(selected_rows)
                )
                ** 0.5
            )
            > 0.0
            else 1.0
        )
        for index in range(len(selected))
    )
    return V2DevelopmentPipeline(
        imputation_values=imputation,
        selected_feature_indices=selected,
        feature_mean=mean,
        feature_std=std,
        pu_prior=sum(row.label == "positive" for row in rows) / len(rows),
    )


def train_v2_development_fold(
    fit_rows: Iterable[MultispeciesV2SiteRow],
    validation_rows: Iterable[MultispeciesV2SiteRow],
    features: dict[tuple[str, int], tuple[float | None, ...]],
    *,
    seed: int,
    holdout_fraction: float = 0.2,
) -> V2DevelopmentFoldResult:
    """Fit PU logistic and score one validation fold without a test argument."""
    fit = tuple(fit_rows)
    validation = tuple(validation_rows)
    if not validation or not any(row.label == "positive" for row in validation):
        raise ValueError("validation rows require at least one positive")
    pipeline = fit_v2_development_pipeline(fit, features)
    fit_features = pipeline.transform(fit, features)
    validation_features = pipeline.transform(validation, features)
    scores = pu_logistic_regression_scores(
        fit_features,
        [row.label for row in fit],
        validation_features,
        seed=seed,
        holdout_fraction=holdout_fraction,
    )
    from sklearn.metrics import average_precision_score  # type: ignore[import-untyped]

    return V2DevelopmentFoldResult(
        pipeline=pipeline,
        validation_ap=float(
            average_precision_score(
                [row.label == "positive" for row in validation], scores
            )
        ),
    )


def select_v2_development_hyperparameters(
    fit_rows: Iterable[MultispeciesV2SiteRow],
    validation_rows: Iterable[MultispeciesV2SiteRow],
    features: dict[tuple[str, int], tuple[float | None, ...]],
    *,
    candidates: tuple[float, ...],
    seed: int,
) -> float:
    """Select PU calibration holdout fraction from development data only."""
    if not candidates or any(not 0.0 < value < 1.0 for value in candidates):
        raise ValueError("PU holdout candidates must be non-empty fractions")
    fit = tuple(fit_rows)
    validation = tuple(validation_rows)
    results = [
        (
            train_v2_development_fold(
                fit,
                validation,
                features,
                seed=seed,
                holdout_fraction=candidate,
            ).validation_ap,
            candidate,
        )
        for candidate in candidates
    ]
    return max(results)[1]


def append_training_event(path: Path, event: TrainingEvent) -> None:
    """Append one complete, machine-readable training event."""
    numeric_values = (
        event.loss,
        event.val_ap,
        event.lr,
        event.wall_seconds,
    )
    if any(value is not None and not math.isfinite(value) for value in numeric_values):
        raise ValueError("training event numeric values must be finite")
    if event.gpu_memory is not None and event.gpu_memory < 0:
        raise ValueError("training event GPU memory must be non-negative")
    if (
        not event.timestamp
        or not event.track
        or event.fold < 0
        or event.seed < 0
        or not event.model
        or event.epoch < 0
        or re.fullmatch(r"cpu|cuda(?::\d+)?", event.device) is None
        or event.wall_seconds < 0.0
    ):
        raise ValueError("training event fields are invalid")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(asdict(event), sort_keys=True) + "\n")


def _task_identity(payload: dict[str, object]) -> tuple[str, int, int, str]:
    track = payload.get("track")
    fold = payload.get("fold")
    seed = payload.get("seed")
    model = payload.get("model")
    if (
        not isinstance(track, str)
        or not isinstance(fold, int)
        or not isinstance(seed, int)
        or not isinstance(model, str)
    ):
        raise RuntimeError("task identity is invalid")
    return track, fold, seed, model


def audit_training_event_log(
    path: Path, task_roster: tuple[TaskFingerprint, ...]
) -> None:
    """Require complete JSONL schema and exact per-task release coverage."""
    expected = {
        (entry.track, entry.fold, entry.seed, entry.model) for entry in task_roster
    }
    if not expected or len(expected) != len(task_roster):
        raise RuntimeError("task roster must contain unique tasks")
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as exc:
        raise RuntimeError("training log is unavailable") from exc
    seen: set[tuple[str, int, int, str]] = set()
    required = set(TrainingEvent.__dataclass_fields__)
    for line in lines:
        try:
            event = json.loads(line)
        except json.JSONDecodeError as exc:
            raise RuntimeError("training log contains invalid JSON") from exc
        if not isinstance(event, dict) or set(event) != required:
            raise RuntimeError("training log required fields are missing")
        _audit_training_event_payload(event)
        identity = _task_identity(event)
        if identity in seen:
            raise RuntimeError("training log contains duplicate task event")
        seen.add(identity)
    if seen != expected:
        raise RuntimeError("training log task coverage is incomplete")


def _audit_training_event_payload(event: dict[str, object]) -> None:
    """Apply the producer's semantic validation to untrusted JSONL bytes."""
    nullable_floats = ("loss", "val_ap", "lr")
    if any(
        event[name] is not None
        and (not isinstance(event[name], int | float) or isinstance(event[name], bool))
        for name in nullable_floats
    ):
        raise RuntimeError("training event fields are invalid")
    gpu_memory = event["gpu_memory"]
    if gpu_memory is not None and (
        not isinstance(gpu_memory, int) or isinstance(gpu_memory, bool)
    ):
        raise RuntimeError("training event fields are invalid")
    try:
        parsed = TrainingEvent(**event)  # type: ignore[arg-type]
        if not isinstance(parsed.timestamp, str) or not isinstance(parsed.track, str):
            raise ValueError
        if not isinstance(parsed.model, str) or not isinstance(parsed.device, str):
            raise ValueError
        if any(
            not isinstance(value, int) or isinstance(value, bool)
            for value in (parsed.fold, parsed.seed, parsed.epoch)
        ):
            raise ValueError
        if not isinstance(parsed.wall_seconds, int | float) or isinstance(
            parsed.wall_seconds, bool
        ):
            raise ValueError
        if re.fullmatch(r"cpu|cuda(?::\d+)?", parsed.device) is None:
            raise ValueError
        numeric_values = (parsed.loss, parsed.val_ap, parsed.lr, parsed.wall_seconds)
        if any(
            value is not None and not math.isfinite(value) for value in numeric_values
        ):
            raise ValueError
        if (
            parsed.fold < 0
            or parsed.seed < 0
            or parsed.epoch < 0
            or parsed.wall_seconds < 0
            or (parsed.gpu_memory is not None and parsed.gpu_memory < 0)
            or not parsed.timestamp
            or not parsed.track
            or not parsed.model
        ):
            raise ValueError
    except (TypeError, ValueError) as exc:
        raise RuntimeError("training event fields are invalid") from exc


def _audit_checkpoint_coverage(
    artifacts: dict[str, object], roster: tuple[TaskFingerprint, ...]
) -> None:
    expected = {_task_identity(asdict(item)): asdict(item) for item in roster}
    observed: dict[tuple[str, int, int, str], dict[str, object]] = {}
    for record in artifacts.values():
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            continue
        artifact_path = Path(record["path"])
        try:
            payload = json.loads(artifact_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(payload, dict) or set(payload) != {
            "fingerprint",
            "state",
            "state_sha256",
        }:
            continue
        fingerprint = payload["fingerprint"]
        if not isinstance(fingerprint, dict):
            continue
        identity = _task_identity(fingerprint)
        if identity in observed:
            raise RuntimeError("runtime checkpoint coverage is duplicated")
        observed[identity] = fingerprint
    if set(observed) != set(expected) or any(
        observed[identity] != expected[identity] for identity in expected
    ):
        raise RuntimeError("runtime checkpoint coverage is incomplete")


def _valid_sha256(value: str) -> bool:
    if len(value) != 64:
        return False
    try:
        int(value, 16)
    except ValueError:
        return False
    return True


def validate_reference_proteome_inputs(cfg: dict[str, object]) -> None:
    """Verify every configured reference FASTA before any Cys enumeration."""
    entries = cfg.get("reference_proteomes")
    if not isinstance(entries, list) or not entries:
        raise RuntimeError("reference_proteomes configuration is required")
    seen_species: set[str] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            raise RuntimeError("reference proteome entry must be a mapping")
        species = entry.get("species")
        path_text = entry.get("path")
        expected = entry.get("sha256")
        if (
            not isinstance(species, str)
            or not species
            or species in seen_species
            or not isinstance(path_text, str)
            or not isinstance(expected, str)
            or not _valid_sha256(expected)
        ):
            raise RuntimeError("reference proteome entry is invalid")
        seen_species.add(species)
        path = Path(path_text)
        if not path.is_file() or _sha256_file(path).lower() != expected.lower():
            raise RuntimeError(f"reference proteome SHA256 mismatch: {species}")


def write_run_manifest(
    path: Path,
    *,
    config_sha256: str,
    config_path: Path | None = None,
    split_sha256: str,
    code_revision: str,
    code_sha256: str,
    input_sha256: dict[str, str],
    command: list[str],
    device: str,
    environment: dict[str, object] | None = None,
    gpu_map: tuple[str, ...] = (),
    artifacts: dict[str, Path] | None = None,
    checkpoint_paths: dict[str, Path] | None = None,
    task_roster: tuple[TaskFingerprint, ...] = (),
    external_environments: dict[str, dict[str, object]] | None = None,
    dirty: bool = False,
    dirty_paths: tuple[str, ...] = (),
) -> dict[str, object]:
    """Write an atomic run manifest after verifying every input hash."""
    if not all(
        _valid_sha256(value) for value in (config_sha256, split_sha256, code_sha256)
    ):
        raise ValueError("config_sha256, split_sha256, and code_sha256 are required")
    if not input_sha256 or any(
        not name or not _valid_sha256(value) for name, value in input_sha256.items()
    ):
        raise ValueError("input_sha256 must contain complete SHA256 values")
    if not code_revision or not command or not device:
        raise ValueError("code_revision, command, and device are required")
    if config_path is not None and (
        not config_path.is_file() or _sha256_file(config_path) != config_sha256
    ):
        raise ValueError("config path does not match config_sha256")
    if (
        any(
            flag in command
            for flag in ("--prepare-development", "--run-literature-baselines")
        )
        and config_path is None
    ):
        raise ValueError("production runtime manifests require config_path")
    if any(re.fullmatch(r"cuda:\d+", value) is None for value in gpu_map):
        raise ValueError("gpu_map entries must be cuda:<integer>")
    environment_payload = environment or {}
    if not {
        "python_executable",
        "python_version",
        "platform",
        "dependencies",
        "gpu_probe",
    }.issubset(environment_payload):
        raise ValueError("environment provenance is incomplete")
    if (
        not all(
            isinstance(environment_payload[name], str) and environment_payload[name]
            for name in ("python_executable", "python_version", "platform")
        )
        or not isinstance(environment_payload["dependencies"], dict)
        or not isinstance(environment_payload["gpu_probe"], list)
    ):
        raise ValueError("environment provenance is invalid")
    for fingerprint in task_roster:
        _validate_task_fingerprint(fingerprint)
    external_environment_payload = external_environments or {}
    for name, payload in external_environment_payload.items():
        if not name or not {
            "python_executable",
            "python_version",
            "platform",
            "dependencies",
            "gpu_probe",
        }.issubset(payload):
            raise ValueError("external environment provenance is incomplete")
    declared_artifacts = {**(artifacts or {}), **(checkpoint_paths or {})}
    if any(
        not name or not artifact.is_file()
        for name, artifact in declared_artifacts.items()
    ):
        raise ValueError("manifest artifacts must name existing files")
    artifact_records = {
        name: {
            "path": artifact.resolve().as_posix(),
            "sha256": _sha256_file(artifact),
        }
        for name, artifact in sorted(declared_artifacts.items())
    }
    manifest: dict[str, object] = {
        "schema_version": 4,
        "config_sha256": config_sha256,
        "config_path": config_path.resolve().as_posix() if config_path else None,
        "split_sha256": split_sha256,
        "code_revision": code_revision,
        "code_sha256": code_sha256,
        "input_sha256": dict(sorted(input_sha256.items())),
        "command": command,
        "device": device,
        "environment": dict(sorted(environment_payload.items())),
        "external_environments": external_environment_payload,
        "gpu_map": list(gpu_map),
        "dirty": dirty,
        "dirty_paths": sorted(dirty_paths),
        "artifacts": artifact_records,
        "task_roster": [asdict(entry) for entry in task_roster],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    serialized = json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    if path.exists():
        if path.read_text(encoding="utf-8") != serialized:
            raise RuntimeError(f"refusing to overwrite non-identical manifest: {path}")
        return manifest
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(serialized)
    temporary.replace(path)
    return manifest


def audit_v2_run_manifest(path: Path) -> dict[str, object]:
    """Fail closed if a v2 runtime manifest or declared artifact changed."""
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"runtime manifest is invalid: {path}") from exc
    required = {
        "schema_version",
        "config_sha256",
        "config_path",
        "split_sha256",
        "code_revision",
        "code_sha256",
        "input_sha256",
        "command",
        "device",
        "environment",
        "external_environments",
        "gpu_map",
        "dirty",
        "dirty_paths",
        "artifacts",
        "task_roster",
    }
    if not isinstance(manifest, dict) or set(manifest) != required:
        raise RuntimeError("runtime manifest required fields are missing")
    if manifest["schema_version"] != 4:
        raise RuntimeError("runtime manifest schema is invalid")
    if not all(
        isinstance(manifest.get(name), str) and _valid_sha256(manifest[name])
        for name in ("config_sha256", "split_sha256", "code_sha256")
    ):
        raise RuntimeError("runtime manifest hashes are invalid")
    artifacts = manifest["artifacts"]
    if not isinstance(artifacts, dict):
        raise RuntimeError("runtime manifest artifacts are invalid")
    for name, record in artifacts.items():
        if not isinstance(name, str) or not isinstance(record, dict):
            raise RuntimeError("runtime manifest artifacts are invalid")
        artifact_path = record.get("path")
        artifact_sha256 = record.get("sha256")
        if (
            not isinstance(artifact_path, str)
            or not isinstance(artifact_sha256, str)
            or not _valid_sha256(artifact_sha256)
            or not Path(artifact_path).is_file()
            or _sha256_file(Path(artifact_path)) != artifact_sha256
        ):
            raise RuntimeError("artifact SHA256 mismatch")
    environment = manifest["environment"]
    if not isinstance(environment, dict) or not {
        "python_executable",
        "python_version",
        "platform",
        "dependencies",
        "gpu_probe",
    }.issubset(environment):
        raise RuntimeError("runtime environment provenance is incomplete")
    if not isinstance(environment["dependencies"], dict) or not isinstance(
        environment["gpu_probe"], list
    ):
        raise RuntimeError("runtime environment provenance is invalid")
    external_environments = manifest["external_environments"]
    if not isinstance(external_environments, dict):
        raise RuntimeError("external runtime environment provenance is invalid")
    for payload in external_environments.values():
        if not isinstance(payload, dict) or not {
            "python_executable",
            "python_version",
            "platform",
            "dependencies",
            "gpu_probe",
        }.issubset(payload):
            raise RuntimeError("external runtime environment provenance is incomplete")
    roster_payload = manifest["task_roster"]
    if not isinstance(roster_payload, list):
        raise RuntimeError("runtime task roster is invalid")
    try:
        roster = tuple(TaskFingerprint(**entry) for entry in roster_payload)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("runtime task roster is invalid") from exc
    for fingerprint in roster:
        _validate_task_fingerprint(fingerprint)
    _audit_configured_task_roster(manifest, roster)
    training_log = artifacts.get("training_log")
    if not isinstance(training_log, dict) or not isinstance(
        training_log.get("path"), str
    ):
        raise RuntimeError("runtime training log is missing")
    audit_training_event_log(Path(training_log["path"]), roster)
    logged_events = [
        json.loads(line)
        for line in Path(training_log["path"]).read_text(encoding="utf-8").splitlines()
    ]
    if any(str(event["device"]).startswith("cuda") for event in logged_events):
        probes = [environment.get("gpu_probe")]
        probes.extend(
            payload.get("gpu_probe")
            for payload in external_environments.values()
            if isinstance(payload, dict)
        )
        if not any(isinstance(probe, list) and probe for probe in probes):
            raise RuntimeError("CUDA task lacks GPU provenance")
    _audit_checkpoint_coverage(artifacts, roster)
    return manifest


_COMPARISON_INPUT_REGISTRIES = (
    Path("data/registry/cross_crop_target_label_free_inputs_v1.tsv"),
    Path("data/registry/supplementary_sources.tsv"),
    Path("data/registry/model_inputs.tsv"),
)


def _is_registered_comparison_input(path: Path) -> bool:
    """Mirror the runtime's registered-input gate so the audit reflects the
    same ready/blocked decision the run itself made, not config presence."""
    if not path.is_file():
        return False
    for registry in _COMPARISON_INPUT_REGISTRIES:
        try:
            assert_registered_input(path, registry)
            return True
        except RuntimeError:
            continue
    return False


def _audit_configured_task_roster(
    manifest: dict[str, object], roster: tuple[TaskFingerprint, ...]
) -> None:
    """Derive production task identities from the immutable experiment config."""
    command = manifest.get("command")
    if not isinstance(command, list) or not all(
        isinstance(item, str) for item in command
    ):
        raise RuntimeError("runtime command provenance is invalid")
    config_text = manifest.get("config_path")
    if not isinstance(config_text, str):
        raise RuntimeError("runtime configuration provenance is missing")
    config_path = Path(config_text)
    if not config_path.is_file() or _sha256_file(config_path) != manifest.get(
        "config_sha256"
    ):
        raise RuntimeError("runtime configuration SHA256 mismatch")
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise RuntimeError("runtime configuration is invalid")

    tracks = {item.track for item in roster}
    if not tracks:
        return
    expected: set[tuple[str, int, int, str]] = set()
    if "strict_cluster_holdout" in tracks:
        strict = config.get("strict_cluster_holdout")
        folds = strict.get("development_folds") if isinstance(strict, dict) else None
        if not isinstance(folds, int) or folds < 1:
            raise RuntimeError("runtime strict task configuration is invalid")
        expected.update(
            ("strict_cluster_holdout", fold, fold, "pu_logistic")
            for fold in range(folds)
        )
    if "literature_random_protein" in tracks:
        random_track = config.get("literature_random_protein")
        seeds = random_track.get("seeds") if isinstance(random_track, dict) else None
        models = config.get("models")
        comparison_inputs = config.get("comparison_inputs")
        if (
            not isinstance(seeds, list)
            or not all(isinstance(seed, int) for seed in seeds)
            or not isinstance(models, list)
            or not all(isinstance(model, str) and model for model in models)
            or not isinstance(comparison_inputs, dict)
        ):
            raise RuntimeError("runtime literature task configuration is invalid")
        configured_models = list(models)
        sul_manifest_value = comparison_inputs.get("sul_environment_manifest")
        external_environments = manifest.get("external_environments") or {}
        if sul_manifest_value and (
            "sul_bertgru" in external_environments
            or _is_registered_comparison_input(Path(str(sul_manifest_value)))
        ):
            configured_models.append("sul_bertgru")
        expected.update(
            ("literature_random_protein", 0, seed, model)
            for seed in seeds
            for model in configured_models
        )
    observed = {(item.track, item.fold, item.seed, item.model) for item in roster}
    if observed != expected:
        raise RuntimeError("runtime configured task roster coverage is incomplete")


def random_protein_train_validation_test(
    rows: Iterable[MultispeciesSiteRow],
    *,
    seed: int,
    test_fraction: float = 0.2,
    validation_fraction_of_remaining: float = 0.2,
) -> dict[str, list[MultispeciesSiteRow]]:
    """Make the Sul-BertGRU-comparable random *protein* split.

    The inner validation partition is also grouped by protein.  This avoids
    the site-row leakage that an ordinary shuffled validation set permits.
    """
    materialized = list(rows)
    proteins = sorted({row.global_protein_id for row in materialized})
    if len(proteins) < 3:
        raise RuntimeError("random protein split requires at least three proteins")
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test_fraction must be in (0, 1)")
    if not 0.0 < validation_fraction_of_remaining < 1.0:
        raise ValueError("validation_fraction_of_remaining must be in (0, 1)")
    rng = random.Random(seed)
    rng.shuffle(proteins)
    n_test = max(1, round(len(proteins) * test_fraction))
    if n_test >= len(proteins):
        raise RuntimeError("random protein split leaves no training proteins")
    test_ids = set(proteins[:n_test])
    remaining = proteins[n_test:]
    n_validation = max(1, round(len(remaining) * validation_fraction_of_remaining))
    if n_validation >= len(remaining):
        raise RuntimeError("random protein split leaves no training proteins")
    validation_ids = set(remaining[:n_validation])
    train_ids = set(remaining[n_validation:])
    result: dict[str, list[MultispeciesSiteRow]] = {
        "train": [],
        "validation": [],
        "test": [],
    }
    for row in materialized:
        if row.global_protein_id in test_ids:
            result["test"].append(row)
        elif row.global_protein_id in validation_ids:
            result["validation"].append(row)
        elif row.global_protein_id in train_ids:
            result["train"].append(row)
        else:  # pragma: no cover - protects future changes to the branches above
            raise RuntimeError(f"unassigned protein: {row.global_protein_id}")
    return result


def subsample_unlabeled_by_partition(
    partitions: dict[str, list[MultispeciesSiteRow]],
    *,
    per_positive: int,
    seed: int,
) -> dict[str, list[MultispeciesSiteRow]]:
    """Sample PU panels only after proving the input partitions are disjoint."""
    if per_positive < 1:
        raise ValueError("per_positive must be >= 1")
    seen: set[tuple[str, int]] = set()
    for name, rows in partitions.items():
        for row in rows:
            key = (row.global_protein_id, row.cys_position)
            if key in seen:
                raise RuntimeError(
                    f"site row appears in more than one partition: {name}:{key}"
                )
            seen.add(key)
    sampled: dict[str, list[MultispeciesSiteRow]] = {}
    for index, (name, rows) in enumerate(sorted(partitions.items())):
        positives = [row for row in rows if row.label == "positive"]
        unlabeled = [row for row in rows if row.label != "positive"]
        keep = min(len(unlabeled), per_positive * len(positives))
        rng = random.Random(seed + index)
        chosen = rng.sample(unlabeled, keep)
        sampled[name] = sorted(
            [*positives, *chosen],
            key=lambda row: (row.global_protein_id, row.cys_position),
        )
    return sampled


def development_fold_rows(
    rows: Iterable[MultispeciesSiteRow],
    frozen_split: FrozenMultispeciesSplit,
    *,
    validation_fold: int,
) -> tuple[list[MultispeciesSiteRow], list[MultispeciesSiteRow]]:
    """Return one train/validation fold without ever returning test rows."""
    fold_by_protein = {
        row.global_protein_id: row.development_fold
        for row in frozen_split.rows
        if row.split == DEVELOPMENT_SPLIT
    }
    test_ids = {
        row.global_protein_id for row in frozen_split.rows if row.split == TEST_SPLIT
    }
    if not 0 <= validation_fold < frozen_split.n_development_folds:
        raise ValueError("validation_fold outside configured development folds")
    train: list[MultispeciesSiteRow] = []
    validation: list[MultispeciesSiteRow] = []
    for row in rows:
        protein_id = row.global_protein_id
        if protein_id in test_ids:
            continue
        fold = fold_by_protein.get(protein_id)
        if fold is None:
            raise RuntimeError(f"site protein missing from frozen split: {protein_id}")
        if fold == validation_fold:
            validation.append(row)
        else:
            train.append(row)
    return train, validation


def assert_test_unlocked(
    config_path: Path,
    code_revision: str,
    split_sha256: str,
    unlock_path: Path | None,
) -> None:
    """Permit test scoring only with an exact pre-approved unlock record."""
    if unlock_path is None or not unlock_path.is_file():
        raise RuntimeError("frozen test is locked")
    try:
        payload = json.loads(unlock_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise RuntimeError("test unlock is not valid JSON") from exc
    expected = {
        "config_sha256": _sha256_file(config_path),
        "code_revision": code_revision,
        "split_sha256": split_sha256,
    }
    if not isinstance(payload, dict) or any(
        payload.get(key) != value for key, value in expected.items()
    ):
        raise RuntimeError("test unlock mismatch")


def run_multispecies_experiment(
    config_path: Path,
    *,
    score_frozen_test: bool = False,
    test_unlock_path: Path | None = None,
    code_revision: str,
) -> MultispeciesExperimentPreparation:
    """Validate the immutable v2 experiment boundary before model execution.

    This preflight deliberately does not materialize or score test rows unless
    a matching unlock record is supplied. Model-specific execution is layered
    on top of this boundary so no adapter can bypass it.
    """
    try:
        cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(
            f"invalid multispecies experiment config: {config_path}"
        ) from exc
    if not isinstance(cfg, dict) or cfg.get("version") != 2:
        raise RuntimeError("multispecies v2 config is required")
    strict = cfg.get("strict_cluster_holdout")
    if not isinstance(strict, dict):
        raise RuntimeError("strict_cluster_holdout configuration is required")
    data_isolation = cfg.get("data_isolation")
    required_isolation = {
        "split_before_unlabeled_sampling": True,
        "development_fit_scope": "training_fold_only",
        "panther_role": "unlabelled_evolutionary_feature_only",
    }
    if data_isolation != required_isolation:
        raise RuntimeError("data_isolation configuration is required")
    validate_reference_proteome_inputs(cfg)
    cluster_key = strict.get("cluster_table")
    if not cluster_key:
        global_cfg = cfg.get("global_mmseqs2")
        if not isinstance(global_cfg, dict) or not global_cfg.get("cluster_table"):
            raise RuntimeError("strict cluster table is required")
        cluster_key = global_cfg["cluster_table"]
    clusters = load_global_cluster_table(Path(str(cluster_key)))
    split = load_frozen_multispecies_split(Path(str(strict["split_path"])))
    split.require_trainable()
    cluster_ids = {row.global_protein_id for row in clusters}
    split_ids = {row.global_protein_id for row in split.rows}
    if cluster_ids != split_ids:
        raise RuntimeError("global cluster table and frozen split disagree")
    if score_frozen_test:
        if bool(strict.get("test_unlock_required", True)):
            assert_test_unlocked(
                config_path, code_revision, split.sha256, test_unlock_path
            )
    return MultispeciesExperimentPreparation(
        config_path=config_path,
        config_sha256=_sha256_file(config_path),
        split=split,
        cluster_count=len({row.cluster_id for row in clusters}),
        test_scoring_enabled=score_frozen_test,
    )
