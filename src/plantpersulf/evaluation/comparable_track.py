"""Eligibility and reporting contracts for the Task 9.4 comparison track."""

from __future__ import annotations

import csv
import hashlib
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.benchmark.literature_random_track import (
    ComparisonModelInput,
    ComparisonModelScores,
    ComparisonSite,
    assert_complete_model_scores,
    prepare_sul_bertgru_adapter_rows,
)


@dataclass(frozen=True)
class ComparatorStatus:
    model: str
    role: str
    status: str
    reason: str


@dataclass(frozen=True)
class SulEnvironment:
    python_executable: Path
    environment_lock: Path
    adapter_path: Path
    repo_commit: str
    parameters: tuple[tuple[str, str], ...] = ()
    input_artifacts: tuple[tuple[Path, str], ...] = ()


def build_registered_structure_features(
    site_keys: set[tuple[str, int]],
    *,
    registry_path: Path,
    registry_base: Path,
) -> dict[tuple[str, int], tuple[tuple[float, float], bool]]:
    """Extract registered AlphaFold features and explicitly mask absence."""
    from plantpersulf.download.alphafold import audit_alphafold_structures
    from plantpersulf.features.structure import extract_cys_structure_features

    sources = {
        source.accession: source
        for source in audit_alphafold_structures(
            registry_path, base_directory=registry_base
        )
    }
    keys_by_accession: dict[str, list[tuple[str, int]]] = {}
    for key in sorted(site_keys):
        global_protein_id, _ = key
        try:
            _, accession = global_protein_id.split("|", 1)
        except ValueError as exc:
            raise RuntimeError(
                f"invalid global protein key for structure feature: {key}"
            ) from exc
        keys_by_accession.setdefault(accession, []).append(key)

    result: dict[tuple[str, int], tuple[tuple[float, float], bool]] = {}
    for accession, keys in keys_by_accession.items():
        source = sources.get(accession)
        if source is None:
            result.update({key: ((0.0, 0.0), False) for key in keys})
            continue
        positions = sorted({position for _, position in keys})
        extracted = {
            row.cys_position: row
            for row in extract_cys_structure_features(
                source.local_path.read_text(encoding="utf-8"),
                accession,
                positions,
            )
        }
        for key in keys:
            feature = extracted[key[1]]
            if (
                not feature.has_structure
                or feature.contact_number_proxy is None
                or feature.plddt is None
            ):
                result[key] = ((0.0, 0.0), False)
            else:
                result[key] = (
                    (feature.contact_number_proxy, feature.plddt),
                    True,
                )
    return result


def _registered_ready(path: Path | None, registered: frozenset[Path]) -> bool:
    return path is not None and path.is_file() and path.resolve() in registered


def comparator_statuses(
    *,
    esm_features: Path | None,
    structure_features: Path | None,
    sul_environment_manifest: Path | None,
    pcysmod_scores: Path | None,
    registered_input_paths: frozenset[Path] = frozenset(),
) -> tuple[ComparatorStatus, ...]:
    """Return the complete comparator roster without fallback substitution."""
    registered = frozenset(path.resolve() for path in registered_input_paths)
    result = [
        ComparatorStatus(model, "direct_baseline", "ready", "in_process_features")
        for model in ("pu_logistic", "random_forest", "xgboost")
    ]
    result.append(
        ComparatorStatus(
            "esm_linear_head",
            "direct_baseline",
            "ready" if _registered_ready(esm_features, registered) else "blocked",
            (
                "registered_esm_features"
                if _registered_ready(esm_features, registered)
                else "missing_registered_esm_features"
            ),
        )
    )
    result.append(
        ComparatorStatus(
            "structure_ranker",
            "direct_baseline",
            (
                "ready"
                if _registered_ready(structure_features, registered)
                else "blocked"
            ),
            (
                "registered_structure_features"
                if _registered_ready(structure_features, registered)
                else "missing_registered_structure_features"
            ),
        )
    )
    result.append(
        ComparatorStatus(
            "sul_bertgru",
            "external_comparator",
            (
                "ready"
                if _registered_ready(sul_environment_manifest, registered)
                else "blocked"
            ),
            (
                "registered_isolated_environment"
                if _registered_ready(sul_environment_manifest, registered)
                else "missing_registered_isolated_environment"
            ),
        )
    )
    result.append(
        ComparatorStatus(
            "pcysmod",
            "external_comparator",
            (
                "ready"
                if _registered_ready(pcysmod_scores, registered)
                else "qualitative_only"
            ),
            (
                "registered_complete_scores_pending_coverage_audit"
                if _registered_ready(pcysmod_scores, registered)
                else "complete_frozen_input_scores_unavailable"
            ),
        )
    )
    result.extend(
        ComparatorStatus(
            model,
            "architecture_reference",
            "architecture_reference_only",
            "not_a_plant_persulfidation_competitor",
        )
        for model in ("tree", "graft")
    )
    return tuple(result)


def _species_ap(
    rows: tuple[ComparisonSite, ...],
    scores: dict[tuple[str, int], float],
    species: str,
) -> dict[str, float | int]:
    selected = [row for row in rows if row.species == species]
    labels = [row.label == "positive" for row in selected]
    if not selected or not any(labels):
        raise RuntimeError(f"comparison partition lacks positives for {species}")
    from sklearn.metrics import average_precision_score  # type: ignore[import-untyped]

    values = [scores[row.site_key] for row in selected]
    return {
        "average_precision": float(average_precision_score(labels, values)),
        "base_rate": sum(labels) / len(labels),
        "site_count": len(selected),
        "positive_count": sum(labels),
    }


def summarize_comparable_scores(
    model_input: ComparisonModelInput,
    model_scores: ComparisonModelScores,
    *,
    primary_species: tuple[str, ...],
    pressure_species: tuple[str, ...],
) -> dict[str, object]:
    """Summarize one shared-panel run without mixing fungal pressure data."""
    if (
        model_scores.model != model_input.model
        or model_scores.seed != model_input.seed
        or model_scores.panel_sha256 != model_input.panel_sha256
    ):
        raise RuntimeError("comparison scores do not match bound panel")
    if set(primary_species) & set(pressure_species):
        raise ValueError("primary and pressure species must be disjoint")
    partitions = dict(model_input.partition_rows)
    scored = dict(model_scores.partition_scores)
    for name in ("validation", "test"):
        assert_complete_model_scores(model_input, name, scored[name])
    test_rows = partitions["test"]
    test_scores = scored["test"]
    validation_rows = partitions["validation"]
    validation_scores = scored["validation"]
    validation_primary = {
        species: _species_ap(validation_rows, validation_scores, species)
        for species in primary_species
    }
    primary = {
        species: _species_ap(test_rows, test_scores, species)
        for species in primary_species
    }
    pressure = {
        species: _species_ap(test_rows, test_scores, species)
        for species in pressure_species
    }
    macro_ap = sum(
        float(primary[species]["average_precision"])
        for species in primary_species
    ) / len(primary_species)
    validation_macro_ap = sum(
        float(validation_primary[species]["average_precision"])
        for species in primary_species
    ) / len(primary_species)
    return {
        "track": "literature_random_protein_development_zone",
        "model": model_input.model,
        "seed": model_input.seed,
        "panel_sha256": model_input.panel_sha256,
        "three_crop_primary_validation": validation_primary,
        "three_crop_validation_macro_average_precision": validation_macro_ap,
        "three_crop_primary_test": primary,
        "three_crop_macro_average_precision": macro_ap,
        "magnaporthe_pressure_test": pressure,
        "limitation": "within_dataset_literature_comparable_not_gate2",
    }


def load_complete_external_scores(
    path: Path,
    model_input: ComparisonModelInput,
    *,
    partition: str,
) -> dict[tuple[str, int], float]:
    """Load an external comparator only when every frozen site is scored."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_columns = ("global_protein_id", "cys_position", "score")
        if tuple(reader.fieldnames or ()) != expected_columns:
            raise RuntimeError("external score table has invalid columns")
        scores: dict[tuple[str, int], float] = {}
        for row in reader:
            key = (row["global_protein_id"], int(row["cys_position"]))
            if key in scores:
                raise RuntimeError(f"duplicate external score: {key}")
            scores[key] = float(row["score"])
    assert_complete_model_scores(model_input, partition, scores)
    return scores


def _manifest_path(manifest: Path, value: object) -> Path:
    path = Path(str(value))
    if not path.is_absolute():
        path = manifest.parent / path
    return path.resolve()


def _file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _artifact_sha256(path: Path) -> str:
    if path.is_file():
        return _file_sha256(path)
    if not path.is_dir():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        digest.update(child.relative_to(path).as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(_file_sha256(child)))
    return digest.hexdigest()


def validate_sul_environment_manifest(path: Path) -> SulEnvironment:
    """Validate an isolated, hash-bound Sul-BertGRU adapter environment."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError("Sul-BertGRU environment manifest is invalid") from exc
    if not isinstance(payload, dict) or payload.get("isolated") is not True:
        raise RuntimeError("Sul-BertGRU environment must be isolated")
    python_executable = _manifest_path(path, payload.get("python_executable"))
    if python_executable == Path(sys.executable).resolve():
        raise RuntimeError("Sul-BertGRU requires a separate Python environment")
    if not python_executable.is_file():
        raise RuntimeError("Sul-BertGRU isolated Python is missing")
    environment_lock = _manifest_path(path, payload.get("environment_lock"))
    adapter_path = _manifest_path(path, payload.get("adapter_path"))
    for artifact, hash_key in (
        (environment_lock, "environment_lock_sha256"),
        (adapter_path, "adapter_sha256"),
    ):
        if not artifact.is_file() or _file_sha256(artifact) != payload.get(hash_key):
            raise RuntimeError(f"Sul-BertGRU {hash_key} mismatch")
    repo_commit = str(payload.get("repo_commit", ""))
    if re.fullmatch(r"[0-9a-f]{40}", repo_commit) is None:
        raise RuntimeError("Sul-BertGRU repository commit is invalid")
    raw_parameters = payload.get("parameters", {})
    if not isinstance(raw_parameters, dict) or any(
        not isinstance(key, str)
        or not key.replace("_", "").isalnum()
        or isinstance(value, (dict, list))
        for key, value in raw_parameters.items()
    ):
        raise RuntimeError("Sul-BertGRU adapter parameters are invalid")
    raw_artifacts = payload.get("input_artifacts", [])
    if not isinstance(raw_artifacts, list):
        raise RuntimeError("Sul-BertGRU input artifacts are invalid")
    input_artifacts: list[tuple[Path, str]] = []
    for artifact in raw_artifacts:
        if not isinstance(artifact, dict):
            raise RuntimeError("Sul-BertGRU input artifact is invalid")
        artifact_path = _manifest_path(path, artifact.get("path"))
        artifact_sha256 = str(artifact.get("sha256", ""))
        if (
            re.fullmatch(r"[0-9a-f]{64}", artifact_sha256) is None
            or _artifact_sha256(artifact_path) != artifact_sha256
        ):
            raise RuntimeError("Sul-BertGRU input artifact SHA256 mismatch")
        input_artifacts.append((artifact_path, artifact_sha256))
    return SulEnvironment(
        python_executable=python_executable,
        environment_lock=environment_lock,
        adapter_path=adapter_path,
        repo_commit=repo_commit,
        parameters=tuple(
            sorted((key, str(value)) for key, value in raw_parameters.items())
        ),
        input_artifacts=tuple(input_artifacts),
    )


def run_sul_bertgru_adapter(
    model_input: ComparisonModelInput,
    sequences: dict[str, str],
    environment: SulEnvironment,
    *,
    work_directory: Path,
    timeout_seconds: float = 21600.0,
    device: str | None = None,
    oom_deviations: list[dict[str, int]] | None = None,
) -> ComparisonModelScores:
    """Run the isolated Sul-BertGRU adapter on the exact shared panel."""
    if timeout_seconds <= 0.0:
        raise ValueError("Sul-BertGRU timeout must be positive")
    adapter_rows = prepare_sul_bertgru_adapter_rows(model_input, sequences)
    work_directory.mkdir(parents=True, exist_ok=False)
    input_path = work_directory / "shared_panel.tsv"
    with input_path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            (
                "partition",
                "global_protein_id",
                "cys_position",
                "sequence_window",
                "adapter_label",
            )
        )
        for row in adapter_rows:
            global_protein_id, cys_position = row.site_key
            writer.writerow(
                (
                    row.partition,
                    global_protein_id,
                    cys_position,
                    row.sequence_window,
                    row.adapter_label,
                )
            )
    parameters = dict(environment.parameters)
    if device is not None:
        if re.fullmatch(r"cpu|cuda(?::\d+)?", device) is None:
            raise ValueError("Sul-BertGRU device is invalid")
        parameters["device"] = device
    batch_values = [
        int(parameters[name])
        for name in ("bert_batch_size", "train_batch_size")
        if name in parameters
    ]
    initial_batch_size = min(batch_values) if batch_values else 1

    def execute(batch_size: int) -> Path:
        output_path = work_directory / f"scores_batch{batch_size}.tsv"
        command = [
                str(environment.python_executable),
                str(environment.adapter_path),
                "--input",
                str(input_path),
                "--output",
                str(output_path),
                "--seed",
                str(model_input.seed),
            ]
        attempt_parameters = dict(parameters)
        for name in ("bert_batch_size", "train_batch_size"):
            if name in attempt_parameters:
                attempt_parameters[name] = str(batch_size)
        for name, value in sorted(attempt_parameters.items()):
            command.extend((f"--{name.replace('_', '-')}", value))
        try:
            subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=timeout_seconds,
            )
        except subprocess.TimeoutExpired as exc:
            raise RuntimeError("Sul-BertGRU adapter timed out") from exc
        except subprocess.CalledProcessError as exc:
            diagnostic = f"{exc.stdout or ''}\n{exc.stderr or ''}"
            if "out of memory" in diagnostic.lower():
                raise RuntimeError("CUDA out of memory") from exc
            raise RuntimeError("Sul-BertGRU adapter execution failed") from exc
        return output_path

    from plantpersulf.workflows.multispecies_v2 import run_with_oom_batch_retry

    output_path, _, deviations = run_with_oom_batch_retry(
        execute, initial_batch_size=initial_batch_size
    )
    if oom_deviations is not None:
        oom_deviations.extend(deviations)
    if not output_path.is_file():
        raise RuntimeError("Sul-BertGRU adapter did not produce scores")

    scores_by_partition: dict[str, dict[tuple[str, int], float]] = {
        "validation": {},
        "test": {},
    }
    with output_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_columns = (
            "partition",
            "global_protein_id",
            "cys_position",
            "score",
        )
        if tuple(reader.fieldnames or ()) != expected_columns:
            raise RuntimeError("Sul-BertGRU score table has invalid columns")
        for score_row in reader:
            partition = score_row["partition"]
            if partition not in scores_by_partition:
                raise RuntimeError(
                    "Sul-BertGRU may only return validation and test scores"
                )
            key = (
                score_row["global_protein_id"],
                int(score_row["cys_position"]),
            )
            partition_scores = scores_by_partition[partition]
            if key in partition_scores:
                raise RuntimeError(f"duplicate Sul-BertGRU score: {key}")
            partition_scores[key] = float(score_row["score"])
    for partition, scores in scores_by_partition.items():
        assert_complete_model_scores(model_input, partition, scores)
    return ComparisonModelScores(
        model="sul_bertgru",
        seed=model_input.seed,
        panel_sha256=model_input.panel_sha256,
        partition_scores=tuple(scores_by_partition.items()),
    )
