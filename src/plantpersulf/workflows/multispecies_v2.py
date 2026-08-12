"""Leakage-safe split helpers for the multispecies v2 experiment runner."""

from __future__ import annotations

import hashlib
import json
import random
import tempfile
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path

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


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(asdict(event), sort_keys=True) + "\n")


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
    split_sha256: str,
    code_revision: str,
    input_sha256: dict[str, str],
    command: list[str],
    device: str,
) -> dict[str, object]:
    """Write an atomic run manifest after verifying every input hash."""
    if not _valid_sha256(config_sha256) or not _valid_sha256(split_sha256):
        raise ValueError("config_sha256 and split_sha256 must be SHA256 values")
    if not input_sha256 or any(
        not name or not _valid_sha256(value) for name, value in input_sha256.items()
    ):
        raise ValueError("input_sha256 must contain complete SHA256 values")
    if not code_revision or not command or not device:
        raise ValueError("code_revision, command, and device are required")
    manifest: dict[str, object] = {
        "config_sha256": config_sha256,
        "split_sha256": split_sha256,
        "code_revision": code_revision,
        "input_sha256": dict(sorted(input_sha256.items())),
        "command": command,
        "device": device,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)
    return manifest


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
