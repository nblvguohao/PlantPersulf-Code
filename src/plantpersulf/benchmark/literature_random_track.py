"""Development-only random-protein comparison track for Task 9.4."""

from __future__ import annotations

import hashlib
import math
import random
from dataclasses import dataclass

from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow


@dataclass(frozen=True)
class LiteratureRandomRun:
    """One immutable 64/16/20 protein split and its shared PU panel."""

    seed: int
    protein_partitions: tuple[tuple[str, str], ...]
    partitions: dict[str, tuple[MultispeciesV2SiteRow, ...]]
    panel_sha256: str


@dataclass(frozen=True)
class ComparisonSite:
    """Immutable shared site/PU-label record supplied to every model."""

    species: str
    protein_accession: str
    cys_position: int
    label: str
    study_accessions: tuple[str, ...]

    @property
    def site_key(self) -> tuple[str, int]:
        return f"{self.species}|{self.protein_accession}", self.cys_position


@dataclass(frozen=True)
class ComparisonModelInput:
    """Model-specific handle bound to one byte-identical comparison panel."""

    model: str
    seed: int
    panel_sha256: str
    partition_rows: tuple[tuple[str, tuple[ComparisonSite, ...]], ...]


@dataclass(frozen=True)
class SulBertGruAdapterRow:
    """Isolated-adapter view; project PU labels remain untouched."""

    partition: str
    site_key: tuple[str, int]
    sequence_window: str
    adapter_label: str


@dataclass(frozen=True)
class ComparisonModelScores:
    model: str
    seed: int
    panel_sha256: str
    partition_scores: tuple[
        tuple[str, dict[tuple[str, int], float]], ...
    ]


def _site_key(row: MultispeciesV2SiteRow) -> tuple[str, int]:
    return row.global_protein_id, row.cys_position


def _panel_sha256(
    partitions: dict[str, tuple[MultispeciesV2SiteRow, ...]],
) -> str:
    lines: list[str] = []
    for partition_name in ("train", "validation", "test"):
        for row in partitions[partition_name]:
            lines.append(
                "\t".join(
                    (
                        partition_name,
                        row.species,
                        row.protein_accession,
                        str(row.cys_position),
                        row.label,
                        ";".join(row.study_accessions),
                    )
                )
            )
    return hashlib.sha256("\n".join(lines).encode("utf-8")).hexdigest()


def build_literature_random_track(
    rows: tuple[MultispeciesV2SiteRow, ...] | list[MultispeciesV2SiteRow],
    *,
    seeds: tuple[int, ...],
    test_fraction: float,
    validation_fraction_of_remaining: float,
    unlabeled_per_positive: int,
    sampling_seed: int,
) -> tuple[LiteratureRandomRun, ...]:
    """Split development proteins, then sample PU backgrounds per partition."""
    if not 0.0 < test_fraction < 1.0:
        raise ValueError("test_fraction must be in (0, 1)")
    if not 0.0 < validation_fraction_of_remaining < 1.0:
        raise ValueError("validation fraction must be in (0, 1)")
    if unlabeled_per_positive < 1:
        raise ValueError("unlabeled_per_positive must be positive")
    materialized = tuple(rows)
    seen: set[tuple[str, int]] = set()
    for row in materialized:
        if row.split != "development":
            raise RuntimeError(
                f"strict frozen test row entered random track: {row.global_protein_id}"
            )
        key = _site_key(row)
        if key in seen:
            raise RuntimeError(f"duplicate comparison site: {key}")
        seen.add(key)

    proteins = sorted({row.global_protein_id for row in materialized})
    if len(proteins) < 3:
        raise RuntimeError("random protein track requires at least three proteins")
    runs: list[LiteratureRandomRun] = []
    for seed in seeds:
        shuffled = proteins[:]
        random.Random(seed).shuffle(shuffled)
        n_test = max(1, round(len(shuffled) * test_fraction))
        remaining = shuffled[n_test:]
        n_validation = max(
            1, round(len(remaining) * validation_fraction_of_remaining)
        )
        if n_test + n_validation >= len(shuffled):
            raise RuntimeError("random protein split leaves no training proteins")
        protein_partition = {
            **{protein: "test" for protein in shuffled[:n_test]},
            **{
                protein: "validation"
                for protein in remaining[:n_validation]
            },
            **{protein: "train" for protein in remaining[n_validation:]},
        }
        unsampled: dict[str, list[MultispeciesV2SiteRow]] = {
            "train": [],
            "validation": [],
            "test": [],
        }
        for row in materialized:
            unsampled[protein_partition[row.global_protein_id]].append(row)

        sampled: dict[str, tuple[MultispeciesV2SiteRow, ...]] = {}
        for index, partition_name in enumerate(("train", "validation", "test")):
            partition = unsampled[partition_name]
            positives = [row for row in partition if row.label == "positive"]
            unlabeled = [row for row in partition if row.label == "unlabeled"]
            keep = min(len(unlabeled), unlabeled_per_positive * len(positives))
            chosen = random.Random(sampling_seed + seed * 3 + index).sample(
                unlabeled, keep
            )
            sampled[partition_name] = tuple(
                sorted((*positives, *chosen), key=_site_key)
            )
        runs.append(
            LiteratureRandomRun(
                seed=seed,
                protein_partitions=tuple(sorted(protein_partition.items())),
                partitions=sampled,
                panel_sha256=_panel_sha256(sampled),
            )
        )
    return tuple(runs)


def bind_model_to_comparison_panel(
    run: LiteratureRandomRun, model: str
) -> ComparisonModelInput:
    """Bind a model name without allowing it to alter sites or PU labels."""
    if competitor_role(model) == "unknown":
        raise ValueError(f"unknown comparison model: {model}")
    partitions: list[tuple[str, tuple[ComparisonSite, ...]]] = []
    for name in ("train", "validation", "test"):
        partitions.append(
            (
                name,
                tuple(
                    ComparisonSite(
                        species=row.species,
                        protein_accession=row.protein_accession,
                        cys_position=row.cys_position,
                        label=row.label,
                        study_accessions=row.study_accessions,
                    )
                    for row in run.partitions[name]
                ),
            )
        )
    return ComparisonModelInput(
        model=model,
        seed=run.seed,
        panel_sha256=run.panel_sha256,
        partition_rows=tuple(partitions),
    )


def assert_complete_model_scores(
    model_input: ComparisonModelInput,
    partition: str,
    scores: dict[tuple[str, int], float],
) -> None:
    """Require exact score coverage before any comparator enters a table."""
    try:
        rows = dict(model_input.partition_rows)[partition]
    except KeyError as exc:
        raise ValueError(f"unknown comparison partition: {partition}") from exc
    expected = {row.site_key for row in rows}
    has_nonfinite = any(not math.isfinite(value) for value in scores.values())
    if set(scores) != expected or has_nonfinite:
        raise RuntimeError(
            f"{model_input.model} did not score the complete frozen input"
        )


def competitor_role(model: str) -> str:
    """Classify methods without promoting architecture inspirations to SOTA."""
    if model in {
        "pu_logistic",
        "random_forest",
        "xgboost",
        "esm_linear_head",
        "structure_ranker",
    }:
        return "direct_baseline"
    if model == "sul_bertgru":
        return "isolated_external_adapter"
    if model == "pcysmod":
        return "quantitative_only_with_complete_scores"
    if model in {"tree", "graft"}:
        return "architecture_reference_only"
    return "unknown"


def prepare_sul_bertgru_adapter_rows(
    model_input: ComparisonModelInput,
    sequences: dict[str, str],
) -> tuple[SulBertGruAdapterRow, ...]:
    """Create 31-aa adapter windows without mutating project PU labels."""
    if model_input.model != "sul_bertgru":
        raise ValueError("Sul-BertGRU adapter requires its bound model input")
    result: list[SulBertGruAdapterRow] = []
    flank = 15
    for partition, rows in model_input.partition_rows:
        for row in rows:
            protein_id, position = row.site_key
            sequence = sequences.get(protein_id)
            if sequence is None:
                raise RuntimeError(f"missing adapter sequence: {protein_id}")
            invalid_position = position < 1 or position > len(sequence)
            if invalid_position or sequence[position - 1] != "C":
                raise RuntimeError(f"adapter coordinate is not Cys: {row.site_key}")
            start = position - 1 - flank
            end = position + flank
            left_padding = "X" * max(0, -start)
            right_padding = "X" * max(0, end - len(sequence))
            window = left_padding + sequence[max(0, start) : min(len(sequence), end)]
            window += right_padding
            if len(window) != 31 or window[15] != "C":
                raise RuntimeError(f"invalid Sul-BertGRU window: {row.site_key}")
            result.append(
                SulBertGruAdapterRow(
                    partition=partition,
                    site_key=row.site_key,
                    sequence_window=window,
                    adapter_label=(
                        "positive" if row.label == "positive" else "putative_negative"
                    ),
                )
            )
    return tuple(result)


def run_tabular_direct_baseline(
    model_input: ComparisonModelInput,
    features: dict[tuple[str, int], tuple[float, ...]],
    parameters: dict[str, int | float] | None = None,
) -> ComparisonModelScores:
    """Fit one direct baseline once and score shared validation/test rows."""
    partitions = dict(model_input.partition_rows)
    train = partitions["train"]
    predict = (*partitions["validation"], *partitions["test"])
    train_labels = [row.label for row in train]
    parameters = parameters or {}

    if model_input.model == "esm_linear_head":
        import numpy as np

        from plantpersulf.models.esm_baseline import esm_linear_head_scores

        try:
            esm_train_features = np.stack(
                [features[row.site_key] for row in train]
            ).astype(np.float32, copy=False)
            esm_predict_features = np.stack(
                [features[row.site_key] for row in predict]
            ).astype(np.float32, copy=False)
        except KeyError as exc:
            message = f"missing registered comparison feature: {exc.args[0]}"
            raise RuntimeError(message) from exc
        scores = esm_linear_head_scores(
            esm_train_features,
            train_labels,
            esm_predict_features,
            model_input.seed,
        )
    else:
        try:
            tabular_train_features = [
                list(features[row.site_key]) for row in train
            ]
            tabular_predict_features = [
                list(features[row.site_key]) for row in predict
            ]
        except KeyError as exc:
            message = f"missing registered comparison feature: {exc.args[0]}"
            raise RuntimeError(message) from exc
        from plantpersulf.models.traditional import (
            TrainOnlyScaler,
            pu_logistic_regression_scores,
            random_forest_scores,
            xgboost_scores,
        )

        scaler = TrainOnlyScaler.fit(tabular_train_features)
        scaled_train = scaler.transform(tabular_train_features)
        scaled_predict = scaler.transform(tabular_predict_features)
        if model_input.model == "pu_logistic":
            scores = pu_logistic_regression_scores(
                scaled_train,
                train_labels,
                scaled_predict,
                seed=model_input.seed,
                holdout_fraction=float(parameters.get("holdout_fraction", 0.2)),
            )
        elif model_input.model == "random_forest":
            scores = random_forest_scores(
                scaled_train,
                train_labels,
                scaled_predict,
                seed=model_input.seed,
                n_estimators=int(parameters.get("n_estimators", 100)),
            )
        elif model_input.model == "xgboost":
            scores = xgboost_scores(
                scaled_train,
                train_labels,
                scaled_predict,
                seed=model_input.seed,
                n_estimators=int(parameters.get("n_estimators", 100)),
            )
        else:
            raise ValueError(
                f"model requires a non-tabular adapter: {model_input.model}"
            )

    n_validation = len(partitions["validation"])
    output: list[tuple[str, dict[tuple[str, int], float]]] = []
    for name, rows, values in (
        ("validation", partitions["validation"], scores[:n_validation]),
        ("test", partitions["test"], scores[n_validation:]),
    ):
        mapped = {
            row.site_key: float(value) for row, value in zip(rows, values, strict=True)
        }
        assert_complete_model_scores(model_input, name, mapped)
        output.append((name, mapped))
    return ComparisonModelScores(
        model=model_input.model,
        seed=model_input.seed,
        panel_sha256=model_input.panel_sha256,
        partition_scores=tuple(output),
    )


def run_structure_direct_baseline(
    model_input: ComparisonModelInput,
    sequence_features: dict[tuple[str, int], tuple[float, ...]],
    esm_features: dict[tuple[str, int], tuple[float, ...]] | None,
    structure_features: dict[
        tuple[str, int], tuple[tuple[float, float], bool]
    ],
    parameters: dict[str, int | float] | None = None,
    device: str = "cpu",
) -> ComparisonModelScores:
    """Fit the current structure-aware ranker on one shared protein panel."""
    if model_input.model != "structure_ranker":
        raise ValueError("structure baseline requires its bound model input")
    from plantpersulf.models.structure_ranker import (
        AblationConfig,
        BranchFeatures,
        structure_ranker_scores,
    )

    partitions = dict(model_input.partition_rows)
    train_rows = partitions["train"]
    predict_rows = (*partitions["validation"], *partitions["test"])
    parameters = parameters or {}
    use_esm = bool(parameters.get("use_esm", 0))
    if use_esm and esm_features is None:
        raise RuntimeError("structure ranker requires registered ESM features")

    def branches(rows: tuple[ComparisonSite, ...]) -> BranchFeatures:
        try:
            sequence = [list(sequence_features[row.site_key]) for row in rows]
            structure_pairs = [structure_features[row.site_key] for row in rows]
        except KeyError as exc:
            raise RuntimeError(
                f"missing registered structure-ranker feature: {exc.args[0]}"
            ) from exc
        esm = (
            [list(esm_features[row.site_key]) for row in rows]
            if use_esm and esm_features is not None
            else [[0.0] for _ in rows]
        )
        return BranchFeatures(
            sequence=sequence,
            esm=esm,
            structure=[list(values) for values, _ in structure_pairs],
            structure_mask=[available for _, available in structure_pairs],
            study_ids=None,
        )

    output = structure_ranker_scores(
        branches(train_rows),
        [row.label for row in train_rows],
        branches(predict_rows),
        seed=model_input.seed,
        ablation=AblationConfig(
            use_esm=use_esm,
            use_study_context=False,
        ),
        hidden=int(parameters.get("hidden", 16)),
        dropout=float(parameters.get("dropout", 0.2)),
        epochs=int(parameters.get("epochs", 200)),
        lr=float(parameters.get("lr", 0.05)),
        holdout_fraction=float(parameters.get("holdout_fraction", 0.2)),
        n_mc_dropout=int(parameters.get("n_mc_dropout", 16)),
        device_name=device,
    )
    n_validation = len(partitions["validation"])
    partition_scores: list[
        tuple[str, dict[tuple[str, int], float]]
    ] = []
    for name, rows, values in (
        (
            "validation",
            partitions["validation"],
            output.scores[:n_validation],
        ),
        ("test", partitions["test"], output.scores[n_validation:]),
    ):
        scores = {
            row.site_key: float(value)
            for row, value in zip(rows, values, strict=True)
        }
        assert_complete_model_scores(model_input, name, scores)
        partition_scores.append((name, scores))
    return ComparisonModelScores(
        model=model_input.model,
        seed=model_input.seed,
        panel_sha256=model_input.panel_sha256,
        partition_scores=tuple(partition_scores),
    )


def run_direct_comparison_roster(
    run: LiteratureRandomRun,
    *,
    sequence_features: dict[tuple[str, int], tuple[float, ...]],
    esm_features: dict[tuple[str, int], tuple[float, ...]],
    structure_features: dict[
        tuple[str, int], tuple[tuple[float, float], bool]
    ],
    parameters: dict[str, dict[str, int | float]],
    device: str = "cpu",
) -> tuple[ComparisonModelScores, ...]:
    """Execute all five direct baselines on one byte-identical panel."""
    results: list[ComparisonModelScores] = []
    for model in ("pu_logistic", "random_forest", "xgboost"):
        results.append(
            run_tabular_direct_baseline(
                bind_model_to_comparison_panel(run, model),
                sequence_features,
                parameters=parameters.get(model),
            )
        )
    results.append(
        run_tabular_direct_baseline(
            bind_model_to_comparison_panel(run, "esm_linear_head"),
            esm_features,
            parameters=parameters.get("esm_linear_head"),
        )
    )
    results.append(
        run_structure_direct_baseline(
            bind_model_to_comparison_panel(run, "structure_ranker"),
            sequence_features,
            esm_features,
            structure_features,
            parameters=parameters.get("structure_ranker"),
            device=device,
        )
    )
    if {result.panel_sha256 for result in results} != {run.panel_sha256}:
        raise RuntimeError("direct roster escaped its bound comparison panel")
    return tuple(results)
