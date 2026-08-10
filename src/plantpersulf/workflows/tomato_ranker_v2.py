"""Frozen repeated-cluster tomato PU workflow.

The runner is intentionally inert until its caller supplies an empty output
directory. It produces a development-only manifest; its fold names never
match the Gate 2 leave-study-out grammar.
"""

from __future__ import annotations

import csv
import json
import logging
import statistics
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.evaluation.metrics import (
    average_precision,
    mean_reciprocal_rank,
    recall_at_k,
)
from plantpersulf.evaluation.model_admission import admit_candidate_model
from plantpersulf.evaluation.rank_ensemble import (
    ApplicabilityEnvelope,
    summarize_rank_runs,
)
from plantpersulf.features.sequence import _load_proteome
from plantpersulf.features.site_biology import (
    BIOLOGY_FEATURE_NAMES,
    build_site_biology_vector_from_sequence,
)
from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)
from plantpersulf.models.traditional import pu_logistic_regression_scores
from plantpersulf.proteomics.tomato_candidate_registry import (
    build_tomato_candidate_registry,
)
from plantpersulf.proteomics.tomato_local_dataset import (
    ARENA_PANEL,
    ARENA_PROTEOME,
    assign_grouped_folds,
    build_tomato_pu_rows,
)
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file


def sensitivity_priors(
    n_positive: int, n_total: int, multipliers: tuple[float, ...], cap: float
) -> tuple[float, ...]:
    if n_positive <= 0 or n_total <= 0 or n_positive > n_total:
        raise ValueError("positive count must be within the training rows")
    if not 0.0 < cap <= 1.0 or not multipliers:
        raise ValueError("prior policy is invalid")
    observed = n_positive / n_total
    return tuple(
        sorted(
            {round(min(cap, observed * multiplier), 12) for multiplier in multipliers}
        )
    )


def select_release_model(candidate_admitted: bool) -> str:
    return "additive_pu" if candidate_admitted else "pu_logistic"


def _load_clusters(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("protein_accession", "cluster_id"):
            raise RuntimeError("tomato cluster file has invalid columns")
        return {row["protein_accession"]: row["cluster_id"] for row in reader}


def _matrix(rows: list[Any], proteome: dict[str, str]) -> list[list[float]]:
    return [
        list(
            build_site_biology_vector_from_sequence(
                row.protein_accession, row.cys_position, proteome[row.protein_accession]
            ).values
        )
        for row in rows
    ]


def _metrics(
    scores: list[float], labels: list[str], k: int
) -> tuple[float, float, float]:
    scored = list(zip(scores, labels, strict=True))
    ap, recall, mrr = (
        average_precision(scored),
        recall_at_k(scored, k),
        mean_reciprocal_rank(scored),
    )
    if ap is None or recall is None or mrr is None:
        raise RuntimeError("held-out fold lacks a positive")
    return recall, ap / (labels.count("positive") / len(labels)), mrr


def _lower_quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    return ordered[int((len(ordered) - 1) * fraction)]


def _percentiles(scores: tuple[float, ...]) -> tuple[float, ...]:
    """Map a score column to descending percentiles with stable ties."""
    order = sorted(range(len(scores)), key=lambda index: (-scores[index], index))
    denominator = max(1, len(scores) - 1)
    percentiles = [0.0] * len(scores)
    for rank, index in enumerate(order):
        percentiles[index] = 1.0 - rank / denominator
    return tuple(percentiles)


def configure_run_logger(log_path: Path) -> logging.Logger:
    """Create a per-run progress logger without mutating global logging."""
    logger = logging.getLogger(f"tomato_ranker_v2.{log_path.resolve()}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    if logger.handlers:
        return logger
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for handler in (
        logging.FileHandler(log_path, encoding="utf-8"),
        logging.StreamHandler(),
    ):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def _validate_config(cfg: dict[str, Any]) -> None:
    if cfg.get("claim_class") != "development_candidate_ranking_not_gate2":
        raise RuntimeError("tomato v2 cannot emit a Gate 2 claim")
    evaluation = cfg["evaluation"]
    if (
        evaluation["primary_arena"],
        evaluation["split"],
        int(evaluation["folds"]),
        int(evaluation["repetitions"]),
    ) != ("panel", "repeated_homology_cluster_cv", 5, 5):
        raise RuntimeError("primary arena or split differs from the frozen policy")
    if (
        cfg["structure"]["delta_default"] != 0
        or cfg["observation_propensity"]["enabled"] is not False
    ):
        raise RuntimeError("unadmitted feature branch enabled")
    if cfg["models"] != {
        "candidate": "additive_pu",
        "baseline": "pu_logistic",
        "fallback_on_admission_failure": True,
    }:
        raise RuntimeError(
            "model comparison or fallback policy differs from frozen config"
        )
    if tuple(cfg["features"]["core"]) != BIOLOGY_FEATURE_NAMES:
        raise RuntimeError(
            "core feature order differs from the frozen biological contract"
        )
    if cfg["compute"] != {"device": "auto", "score_batch_size": 16384}:
        raise RuntimeError("compute policy differs from the frozen config")


def run_tomato_ranker_v2(config_path: Path, output_dir: Path) -> dict[str, object]:
    """Run the pre-registered development workflow once into a new directory."""
    cfg: dict[str, Any] = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    _validate_config(cfg)
    if output_dir.exists():
        raise FileExistsError(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    logger = configure_run_logger(output_dir / "run.log")
    logger.info("phase=run_start config=%s", config_path)
    data = cfg["data"]
    xlsx, fasta, clusters = (
        Path(data["kiae271_xlsx"]),
        Path(data["reference_proteome"]),
        Path(data["cluster_file"]),
    )
    supplementary, model_inputs = (
        Path(data["supplementary_registry"]),
        Path(data["model_input_registry"]),
    )
    assert_registered_input(xlsx, supplementary)
    assert_registered_input(fasta, model_inputs)
    assert_registered_input(clusters, model_inputs)
    proteome, cluster_map = _load_proteome(fasta), _load_clusters(clusters)
    if set(proteome) - set(cluster_map):
        raise RuntimeError("registered cluster table misses proteins")
    seed = int(data["subsample_seed"])
    device = str(cfg["compute"]["device"])
    score_batch_size = int(cfg["compute"]["score_batch_size"])
    folds, repetitions = (
        int(cfg["evaluation"]["folds"]),
        int(cfg["evaluation"]["repetitions"]),
    )
    candidate_runs: dict[str, tuple[float, float, float]] = {}
    baseline_runs: dict[str, tuple[float, float, float]] = {}
    oof_rows: list[dict[str, object]] = []
    for arena in (ARENA_PANEL, ARENA_PROTEOME):
        logger.info("phase=arena_start arena=%s", arena)
        rows = build_tomato_pu_rows(
            xlsx, proteome, arena, int(data["arena_ratio"]), seed
        )
        features = _matrix(rows, proteome)
        for repetition in range(repetitions):
            logger.info(
                "phase=repetition_start arena=%s repetition=%s", arena, repetition
            )
            folded = assign_grouped_folds(rows, cluster_map, folds, seed + repetition)
            for fold in range(folds):
                logger.info(
                    "phase=fold_start arena=%s repetition=%s fold=%s",
                    arena,
                    repetition,
                    fold,
                )
                train = [i for i, row in enumerate(folded) if row.fold != fold]
                test = [i for i, row in enumerate(folded) if row.fold == fold]
                train_x, test_x = (
                    [features[i] for i in train],
                    [features[i] for i in test],
                )
                train_y, test_y = (
                    [folded[i].label for i in train],
                    [folded[i].label for i in test],
                )
                priors = sensitivity_priors(
                    train_y.count("positive"),
                    len(train_y),
                    tuple(float(x) for x in cfg["pu"]["prior_multipliers"]),
                    float(cfg["pu"]["prior_cap"]),
                )
                candidate_columns, baseline_columns = [], []
                proteins = [folded[i].protein_accession for i in train]
                for model_seed in cfg["evaluation"]["model_seeds"]:
                    baseline_columns.append(
                        pu_logistic_regression_scores(
                            train_x, train_y, test_x, int(model_seed)
                        )
                    )
                    for prior in priors:
                        candidate_columns.append(
                            fit_additive_pu_ranker(
                                train_x,
                                train_y,
                                proteins,
                                BIOLOGY_FEATURE_NAMES,
                                AdditivePuConfig(
                                    class_prior=prior,
                                    seed=int(model_seed),
                                    device=device,
                                ),
                            ).score(
                                test_x,
                                device=device,
                                batch_size=score_batch_size,
                            )
                        )
                candidate_scores = [
                    statistics.median(column[i] for column in candidate_columns)
                    for i in range(len(test))
                ]
                baseline_scores = [
                    statistics.median(column[i] for column in baseline_columns)
                    for i in range(len(test))
                ]
                key = f"{arena}_r{repetition}f{fold}"
                candidate_runs[key] = _metrics(
                    candidate_scores, test_y, int(cfg["evaluation"]["development_k"])
                )
                baseline_runs[key] = _metrics(
                    baseline_scores, test_y, int(cfg["evaluation"]["development_k"])
                )
                envelope = ApplicabilityEnvelope.fit(train_x)
                oof_rows.extend(
                    {
                        "arena": arena,
                        "run": key,
                        "protein_accession": folded[index].protein_accession,
                        "cys_position": folded[index].cys_position,
                        "label": folded[index].label,
                        "candidate_score": candidate_scores[local],
                        "baseline_score": baseline_scores[local],
                        "in_domain": envelope.contains(test_x[local]),
                    }
                    for local, index in enumerate(test)
                )
                logger.info(
                    "phase=fold_complete arena=%s repetition=%s fold=%s rows=%s",
                    arena,
                    repetition,
                    fold,
                    len(test),
                )
    panel_candidate = {
        key: value for key, value in candidate_runs.items() if key.startswith("panel_")
    }
    panel_baseline = {
        key: value for key, value in baseline_runs.items() if key.startswith("panel_")
    }
    blocks = {key: key.rsplit("f", maxsplit=1)[0] for key in panel_candidate}
    decision = admit_candidate_model(panel_candidate, panel_baseline, blocks, seed)
    rank_summaries: dict[str, list[dict[str, object]]] = {}
    for arena in (ARENA_PANEL, ARENA_PROTEOME):
        site_keys = tuple(
            sorted(
                {
                    f"{row['protein_accession']}:C{row['cys_position']}"
                    for row in oof_rows
                    if row["arena"] == arena
                }
            )
        )
        run_scores: list[tuple[float, ...]] = []
        for repetition in range(repetitions):
            by_site = {
                f"{row['protein_accession']}:C{row['cys_position']}": float(
                    cast(float, row["candidate_score"])
                )
                for row in oof_rows
                if row["arena"] == arena
                and str(row["run"]).startswith(f"{arena}_r{repetition}f")
            }
            run_scores.append(tuple(by_site[key] for key in site_keys))
        rank_summaries[arena] = [
            {
                "site_key": item.site_key,
                "median_score": item.median_score,
                "median_rank": item.median_rank,
                "best_rank": item.best_rank,
                "worst_rank": item.worst_rank,
            }
            for item in summarize_rank_runs(site_keys, tuple(run_scores))
        ]

    selected_model = select_release_model(decision.admitted)
    logger.info(
        "phase=admission_complete selected_model=%s admitted=%s",
        selected_model,
        decision.admitted,
    )
    full_rows = build_tomato_pu_rows(
        xlsx, proteome, ARENA_PANEL, int(data["arena_ratio"]), seed
    )
    full_x = _matrix(full_rows, proteome)
    full_y = [row.label for row in full_rows]
    full_proteins = [row.protein_accession for row in full_rows]
    full_priors = sensitivity_priors(
        full_y.count("positive"),
        len(full_y),
        tuple(float(value) for value in cfg["pu"]["prior_multipliers"]),
        float(cfg["pu"]["prior_cap"]),
    )
    registry = build_tomato_candidate_registry(xlsx, fasta, supplementary, model_inputs)
    logger.info("phase=candidate_registry_complete count=%s", len(registry.candidates))
    candidate_x = _matrix(list(registry.candidates), proteome)
    release_models: list[dict[str, object]] = []
    percentile_columns: list[tuple[float, ...]] = []
    if selected_model == "additive_pu":
        for model_seed in cfg["evaluation"]["model_seeds"]:
            for prior in full_priors:
                logger.info(
                    "phase=release_refit_start model=additive_pu seed=%s prior=%s",
                    model_seed,
                    prior,
                )
                model = fit_additive_pu_ranker(
                    full_x,
                    full_y,
                    full_proteins,
                    BIOLOGY_FEATURE_NAMES,
                    AdditivePuConfig(
                        class_prior=prior,
                        seed=int(model_seed),
                        device=device,
                    ),
                )
                release_models.append(
                    {
                        "seed": int(model_seed),
                        "class_prior": prior,
                        "model": model.to_dict(),
                    }
                )
                percentile_columns.append(
                    _percentiles(
                        model.score(
                            candidate_x,
                            device=device,
                            batch_size=score_batch_size,
                        )
                    )
                )
                logger.info(
                    "phase=release_refit_complete model=additive_pu seed=%s prior=%s",
                    model_seed,
                    prior,
                )
    else:
        for model_seed in cfg["evaluation"]["model_seeds"]:
            logger.info(
                "phase=release_refit_start model=pu_logistic seed=%s", model_seed
            )
            scores = pu_logistic_regression_scores(
                full_x, full_y, candidate_x, int(model_seed)
            )
            release_models.append(
                {
                    "seed": int(model_seed),
                    "model_type": "pu_logistic",
                    "refit_recipe": "traditional.pu_logistic_regression_scores",
                }
            )
            percentile_columns.append(_percentiles(tuple(scores)))
            logger.info(
                "phase=release_refit_complete model=pu_logistic seed=%s", model_seed
            )
    release_envelope = ApplicabilityEnvelope.fit(full_x)
    candidate_score_rows = [
        {
            "site_key": candidate.site_key,
            "protein_accession": candidate.protein_accession,
            "cys_position": candidate.cys_position,
            "median_percentile": statistics.median(
                column[index] for column in percentile_columns
            ),
            "conservative_percentile": _lower_quantile(
                [column[index] for column in percentile_columns], 0.2
            ),
            "in_domain": release_envelope.contains(candidate_x[index]),
        }
        for index, candidate in enumerate(registry.candidates)
    ]
    model_release_created = bool(percentile_columns and candidate_score_rows)
    if not model_release_created:
        raise RuntimeError("selected model produced no frozen candidate scores")
    summary: dict[str, object] = {
        "experiment": cfg["experiment"]["name"],
        "claim_class": cfg["claim_class"],
        "candidate_model_admitted": decision.admitted,
        "admission_reason": decision.reason,
        "admission_paired_intervals": decision.paired_intervals,
        "candidate_runs": candidate_runs,
        "baseline_runs": baseline_runs,
        "structure_delta": 0,
        "selected_model": selected_model,
        "model_release_created": model_release_created,
        "candidate_registry_count": len(registry.candidates),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "oof_scores.json").write_text(
        json.dumps(oof_rows, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "rank_summaries.json").write_text(
        json.dumps(rank_summaries, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "frozen_models.json").write_text(
        json.dumps(release_models, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "candidate_scores.json").write_text(
        json.dumps(candidate_score_rows, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "complete": True,
        "claim_class": cfg["claim_class"],
        "config_sha256": hash_file(config_path, "sha256"),
        "kiae271_sha256": hash_file(xlsx, "sha256"),
        "proteome_sha256": hash_file(fasta, "sha256"),
        "clusters_sha256": hash_file(clusters, "sha256"),
        "selected_model": selected_model,
        "model_release_created": model_release_created,
        "frozen_models_sha256": hash_file(output_dir / "frozen_models.json", "sha256"),
        "candidate_scores_sha256": hash_file(
            output_dir / "candidate_scores.json", "sha256"
        ),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    logger.info("phase=run_complete output_dir=%s", output_dir)
    return summary
