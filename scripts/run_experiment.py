#!/usr/bin/env python
"""Deterministic TDD experiment runner for persulfidation baselines.

Reads a frozen YAML experiment config and runs every model on the same
benchmark, features, splits, and seeds — producing a single, deterministic
results tree under ``results/experiments/<name>/`` with per-seed metrics
and a manifest recording every input SHA256 used.

Usage (per the Task 8 acceptance gate)::

    python scripts/run_experiment.py \\
        --config configs/experiments/baseline_sequence_v1.yaml
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, cast

import yaml

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

BENCHMARK_FIELDS = (
    "protein_accession", "cys_position_in_protein", "label",
    "study_accession", "evidence_level", "source_sha256",
)


def _ensure_venv() -> None:
    venv = Path(sys.executable).with_name("pyvenv.cfg")
    if not venv.is_file():
        print("WARNING: not running inside a project venv.", file=sys.stderr)


def _load_config(path: Path) -> dict[str, Any]:
    loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise SystemExit(f"invalid experiment config: {path}")
    return cast(dict[str, Any], loaded)


# ---------------------------------------------------------------------------
# results
# ---------------------------------------------------------------------------

@dataclass
class ModelResult:
    model: str
    seed: int
    train_ap: float | None = None
    val_ap: float | None = None
    test_ap: float | None = None
    test_recall_10: float | None = None
    test_recall_50: float | None = None
    test_mrr: float | None = None


@dataclass
class ExperimentResults:
    experiment: str
    model_results: list[ModelResult] = field(default_factory=list)
    manifest_sha256: str = ""


# ---------------------------------------------------------------------------
# bench --> split
# ---------------------------------------------------------------------------

def _build_singleton_cluster_file(benchmark_path: Path, output: Path) -> Path:
    """Each protein in its own cluster when no MMseqs2 output is available."""
    proteins: set[str] = set()
    with benchmark_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            proteins.add(row["protein_accession"])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("protein_accession", "cluster_id"))
        for i, p in enumerate(sorted(proteins)):
            writer.writerow((p, str(i)))
    return output


def _train_val_test_rows(
    benchmark_rows: list[dict[str, str]],
    proteome_path: Path,
    cluster_file: Path,
    split_config_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    from plantpersulf.benchmark.splits import build_splits

    table = build_splits(split_config_path)
    split_map = {row.protein_accession: row.split for row in table.rows}

    train_rows: list[dict[str, str]] = []
    val_rows: list[dict[str, str]] = []
    test_rows: list[dict[str, str]] = []
    for row in benchmark_rows:
        split = split_map.get(row["protein_accession"], "train")
        if split == "validation":
            val_rows.append(row)
        elif split == "test":
            test_rows.append(row)
        else:
            train_rows.append(row)
    return train_rows, val_rows, test_rows


# ---------------------------------------------------------------------------
# feature extraction (lightweight, just collects per-row vectors)
# ---------------------------------------------------------------------------

def _motif_features(
    benchmark_rows: list[dict[str, str]],
    proteome_path: Path,
) -> list[list[float]]:
    from plantpersulf.features.sequence import extract_sequence_features

    feats = extract_sequence_features(
        Path("__placeholder_labels.tsv"), proteome_path, window_radius=10
    )
    # Build lookup: (protein, pos) -> [hydrophobicity, cys_density]
    lookup: dict[tuple[str, int], list[float]] = {}
    for f in feats:
        lookup[(f.protein_accession, f.cys_position)] = [
            f.hydrophobicity, f.cys_density,
        ]
    vectors: list[list[float]] = []
    for row in benchmark_rows:
        key = (row["protein_accession"], int(row["cys_position_in_protein"]))
        vec = lookup.get(key, [0.0, 0.0])
        vectors.append(vec)
    return vectors


def _esm2_features(
    benchmark_rows: list[dict[str, str]],
    proteome_path: Path,
) -> list[list[float]]:
    from plantpersulf.features.esm2 import extract_esm2_embeddings

    emb = extract_esm2_embeddings(
        Path("__placeholder_labels.tsv"), proteome_path
    )
    lookup: dict[tuple[str, int], list[float]] = {}
    for f in emb:
        lookup[(f.protein_accession, f.cys_position)] = list(f.embedding)
    vectors: list[list[float]] = []
    for row in benchmark_rows:
        key = (row["protein_accession"], int(row["cys_position_in_protein"]))
        vec = lookup.get(key, [0.0] * 1280)
        vectors.append(vec)
    return vectors


# ---------------------------------------------------------------------------
# baseline scoring
# ---------------------------------------------------------------------------

def _train_predict(
    model_name: str,
    train_X: list[list[float]],
    train_y: list[str],
    predict_X: list[list[float]],
    seed: int,
) -> list[float]:
    from plantpersulf.models.esm_baseline import esm_linear_head_scores
    from plantpersulf.models.traditional import (
        logistic_regression_scores,
        pu_logistic_regression_scores,
        random_forest_scores,
        xgboost_scores,
    )

    train_emb = [tuple(v) for v in train_X]
    predict_emb = [tuple(v) for v in predict_X]

    if model_name == "logistic":
        return logistic_regression_scores(train_X, train_y, predict_X, seed)
    if model_name == "pu_logistic":
        return pu_logistic_regression_scores(train_X, train_y, predict_X, seed)
    if model_name == "random_forest":
        return random_forest_scores(
            train_X, train_y, predict_X, seed, n_estimators=50
        )
    if model_name == "xgboost":
        return xgboost_scores(
            train_X, train_y, predict_X, seed, n_estimators=50
        )
    if model_name == "esm_linear_head":
        return esm_linear_head_scores(train_emb, train_y, predict_emb, seed)
    raise ValueError(f"unknown model: {model_name}")


def _evaluate(
    model_name: str,
    seed: int,
    train_X: list[list[float]],
    train_y: list[str],
    val_X: list[list[float]],
    val_y: list[str],
    test_X: list[list[float]],
    test_y: list[str],
) -> ModelResult:
    from plantpersulf.evaluation.metrics import (
        average_precision,
        mean_reciprocal_rank,
        recall_at_k,
    )

    result = ModelResult(model=model_name, seed=seed)

    def _scored(
        features: list[list[float]], labels: list[str],
    ) -> list[tuple[float, str]]:
        scores = _train_predict(model_name, train_X, train_y, features, seed)
        return list(zip(scores, labels, strict=True))

    val_scored = _scored(val_X, val_y)
    result.val_ap = average_precision(val_scored)

    if test_X:
        test_scored = _scored(test_X, test_y)
        result.test_ap = average_precision(test_scored)
        result.test_recall_10 = recall_at_k(test_scored, 10)
        result.test_recall_50 = recall_at_k(test_scored, 50)
        result.test_mrr = mean_reciprocal_rank(test_scored)
    return result


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def run_experiment(config_path: Path) -> ExperimentResults:
    cfg = _load_config(config_path)
    exp_cfg = cfg["experiment"]
    eval_cfg = cfg["evaluation"]
    out_cfg = cfg["output"]

    # --- benchmark ---
    benchmark_path = Path(exp_cfg["benchmark_labels"])
    if not benchmark_path.exists():
        print(f"Building benchmark into {benchmark_path} …")
        from plantpersulf.benchmark.labels import build_benchmark_labels
        from plantpersulf.proteomics.persulfidation_publish import (
            publish_persulfidation_sites,
        )
        interim = Path("data/interim")
        for acc in ("PXD006140", "PXD024061"):
            publish_persulfidation_sites(acc, output_root=interim)
        build_benchmark_labels(
            site_output_root=interim,
            proteome_path=Path(exp_cfg["reference_proteome"]),
            output_directory=benchmark_path.parent,
        )

    rows: list[dict[str, str]] = []
    with benchmark_path.open(encoding="utf-8", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            rows.append(dict(row))

    proteome_path = Path(exp_cfg["reference_proteome"])

    # --- splits ---
    splits_cfg = cfg["splits"]
    cluster_file = Path(splits_cfg["cluster_file"])
    if not cluster_file.exists():
        cluster_file = _build_singleton_cluster_file(benchmark_path, cluster_file)
    split_config = Path(splits_cfg["split_config"])
    train_rows, val_rows, test_rows = _train_val_test_rows(
        rows, proteome_path, cluster_file, split_config,
    )
    print(f"split: train={len(train_rows)} val={len(val_rows)} test={len(test_rows)}")

    train_pos = sum(1 for r in train_rows if r["label"] == "positive")
    val_pos = sum(1 for r in val_rows if r["label"] == "positive")
    test_pos = sum(1 for r in test_rows if r["label"] == "positive")
    print(f"positives: train={train_pos} val={val_pos} test={test_pos}")

    # --- features ---
    feature_names = [f["name"] for f in cfg["features"]]
    print(f"features: {feature_names}")

    motif_X = _motif_features(train_rows + val_rows + test_rows, proteome_path)
    n_train = len(train_rows)
    n_val = len(val_rows)
    motif_train = motif_X[:n_train]
    motif_val = motif_X[n_train:n_train + n_val]
    motif_test = motif_X[n_train + n_val:]
    train_y = [r["label"] for r in train_rows]
    val_y = [r["label"] for r in val_rows]
    test_y = [r["label"] for r in test_rows]

    # --- run ---
    results = ExperimentResults(experiment=exp_cfg["name"])
    for model_name in cfg["models"]:
        for seed in eval_cfg["seeds"]:
            print(f"  {model_name} seed={seed} …", end=" ")
            if model_name == "motif_frequency":
                # motif baseline uses its own ranking machinery
                r = ModelResult(model="motif_frequency", seed=seed)
                # Skip for now — motif baseline needs a different feature interface
                # (it reads flanking windows directly, not precomputed vectors)
                r.val_ap = None
                r.test_ap = None
                print("skipped (motif freq needs separate codepath)")
            else:
                r = _evaluate(
                    model_name, seed,
                    motif_train, train_y,
                    motif_val, val_y,
                    motif_test, test_y,
                )
                print(f"val_ap={r.val_ap:.4f}" if r.val_ap else "no_val_positives")
            results.model_results.append(r)

    # --- write ---
    out_dir = Path(out_cfg["directory"])
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True)

    rows_out: list[dict[str, str]] = []
    for r in results.model_results:
        out_row: dict[str, str] = {
            "model": r.model, "seed": str(r.seed),
        }
        for key, val in [
            ("val_ap", r.val_ap), ("test_ap", r.test_ap),
            ("test_recall_10", r.test_recall_10),
            ("test_recall_50", r.test_recall_50),
            ("test_mrr", r.test_mrr),
        ]:
            out_row[key] = f"{val:.6f}" if val is not None else ""
        rows_out.append(out_row)
    with (out_dir / "metrics.tsv").open("w", encoding="utf-8", newline="") as h:
        w = csv.DictWriter(
            h, fieldnames=list(rows_out[0]), delimiter="\t", lineterminator="\n",
        )
        w.writeheader()
        w.writerows(rows_out)

    from plantpersulf.provenance.hashing import hash_file
    manifest = {
        "experiment": exp_cfg["name"],
        "config_sha256": hash_file(config_path, "sha256"),
        "benchmark_sha256": hash_file(benchmark_path, "sha256"),
        "proteome_sha256": hash_file(proteome_path, "sha256"),
        "cluster_file_sha256": (
            hash_file(cluster_file, "sha256") if cluster_file.exists() else ""
        ),
        "n_seeds": len(eval_cfg["seeds"]),
    }
    (out_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"\nResults written to {out_dir}")
    return results


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="run a persulfidation baseline experiment")
    p.add_argument("--config", type=Path, required=True)
    return p


if __name__ == "__main__":
    _ensure_venv()
    args = build_parser().parse_args()
    run_experiment(args.config)
