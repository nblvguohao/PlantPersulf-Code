#!/usr/bin/env python
"""P2 — regenerate per-site scored predictions for a frozen release.

Deterministic (CPU) re-scoring of the frozen ``pu_ranker_v1`` release for the
Gate-2 statistical wiring: for each leave-study-out fold, trains the release
ablation, the structure-ablated ablation, and the PU baseline on the fold's
training rows — using exactly the same partition, subsample, and train/val
split as ``scripts/run_experiment.py`` — and writes per-site test scores for
every seed. Seed-ensembling and all paired statistics live in
``scripts/validate_external.py``; this script only produces aligned scores.

The default arms (``seq_structure`` / ``sequence_only`` / ``pu_logistic``)
are ESM-free, so no GPU or ESM weights are needed and the run finishes in
minutes on CPU. Scores are bit-reproducible for a fixed device class (the
ranker defaults to CPU); they are a fresh, self-consistent measurement whose
point APs should match ``metrics.tsv`` within device floating-point tolerance.

Usage::

    python scripts/score_release.py \
        --config configs/experiments/pu_ranker_v1.yaml \
        --output-dir results/external_validation/pu_ranker_v1/scored
"""

from __future__ import annotations

import argparse
import csv
import random
import shutil
import tempfile
from pathlib import Path
from typing import Any

BENCHMARK_FIELDS = (
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "evidence_level",
    "source_sha256",
)

SCORED_FIELDS = (
    "fold",
    "seed",
    "protein_accession",
    "cys_position",
    "label",
    "cluster_id",
    "has_structure",
    "score",
)


def _read_benchmark(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != BENCHMARK_FIELDS:
            raise RuntimeError(f"benchmark labels have invalid columns: {path}")
        return [dict(r) for r in reader]


def _read_clusters(path: Path) -> dict[str, str]:
    """protein_accession -> cluster_id; missing proteins become singletons
    downstream so they are never lumped into one fake cluster."""
    mapping: dict[str, str] = {}
    if not path.is_file():
        print(f"NOTE: cluster file {path} not found — singleton clusters.")
        return mapping
    with path.open(encoding="utf-8", newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            mapping[r["protein_accession"]] = r["cluster_id"]
    return mapping


def _cluster_of(mapping: dict[str, str], accession: str) -> str:
    return mapping.get(accession, f"__singleton__{accession}")


def _ablation_named(cfg: dict[str, Any], name: str) -> Any:
    from run_experiment import _ablation_from_spec  # type: ignore[import-not-found]

    for spec in cfg.get("ablations") or []:
        if spec.get("name") == name:
            return _ablation_from_spec(spec)
    raise ValueError(f"ablation {name!r} not found in config")


def _write_scored(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(SCORED_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def score_release(
    config_path: Path,
    clusters_path: Path,
    release_ablation: str,
    ablated_ablation: str,
    baseline_model: str,
    output_dir: Path,
) -> None:
    from run_experiment import (  # type: ignore[import-not-found]
        _build_branch_features,
        _load_config,
        _sequence_feature_vectors,
        _study_fold_rows,
        _subsample_unlabeled,
        _train_predict,
    )

    from plantpersulf.evaluation.metrics import average_precision
    from plantpersulf.models.structure_ranker import structure_ranker_scores

    cfg = _load_config(config_path)
    exp = cfg["experiment"]
    all_rows = _read_benchmark(Path(exp["benchmark_labels"]))
    proteome_path = Path(exp["reference_proteome"])
    cluster_map = _read_clusters(clusters_path)

    ratio = int(cfg.get("subsample", {}).get("unlabeled_per_positive", 50))
    sub_seed = int(cfg.get("subsample", {}).get("seed", 12345))
    seeds = [int(s) for s in cfg["evaluation"]["seeds"]]
    ranker_params = cfg.get("ranker", {})

    ab_release = _ablation_named(cfg, release_ablation)
    ab_ablated = _ablation_named(cfg, ablated_ablation)
    need_esm = bool(ab_release.use_esm or ab_ablated.use_esm)
    if need_esm:
        print("NOTE: an arm uses the ESM branch — extraction will be slow.")

    arms: dict[str, list[dict[str, Any]]] = {
        "model": [],
        "ablated": [],
        "baseline": [],
    }

    for study in cfg["splits"]["studies"]:
        print(f"\n===== Fold: leave_{study}_out =====")
        train_rows, test_rows = _study_fold_rows(all_rows, study)
        train_rows = _subsample_unlabeled(train_rows, ratio, sub_seed)
        test_rows = _subsample_unlabeled(test_rows, ratio, sub_seed + 1)

        # Identical train/val split to the runner (val rows are excluded
        # from training; val scores themselves are not needed here).
        rng = random.Random(sub_seed)
        train_dedup = train_rows[:]
        rng.shuffle(train_dedup)
        n_val = max(1, int(len(train_dedup) * 0.2))
        train_rows_fold = train_dedup[n_val:]
        train_y = [r["label"] for r in train_rows_fold]
        print(f"  train={len(train_rows_fold)} test={len(test_rows)}")

        scratch = Path(tempfile.mkdtemp(prefix="score_release_"))
        try:
            branch_train = _build_branch_features(
                train_rows_fold, proteome_path, scratch, "train", need_esm
            )
            branch_test = _build_branch_features(
                test_rows, proteome_path, scratch, "test", need_esm
            )
            feat_train = _sequence_feature_vectors(
                train_rows_fold, proteome_path, scratch, "train_flat"
            )
            feat_test = _sequence_feature_vectors(
                test_rows, proteome_path, scratch, "test_flat"
            )
        finally:
            shutil.rmtree(scratch, ignore_errors=True)

        def _emit(
            arm: str,
            seed: int,
            scores: list[float],
            rows: list[dict[str, str]] = test_rows,
            mask: list[bool] = branch_test.structure_mask,
            fold: str = study,
        ) -> None:
            for row, score, has_struct in zip(rows, scores, mask, strict=True):
                acc = row["protein_accession"]
                arms[arm].append(
                    {
                        "fold": fold,
                        "seed": seed,
                        "protein_accession": acc,
                        "cys_position": int(row["cys_position_in_protein"]),
                        "label": row["label"],
                        "cluster_id": _cluster_of(cluster_map, acc),
                        "has_structure": "1" if has_struct else "0",
                        "score": f"{score:.10g}",
                    }
                )

        for seed in seeds:
            for arm, ab in (("model", ab_release), ("ablated", ab_ablated)):
                out = structure_ranker_scores(
                    branch_train,
                    train_y,
                    branch_test,
                    seed=seed,
                    ablation=ab,
                    hidden=int(ranker_params.get("hidden", 16)),
                    dropout=float(ranker_params.get("dropout", 0.2)),
                    epochs=int(ranker_params.get("epochs", 200)),
                    lr=float(ranker_params.get("lr", 0.05)),
                    n_mc_dropout=int(ranker_params.get("n_mc_dropout", 16)),
                )
                _emit(arm, seed, list(out.scores))
                ap = average_precision(
                    list(zip(out.scores, [r["label"] for r in test_rows], strict=True))
                )
                print(f"  {arm} seed={seed} test_ap={ap:.4f}")

            base_scores = _train_predict(
                baseline_model, feat_train, train_y, feat_test, seed
            )
            _emit("baseline", seed, base_scores)
            ap_b = average_precision(
                list(zip(base_scores, [r["label"] for r in test_rows], strict=True))
            )
            print(f"  baseline[{baseline_model}] seed={seed} test_ap={ap_b:.4f}")

    names = {
        "model": release_ablation,
        "ablated": ablated_ablation,
        "baseline": baseline_model,
    }
    for arm, rows in arms.items():
        _write_scored(output_dir / f"{arm}.tsv", rows)
        print(f"wrote {len(rows)} rows -> {output_dir / f'{arm}.tsv'} ({names[arm]})")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="per-site score regeneration for Gate-2 statistics"
    )
    p.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/pu_ranker_v1.yaml"),
    )
    p.add_argument(
        "--clusters",
        type=Path,
        default=Path("data/processed/clusters/protein_clusters_v1.tsv"),
    )
    p.add_argument("--release-ablation", default="seq_structure")
    p.add_argument("--ablated-ablation", default="sequence_only")
    p.add_argument("--baseline", default="pu_logistic")
    p.add_argument("--output-dir", type=Path, required=True)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    score_release(
        config_path=args.config,
        clusters_path=args.clusters,
        release_ablation=str(args.release_ablation),
        ablated_ablation=str(args.ablated_ablation),
        baseline_model=str(args.baseline),
        output_dir=args.output_dir,
    )
