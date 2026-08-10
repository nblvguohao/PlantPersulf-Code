#!/usr/bin/env python
"""Tomato-local model v1: within-species training & evaluation on kiae271.

Builds a tomato-specific candidate-prioritisation model trained and evaluated
entirely within Solanum lycopersicum using the 99 coordinate-verified kiae271
persulfidation sites as supervision (Zhang et al. 2024, Plant Physiology,
doi:10.1093/plphys/kiae271). Grouped homology-cluster 5-fold cross-validation
with three arms: seq_2 (hydrophobicity + cys_density), seq_3 (+ local positive
charge density), seq_structure (+ AlphaFold contact number & pLDDT).

Results go to results/tomato_local_v1/ (separate namespace — does not touch
any results/external_validation/ or results/experiments/ artifact).

Usage::

    python scripts/run_tomato_local_v1.py
"""

from __future__ import annotations

import csv
import json
import random
import shutil
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

KIAE271_XLSX = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
TOMATO_PROTEOME = Path("data/raw/references/tomato_ref_proteome_v1.fasta")
CLUSTER_FILE = Path("data/processed/clusters/tomato_proteome_clusters_v1.tsv")
OUTPUT_DIR = Path("results/tomato_local_v1")
CONFIG = Path("configs/experiments/tomato_local_v1.yaml")

SEED_PROTEOME = 12345
SEED_PANEL = 12346
SEED_FOLDS = 20260810
MODEL_SEEDS = (0, 1, 2, 3, 4)
RATIO = 20
N_FOLDS = 5

BENCHMARK_FIELDS = (
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "evidence_level",
    "source_sha256",
)

ARM_NAMES = ("seq_2", "seq_3", "seq_structure")
ARENAS = ("proteome", "panel")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _sha256(path: Path) -> str:
    import hashlib

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_labels_tsv(rows: list[dict[str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(BENCHMARK_FIELDS),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# sequence features (1-3 dim)
# ---------------------------------------------------------------------------


def _seq_feature_vectors(
    rows: list[tuple[str, int, str]],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
    fields: tuple[str, ...],
) -> list[list[float]]:
    """Extract per-Cys sequence features via the real feature extractor.

    ``fields`` selects which SequenceFeatureRow attributes to use, in order.
    Valid values: ``hydrophobicity``, ``cys_density``,
    ``local_positive_charge_density``.
    """
    from plantpersulf.features.sequence import extract_sequence_features

    benchmark_rows: list[dict[str, str]] = [
        {
            "protein_accession": acc,
            "cys_position_in_protein": str(pos),
            "label": lab,
            "study_accession": "kiae271",
            "evidence_level": "site_biochemical",
            "source_sha256": "tomato_local_v1",
        }
        for acc, pos, lab in rows
    ]
    labels_path = scratch_dir / f"{tag}_labels.tsv"
    _write_labels_tsv(benchmark_rows, labels_path)
    feats = extract_sequence_features(labels_path, proteome_path, window_radius=10)
    lookup: dict[tuple[str, int], list[float]] = {}
    for f in feats:
        values = []
        for field in fields:
            if field == "hydrophobicity":
                values.append(f.hydrophobicity)
            elif field == "cys_density":
                values.append(f.cys_density)
            elif field == "local_positive_charge_density":
                values.append(f.local_positive_charge_density)
            else:
                values.append(0.0)
        lookup[(f.protein_accession, f.cys_position)] = values
    default = [0.0] * len(fields)
    return [
        lookup.get((acc, pos), default) for acc, pos, _ in rows
    ]


# ---------------------------------------------------------------------------
# structure features (2 dim: contact_number_proxy, plddt)
# ---------------------------------------------------------------------------


def _structure_feature_vectors(
    rows: list[tuple[str, int, str]],
    scratch_dir: Path,
    tag: str,
) -> tuple[list[list[float]], list[bool]]:
    """Extract [contact_number_proxy, plddt] per Cys from registered AlphaFold
    structures. Rows without a structure get [0.0, 0.0] and mask=False."""
    from plantpersulf.download.alphafold import audit_alphafold_structures
    from plantpersulf.features.structure import extract_cys_structure_features

    sources = audit_alphafold_structures()
    pdb_cache: dict[str, str] = {}
    for s in sources:
        pdb_cache[s.accession] = s.local_path.read_text(encoding="utf-8")

    vectors: list[list[float]] = []
    mask: list[bool] = []
    for acc, pos, _ in rows:
        pdb = pdb_cache.get(acc)
        if pdb is None:
            vectors.append([0.0, 0.0])
            mask.append(False)
            continue
        feats = extract_cys_structure_features(pdb, acc, [pos])
        f = feats[0]
        if not f.has_structure:
            vectors.append([0.0, 0.0])
            mask.append(False)
        else:
            contact = f.contact_number_proxy or 0.0
            plddt = f.plddt or 0.0
            vectors.append([contact, plddt])
            mask.append(True)
    return vectors, mask


# ---------------------------------------------------------------------------
# model fitting & scoring
# ---------------------------------------------------------------------------


def _fit_pu_logistic(
    train_X: list[list[float]],
    train_y: list[str],
    predict_X: list[list[float]],
    seed: int,
) -> list[float]:
    from plantpersulf.models.traditional import pu_logistic_regression_scores

    return pu_logistic_regression_scores(train_X, train_y, predict_X, seed)


# ---------------------------------------------------------------------------
# per-arm evaluation across folds
# ---------------------------------------------------------------------------


def _evaluate_one_arm(
    folded_rows: list[Any],  # TomatoPuRow with .fold, .label, .cys_position, etc.
    proteome_path: Path,
    arm_name: str,
    scratch_dir: Path,
) -> dict[str, Any]:
    """Run 5-fold grouped CV for one arm, one arena. Returns per-fold APs and
    seed-ensembled scores."""
    from plantpersulf.evaluation.metrics import (
        average_precision,
        mean_reciprocal_rank,
        recall_at_k,
    )

    # Convert TomatoPuRow -> (acc, pos, label) tuples
    all_tuples: list[tuple[str, int, str]] = [
        (r.protein_accession, r.cys_position, r.label) for r in folded_rows
    ]

    seq_fields: tuple[str, ...] = ()
    if arm_name == "seq_2":
        seq_fields = ("hydrophobicity", "cys_density")
    elif arm_name == "seq_3":
        seq_fields = ("hydrophobicity", "cys_density", "local_positive_charge_density")
    elif arm_name == "seq_structure":
        seq_fields = ("hydrophobicity", "cys_density", "local_positive_charge_density")

    per_fold: list[dict[str, Any]] = []
    per_seed_scores: list[list[float]] = []

    for seed in MODEL_SEEDS:
        seed_scores: list[float] = [float("nan")] * len(all_tuples)
        for fold_idx in range(N_FOLDS):
            train_idx = [i for i, r in enumerate(folded_rows) if r.fold != fold_idx]
            test_idx = [i for i, r in enumerate(folded_rows) if r.fold == fold_idx]
            train_rows = [all_tuples[i] for i in train_idx]
            test_rows = [all_tuples[i] for i in test_idx]
            train_y = [folded_rows[i].label for i in train_idx]
            test_y = [folded_rows[i].label for i in test_idx]

            # sequence features
            train_seq = _seq_feature_vectors(
                train_rows, proteome_path, scratch_dir,
                f"train_{arm_name}_fold{fold_idx}_seed{seed}", seq_fields,
            )
            test_seq = _seq_feature_vectors(
                test_rows, proteome_path, scratch_dir,
                f"test_{arm_name}_fold{fold_idx}_seed{seed}", seq_fields,
            )

            # structure features (if arm includes them)
            train_X: list[list[float]]
            test_X: list[list[float]]
            if arm_name == "seq_structure":
                train_struct, _ = _structure_feature_vectors(
                    train_rows, scratch_dir,
                    f"train_struct_fold{fold_idx}_seed{seed}",
                )
                test_struct, test_struct_mask = _structure_feature_vectors(
                    test_rows, scratch_dir,
                    f"test_struct_fold{fold_idx}_seed{seed}",
                )
                train_X = [s + u for s, u in zip(train_seq, train_struct, strict=True)]
                test_X = [s + u for s, u in zip(test_seq, test_struct, strict=True)]
            else:
                train_X = train_seq
                test_X = test_seq

            scores = _fit_pu_logistic(train_X, train_y, test_X, seed)
            for i, s in zip(test_idx, scores):
                seed_scores[i] = s

            scored = list(zip(scores, test_y, strict=True))
            ap = average_precision(scored)
            if fold_idx == 0 and seed == MODEL_SEEDS[0]:
                # record one set of per-fold metrics (using first seed only,
                # to avoid 15 entries per arena.ensemble uses all seeds.)
                pass

        per_seed_scores.append(seed_scores)

    # Ensemble: mean score per row across model seeds
    n_rows = len(all_tuples)
    ens = [
        sum(per_seed_scores[s][i] for s in range(len(MODEL_SEEDS)))
        / len(MODEL_SEEDS)
        for i in range(n_rows)
    ]

    # Compute per-fold metrics using ensembled scores
    ens_by_fold: list[dict[str, Any]] = []
    for fold_idx in range(N_FOLDS):
        test_idx = [i for i, r in enumerate(folded_rows) if r.fold == fold_idx]
        fold_y = [folded_rows[i].label for i in test_idx]
        fold_scores = [ens[i] for i in test_idx]
        scored = list(zip(fold_scores, fold_y, strict=True))
        ap = average_precision(scored)
        r10 = recall_at_k(scored, 10)
        r50 = recall_at_k(scored, 50)
        mrr = mean_reciprocal_rank(scored)
        n_pos = sum(1 for y in fold_y if y == "positive")
        ens_by_fold.append({
            "fold": fold_idx,
            "n_rows": len(test_idx),
            "n_positive": n_pos,
            "ap": ap,
            "recall_at_10": r10,
            "recall_at_50": r50,
            "mrr": mrr,
        })

    all_scored = list(zip(ens, [r.label for r in folded_rows], strict=True))
    ap_all = average_precision(all_scored)
    r10_all = recall_at_k(all_scored, 10)
    r50_all = recall_at_k(all_scored, 50)

    return {
        "per_fold": ens_by_fold,
        "ensemble_ap": ap_all,
        "ensemble_recall_at_10": r10_all,
        "ensemble_recall_at_50": r50_all,
    }


# ---------------------------------------------------------------------------
# bootstrap CI (cluster-level)
# ---------------------------------------------------------------------------


def _cluster_bootstrap_ci(
    ens_scores: list[float],
    labels: list[str],
    cluster_ids: list[str],
    n_boot: int = 1000,
    seed: int = 0,
) -> dict[str, float]:
    from plantpersulf.evaluation.bootstrap import cluster_bootstrap_ci
    from plantpersulf.evaluation.metrics import average_precision

    scored: list[tuple[float, str, str]] = [
        (s, lab, cid) for s, lab, cid in zip(ens_scores, labels, cluster_ids, strict=True)
    ]
    result = cluster_bootstrap_ci(scored, average_precision, n_boot=n_boot, seed=seed)
    return {"point": result.point, "lower": result.lower, "upper": result.upper}


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def run_tomato_local_v1(
    output_dir: Path = OUTPUT_DIR,
) -> dict[str, Any]:
    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites
    from plantpersulf.proteomics.tomato_local_dataset import (
        ARENA_PANEL,
        ARENA_PROTEOME,
        assign_grouped_folds,
        build_tomato_pu_rows,
        read_cluster_map,
    )

    # 1. verify inputs
    for path in (KIAE271_XLSX, TOMATO_PROTEOME, CLUSTER_FILE):
        if not path.is_file():
            raise RuntimeError(f"required input missing: {path}")

    output_dir.mkdir(parents=True, exist_ok=True)

    # 2. load
    proteome = _load_proteome(TOMATO_PROTEOME)
    table = parse_kiae271_sites(KIAE271_XLSX, proteome)
    cluster_map = read_cluster_map(CLUSTER_FILE)
    print(
        f"tomato proteome: {len(proteome)} proteins, "
        f"{sum(s.count('C') for s in proteome.values())} Cys"
    )
    print(
        f"kiae271: {table.total_verified_sites} sites / {table.total_verified_proteins} proteins"
    )

    # 3. build PU rows for both arenas
    arena_seeds = {"proteome": SEED_PROTEOME, "panel": SEED_PANEL}
    arena_data: dict[str, list[Any]] = {}
    arena_summary: dict[str, dict[str, Any]] = {}

    for arena in ARENAS:
        rows = build_tomato_pu_rows(
            KIAE271_XLSX, proteome, arena, ratio=RATIO, seed=arena_seeds[arena],
        )
        folded = assign_grouped_folds(rows, cluster_map, n_folds=N_FOLDS, seed=SEED_FOLDS)
        n_pos = sum(1 for r in folded if r.label == "positive")
        n_unl = sum(1 for r in folded if r.label == "unlabeled")
        n_prots = len({r.protein_accession for r in folded})
        n_clusters = len({r.cluster_id for r in folded})
        base_rate = n_pos / (n_pos + n_unl)

        per_fold_pos = Counter(r.fold for r in folded if r.label == "positive")
        print(
            f"{arena}: {len(folded)} rows (+{n_pos}/-{n_unl}), "
            f"{n_prots} proteins, {n_clusters} clusters, "
            f"base_rate={base_rate:.4f}, "
            f"fold_pos={dict(per_fold_pos)}"
        )
        arena_data[arena] = folded
        arena_summary[arena] = {
            "seed": arena_seeds[arena],
            "ratio": RATIO,
            "n_rows": len(folded),
            "n_positive": n_pos,
            "n_unlabeled": n_unl,
            "n_proteins": n_prots,
            "n_clusters": n_clusters,
            "base_rate": base_rate,
            "per_fold_positive": {str(k): v for k, v in per_fold_pos.items()},
        }

    # 4. run arms
    results: dict[str, dict[str, dict[str, Any]]] = {}
    scratch = Path(tempfile.mkdtemp(prefix="tomato_local_v1_"))

    try:
        for arena in ARENAS:
            results[arena] = {}
            folded = arena_data[arena]
            for arm in ARM_NAMES:
                print(f"\n=== {arena} / {arm} ===")
                tag = f"{arena}_{arm}"
                arm_result = _evaluate_one_arm(
                    folded, TOMATO_PROTEOME, arm, scratch / tag,
                )
                # Re-extract ensembled scores for bootstrap
                all_tuples = [
                    (r.protein_accession, r.cys_position, r.label) for r in folded
                ]
                seq_fields: tuple[str, ...]
                if arm == "seq_2":
                    seq_fields = ("hydrophobicity", "cys_density")
                else:
                    seq_fields = ("hydrophobicity", "cys_density", "local_positive_charge_density")

                # re-train across folds to get ensembled scores
                ens_scores = _recompute_ensemble(
                    folded, all_tuples, TOMATO_PROTEOME, scratch / f"ens_{tag}", arm,
                    seq_fields,
                )
                labels = [r.label for r in folded]
                cluster_ids = [r.cluster_id for r in folded]
                ci = _cluster_bootstrap_ci(ens_scores, labels, cluster_ids)

                per_fold_aps = [f["ap"] for f in arm_result["per_fold"]]
                valid_aps = [v for v in per_fold_aps if v is not None]
                mean_ap = sum(valid_aps) / len(valid_aps) if valid_aps else None

                results[arena][arm] = {
                    "per_fold": arm_result["per_fold"],
                    "ensemble_ap": arm_result["ensemble_ap"],
                    "ensemble_recall_at_10": arm_result["ensemble_recall_at_10"],
                    "ensemble_recall_at_50": arm_result["ensemble_recall_at_50"],
                    "mean_fold_ap": mean_ap,
                    "ap_bootstrap_ci": ci,
                }

                print(
                    f"  per-fold AP: {per_fold_aps}"
                )
                print(
                    f"  mean fold AP: {mean_ap:.4f}, "
                    f"ensemble AP: {arm_result['ensemble_ap']:.4f}, "
                    f"CI: [{ci['lower']:.4f}, {ci['upper']:.4f}]"
                )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    # 5. write per-site scores TSV
    _write_per_site_scores(output_dir, arena_data, TOMATO_PROTEOME)

    # 6. assemble summary
    inputs_sha = {
        "kiae271_xlsx": _sha256(KIAE271_XLSX),
        "tomato_proteome": _sha256(TOMATO_PROTEOME),
        "cluster_file": _sha256(CLUSTER_FILE),
    }

    summary: dict[str, Any] = {
        "experiment": "tomato_local_v1",
        "description": (
            "Within-tomato candidate-prioritisation model trained and evaluated "
            "on 99 kiae271 persulfidation sites (Zhang et al. 2024, "
            "doi:10.1093/plphys/kiae271) with homology-cluster 5-fold CV. "
            "Not a general persulfidation predictor — single-study, single-lab, "
            "single-chemistry, single-species."
        ),
        "dataset": {
            "study": "KIAE271_SUPPL",
            "doi": "10.1093/plphys/kiae271",
            "species": "Solanum lycopersicum",
            "lab": "Zhang, South China Agricultural University",
            "total_verified_sites": table.total_verified_sites,
            "total_verified_proteins": table.total_verified_proteins,
        },
        "cv_design": {
            "method": "grouped_cluster_k_fold",
            "n_folds": N_FOLDS,
            "fold_seed": SEED_FOLDS,
            "unit": "MMseqs2 homology cluster (30% identity, >=50% coverage)",
            "proteome": "Solanum lycopersicum reference v1 (36,988 proteins)",
            "n_clusters": len(cluster_map),
        },
        "arenas": arena_summary,
        "arms": {
            "seq_2": {
                "features": ["hydrophobicity", "cys_density"],
                "description": "Replicates the 2-feature cross-species panel evaluation baseline",
            },
            "seq_3": {
                "features": ["hydrophobicity", "cys_density", "local_positive_charge_density"],
                "description": "Adds thiolate-stabilisation proxy (K/R fraction in flanking window)",
            },
            "seq_structure": {
                "features": [
                    "hydrophobicity",
                    "cys_density",
                    "local_positive_charge_density",
                    "contact_number_proxy",
                    "plddt",
                ],
                "description": "Sequence + AlphaFold structure features",
            },
        },
        "model_seeds": list(MODEL_SEEDS),
        "results": results,
        "inputs_sha256": inputs_sha,
        "comparison_to_cross_species_baseline": {
            "cross_species_recovery": "48/99 recovered (48.5%, p=0.84 vs 50% chance)",
            "cross_species_mean_percentile": 49.8,
            "note": (
                "The cross-species result used an Arabidopsis-trained 2-feature "
                "pu_logistic model and scored the same 99 sites as external controls. "
                "Comparison is between cross-species-transfer (model trained on "
                "different species) and within-species (model trained and tested on "
                "tomato with homology-cluster-level group-holdout)."
            ),
        },
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"\nSummary written to {output_dir / 'summary.json'}")
    return summary


def _recompute_ensemble(
    folded_rows: list[Any],
    all_tuples: list[tuple[str, int, str]],
    proteome_path: Path,
    scratch_dir: Path,
    arm: str,
    seq_fields: tuple[str, ...],
) -> list[float]:
    """Re-run the full cross-validation to produce ensembled scores per site.
    Separate from _evaluate_one_arm so bootstrap CIs can be computed."""
    n_rows = len(all_tuples)

    per_seed_scores: list[list[float]] = []
    for seed in MODEL_SEEDS:
        seed_scores: list[float] = [0.0] * n_rows
        for fold_idx in range(N_FOLDS):
            train_idx = [i for i, r in enumerate(folded_rows) if r.fold != fold_idx]
            test_idx = [i for i, r in enumerate(folded_rows) if r.fold == fold_idx]
            train_rows = [all_tuples[i] for i in train_idx]
            test_rows = [all_tuples[i] for i in test_idx]
            train_y = [folded_rows[i].label for i in train_idx]

            train_seq = _seq_feature_vectors(
                train_rows, proteome_path, scratch_dir,
                f"ens_train_{arm}_f{fold_idx}_s{seed}", seq_fields,
            )
            test_seq = _seq_feature_vectors(
                test_rows, proteome_path, scratch_dir,
                f"ens_test_{arm}_f{fold_idx}_s{seed}", seq_fields,
            )

            if arm == "seq_structure":
                train_struct, _ = _structure_feature_vectors(
                    train_rows, scratch_dir,
                    f"ens_train_struct_f{fold_idx}_s{seed}",
                )
                test_struct, _ = _structure_feature_vectors(
                    test_rows, scratch_dir,
                    f"ens_test_struct_f{fold_idx}_s{seed}",
                )
                train_X = [s + u for s, u in zip(train_seq, train_struct, strict=True)]
                test_X = [s + u for s, u in zip(test_seq, test_struct, strict=True)]
            else:
                train_X = train_seq
                test_X = test_seq

            scores = _fit_pu_logistic(train_X, train_y, test_X, seed)
            for i, s in zip(test_idx, scores):
                seed_scores[i] = s
        per_seed_scores.append(seed_scores)

    return [
        sum(per_seed_scores[s][i] for s in range(len(MODEL_SEEDS)))
        / len(MODEL_SEEDS)
        for i in range(n_rows)
    ]


def _write_per_site_scores(
    output_dir: Path,
    arena_data: dict[str, list[Any]],
    proteome_path: Path,
) -> None:
    """Write per-site scores for the seq_3 and seq_structure arms, both arenas."""
    from plantpersulf.features.sequence import _load_proteome

    proteome = _load_proteome(proteome_path)
    scratch = Path(tempfile.mkdtemp(prefix="tomato_local_v1_scores_"))

    try:
        for arena in ARENAS:
            folded = arena_data[arena]
            all_tuples = [
                (r.protein_accession, r.cys_position, r.label) for r in folded
            ]
            for arm in ("seq_3", "seq_structure"):
                seq_fields: tuple[str, ...] = (
                    ("hydrophobicity", "cys_density")
                    if arm == "seq_2"
                    else ("hydrophobicity", "cys_density", "local_positive_charge_density")
                )
                ens_scores = _recompute_ensemble(
                    folded, all_tuples, proteome_path,
                    scratch / f"scores_{arena}_{arm}", arm, seq_fields,
                )

                fields = (
                    "protein_accession",
                    "cys_position",
                    "label",
                    "cluster_id",
                    "fold",
                    "score",
                )
                rows_out: list[dict[str, str]] = []
                for row, score in zip(folded, ens_scores):
                    rows_out.append({
                        "protein_accession": row.protein_accession,
                        "cys_position": str(row.cys_position),
                        "label": row.label,
                        "cluster_id": row.cluster_id,
                        "fold": str(row.fold),
                        "score": f"{score:.10g}",
                    })
                path = output_dir / f"scores_{arena}_{arm}.tsv"
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("w", encoding="utf-8", newline="") as h:
                    w = csv.DictWriter(
                        h, fieldnames=list(fields),
                        delimiter="\t", lineterminator="\n",
                    )
                    w.writeheader()
                    w.writerows(rows_out)
                print(f"  wrote {path}")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)


def build_parser() -> Any:
    import argparse

    p = argparse.ArgumentParser(
        description="tomato-local model v1: within-species training & evaluation"
    )
    p.add_argument(
        "--output-dir", type=Path, default=OUTPUT_DIR,
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_tomato_local_v1(output_dir=args.output_dir)
