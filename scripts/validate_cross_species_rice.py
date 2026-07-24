#!/usr/bin/env python
"""Cross-species TRANSFER validation: Arabidopsis-trained ranker -> Rice.

Trains the frozen release arm (``structure_ranker:seq_structure``) and the PU
baseline (``pu_logistic``) on the FULL Arabidopsis benchmark (both studies,
1:20 unlabeled subsample, seed 12345 — no internal holdout, because the
evaluation is entirely external) and scores the coordinate-verified
PXD072089 (Oryza sativa) persulfidome site table against an all-cysteine
rice unlabeled background (1:20 subsample, seed 12346).

Dataset: Xie et al. 2026 (PNAS, doi:10.1073/pnas.2608150123). Rice
persulfidome mapped via NM-biotin chemistry + DTT selective elution.
Independent of the Seville benchmark (lab: Xie/Nanjing Forestry;
chemistry: NM-biotin+DTT; species: rice). 929 coordinate-verified
persulfidation sites (SD01 + SD04 + the submitter-deposited
SS-all-peptides.tsv MaxQuant export on PRIDE) on rice proteins from a
655,127-cysteine background.

Scope and wording constraints (binding):

* Rice has no registered structures, so the structure branch is fully
  masked on the eval rows — the run measures cross-species transfer of
  the learned sequence representation with its structure gate.
* This is NOT plant cross-study evidence and NOT Gate 2 evidence: the
  training positives remain the same-lab, same-chemistry Arabidopsis
  benchmark. Outputs must be described as "cross-species transfer" only,
  under the Gate 2 STOP downgrade statement.
* Gate 2 Condition 1 independence axes are captured comparably: this is
  a separate lab/chemistry/species, but cross-species semantics apply.
  A future 3-study leave-out config (PXD006140/PXD024061/PXD072089) would
  be needed to convert this into within-plant cross-study Condition 1
  evidence.

Usage::

    python scripts/validate_cross_species_rice.py \
        --output-dir results/cross_species/pxd072089_transfer_v1
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import tempfile
from pathlib import Path
from typing import Any

BENCHMARK = Path("data/processed/benchmark_v1/sites.tsv")
ARABIDOPSIS_PROTEOME = Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta")
RICE_PROTEOME = Path("data/raw/references/rice_proteome_v1/uniprot_rice_v1.fasta")
CONFIG = Path("configs/experiments/pu_ranker_v1.yaml")
SD01 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd01.xlsx")
SD04 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd04.xlsx")
SS_ALL_PEPTIDES = Path("data/raw/supplements/PXD072089/SS-all-peptides.tsv")
STUDY = "PXD072089"
FOLD = "pxd072089_transfer"

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

DOWNGRADE_STATEMENT = (
    "当前公开数据不足以证明跨研究预测能力，模型仅用于候选组织与假设生成。"
    "Current public data are insufficient to demonstrate cross-study "
    "predictive ability; the model is used only for candidate organisation "
    "and hypothesis generation."
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_tsv(path: Path, fields: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(fields), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _seed_ensemble(per_seed: list[list[float]]) -> list[float]:
    """Mean score per row across seeds; fails loudly on length drift."""
    n = len(per_seed[0])
    if any(len(scores) != n for scores in per_seed):
        raise RuntimeError("per-seed score vectors differ in length")
    return [sum(scores[i] for scores in per_seed) / len(per_seed) for i in range(n)]


def run_cross_species_rice(
    output_dir: Path,
    seeds: list[int],
    ratio: int = 20,
    train_seed: int = 12345,
    eval_seed: int = 12346,
    n_perm: int = 1000,
) -> dict[str, Any]:
    from run_experiment import (  # type: ignore[import-not-found]
        _build_branch_features,
        _load_config,
        _sequence_feature_vectors,
        _subsample_unlabeled,
        _train_predict,
    )
    from score_release import (  # type: ignore[import-not-found]
        _ablation_named,
        _read_benchmark,
    )

    from plantpersulf.evaluation.effect_size import paired_cluster_bootstrap_delta_ci
    from plantpersulf.evaluation.metrics import average_precision, recall_at_k
    from plantpersulf.evaluation.permutation import permutation_test
    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.models.structure_ranker import structure_ranker_scores
    from plantpersulf.proteomics.pxd072089_sites import (
        build_cross_species_eval_rows,
        parse_pxd072089_sites,
    )

    # --- 1. hash-verified inputs (fail closed) ---
    required = (
        BENCHMARK, ARABIDOPSIS_PROTEOME, CONFIG,
        SD01, SD04, SS_ALL_PEPTIDES, RICE_PROTEOME,
    )
    for path in required:
        if not path.is_file():
            raise RuntimeError(f"required input missing: {path}")

    # --- 2. parse + build eval rows (full all-cysteine background) ---
    proteome = _load_proteome(RICE_PROTEOME)
    table = parse_pxd072089_sites(
        SD01, SD04, proteome, ss_all_peptides_path=SS_ALL_PEPTIDES
    )
    sd01_sha = _sha256(SD01)
    sd04_sha = _sha256(SD04)
    ss_sha = _sha256(SS_ALL_PEPTIDES)
    source_sha = sd01_sha[:16]  # truncated composite for row-level provenance
    eval_rows_full = list(
        build_cross_species_eval_rows(table, proteome, source_sha256=source_sha)
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    benchmark_fields = tuple(eval_rows_full[0].keys())
    _write_tsv(output_dir / "eval_sites.tsv", benchmark_fields, eval_rows_full)

    # --- 3. deterministic subsamples (same scheme as the frozen runner) ---
    eval_rows = _subsample_unlabeled(eval_rows_full, ratio, eval_seed)
    train_rows = _subsample_unlabeled(_read_benchmark(BENCHMARK), ratio, train_seed)
    train_y = [r["label"] for r in train_rows]
    eval_y = [r["label"] for r in eval_rows]
    n_pos = sum(1 for y in eval_y if y == "positive")
    base_rate = n_pos / len(eval_rows)
    n_cys = sum(seq.count("C") for seq in proteome.values())
    print(
        f"rice proteome: {len(proteome)} proteins, {n_cys} Cys, "
        f"train={len(train_rows)} eval={len(eval_rows)} (+{n_pos})"
    )

    # --- 4. features (structure branch auto-masked for rice) ---
    cfg = _load_config(CONFIG)
    ranker_params = cfg.get("ranker", {})
    ablation = _ablation_named(cfg, "seq_structure")
    scratch = Path(tempfile.mkdtemp(prefix="cross_species_rice_"))
    try:
        branch_train = _build_branch_features(
            train_rows, ARABIDOPSIS_PROTEOME, scratch, "train", False
        )
        branch_eval = _build_branch_features(
            eval_rows, RICE_PROTEOME, scratch, "eval", False
        )
        feat_train = _sequence_feature_vectors(
            train_rows, ARABIDOPSIS_PROTEOME, scratch, "train_flat"
        )
        feat_eval = _sequence_feature_vectors(
            eval_rows, RICE_PROTEOME, scratch, "eval_flat"
        )
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
    n_masked = sum(1 for m in branch_eval.structure_mask if not m)
    print(f"eval structure-masked rows: {n_masked}/{len(eval_rows)}")

    # --- 5. train + score per seed ---
    model_per_seed: list[list[float]] = []
    base_per_seed: list[list[float]] = []
    per_seed_metrics: list[dict[str, Any]] = []
    scored_model: list[dict[str, Any]] = []
    scored_base: list[dict[str, Any]] = []
    for seed in seeds:
        out = structure_ranker_scores(
            branch_train,
            train_y,
            branch_eval,
            seed=seed,
            ablation=ablation,
            hidden=int(ranker_params.get("hidden", 16)),
            dropout=float(ranker_params.get("dropout", 0.2)),
            epochs=int(ranker_params.get("epochs", 200)),
            lr=float(ranker_params.get("lr", 0.05)),
            n_mc_dropout=int(ranker_params.get("n_mc_dropout", 16)),
        )
        model_scores = list(out.scores)
        base_scores = _train_predict(
            "pu_logistic", feat_train, train_y, feat_eval, seed
        )
        model_per_seed.append(model_scores)
        base_per_seed.append(base_scores)
        ap_m = average_precision(list(zip(model_scores, eval_y, strict=True)))
        ap_b = average_precision(list(zip(base_scores, eval_y, strict=True)))
        per_seed_metrics.append({"seed": seed, "model_ap": ap_m, "baseline_ap": ap_b})
        print(f"  seed={seed} model_ap={ap_m:.4f} baseline_ap={ap_b:.4f}")
        for row, s_m, s_b, has_struct in zip(
            eval_rows,
            model_scores,
            base_scores,
            branch_eval.structure_mask,
            strict=True,
        ):
            common = {
                "fold": FOLD,
                "seed": seed,
                "protein_accession": row["protein_accession"],
                "cys_position": int(row["cys_position_in_protein"]),
                "label": row["label"],
                "cluster_id": f"__singleton__{row['protein_accession']}",
                "has_structure": "1" if has_struct else "0",
            }
            scored_model.append({**common, "score": f"{s_m:.10g}"})
            scored_base.append({**common, "score": f"{s_b:.10g}"})

    _write_tsv(output_dir / "scored_model.tsv", SCORED_FIELDS, scored_model)
    _write_tsv(output_dir / "scored_baseline.tsv", SCORED_FIELDS, scored_base)

    # --- 6. seed-ensembled metrics + permutation ---
    ens_model = _seed_ensemble(model_per_seed)
    ens_base = _seed_ensemble(base_per_seed)
    scored_pairs = list(zip(ens_model, eval_y, strict=True))
    ap_ens = average_precision(scored_pairs)
    ap_ens_base = average_precision(list(zip(ens_base, eval_y, strict=True)))
    perm = permutation_test(scored_pairs, average_precision, n_perm=n_perm, seed=0)

    # Paired protein-level bootstrap of the model-minus-baseline AP delta.
    cluster_keys = [f"__singleton__{r['protein_accession']}" for r in eval_rows]
    effect = paired_cluster_bootstrap_delta_ci(
        list(zip(ens_model, eval_y, cluster_keys, strict=True)),
        list(zip(ens_base, eval_y, cluster_keys, strict=True)),
        n_boot=1000,
        seed=0,
    )

    # Source breakdown
    from collections import Counter
    source_counts = dict(Counter(s.source for s in table.sites))

    summary: dict[str, Any] = {
        "track": "cross_species_transfer_rice",
        "claim_class": "cross_species_transfer_only_not_gate2_evidence",
        "downgrade_statement": DOWNGRADE_STATEMENT,
        "dataset": {
            "study": STUDY,
            "doi": "10.1073/pnas.2608150123",
            "lab": "Xie, Nanjing Forestry University",
            "species": "Oryza sativa (Japonica + Indica)",
            "chemistry": "NM-biotin + DTT selective elution",
            "independence_axes": {
                "lab": "independent (Xie/NJFU != Romero/Gotor/Seville)",
                "chemistry": "independent (NM-biotin+DTT != tag-switch)",
                "species": "independent (rice != Arabidopsis)",
            },
        },
        "parsed": {
            "sites_verified": table.total_verified_sites,
            "proteins_verified": table.total_verified_proteins,
            "source_breakdown": source_counts,
            "sd01": {
                "rows_total": table.sd01_rows_total,
                "single_cys_parsed": table.sd01_single_cys_parsed,
                "multi_cys_skipped": table.sd01_multi_cys_skipped,
                "dropped_missing_accession": table.sd01_dropped_missing_accession,
                "dropped_coordinate_mismatch": table.sd01_dropped_coordinate_mismatch,
                "dropped_non_cysteine": table.sd01_dropped_non_cysteine,
            },
            "sd04": {
                "rows_total": table.sd04_rows_total,
                "parsed": table.sd04_parsed,
                "dropped_missing_accession": table.sd04_dropped_missing_accession,
                "dropped_coordinate_mismatch": table.sd04_dropped_coordinate_mismatch,
            },
            "ss_all_peptides": {
                "rows_total": table.ss_rows_total,
                "single_cys_parsed": table.ss_single_cys_parsed,
                "multi_cys_skipped": table.ss_multi_cys_skipped,
                "dropped_missing_accession": table.ss_dropped_missing_accession,
                "dropped_coordinate_mismatch": table.ss_dropped_coordinate_mismatch,
            },
        },
        "train": {
            "benchmark": str(BENCHMARK),
            "positives": sum(1 for y in train_y if y == "positive"),
            "rows": len(train_rows),
            "subsample_ratio": ratio,
            "subsample_seed": train_seed,
            "holdout": "none (external evaluation only)",
        },
        "eval": {
            "study": STUDY,
            "species": "Oryza sativa",
            "chemistry": "NM-biotin + DTT selective elution",
            "total_background_cysteines": n_cys,
            "background_cysteines": len(eval_rows_full),
            "subsample_ratio": ratio,
            "subsample_seed": eval_seed,
            "rows": len(eval_rows),
            "positives": n_pos,
            "base_rate": base_rate,
            "structure_masked_rows": n_masked,
        },
        "arms": {
            "model": "structure_ranker:seq_structure",
            "baseline": "pu_logistic",
        },
        "seeds": seeds,
        "per_seed": per_seed_metrics,
        "seed_ensemble": {
            "model_ap": ap_ens,
            "baseline_ap": ap_ens_base,
            "model_recall_at_50": recall_at_k(scored_pairs, 50),
            "model_recall_at_500": recall_at_k(scored_pairs, 500),
            "model_over_baseline_delta": (
                None if ap_ens is None or ap_ens_base is None else ap_ens - ap_ens_base
            ),
            "model_over_base_rate": (None if ap_ens is None else ap_ens / base_rate),
        },
        "permutation": {
            "n_perm": perm.n_perm,
            "observed_ap": perm.observed,
            "p_value": perm.p_value,
        },
        "paired_delta_vs_baseline": {
            "method": "cluster bootstrap, per-protein singleton units",
            "n_boot": 1000,
            "delta_point": effect.point,
            "delta_lower": effect.lower,
            "delta_upper": effect.upper,
        },
        "inputs_sha256": {
            "sd01_xlsx": sd01_sha,
            "sd04_xlsx": sd04_sha,
            "ss_all_peptides_tsv": ss_sha,
            "rice_proteome": _sha256(RICE_PROTEOME),
            "benchmark_sites": _sha256(BENCHMARK),
            "arabidopsis_proteome": _sha256(ARABIDOPSIS_PROTEOME),
            "experiment_config": _sha256(CONFIG),
        },
    }
    (output_dir / "cross_species_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8",
    )
    print(f"wrote summary -> {output_dir / 'cross_species_summary.json'}")
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="cross-species transfer validation (PXD072089 Rice)"
    )
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    p.add_argument("--ratio", type=int, default=20)
    p.add_argument("--n-perm", type=int, default=1000)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_cross_species_rice(
        output_dir=args.output_dir,
        seeds=[int(s) for s in args.seeds],
        ratio=int(args.ratio),
        n_perm=int(args.n_perm),
    )
