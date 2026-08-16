#!/usr/bin/env python
"""Cross-species TRANSFER validation: Arabidopsis-trained ranker -> Magnaporthe.

Trains the frozen release arm (``structure_ranker:seq_structure``) and the PU
baseline (``pu_logistic``) on the FULL Arabidopsis benchmark (both studies,
1:20 unlabeled subsample, seed 12345 — no internal holdout, because the
evaluation is entirely external) and scores the coordinate-verified
PXD063170 (Magnaporthe oryzae) S-sulfhydration site table against an
all-cysteine MG8 unlabeled background (1:20 subsample, seed 12346).

Scope and wording constraints (binding):

* Magnaporthe has no registered structures, so the structure branch is
  fully masked on the eval rows — the run measures cross-species transfer
  of the learned sequence representation with its structure gate.
* This is NOT plant cross-study evidence and NOT Gate 2 evidence: the
  training positives remain the same-lab, same-chemistry Arabidopsis
  benchmark. Outputs must be described as "cross-species transfer" only,
  under the Gate 2 STOP downgrade statement.

All scientific inputs are hash-verified through the supplementary-source
registry before use; the xlsx is the registered source of truth and the TSV
it is converted from must sit next to it.

Usage::

    python scripts/validate_cross_species.py \
        --output-dir results/cross_species/pxd063170_transfer_v1
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
CONFIG = Path("configs/experiments/pu_ranker_v1.yaml")
SITE_TSV_NAME = "PXD063170_sites_moesm3.tsv"
XLSX_MEMBER = "41467_2025_61582_MOESM3_ESM.xlsx"
PROTEOME_MEMBER = "Magnaporthe_oryzae.MG8.pep.all.fa"
STUDY = "PXD063170"
FOLD = "pxd063170_transfer"

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


def run_cross_species(
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
    from plantpersulf.models.structure_ranker import structure_ranker_scores
    from plantpersulf.proteomics.pxd063170_sites import (
        build_cross_species_eval_rows,
        load_ensembl_fungi_proteome,
        parse_pxd063170_sites,
    )
    from plantpersulf.provenance.supplementary import audit_supplementary_sources

    # --- 1. hash-verified inputs (fail closed) ---
    sources = {(s.study_accession, s.member): s for s in audit_supplementary_sources()}
    xlsx = sources[(STUDY, XLSX_MEMBER)].local_path
    proteome_path = sources[(STUDY, PROTEOME_MEMBER)].local_path
    site_tsv = xlsx.with_name(SITE_TSV_NAME)
    if not site_tsv.is_file():
        raise RuntimeError(
            f"derived site TSV missing (convert the registered xlsx first): {site_tsv}"
        )
    for required in (BENCHMARK, ARABIDOPSIS_PROTEOME, CONFIG):
        if not required.is_file():
            raise RuntimeError(f"required input missing: {required}")

    # --- 2. parse + build eval rows (full all-cysteine background) ---
    proteome = load_ensembl_fungi_proteome(proteome_path)
    table = parse_pxd063170_sites(site_tsv, proteome)
    tsv_sha = _sha256(site_tsv)
    eval_rows_full = list(
        build_cross_species_eval_rows(table, proteome, source_sha256=tsv_sha)
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
    print(f"train={len(train_rows)} eval={len(eval_rows)} (+{n_pos})")

    # --- 4. features (structure branch auto-masked for Magnaporthe) ---
    cfg = _load_config(CONFIG)
    ranker_params = cfg.get("ranker", {})
    ablation = _ablation_named(cfg, "seq_structure")
    scratch = Path(tempfile.mkdtemp(prefix="cross_species_"))
    try:
        branch_train = _build_branch_features(
            train_rows, ARABIDOPSIS_PROTEOME, scratch, "train", False
        )
        branch_eval = _build_branch_features(
            eval_rows, proteome_path, scratch, "eval", False
        )
        feat_train = _sequence_feature_vectors(
            train_rows, ARABIDOPSIS_PROTEOME, scratch, "train_flat"
        )
        feat_eval = _sequence_feature_vectors(
            eval_rows, proteome_path, scratch, "eval_flat"
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
                # No cluster file exists for MG8; per-protein singletons keep
                # bootstrap units honest (same caveat as benchmark v1).
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

    # Paired protein-level bootstrap of the model-minus-baseline AP delta on
    # the seed-ensembled scores (same protein-granularity caveat as the
    # singleton cluster file in benchmark v1: whole proteins resampled, no
    # homology families).
    cluster_keys = [f"__singleton__{r['protein_accession']}" for r in eval_rows]
    effect = paired_cluster_bootstrap_delta_ci(
        list(zip(ens_model, eval_y, cluster_keys, strict=True)),
        list(zip(ens_base, eval_y, cluster_keys, strict=True)),
        n_boot=1000,
        seed=0,
    )

    summary: dict[str, Any] = {
        "track": "cross_species_transfer",
        "claim_class": "cross_species_transfer_only_not_gate2_evidence",
        "downgrade_statement": DOWNGRADE_STATEMENT,
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
            "species": "Magnaporthe oryzae",
            "chemistry": "IAA-PEO-biotin (CSE_OE vs WT)",
            "sites_parsed": len(table.sites),
            "dropped_low_localization": table.dropped_low_localization,
            "dropped_coordinate_mismatch": table.dropped_coordinate_mismatch,
            "dropped_missing_accession": table.dropped_missing_accession,
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
            "site_table_xlsx": sources[(STUDY, XLSX_MEMBER)].sha256,
            "site_table_tsv_derived": tsv_sha,
            "mg8_proteome": sources[(STUDY, PROTEOME_MEMBER)].sha256,
            "benchmark_sites": _sha256(BENCHMARK),
            "arabidopsis_proteome": _sha256(ARABIDOPSIS_PROTEOME),
            "experiment_config": _sha256(CONFIG),
        },
    }
    (output_dir / "cross_species_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"wrote summary -> {output_dir / 'cross_species_summary.json'}")
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="cross-species transfer validation (PXD063170 Magnaporthe)"
    )
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    p.add_argument("--ratio", type=int, default=20)
    p.add_argument("--n-perm", type=int, default=1000)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_cross_species(
        output_dir=args.output_dir,
        seeds=[int(s) for s in args.seeds],
        ratio=int(args.ratio),
        n_perm=int(args.n_perm),
    )
