#!/usr/bin/env python
"""Multi-species joint persulfidation-site ranker (v1).

Trains the shared Elkan-Noto PU logistic scorer on the merged four-species
positive set (~2,900 sites: Arabidopsis benchmark, tomato kiae271, rice
PXD072089, Magnaporthe PXD063170) with:

* sequence features: hydrophobicity, cys density, local positive charge
  density (thiolate-stabilisation prior)
* family conservation proxies (v1): PANTHER family species breadth and
  member count (annotation-layer only, no label leakage)
* species one-hot

Evaluation: 5-fold family-grouped CV (PANTHER family level, so homology
families never span folds), with the Arabidopsis-only model trained on the
same folds as the direct comparison baseline. Species-stratified AP and
recall@k reported per fold. The unlabeled background is a 1:20 subsample
of all other cysteines in each species' reference proteome.

Usage::

    python scripts/train_multispecies_v1.py \
        --output-dir results/multispecies_v1
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from plantpersulf.models.traditional import pu_logistic_regression_scores
from plantpersulf.proteomics.multispecies_dataset import (
    MultispeciesSite,
    build_multispecies_sites,
    family_grouped_folds,
    family_stats_from_panther,
    load_panther_family_ids,
)

BENCHMARK = Path("data/processed/benchmark_v1/sites.tsv")
PROTEOMES = {
    "arabidopsis": Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta"),
    "tomato": Path("data/raw/references/tomato_ref_proteome_v1.fasta"),
    "rice": Path("data/raw/references/rice_proteome_v1/uniprot_rice_v1.fasta"),
    "magnaporthe": Path(
        "data/raw/supplements/PXD063170/Magnaporthe_oryzae.MG8.pep.all.fa"
    ),
}
PANTHER_FILES = {
    "arabidopsis": Path(
        "data/raw/references/panther_annotations_v1/panther_arabidopsis_taxon3702.tsv"
    ),
    "tomato": Path(
        "data/raw/references/panther_annotations_v1/panther_tomato_taxon4081.tsv"
    ),
    "rice": Path(
        "data/raw/references/panther_annotations_v1/panther_rice_taxon4530.tsv"
    ),
    "magnaporthe": Path(
        "data/raw/references/panther_annotations_v1/panther_magnaporthe_taxon242507.tsv"
    ),
}
KIAE271_XLSX = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
SD01 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd01.xlsx")
SD04 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd04.xlsx")
SS_ALL = Path("data/raw/supplements/PXD072089/SS-all-peptides.tsv")
MAGNAPORTHE_TSV = Path("data/raw/supplements/PXD063170/PXD063170_sites_moesm3.tsv")

NO_FAMILY = "__nofamily__"
SUBSAMPLE_RATIO = 20
N_FOLDS = 5
FOLD_SEED = 20260810


def _load_proteome(path: Path) -> dict[str, str]:
    from plantpersulf.features.sequence import _load_proteome as _lp

    return _lp(path)


def _load_magnaporthe_proteome(path: Path) -> dict[str, str]:
    from plantpersulf.proteomics.pxd063170_sites import load_ensembl_fungi_proteome

    return load_ensembl_fungi_proteome(path)


def _arabidopsis_positives(benchmark_path: Path) -> list[MultispeciesSite]:
    sites: list[MultispeciesSite] = []
    with benchmark_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("label") == "positive":
                sites.append(
                    MultispeciesSite(
                        protein_accession=row["protein_accession"],
                        cys_position=int(row["cys_position_in_protein"]),
                        species="arabidopsis",
                        study_accession=row.get("study_accession", ""),
                        panther_family="",  # filled after merge
                    )
                )
    return sites


def _tomato_sites(xlsx_path: Path, proteome: dict[str, str]) -> list[MultispeciesSite]:
    from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites

    table = parse_kiae271_sites(xlsx_path, proteome)
    return [
        MultispeciesSite(
            protein_accession=s.protein_accession,
            cys_position=s.cys_position,
            species="tomato",
            study_accession="KIAE271_SUPPL",
            panther_family="",
        )
        for s in table.sites
    ]


def _rice_sites(
    sd01: Path, sd04: Path, proteome: dict[str, str], ss: Path
) -> list[MultispeciesSite]:
    from plantpersulf.proteomics.pxd072089_sites import parse_pxd072089_sites

    table = parse_pxd072089_sites(sd01, sd04, proteome, ss_all_peptides_path=ss)
    return [
        MultispeciesSite(
            protein_accession=s.protein_accession,
            cys_position=s.cys_position,
            species="rice",
            study_accession="PXD072089",
            panther_family="",
        )
        for s in table.sites
    ]


def _magnaporthe_sites(
    tsv_path: Path, proteome: dict[str, str]
) -> list[MultispeciesSite]:
    from plantpersulf.proteomics.pxd063170_sites import parse_pxd063170_sites

    table = parse_pxd063170_sites(tsv_path, proteome)
    return [
        MultispeciesSite(
            protein_accession=s.protein_accession,
            cys_position=s.cys_position,
            species="magnaporthe",
            study_accession="PXD063170",
            panther_family="",
        )
        for s in table.sites
    ]


def _sequence_features(
    proteomes: dict[str, dict[str, str]],
    sites: list[MultispeciesSite],
    window_radius: int = 10,
) -> list[list[float]]:
    from plantpersulf.features.sequence import (
        _flanking_window,
        _hydrophobicity,
        _local_positive_charge_density,
    )

    rows: list[list[float]] = []
    for site in sites:
        seq = proteomes[site.species].get(site.protein_accession)
        if seq is None or site.cys_position < 1 or site.cys_position > len(seq):
            rows.append([0.0, 0.0, 0.0])
            continue
        if seq[site.cys_position - 1] != "C":
            rows.append([0.0, 0.0, 0.0])
            continue
        flank = _flanking_window(seq, site.cys_position, window_radius)
        rows.append(
            [
                _hydrophobicity(flank),
                seq.count("C") / len(seq),
                _local_positive_charge_density(flank),
            ]
        )
    return rows


def _species_onehot(species: str) -> list[float]:
    return [
        1.0 if species == s else 0.0
        for s in ("arabidopsis", "tomato", "rice", "magnaporthe")
    ]


def _feature_matrix(
    proteomes: dict[str, dict[str, str]],
    sites: list[MultispeciesSite],
    family_stats: dict[str, tuple[int, int]],
) -> list[list[float]]:
    seq_feats = _sequence_features(proteomes, sites)
    matrix: list[list[float]] = []
    for site, seq in zip(sites, seq_feats, strict=True):
        breadth, members = family_stats.get(site.panther_family, (0, 0))
        matrix.append(
            [*seq, float(breadth), float(members), *_species_onehot(site.species)]
        )
    return matrix


def _unlabeled_background(
    proteomes: dict[str, dict[str, str]],
    positive_keys: dict[str, set[tuple[str, int]]],
    family_by_species: dict[str, dict[str, str]],
    family_stats: dict[str, tuple[int, int]],
    ratio: int,
    seed: int,
    sample_cap: int = 20000,
) -> tuple[list[dict[str, str]], list[list[float]]]:
    """All-other-Cys background rows (PU semantics: unlabeled, never
    negatives), subsampled per species, with the same feature vector
    layout as the positives — including the REAL PANTHER family stats of
    the background protein (a zeroed family feature would let the model
    trivially separate "has a family annotation" from "is a positive",
    which is an artefact, not signal)."""
    from plantpersulf.features.sequence import (
        _flanking_window,
        _hydrophobicity,
        _local_positive_charge_density,
    )

    rng = random.Random(seed)
    rows: list[dict[str, str]] = []
    feats: list[list[float]] = []
    for species in ("arabidopsis", "tomato", "rice", "magnaporthe"):
        proteome = proteomes[species]
        keys = positive_keys[species]
        candidates: list[tuple[str, int]] = []
        for acc, seq in proteome.items():
            for pos, res in enumerate(seq, start=1):
                if res == "C" and (acc, pos) not in keys:
                    candidates.append((acc, pos))
        rng.shuffle(candidates)
        sampled = candidates[: min(len(candidates) // ratio, sample_cap)]
        for acc, pos in sampled:
            seq = proteome[acc]
            flank = _flanking_window(seq, pos, 10)
            family = family_by_species[species].get(acc, NO_FAMILY)
            breadth, members = family_stats.get(family, (0, 0))
            rows.append(
                {
                    "protein_accession": acc,
                    "cys_position_in_protein": str(pos),
                    "species": species,
                    "label": "unlabeled",
                }
            )
            feats.append(
                [
                    _hydrophobicity(flank),
                    seq.count("C") / len(seq),
                    _local_positive_charge_density(flank),
                    float(breadth),
                    float(members),
                    *_species_onehot(species),
                ]
            )
    return rows, feats


def _run_fold(
    fold_idx: int,
    train_idx: list[int],
    test_idx: list[int],
    sites: list[MultispeciesSite],
    family_stats: dict[str, tuple[int, int]],
    proteomes: dict[str, dict[str, str]],
    unlabeled_rows: list[dict[str, str]],
    unlabeled_feats: list[list[float]],
    seeds: list[int],
) -> dict[str, Any]:
    train_sites = [sites[i] for i in train_idx]
    test_sites = [sites[i] for i in test_idx]

    # Multi-species arm: train on ALL species' positives
    train_feats = _feature_matrix(proteomes, train_sites, family_stats)
    train_y = ["positive"] * len(train_sites)
    combined_X = [*train_feats, *unlabeled_feats]
    combined_y = [*train_y, *["unlabeled"] * len(unlabeled_feats)]

    test_feats = _feature_matrix(proteomes, test_sites, family_stats)
    test_y = ["positive"] * len(test_sites)
    # evaluation mixes test positives with the unlabeled background —
    # scoring every test positive against the background's own rows
    # reproduces the real "rank all candidates in the proteome" use case.
    eval_feats = [*test_feats, *unlabeled_feats]
    eval_y = [*test_y, *["unlabeled"] * len(unlabeled_feats)]

    per_seed_scores: list[list[float]] = []
    for seed in seeds:
        scores = pu_logistic_regression_scores(
            combined_X, combined_y, eval_feats, seed=seed
        )
        per_seed_scores.append(scores)
    n = len(per_seed_scores[0])
    ens_scores = [
        sum(scores[i] for scores in per_seed_scores) / len(per_seed_scores)
        for i in range(n)
    ]

    # Arabidopsis-only baseline arm: train on arabidopsis positives only
    ath_idx = [i for i in train_idx if sites[i].species == "arabidopsis"]
    ath_sites = [sites[i] for i in ath_idx]
    ath_feats = _feature_matrix(proteomes, ath_sites, family_stats)
    ath_y = ["positive"] * len(ath_sites)
    ath_combined_X = [*ath_feats, *unlabeled_feats]
    ath_combined_y = [*ath_y, *["unlabeled"] * len(unlabeled_feats)]
    ath_scores = pu_logistic_regression_scores(
        ath_combined_X, ath_combined_y, eval_feats, seed=seeds[0]
    )

    # Metrics
    from plantpersulf.evaluation.metrics import average_precision, recall_at_k

    pairs = list(zip(ens_scores, eval_y, strict=True))
    ap = average_precision(pairs)
    ath_pairs = list(zip(ath_scores, eval_y, strict=True))
    ap_ath = average_precision(ath_pairs)

    # per-species AP: that species' test positives mixed with THAT
    # species' unlabeled background rows (not the whole cross-species
    # background, which would dilute per-species numbers asymmetrically).
    unlabeled_by_species: dict[str, list[int]] = defaultdict(list)
    for i, row in enumerate(unlabeled_rows):
        unlabeled_by_species[row["species"]].append(i)

    per_species_ap: dict[str, float] = {}
    for sp in ("arabidopsis", "tomato", "rice", "magnaporthe"):
        sp_test_idx = [i for i, site in enumerate(test_sites) if site.species == sp]
        sp_unl_idx = unlabeled_by_species.get(sp, [])
        if not sp_test_idx or not sp_unl_idx:
            continue
        sp_scores = [ens_scores[i] for i in sp_test_idx] + [
            ens_scores[len(test_sites) + j] for j in sp_unl_idx
        ]
        sp_labels = ["positive"] * len(sp_test_idx) + ["unlabeled"] * len(sp_unl_idx)
        per_species_ap[sp] = average_precision(
            list(zip(sp_scores, sp_labels, strict=True))
        )

    return {
        "fold": fold_idx,
        "n_test_positives": len(test_sites),
        "test_positives_by_species": dict(Counter(s.species for s in test_sites)),
        "multispecies_ap": ap,
        "multispecies_recall_at_50": recall_at_k(pairs, 50),
        "multispecies_recall_at_200": recall_at_k(pairs, 200),
        "arabidopsis_only_ap": ap_ath,
        "per_species_ap": per_species_ap,
        "base_rate": (
            sum(1 for y in eval_y if y == "positive") / len(eval_y) if eval_y else 0.0
        ),
    }


def _leave_one_species_out(
    sites: list[MultispeciesSite],
    held_out: str,
    proteomes: dict[str, dict[str, str]],
    unlabeled_rows: list[dict[str, str]],
    unlabeled_feats: list[list[float]],
    family_stats: dict[str, tuple[int, int]],
    seeds: list[int],
) -> dict[str, Any]:
    """Train on the other 3 species, evaluate on the held-out species
    (its positives + its own unlabeled background). Directly answers
    whether a learned signal transfers to an UNSEEN species — the
    species-specificity question."""
    from plantpersulf.evaluation.metrics import average_precision

    train_sites = [s for s in sites if s.species != held_out]
    test_sites = [s for s in sites if s.species == held_out]
    train_feats = _feature_matrix(proteomes, train_sites, family_stats)
    train_y = ["positive"] * len(train_sites)

    held_unl_idx = [i for i, r in enumerate(unlabeled_rows) if r["species"] == held_out]
    held_unl_feats = [unlabeled_feats[i] for i in held_unl_idx]

    combined_X = [*train_feats, *held_unl_feats]
    combined_y = [*train_y, *["unlabeled"] * len(held_unl_feats)]
    eval_feats = [
        *_feature_matrix(proteomes, test_sites, family_stats),
        *held_unl_feats,
    ]
    eval_y = ["positive"] * len(test_sites) + ["unlabeled"] * len(held_unl_feats)

    per_seed = [
        pu_logistic_regression_scores(combined_X, combined_y, eval_feats, seed=s)
        for s in seeds
    ]
    n = len(per_seed[0])
    ens = [sum(scores[i] for scores in per_seed) / len(per_seed) for i in range(n)]
    ap = average_precision(list(zip(ens, eval_y, strict=True)))
    base_rate = len(test_sites) / len(eval_y) if eval_y else 0.0
    return {
        "held_out_species": held_out,
        "n_train_positives": len(train_sites),
        "n_test_positives": len(test_sites),
        "n_unlabeled": len(held_unl_feats),
        "ap": ap,
        "enrichment_over_base_rate": ap / base_rate if base_rate else None,
        "base_rate": base_rate,
    }


def run_multispecies_v1(
    output_dir: Path,
    seeds: list[int] | None = None,
    ratio: int = SUBSAMPLE_RATIO,
    n_folds: int = N_FOLDS,
    fold_seed: int = FOLD_SEED,
) -> dict[str, Any]:
    seeds = seeds or [0, 1, 2, 3, 4]
    for path in (
        *PROTEOMES.values(),
        *PANTHER_FILES.values(),
        BENCHMARK,
        KIAE271_XLSX,
        SD01,
        SD04,
        SS_ALL,
        MAGNAPORTHE_TSV,
    ):
        if not path.is_file():
            raise RuntimeError(f"required input missing: {path}")

    proteomes = {
        "arabidopsis": _load_proteome(PROTEOMES["arabidopsis"]),
        "tomato": _load_proteome(PROTEOMES["tomato"]),
        "rice": _load_proteome(PROTEOMES["rice"]),
        "magnaporthe": _load_magnaporthe_proteome(PROTEOMES["magnaporthe"]),
    }
    panther_by_species = {
        sp: load_panther_family_ids(path) for sp, path in PANTHER_FILES.items()
    }
    family_stats = family_stats_from_panther(panther_by_species)

    merged = build_multispecies_sites(
        arabidopsis=_arabidopsis_positives(BENCHMARK),
        tomato=_tomato_sites(KIAE271_XLSX, proteomes["tomato"]),
        rice=_rice_sites(SD01, SD04, proteomes["rice"], SS_ALL),
        magnaporthe=_magnaporthe_sites(MAGNAPORTHE_TSV, proteomes["magnaporthe"]),
    )
    # assign PANTHER family ids (family level, from each species' annotation file)
    sites = [
        MultispeciesSite(
            protein_accession=s.protein_accession,
            cys_position=s.cys_position,
            species=s.species,
            study_accession=s.study_accession,
            panther_family=panther_by_species[s.species].get(
                s.protein_accession, NO_FAMILY
            ),
        )
        for s in merged
    ]

    n_by_species = dict(Counter(s.species for s in sites))
    print(f"merged positive sites: {len(sites)} {n_by_species}")

    positive_keys: dict[str, set[tuple[str, int]]] = {
        sp: {(s.protein_accession, s.cys_position) for s in sites if s.species == sp}
        for sp in ("arabidopsis", "tomato", "rice", "magnaporthe")
    }
    unlabeled_rows, unlabeled_feats = _unlabeled_background(
        proteomes, positive_keys, panther_by_species, family_stats, ratio, seed=12345
    )
    print(f"unlabeled background rows: {len(unlabeled_rows)}")

    folds = family_grouped_folds(sites, n_folds=n_folds, seed=fold_seed)
    fold_results = [
        _run_fold(
            i,
            tr,
            te,
            sites,
            family_stats,
            proteomes,
            unlabeled_rows,
            unlabeled_feats,
            seeds,
        )
        for i, (tr, te) in enumerate(folds)
    ]

    ap_ms = [f["multispecies_ap"] for f in fold_results]
    ap_ath = [f["arabidopsis_only_ap"] for f in fold_results]

    # Leave-one-species-out: does the learned signal transfer to an
    # UNSEEN species? (the species-specificity question)
    loso_results = [
        _leave_one_species_out(
            sites, sp, proteomes, unlabeled_rows, unlabeled_feats, family_stats, seeds
        )
        for sp in ("arabidopsis", "tomato", "rice", "magnaporthe")
    ]

    summary: dict[str, Any] = {
        "track": "multispecies_joint_training_v1",
        "claim_class": "multispecies_candidate_ranking_not_gate2",
        "design": {
            "species": list(n_by_species),
            "positive_counts_by_species": n_by_species,
            "total_positives": len(sites),
            "unlabeled_background_rows": len(unlabeled_rows),
            "subsample_ratio": ratio,
            "folds": n_folds,
            "fold_seed": fold_seed,
            "grouping": "PANTHER family-level homology clusters",
            "features": [
                "hydrophobicity",
                "cys_density",
                "local_positive_charge_density",
                "family_species_breadth",
                "family_member_count",
                "species_onehot_x4",
            ],
            "model": "Elkan-Noto PU logistic (shared across species)",
            "leakage_control": "family members never span CV folds; "
            "conservation features from PANTHER annotation layer only",
        },
        "per_fold": fold_results,
        "seed_ensembled_pooled": {
            "multispecies_ap_mean": sum(ap_ms) / len(ap_ms),
            "arabidopsis_only_ap_mean": sum(ap_ath) / len(ap_ath),
        },
        "leave_one_species_out": loso_results,
        "note": (
            "Same-lab Arabidopsis-only baseline vs multi-species joint "
            "training, evaluated on identical family-grouped folds. PU "
            "semantics: unlabeled background is never a hard negative. "
            "This is candidate-ranking evidence, NOT cross-study "
            "predictive validation (training positives remain dominated "
            "by two Seville studies; the Gate 2 STOP decision is "
            "unchanged)."
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "multispecies_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"wrote summary -> {output_dir / 'multispecies_summary.json'}")
    print(
        f"mean AP: multispecies={sum(ap_ms) / len(ap_ms):.4f}  "
        f"arabidopsis-only={sum(ap_ath) / len(ap_ath):.4f}"
    )
    for f in fold_results:
        print(
            f"  fold {f['fold']}: n={f['n_test_positives']} "
            f"ms_ap={f['multispecies_ap']:.4f} ath_ap={f['arabidopsis_only_ap']:.4f} "
            f"by_species={ {k: round(v, 3) for k, v in f['per_species_ap'].items()} }"
        )
    print("\nleave-one-species-out (unseen-species transfer):")
    for r in loso_results:
        print(
            f"  held-out {r['held_out_species']}: ap={r['ap']:.4f} "
            f"base_rate={r['base_rate']:.4f} "
            f"enrichment={r['enrichment_over_base_rate']:.2f}x "
            f"(train n={r['n_train_positives']}, test n={r['n_test_positives']})"
        )
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="multi-species joint persulfidation ranker v1"
    )
    p.add_argument("--output-dir", type=Path, default=Path("results/multispecies_v1"))
    p.add_argument("--seeds", type=int, nargs="+", default=[0, 1, 2, 3, 4])
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_multispecies_v1(output_dir=args.output_dir, seeds=list(args.seeds))
