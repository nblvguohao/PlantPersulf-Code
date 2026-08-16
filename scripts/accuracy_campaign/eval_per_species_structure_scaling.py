#!/usr/bin/env python
"""Diagnostic B: does per-species structure scaling beat the global scaler?

The W1 evaluation (``w1_tomato_structures/report.md``) showed that swapping in
tomato structure coverage regresses arabidopsis/rice badly (macro
0.3863->0.2117, strict 0.5383->0.3257). Mechanism lead #1: the frozen
``_BranchScalers.struct`` is fit on ALL present-structure rows pooled, and
after release v3 tomato is ~78% of structure rows, so the global scaler is
tomato-dominated and arabidopsis/rice structure inputs get re-scaled against
tomato statistics.

This diagnostic tests that lead on the literature track, same v11 frozen
hyperparameters, same 10 seeds, and the SAME v3 structure registry as W1 —
differing ONLY in the structure scaling:

- **arm ``global``**: structure features standardized exactly as the frozen
  model does (one global ``TrainOnlyScaler`` over present rows) — this is the
  W1 v3 reproduction (reproduction_report must match W1's 0.2117 macro).
- **arm ``per_species``**: structure features first z-scored WITHIN each
  species (scaler fit on train present rows only), then the global scaler
  runs on the aligned values. After per-species z-scoring, columns have
  mean 0 / variance 1 per species, so the subsequent global fit is ~identity
  and the per-species alignment survives.

Both arms share the identical panel (``panel_sha256`` gate against the frozen
v11 summary), features, seeds, and scoring path. Only the structure values
entering the model differ. Claim class ``diagnostic_only``; touches no frozen
artifact; per-species scaling itself would be a new-release candidate change
if positive.

Output: per-arm per-seed TSVs + ``reproduction_report.json`` (global arm APs
vs frozen v11 summary) + ``comparison.json`` (per-species AP@200 / recall@200
for both arms, paired by seed).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import yaml

from plantpersulf.benchmark.literature_random_track import (
    build_literature_random_track,
)
from plantpersulf.evaluation.comparable_track import (
    build_registered_structure_features,
)
from plantpersulf.evaluation.species_structure_scaling import (
    fit_species_struct_scalers,
    transform_species_struct,
)
from plantpersulf.models.structure_ranker import (
    AblationConfig,
    BranchFeatures,
    structure_ranker_scores,
)
from plantpersulf.proteomics.multispecies_v2_sources import (
    sequence_feature_map,
)

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

STRUCTURE_RANKER_PARAMETERS = {
    "epochs": 200,
    "hidden": 16,
    "dropout": 0.2,
    "lr": 0.05,
    "holdout_fraction": 0.2,
    "n_mc_dropout": 16,
}

# K used for the paired AP@K / recall@K comparison (claim-gate style).
K = 200


def load_expected_panel_sha256s(summary_path: Path) -> dict[int, str]:
    raw = json.loads(summary_path.read_text(encoding="utf-8"))
    runs = raw["runs"]
    by_seed: dict[int, dict[str, str]] = {}
    for run in runs:
        seed = int(run["seed"])
        sha256 = str(run["panel_sha256"])
        by_seed.setdefault(seed, {})[run["model"]] = sha256
    if not by_seed:
        raise RuntimeError(f"summary contains no runs: {summary_path}")
    out: dict[int, str] = {}
    for seed, models in by_seed.items():
        if len(set(models.values())) != 1:
            raise RuntimeError(
                f"summary panel_sha256 disagree across models for seed {seed}"
            )
        out[seed] = next(iter(models.values()))
    return out


def write_ranker_scores_tsv(
    path: Path,
    *,
    rows: tuple[tuple[str, int, str], ...],
    scores: list[float],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    lines = ["global_protein_id\tcys_position\tlabel\tscore"]
    for (protein_id, position, label), score in zip(
        rows, scores, strict=True
    ):
        lines.append(f"{protein_id}\t{position}\t{label}\t{score}")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(path)


def _species_of(protein_id: str) -> str:
    return protein_id.split("|", 1)[0]


def per_species_apk_and_recall(
    rows: list[tuple[str, int, str]], scores: list[float], k: int = K
) -> dict[str, dict[str, float]]:
    """Per-species average precision @K and recall @K (positive reach)."""
    import sklearn.metrics

    by_species: dict[str, list[tuple[float, bool]]] = {}
    for (protein_id, _position, label), score in zip(rows, scores, strict=True):
        species = _species_of(protein_id)
        by_species.setdefault(species, []).append((score, label == "positive"))
    out: dict[str, dict[str, float]] = {}
    for species, pairs in by_species.items():
        positives = sum(1 for _, positive in pairs if positive)
        if positives == 0:
            continue
        ranked = sorted(pairs, key=lambda pair: pair[0], reverse=True)
        top_k = ranked[:k]
        recall = sum(1 for _, positive in top_k if positive) / positives
        out[species] = {
            "ap_k": float(
                sklearn.metrics.average_precision_score(
                    [positive for _, positive in pairs],
                    [score for score, _ in pairs],
                )
            ),
            "recall_k": recall,
            "n_positive": positives,
        }
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(
            "configs/experiments/multispecies_v2_global_clusters_v11.yaml"
        ),
    )
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--structure-registry",
        type=Path,
        required=True,
        help="v3 tomato structure release registry (same as W1).",
    )
    parser.add_argument("--seeds", type=str, default="0,1,2,3,4,5,6,7,8,9")
    parser.add_argument("--device", type=str, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument(
        "--w1-reproduction",
        type=Path,
        default=None,
        help="W1 v3 reproduction_report.json; when given, the global arm's "
        "test-only per-species AP must match it (reproduction gate).",
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        raise RuntimeError("multispecies v2 config is required")

    compute_cfg = cfg.get("compute")
    if not isinstance(compute_cfg, dict):
        raise RuntimeError("compute configuration is required")
    device = args.device or str(compute_cfg.get("device", "cpu"))
    batch_size = args.batch_size or int(compute_cfg.get("score_batch_size", 16384))

    # --- identical panel construction to train_multispecies_v2.py ---
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from train_multispecies_v2 import _development_rows

    rows, _frozen, proteomes = _development_rows(args.config, sample_unlabeled=False)
    random_cfg = cfg["literature_random_protein"]
    seeds = tuple(int(item) for item in args.seeds.split(",") if item)
    runs = build_literature_random_track(
        rows,
        seeds=seeds,
        test_fraction=random_cfg["test_fraction"],
        validation_fraction_of_remaining=random_cfg[
            "validation_fraction_of_remaining"
        ],
        unlabeled_per_positive=cfg["unlabeled_panel"]["per_positive"],
        sampling_seed=cfg["unlabeled_panel"]["seed"],
    )
    expected = load_expected_panel_sha256s(args.summary)
    missing = [seed for seed in seeds if seed not in expected]
    if missing:
        raise RuntimeError(f"summary lacks panel_sha256 for seeds {missing}")
    for run in runs:
        if run.panel_sha256 != expected[run.seed]:
            raise RuntimeError(
                f"rebuilt panel for seed {run.seed} does not match the frozen "
                f"v11 summary (built={run.panel_sha256[:16]}... "
                f"expected={expected[run.seed][:16]}...)"
            )

    # --- features over the union of all selected panel rows ---
    selected_rows = {
        (row.global_protein_id, row.cys_position): row
        for run in runs
        for partition in run.partitions.values()
        for row in partition
    }
    raw_features = sequence_feature_map(selected_rows.values(), proteomes)
    if any(any(value is None for value in values) for values in raw_features.values()):
        raise RuntimeError("sequence comparison features contain missing values")
    sequence_features = {
        key: tuple(float(value) for value in values)
        for key, values in raw_features.items()
    }

    build_cfg = cfg.get("comparison_feature_build")
    if not isinstance(build_cfg, dict):
        raise RuntimeError("comparison_feature_build configuration is required")
    structure_features = build_registered_structure_features(
        set(selected_rows),
        registry_path=args.structure_registry,
        registry_base=Path(str(build_cfg["structure_registry_base"])),
    )

    def branches(
        partition_rows: tuple, *, structure_override: list[list[float]] | None
    ) -> BranchFeatures:
        sequence = [
            list(sequence_features[(row.global_protein_id, row.cys_position)])
            for row in partition_rows
        ]
        structure_pairs = [
            structure_features[(row.global_protein_id, row.cys_position)]
            for row in partition_rows
        ]
        raw_structure = [
            [float(values) for values in raw]
            for raw, _available in structure_pairs
        ]
        structure = (
            structure_override
            if structure_override is not None
            else raw_structure
        )
        return BranchFeatures(
            sequence=sequence,
            esm=[[0.0] for _ in partition_rows],
            structure=structure,
            structure_mask=[available for _, available in structure_pairs],
            study_ids=None,
        )

    # --- per-seed per-arm scoring -------------------------------------------
    report: dict[str, object] = {
        "config": str(args.config),
        "summary": str(args.summary),
        "structure_registry": str(args.structure_registry),
        "device": device,
        "seeds": {seed: {"global": {}, "per_species": {}} for seed in seeds},
    }
    comparison: dict[int, dict[str, object]] = {}
    seeds_report = report["seeds"]
    assert isinstance(seeds_report, dict)

    for run in runs:
        partitions = dict(run.partitions)
        train_rows = partitions["train"]
        predict_rows = (*partitions["validation"], *partitions["test"])
        train_structure_pairs = [
            structure_features[(row.global_protein_id, row.cys_position)]
            for row in train_rows
        ]
        train_raw = [
            [float(values) for values in raw]
            for raw, _available in train_structure_pairs
        ]
        train_masks = [available for _, available in train_structure_pairs]
        train_species = [_species_of(row.global_protein_id) for row in train_rows]

        # Arm A: global scaling (W1 v3 reproduction — frozen behavior).
        global_branches = branches(train_rows, structure_override=None)
        predict_global = branches(predict_rows, structure_override=None)
        out_global = structure_ranker_scores(
            global_branches,
            [row.label for row in train_rows],
            predict_global,
            seed=run.seed,
            ablation=AblationConfig(use_esm=False, use_study_context=False),
            **STRUCTURE_RANKER_PARAMETERS,
            device_name=device,
            batch_size=batch_size,
        )

        # Arm B: per-species structure scaling (scalers fit on TRAIN present
        # rows only, then applied to train + predict).
        species_scalers = fit_species_struct_scalers(
            train_raw, train_masks, train_species
        )
        train_scaled = transform_species_struct(
            train_raw, train_masks, train_species, species_scalers
        )
        predict_structure_pairs = [
            structure_features[(row.global_protein_id, row.cys_position)]
            for row in predict_rows
        ]
        predict_raw = [
            [float(values) for values in raw]
            for raw, _available in predict_structure_pairs
        ]
        predict_masks = [available for _, available in predict_structure_pairs]
        predict_species = [
            _species_of(row.global_protein_id) for row in predict_rows
        ]
        predict_scaled = transform_species_struct(
            predict_raw, predict_masks, predict_species, species_scalers
        )
        scaled_branches = branches(train_rows, structure_override=train_scaled)
        predict_scaled_branches = branches(
            predict_rows, structure_override=predict_scaled
        )
        out_scaled = structure_ranker_scores(
            scaled_branches,
            [row.label for row in train_rows],
            predict_scaled_branches,
            seed=run.seed,
            ablation=AblationConfig(use_esm=False, use_study_context=False),
            **STRUCTURE_RANKER_PARAMETERS,
            device_name=device,
            batch_size=batch_size,
        )

        # --- partition-aligned metrics (validation then test) ---------------
        # ``structure_ranker_scores`` scores ``predict_rows`` in order
        # (validation, test); split the score vector at ``n_validation`` so
        # the test partition is scored alone — the W1-v3-comparable number.
        n_validation = len(partitions["validation"])
        test_rows = tuple(
            (row.global_protein_id, row.cys_position, row.label)
            for row in partitions["test"]
        )
        val_rows = tuple(
            (row.global_protein_id, row.cys_position, row.label)
            for row in partitions["validation"]
        )
        per_seed: dict[str, object] = {}
        for arm, output, name in (
            ("global", out_global, "global"),
            ("per_species", out_scaled, "per_species"),
        ):
            test_scores = list(output.scores[n_validation:])
            val_scores = list(output.scores[:n_validation])
            write_ranker_scores_tsv(
                args.output / f"seed{run.seed}_{name}_test.tsv",
                rows=test_rows,
                scores=test_scores,
            )
            write_ranker_scores_tsv(
                args.output / f"seed{run.seed}_{name}_validation.tsv",
                rows=val_rows,
                scores=val_scores,
            )
            test_metrics = per_species_apk_and_recall(
                list(test_rows), test_scores
            )
            val_metrics = per_species_apk_and_recall(list(val_rows), val_scores)
            seeds_report[run.seed][arm] = {
                "test": test_metrics,
                "validation": val_metrics,
            }
            per_seed[arm] = {
                "test_macro_ap_k": round(
                    sum(v["ap_k"] for v in test_metrics.values())
                    / len(test_metrics)
                    if test_metrics
                    else 0.0,
                    6,
                ),
                "test_species": test_metrics,
                "validation_macro_ap_k": round(
                    sum(v["ap_k"] for v in val_metrics.values())
                    / len(val_metrics)
                    if val_metrics
                    else 0.0,
                    6,
                ),
            }
            print(
                f"seed {run.seed} {arm}: test macro_ap@{K} "
                f"{per_seed[arm]['test_macro_ap_k']:.4f} "
                f"{ {s: round(v['ap_k'], 4) for s, v in test_metrics.items()} }"
            )

            # --- W1 v3 reproduction gate (global arm only) -------------------
            if arm == "global" and args.w1_reproduction is not None:
                w1 = json.loads(
                    args.w1_reproduction.read_text(encoding="utf-8")
                )
                expected = w1["seeds"][str(run.seed)]["test"]
                mismatches = {
                    species: (test_metrics[species]["ap_k"], expected[species])
                    for species in expected
                    if abs(test_metrics[species]["ap_k"] - expected[species]) > 1e-6
                }
                if mismatches:
                    raise RuntimeError(
                        f"seed {run.seed} global test AP diverges from W1 v3: "
                        f"{mismatches}"
                    )
                print(f"seed {run.seed} global reproduces W1 v3 (test AP)")

        comparison[run.seed] = per_seed

    # --- paired per-species deltas (test partition, per_species minus global)
    paired: dict[str, dict[str, object]] = {}
    for species in sorted(
        {
            s
            for seed in comparison.values()
            for arm in seed.values()
            for s in arm["test_species"]
        }
    ):
        deltas = [
            seed["per_species"]["test_species"].get(species, {}).get("ap_k", 0.0)
            - seed["global"]["test_species"].get(species, {}).get("ap_k", 0.0)
            for seed in comparison.values()
        ]
        paired[species] = {
            "mean_delta_ap_k": round(sum(deltas) / len(deltas), 6) if deltas else 0.0,
            "n_seeds_with_species": sum(
                1
                for seed in comparison.values()
                if species in seed["global"]["test_species"]
            ),
            "positive_delta_seeds": sum(1 for d in deltas if d > 0.0),
        }

    document = {
        "track": "w4_per_species_structure_scaling",
        "claim_class": "diagnostic_only",
        "note": (
            "Diagnostic B: per-species structure scaling vs the frozen global "
            "scaler, same v11 hyperparameters, same v3 tomato structure "
            "registry as W1, same 10 seeds and panel (panel_sha256 gate). "
            "Per-species scalers fit on train present rows only. Metrics are "
            "test-partition AP@K / recall@K (claim-gate style), so the global "
            "arm is directly comparable to W1 v3's test AP (reproduction gate "
            "when --w1-reproduction is given). This tests W1 mechanism lead "
            "#1 (global struct scaler is tomato-dominated)."
        ),
        "k": K,
        "seeds": seeds,
        "report": report,
        "paired_per_species_ap_delta": paired,
        "summary": {
            "headline": (
                "per-species test macro AP@200 minus global (positive = "
                "per-species scaling helps)"
            ),
            "per_species": {
                species: {k2: v for k2, v in values.items()}
                for species, values in paired.items()
            },
        },
    }
    args.output.mkdir(parents=True, exist_ok=True)
    comparison_path = args.output / "comparison.json"
    comparison_path.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {comparison_path}")


if __name__ == "__main__":
    main()
