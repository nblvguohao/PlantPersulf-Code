#!/usr/bin/env python
"""W1 strict-track arm: structure_ranker 5-fold development CV, per species.

Fits ``structure_ranker`` with frozen v11 hyperparameters on each of the five
strict-cluster development folds (identical fold construction to the
production ``prepare_v2_development_fold`` path) and reports per-species
validation AP per fold. Running with two structure registries — the frozen v2
release (tomato masked) and the campaign v3 release (tomato structures) —
quantifies the tomato structure-coverage lever on the strict track.

The frozen 20% test partition is never touched here (the one-shot unlock is
consumed); this is development-CV evidence only.

Usage::

    PYTHONPATH='src;scripts' python \\
        scripts/accuracy_campaign/eval_strict_folds_ranker.py \\
        --config configs/experiments/multispecies_v2_global_clusters_v11.yaml \\
        --output results/.../w1_tomato_structures/strict_folds_v2 \\
        [--structure-registry <registered-release-tsv>]
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import yaml

from plantpersulf.evaluation.comparable_track import (  # noqa: E402
    build_registered_structure_features,
)
from plantpersulf.models.structure_ranker import (  # noqa: E402
    AblationConfig,
    BranchFeatures,
    structure_ranker_scores,
)
from plantpersulf.proteomics.multispecies_v2_dataset import (  # noqa: E402
    MultispeciesV2SiteRow,
    prepare_v2_development_fold,
)
from plantpersulf.proteomics.multispecies_v2_sources import (  # noqa: E402
    sequence_feature_map,
)

PRIMARY_SPECIES = ("arabidopsis", "rice", "tomato")


def _branch_features(
    rows: tuple[MultispeciesV2SiteRow, ...],
    sequence_features: dict[tuple[str, int], tuple[float, ...]],
    structure_features: dict[tuple[str, int], tuple[tuple[float, ...], bool]],
) -> BranchFeatures:
    keys = [(row.global_protein_id, row.cys_position) for row in rows]
    return BranchFeatures(
        sequence=[list(sequence_features[key]) for key in keys],
        esm=[[0.0] for _ in keys],
        structure=[list(structure_features[key][0]) for key in keys],
        structure_mask=[structure_features[key][1] for key in keys],
        study_ids=None,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/multispecies_v2_global_clusters_v11.yaml"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--structure-registry", type=Path, default=None)
    parser.add_argument("--device", type=str, default="cuda:0")
    parser.add_argument("--batch-size", type=int, default=16384)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    from train_multispecies_v2 import _development_rows

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    if not isinstance(cfg, dict):
        raise RuntimeError("multispecies v2 config is required")

    panel, _frozen, proteomes = _development_rows(args.config, sample_unlabeled=True)
    raw_features = sequence_feature_map(panel, proteomes)
    if any(any(value is None for value in values) for values in raw_features.values()):
        raise RuntimeError("sequence features contain missing values")
    sequence_features = {
        key: tuple(float(value) for value in values if value is not None)
        for key, values in raw_features.items()
    }
    inputs_cfg = cfg.get("comparison_inputs", {})
    build_cfg = cfg.get("comparison_feature_build")
    if not isinstance(build_cfg, dict):
        raise RuntimeError("comparison_feature_build configuration is required")
    registry_path = args.structure_registry or Path(
        str(inputs_cfg["structure_features"])
    )
    if args.structure_registry is not None:
        from plantpersulf.provenance.audit import assert_registered_input

        for registry in (
            Path("data/registry/cross_crop_target_label_free_inputs_v1.tsv"),
            Path("data/registry/supplementary_sources.tsv"),
            Path("data/registry/model_inputs.tsv"),
        ):
            try:
                assert_registered_input(args.structure_registry, registry)
                break
            except RuntimeError:
                pass
        else:
            raise RuntimeError(
                f"structure registry override is absent from every approved "
                f"registry: {args.structure_registry}"
            )
    structure_features = build_registered_structure_features(
        set(sequence_features),
        registry_path=registry_path,
        registry_base=Path(str(build_cfg["structure_registry_base"])),
    )

    params = cfg["literature_baseline_parameters"]["structure_ranker"]
    results: list[dict[str, object]] = []
    for fold in range(5):
        dev_fold = prepare_v2_development_fold(
            panel, validation_fold=fold, n_development_folds=5
        )
        output = structure_ranker_scores(
            _branch_features(dev_fold.fit_rows, sequence_features, structure_features),
            [row.label for row in dev_fold.fit_rows],
            _branch_features(
                dev_fold.validation_rows, sequence_features, structure_features
            ),
            seed=fold,
            ablation=AblationConfig(use_esm=False, use_study_context=False),
            hidden=int(params["hidden"]),
            dropout=float(params["dropout"]),
            epochs=int(params["epochs"]),
            lr=float(params["lr"]),
            holdout_fraction=float(params["holdout_fraction"]),
            n_mc_dropout=int(params["n_mc_dropout"]),
            device_name=args.device,
            batch_size=args.batch_size,
        )
        per_species: dict[str, float] = {}
        import sklearn.metrics

        by_species: dict[str, list[tuple[float, bool]]] = {}
        for row, score in zip(
            dev_fold.validation_rows, output.scores, strict=True
        ):
            by_species.setdefault(row.species, []).append(
                (score, row.label == "positive")
            )
        for species, pairs in by_species.items():
            if any(positive for _, positive in pairs):
                per_species[species] = float(
                    sklearn.metrics.average_precision_score(
                        [positive for _, positive in pairs],
                        [score for score, _ in pairs],
                    )
                )
        results.append(
            {
                "fold": fold,
                "fit_rows": len(dev_fold.fit_rows),
                "validation_rows": len(dev_fold.validation_rows),
                "per_species_ap": per_species,
            }
        )
        print(
            f"fold {fold}: fit={len(dev_fold.fit_rows)} "
            f"validation={len(dev_fold.validation_rows)} "
            f"{ {s: round(v, 4) for s, v in per_species.items()} }"
        )

    summary = {
        species: statistics.mean(
            result["per_species_ap"][species]  # type: ignore[index]
            for result in results
            if species in result["per_species_ap"]  # type: ignore[index]
        )
        for species in PRIMARY_SPECIES
    }
    macro = statistics.mean(summary.values())
    print(f"mean fold AP: { {s: round(v, 4) for s, v in summary.items()} }")
    print(f"three-crop macro: {macro:.4f}")

    args.output.mkdir(parents=True, exist_ok=True)
    payload = {
        "structure_registry": str(registry_path),
        "per_fold": results,
        "mean_per_species_ap": summary,
        "three_crop_macro": macro,
    }
    (args.output / "strict_folds_report.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output / 'strict_folds_report.json'}")


if __name__ == "__main__":
    main()
