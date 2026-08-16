#!/usr/bin/env python
"""Accuracy-campaign runner: literature-track structure_ranker per-site scores.

Re-runs the v11 literature random-protein track for ``structure_ranker`` only,
capturing **per-site scores and MC-dropout uncertainty** (the production track
summarizes these away), under a hard reproducibility gate: every rebuilt panel
must match the ``panel_sha256`` recorded in the frozen v11 ``summary.json``.

The runner reuses the exact production call path
(``_development_rows`` + ``build_literature_random_track`` +
``sequence_feature_map`` + ``build_registered_structure_features`` +
``structure_ranker_scores``) and adds no new modeling logic. Campaign-only
variants are selected by flags:

- ``--use-esm 1`` — W3-L2: enable the ESM branch (features already cached).
- ``--structure-registry PATH`` — W1: point at a newer structure registry
  release (e.g. tomato release v3); must be registered in
  ``data/registry/model_inputs.tsv``.
- ``--device`` / ``--batch-size`` — compute overrides (defaults from config).

Outputs one tsv per seed: ``partition, global_protein_id, cys_position,
label, score, uncertainty`` for validation+test rows, plus
``reproduction_report.json`` comparing rebuilt per-species APs against the
frozen v11 summary. Training rows are not scored (the production track does
not score them either; stacking combiners fit on validation rows).

Usage::

    PYTHONPATH='src;scripts' python \\
        scripts/accuracy_campaign/run_literature_ranker_scores.py \\
        --config configs/experiments/multispecies_v2_global_clusters_v11.yaml \\
        --summary results/.../v11/literature_random_protein/summary.json \\
        --output results/experiments/accuracy_campaign/w2_tool_stacking
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
from plantpersulf.evaluation.comparable_track import (  # noqa: E402
    build_registered_structure_features,
)
from plantpersulf.features.esm2_windows import (  # noqa: E402
    load_window_embedding_artifact,
)
from plantpersulf.models.structure_ranker import (  # noqa: E402
    AblationConfig,
    BranchFeatures,
    structure_ranker_scores,
)
from plantpersulf.proteomics.multispecies_v2_sources import (  # noqa: E402
    sequence_feature_map,
)

STRUCTURE_RANKER_PARAMETERS = {
    "epochs": 200,
    "hidden": 16,
    "dropout": 0.2,
    "lr": 0.05,
    "holdout_fraction": 0.2,
    "n_mc_dropout": 16,
}


def load_expected_panel_sha256s(summary_path: Path) -> dict[int, str]:
    """Read the frozen v11 summary and return ``{seed: panel_sha256}``.

    All six models share one panel per seed; the function verifies that
    cross-model consistency before returning the map.
    """
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
    partition: str,
    rows: tuple[tuple[str, int, str], ...],
    scores: list[float],
    uncertainty: list[float],
) -> None:
    """Write ``(global_protein_id, cys_position, label)`` rows with scores.

    ``rows`` elements are ``(global_protein_id, cys_position, label)`` aligned
    with ``scores``/``uncertainty``. Writes atomically (temp file + replace).
    """
    if not (len(rows) == len(scores) == len(uncertainty)):
        raise ValueError("rows, scores, and uncertainty must align")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    lines = ["partition\tglobal_protein_id\tcys_position\tlabel\tscore\tuncertainty"]
    for (protein_id, position, label), score, sigma in zip(
        rows, scores, uncertainty, strict=True
    ):
        lines.append(
            f"{partition}\t{protein_id}\t{position}\t{label}\t{score}\t{sigma}"
        )
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    temporary.replace(path)


def per_species_average_precision(
    rows: list[tuple[str, int, str]], scores: list[float]
) -> dict[str, float]:
    """sklearn AP per species prefix of ``global_protein_id`` (macro track)."""
    import sklearn.metrics

    by_species: dict[str, list[tuple[float, bool]]] = {}
    for (protein_id, _position, label), score in zip(rows, scores, strict=True):
        species = protein_id.split("|", 1)[0]
        by_species.setdefault(species, []).append((score, label == "positive"))
    return {
        species: float(
            sklearn.metrics.average_precision_score(
                [positive for _, positive in pairs],
                [score for score, _ in pairs],
            )
        )
        for species, pairs in by_species.items()
        if any(positive for _, positive in pairs)
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/multispecies_v2_global_clusters_v11.yaml"),
    )
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--use-esm", type=int, default=0, choices=(0, 1))
    parser.add_argument("--structure-registry", type=Path, default=None)
    parser.add_argument("--seeds", type=str, default="0,1,2,3,4,5,6,7,8,9")
    parser.add_argument("--device", type=str, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
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
    # Lazy import: scripts/ is only on sys.path when this file runs as a
    # script (unit tests import the pure helpers without the scripts path).
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

    use_esm = bool(args.use_esm)
    esm_features = None
    if use_esm:
        inputs_cfg = cfg.get("comparison_inputs", {})
        esm_manifest = Path(str(inputs_cfg["esm_features"]))
        esm_features = load_window_embedding_artifact(
            esm_manifest.parent,
            required_site_keys=set(selected_rows),
        )

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
        set(selected_rows),
        registry_path=registry_path,
        registry_base=Path(str(build_cfg["structure_registry_base"])),
    )

    def branches(partition_rows: tuple) -> BranchFeatures:
        sequence = [
            list(sequence_features[(row.global_protein_id, row.cys_position)])
            for row in partition_rows
        ]
        structure_pairs = [
            structure_features[(row.global_protein_id, row.cys_position)]
            for row in partition_rows
        ]
        esm = (
            [
                list(esm_features[(row.global_protein_id, row.cys_position)])
                for row in partition_rows
            ]
            if use_esm and esm_features is not None
            else [[0.0] for _ in partition_rows]
        )
        return BranchFeatures(
            sequence=sequence,
            esm=esm,
            structure=[list(values) for values, _ in structure_pairs],
            structure_mask=[available for _, available in structure_pairs],
            study_ids=None,
        )

    report: dict[str, object] = {
        "config": str(args.config),
        "summary": str(args.summary),
        "use_esm": use_esm,
        "structure_registry": str(registry_path),
        "device": device,
        "seeds": {seed: {} for seed in seeds},
    }
    seeds_report = report["seeds"]
    assert isinstance(seeds_report, dict)
    for run in runs:
        partitions = dict(run.partitions)
        train_rows = partitions["train"]
        predict_rows = (*partitions["validation"], *partitions["test"])
        output = structure_ranker_scores(
            branches(train_rows),
            [row.label for row in train_rows],
            branches(predict_rows),
            seed=run.seed,
            ablation=AblationConfig(use_esm=use_esm, use_study_context=False),
            **STRUCTURE_RANKER_PARAMETERS,
            device_name=device,
            batch_size=batch_size,
        )
        n_validation = len(partitions["validation"])
        for name, rows_slice, start, stop in (
            ("validation", partitions["validation"], 0, n_validation),
            ("test", partitions["test"], n_validation, len(output.scores)),
        ):
            site_rows = tuple(
                (row.global_protein_id, row.cys_position, row.label)
                for row in rows_slice
            )
            write_ranker_scores_tsv(
                args.output / f"structure_ranker_seed{run.seed}_{name}.tsv",
                partition=name,
                rows=site_rows,
                scores=list(output.scores[start:stop]),
                uncertainty=list(output.uncertainty[start:stop]),
            )
            ap = per_species_average_precision(
                list(site_rows), list(output.scores[start:stop])
            )
            seeds_report[run.seed][name] = ap
            print(
                f"seed {run.seed} {name}: "
                f"{ {s: round(v, 4) for s, v in ap.items()} }"
            )

    report_path = args.output / "reproduction_report.json"
    args.output.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {report_path}")


if __name__ == "__main__":
    main()
