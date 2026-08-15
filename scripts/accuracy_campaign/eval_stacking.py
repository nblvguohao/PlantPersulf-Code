#!/usr/bin/env python
"""W2 stacking evaluation: structure_ranker x Sul-BertGRU on the literature track.

Per seed, joins the campaign re-run's per-site structure_ranker scores with the
frozen v11 Sul-BertGRU adapter scores (same shared panel), then evaluates four
schemes on the **test** partition only:

- ``sr_only`` — structure_ranker alone (must reproduce the frozen v11 summary).
- ``sul_only`` — Sul-BertGRU alone.
- ``rank_aggregate`` — training-free mean-rank aggregation of the two models.
- ``val_fit_combiner`` — convex combination weight grid-fit on the validation
  partition (never sees test rows).

The reproduction gate compares ``sr_only`` test APs against the frozen
``summary.json`` values; any divergence is reported per seed, never silenced.

Usage::

    PYTHONPATH='src;scripts' python scripts/accuracy_campaign/eval_stacking.py \\
        --sr-scores results/experiments/accuracy_campaign/w2_tool_stacking \\
        --sul-scores results/.../v11/literature_random_protein \\
        --summary results/.../v11/literature_random_protein/summary.json \\
        --output results/experiments/accuracy_campaign/w2_tool_stacking
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from plantpersulf.evaluation.stacking import (  # noqa: E402
    fit_score_combiner,
    format_site_key,
    rank_aggregate_scores,
)
from scripts.accuracy_campaign.run_literature_ranker_scores import (  # noqa: E402
    per_species_average_precision,
)

PRIMARY_SPECIES = ("arabidopsis", "rice", "tomato")

SiteKey = tuple[str, int]  # (global_protein_id, cys_position)


def read_sul_scores(path: Path) -> dict[str, dict[SiteKey, float]]:
    """Read a Sul-BertGRU adapter tsv into ``{partition: {(pid, pos): score}}``."""
    out: dict[str, dict[SiteKey, float]] = {}
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            key = (row["global_protein_id"], int(row["cys_position"]))
            out.setdefault(row["partition"], {})[key] = float(row["score"])
    return out


def read_sr_scores(
    directory: Path, seed: int
) -> dict[str, dict[SiteKey, tuple[float, float, str]]]:
    """Read campaign re-run tsvs into ``{partition: {(pid, pos): ...}}``.

    sigma is the MC-dropout uncertainty, label the PU label.
    """
    out: dict[str, dict[SiteKey, tuple[float, float, str]]] = {}
    for partition in ("validation", "test"):
        path = directory / f"structure_ranker_seed{seed}_{partition}.tsv"
        table: dict[SiteKey, tuple[float, float, str]] = {}
        with path.open(encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                table[(row["global_protein_id"], int(row["cys_position"]))] = (
                    float(row["score"]),
                    float(row["uncertainty"]),
                    row["label"],
                )
        out[partition] = table
    return out


def macro_ap(per_species: dict[str, float]) -> float:
    values = [
        per_species[species] for species in PRIMARY_SPECIES if species in per_species
    ]
    if not values:
        return float("nan")
    return sum(values) / len(values)


def ordered_partition(
    sr_partition: dict[SiteKey, tuple[float, float, str]],
    sul_partition: dict[SiteKey, float],
    seed: int,
    partition: str,
) -> tuple[list[SiteKey], list[float], list[float], list[str]]:
    """Align sr and sul maps on identical site keys (sorted, fail-loud)."""
    keys = sorted(sr_partition)
    if set(sul_partition) != set(keys):
        raise RuntimeError(
            f"seed {seed} {partition}: sul/sr site-key sets differ "
            f"(sr-only={len(set(keys) - set(sul_partition))}, "
            f"sul-only={len(set(sul_partition) - set(keys))})"
        )
    return (
        keys,
        [sr_partition[key][0] for key in keys],
        [sul_partition[key] for key in keys],
        [sr_partition[key][2] for key in keys],
    )


def evaluate_seed(
    seed: int,
    sr_directory: Path,
    sul_directory: Path,
    summary: dict[str, dict],
) -> dict[str, object]:
    sul = read_sul_scores(
        sul_directory / f"sul_bertgru_seed{seed}" / "scores_batch128.tsv"
    )
    sr = read_sr_scores(sr_directory, seed)

    val_keys, val_sr, val_sul, val_labels = ordered_partition(
        sr["validation"], sul["validation"], seed, "validation"
    )
    alpha = fit_score_combiner(tuple(val_sr), tuple(val_sul), tuple(val_labels))

    test_keys, test_sr_scores, test_sul_scores, test_labels = ordered_partition(
        sr["test"], sul["test"], seed, "test"
    )
    test_rows = [
        (pid, pos, label)
        for (pid, pos), label in zip(test_keys, test_labels, strict=True)
    ]
    sr_test_ap = per_species_average_precision(test_rows, test_sr_scores)
    sul_test_ap = per_species_average_precision(test_rows, test_sul_scores)
    rank_combined = rank_aggregate_scores(
        tuple(format_site_key(pid, pos) for pid, pos in test_keys),
        (tuple(test_sr_scores), tuple(test_sul_scores)),
    )
    rank_test_ap = per_species_average_precision(test_rows, list(rank_combined))
    combined = tuple(
        alpha * a + (1.0 - alpha) * b
        for a, b in zip(test_sr_scores, test_sul_scores, strict=True)
    )
    combo_test_ap = per_species_average_precision(test_rows, list(combined))

    # reproduction gate vs frozen v11 summary
    expected = summary["structure_ranker"][seed]["three_crop_primary_test"]
    mismatches = {
        species: (expected[species]["average_precision"], sr_test_ap.get(species))
        for species in PRIMARY_SPECIES
        if expected[species]["average_precision"] != sr_test_ap.get(species)
    }
    if mismatches:
        print(
            f"  WARNING seed {seed} sr_only test AP diverges "
            f"from frozen v11: {mismatches}"
        )
    else:
        print(f"  seed {seed} sr_only test AP reproduces frozen v11 exactly")

    return {
        "seed": seed,
        "combiner_alpha": alpha,
        "reproduction_mismatches": mismatches,
        "test": {
            "sr_only": {"per_species": sr_test_ap, "macro": macro_ap(sr_test_ap)},
            "sul_only": {"per_species": sul_test_ap, "macro": macro_ap(sul_test_ap)},
            "rank_aggregate": {
                "per_species": rank_test_ap,
                "macro": macro_ap(rank_test_ap),
            },
            "val_fit_combiner": {
                "alpha": alpha,
                "per_species": combo_test_ap,
                "macro": macro_ap(combo_test_ap),
            },
        },
    }


def load_summary_by_seed(summary_path: Path) -> dict[str, dict]:
    raw = json.loads(summary_path.read_text(encoding="utf-8"))
    out: dict[str, dict] = {}
    for run in raw["runs"]:
        out.setdefault(run["model"], {})[int(run["seed"])] = run
    return out


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sr-scores", type=Path, required=True)
    parser.add_argument("--sul-scores", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=str, default="0,1,2,3,4,5,6,7,8,9")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    seeds = tuple(int(item) for item in args.seeds.split(",") if item)
    summary = load_summary_by_seed(args.summary)
    if "structure_ranker" not in summary:
        raise RuntimeError("summary lacks structure_ranker runs")
    results: list[dict[str, object]] = []
    for seed in seeds:
        print(f"evaluating seed {seed}")
        results.append(evaluate_seed(seed, args.sr_scores, args.sul_scores, summary))

    schemes = ("sr_only", "sul_only", "rank_aggregate", "val_fit_combiner")
    table: dict[str, dict[str, object]] = {}
    for scheme in schemes:
        macros = [result["test"][scheme]["macro"] for result in results]  # type: ignore[index]
        table[scheme] = {
            "macro_mean": statistics.mean(macros),
            "macro_std": statistics.stdev(macros) if len(macros) > 1 else 0.0,
            "macro_per_seed": macros,
            "tomato_mean": statistics.mean(
                result["test"][scheme]["per_species"].get("tomato", float("nan"))  # type: ignore[index]
                for result in results
            ),
        }
        print(
            f"{scheme}: macro {table[scheme]['macro_mean']:.4f} "
            f"± {table[scheme]['macro_std']:.4f}"
        )

    args.output.mkdir(parents=True, exist_ok=True)
    payload = {"seeds": list(seeds), "schemes": table, "per_seed": results}
    (args.output / "stacking_report.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output / 'stacking_report.json'}")


if __name__ == "__main__":
    main()
