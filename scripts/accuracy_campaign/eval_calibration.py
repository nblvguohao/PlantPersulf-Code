#!/usr/bin/env python
"""W3-L8 calibration / selection-rule evaluation on literature-track scores.

Evaluates, per seed and on the **test** partition only, whether the kiae271
diagnostic ("signal is mid-list, top is crowded by feature-extreme sites")
can be addressed by selection rules — WITHOUT changing the frozen model:

- ``raw`` — anchor: raw scores, pooled hits@K and recall@K over the FULL
  original positive count (K=50/200, three primary crops pooled).
- ``unc_gated_q{q}`` — drop the top q% highest-MC-dropout-uncertainty rows
  per species (q in 0.1/0.25/0.5), then re-rank. Reports how many positives
  were dropped, and hits/recall over the full original positive count (a
  dropped positive can never be found, so recall_full cannot be inflated).
- ``calib_pooled`` — isotonic calibration fit on validation rows (pooled),
  monotone transform of test scores, pooled hits@K/recall_full@K
  (per-species AP is rank-invariant under monotone transforms, so it is not
  reported for calibrated schemes).
- ``calib_per_species`` — same, calibrator fit per species.

Usage::

    PYTHONPATH='src;scripts' python scripts/accuracy_campaign/eval_calibration.py \\
        --sr-scores results/experiments/accuracy_campaign/w2_tool_stacking \\
        --output results/experiments/accuracy_campaign/w3_calibration_esm
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))



from scripts.accuracy_campaign.run_literature_ranker_scores import (  # noqa: E402
    per_species_average_precision,
)

PRIMARY_SPECIES = ("arabidopsis", "rice", "tomato")
GATE_QUANTILES = (0.1, 0.25, 0.5)
RECALL_KS = (50, 200)


def read_sr_scores(
    directory: Path, seed: int
) -> dict[str, dict[tuple[str, int], tuple[float, float, str]]]:
    out: dict[str, dict[tuple[str, int], tuple[float, float, str]]] = {}
    for partition in ("validation", "test"):
        path = directory / f"structure_ranker_seed{seed}_{partition}.tsv"
        table: dict[tuple[str, int], tuple[float, float, str]] = {}
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


def primary_rows(
    table: dict[tuple[str, int], tuple[float, float, str]]
) -> list[tuple[tuple[str, int], float, float, str]]:
    """Test rows of the three primary crops as ``(key, score, sigma, label)``."""
    return [
        (key, score, sigma, label)
        for key, (score, sigma, label) in table.items()
        if key[0].split("|", 1)[0] in PRIMARY_SPECIES
    ]


def pooled_hits(
    rows: list[tuple[tuple[str, int], float, float, str]], k: int
) -> int:
    """Absolute positives found in the top-k of the pooled ranking."""
    ranked = sorted(rows, key=lambda item: (-item[1], item[0][0], item[0][1]))
    return sum(1 for _key, _score, _sigma, label in ranked[:k] if label == "positive")


def evaluate_seed(seed: int, directory: Path) -> dict[str, object]:
    sr = read_sr_scores(directory, seed)
    val_rows = primary_rows(sr["validation"])
    test_rows = primary_rows(sr["test"])
    total_positives = sum(
        1 for _key, _score, _sigma, label in test_rows if label == "positive"
    )
    result: dict[str, object] = {"seed": seed, "total_positives": total_positives}

    def metric_block(
        rows: list[tuple[tuple[str, int], float, float, str]]
    ) -> dict[str, dict[int, float | None]]:
        """hits@k (absolute) and recall@k over the FULL original positive count."""
        hits = {k: pooled_hits(rows, k) for k in RECALL_KS}
        recall_full = {
            k: (hits[k] / total_positives if total_positives else None)
            for k in RECALL_KS
        }
        return {"hits": hits, "recall_full": recall_full}

    def ap_block(
        rows: list[tuple[tuple[str, int], float, float, str]]
    ) -> dict[str, float]:
        return per_species_average_precision(
            [(key[0], key[1], label) for key, _score, _sigma, label in rows],
            [score for _, score, _sigma, _label in rows],
        )

    result["raw"] = {
        "metrics": metric_block(test_rows),
        "per_species_ap": ap_block(test_rows),
    }

    # uncertainty-gated selection rules: drop the top q% highest-uncertainty
    # rows per species (q = dropped fraction, not the kept fraction).
    by_species: dict[str, list[tuple[tuple[str, int], float, float, str]]] = {}
    for row in test_rows:
        species = row[0][0].split("|", 1)[0]
        by_species.setdefault(species, []).append(row)
    for q in GATE_QUANTILES:
        kept: list[tuple[tuple[str, int], float, float, str]] = []
        dropped_positives = 0
        for rows in by_species.values():
            sigmas = sorted(row[2] for row in rows)
            if not sigmas:
                continue
            keep_fraction = 1.0 - q
            threshold = sigmas[min(len(sigmas) - 1, int(keep_fraction * len(sigmas)))]
            for key, score, sigma, label in rows:
                if sigma <= threshold:
                    kept.append((key, score, sigma, label))
                elif label == "positive":
                    dropped_positives += 1
        result[f"unc_gated_q{q}"] = {
            "metrics": metric_block(kept),
            "per_species_ap": ap_block(kept),
            "dropped_positives": dropped_positives,
            "kept_rows": len(kept),
        }

    # isotonic calibration (fit on validation rows only)
    from sklearn.isotonic import IsotonicRegression

    def calibrate(
        fit_rows: list[tuple[tuple[str, int], float, float, str]],
        apply_rows: list[tuple[tuple[str, int], float, float, str]],
        per_species: bool,
    ) -> list[tuple[tuple[str, int], float, float, str]]:
        if not per_species:
            calibrator = IsotonicRegression(out_of_bounds="clip")
            calibrator.fit(
                [score for _, score, _sigma, _ in fit_rows],
                [label == "positive" for _, _, _, label in fit_rows],
            )
            return [
                (key, float(calibrator.predict([score])[0]), sigma, label)
                for key, score, sigma, label in apply_rows
            ]
        by_species_fit: dict[str, list[tuple[tuple[str, int], float, float, str]]] = {}
        for row in fit_rows:
            by_species_fit.setdefault(row[0][0].split("|", 1)[0], []).append(row)
        calibrated: dict[str, list[tuple[tuple[str, int], float, float, str]]] = {}
        for species, rows in by_species_fit.items():
            calibrator = IsotonicRegression(out_of_bounds="clip")
            calibrator.fit(
                [score for _, score, _sigma, _ in rows],
                [label == "positive" for _, _, _, label in rows],
            )
            calibrated[species] = [
                (key, float(calibrator.predict([score])[0]), sigma, label)
                for key, score, sigma, label in apply_rows
                if key[0].split("|", 1)[0] == species
            ]
        return [row for species in sorted(calibrated) for row in calibrated[species]]

    for name, per_species in (("calib_pooled", False), ("calib_per_species", True)):
        calibrated_rows = calibrate(val_rows, test_rows, per_species)
        # Per-species AP is rank-invariant under these monotone transforms
        # (up to plateau ties), so calibration is evaluated on pooled
        # recall@K only; per-species AP is not reported here.
        result[name] = {"metrics": metric_block(calibrated_rows)}
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sr-scores", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seeds", type=str, default="0,1,2,3,4,5,6,7,8,9")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    seeds = tuple(int(item) for item in args.seeds.split(",") if item)
    results = [evaluate_seed(seed, args.sr_scores) for seed in seeds]

    summary: dict[str, dict[str, object]] = {}
    all_schemes = (
        "raw",
        *(f"unc_gated_q{q}" for q in GATE_QUANTILES),
        "calib_pooled",
        "calib_per_species",
    )
    for scheme in all_schemes:
        means: dict[str, object] = {}
        for k in RECALL_KS:
            values = [
                result[scheme]["metrics"]["recall_full"][k]  # type: ignore[index]
                for result in results
                if result[scheme]["metrics"]["recall_full"][k] is not None  # type: ignore[index]
            ]
            hits = [
                result[scheme]["metrics"]["hits"][k]  # type: ignore[index]
                for result in results
            ]
            missing = len(results) - len(values)
            means[f"recall_full@{k}"] = statistics.mean(values) if values else None
            means[f"hits@{k}"] = statistics.mean(hits)
            if missing:
                means[f"recall_full@{k}_missing_seeds"] = missing
        if scheme.startswith("unc_gated"):
            dropped = [
                result[scheme]["dropped_positives"]  # type: ignore[index]
                for result in results
            ]
            means["dropped_positives_mean"] = statistics.mean(dropped)
        tomato_values = [
            result[scheme]["per_species_ap"].get("tomato", float("nan"))  # type: ignore[index]
            for result in results
            if "per_species_ap" in result[scheme]  # type: ignore[index]
        ]
        tomato_ap_mean = statistics.mean(
            [v for v in tomato_values if v == v]  # drop NaN
        ) if tomato_values else float("nan")
        summary[scheme] = {**means, "tomato_ap_mean": tomato_ap_mean}
        dropped_note = (
            f" dropped_pos {means['dropped_positives_mean']:.0f}"
            if "dropped_positives_mean" in means
            else ""
        )
        print(
            f"{scheme}: recall_full@50 {means['recall_full@50']:.4f} "
            f"hits@50 {means['hits@50']:.1f} "
            f"recall_full@200 {means['recall_full@200']:.4f} "
            f"tomato_ap {tomato_ap_mean:.4f}{dropped_note}"
        )

    args.output.mkdir(parents=True, exist_ok=True)
    payload = {"seeds": list(seeds), "summary": summary, "per_seed": results}
    (args.output / "calibration_report.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output / 'calibration_report.json'}")


if __name__ == "__main__":
    main()
