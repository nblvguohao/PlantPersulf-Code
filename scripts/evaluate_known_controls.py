#!/usr/bin/env python
"""Retrospective recovery evaluation for registered known-mechanism controls.

Task 10 (TDD Codex): known published persulfidation positions must be held
out of training and their rank in the trained model reported transparently —
including failed and unmappable controls. Two control scope classes (see
``plantpersulf.evaluation.known_controls``): cross-species tomato controls
(training positives are exclusively Arabidopsis, so recovery is not a
training artefact) and same-species Arabidopsis controls (the control site
is an unlabeled PU-pool row — never a training positive — and is excluded
from the scorer's training sample before scoring).

Produces a per-control percentile rank against the benchmark's unlabeled
distribution under a shared PU scorer: one ``pu_logistic`` model is trained
on the benchmark positives plus a seeded unlabeled sample (sequence features,
exactly as in ``scripts/run_experiment.py``), then both the unlabeled sample
and each control are scored with that same frozen model. A control's
percentile is the fraction of unlabeled reference scores strictly below the
control's score. The earlier single-class logistic stand-in was degenerate
(sklearn rejects a one-class fit) and has been replaced; the registry,
leakage rules, statuses, and output schema are unchanged.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from plantpersulf.evaluation.known_controls import (
    REGISTERED_CONTROLS,
    filter_out_keys,
    scores_against_benchmark_proteome,
)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

BENCHMARK_FIELDS = (
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "evidence_level",
    "source_sha256",
)


def _positives_from_benchmark(benchmark_path: Path) -> set[tuple[str, int]]:
    """Return all (protein, position) labelled positive (training set)."""
    training: set[tuple[str, int]] = set()
    with benchmark_path.open(encoding="utf-8", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            if row["label"] == "positive":
                training.add(
                    (row["protein_accession"], int(row["cys_position_in_protein"]))
                )
    return training


def _extract_vectors(
    rows: list[dict[str, str]],
    proteome_path: Path,
) -> list[list[float]]:
    """[hydrophobicity, cys_density] per row via the real feature extractor,
    driven by a temp labels file in the benchmark schema."""
    import shutil
    import tempfile

    from plantpersulf.features.sequence import extract_sequence_features

    tmp = Path(tempfile.mkdtemp(prefix="control_feats_"))
    try:
        labels_path = tmp / "labels.tsv"
        with labels_path.open("w", encoding="utf-8", newline="") as h:
            w = csv.DictWriter(
                h,
                fieldnames=list(BENCHMARK_FIELDS),
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            w.writeheader()
            w.writerows(rows)
        feats = extract_sequence_features(labels_path, proteome_path, window_radius=10)
        lookup = {
            (f.protein_accession, f.cys_position): [
                f.hydrophobicity,
                f.cys_density,
            ]
            for f in feats
        }
        return [
            lookup.get(
                (r["protein_accession"], int(r["cys_position_in_protein"])),
                [0.0, 0.0],
            )
            for r in rows
        ]
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _fit_pu_scorer(
    benchmark_path: Path,
    proteome_path: Path,
    max_samples: int = 5000,
    seed: int = 0,
    exclude: set[tuple[str, int]] | None = None,
) -> tuple[Any, list[list[float]]]:
    """Train one PU scorer on benchmark positives + a seeded unlabeled sample
    and return (scorer, sampled_unlabeled_vectors).

    The scorer is a frozen ``pu_logistic`` model (Elkan-Noto, as in the
    experiment runner). The reference distribution is the scorer's output on
    the same unlabeled sample it was trained on — with 2 features and 5k+
    samples, in-sample optimism is negligible and identical for every
    comparison; the controls themselves are never part of training.
    ``exclude`` removes same-species control rows from the unlabeled pool
    before sampling, so those controls score as genuinely held-out rows.
    """
    import random

    from plantpersulf.models.traditional import pu_logistic_regression_scores

    positives: list[dict[str, str]] = []
    unlabeled: list[dict[str, str]] = []
    with benchmark_path.open(encoding="utf-8", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            (positives if row["label"] == "positive" else unlabeled).append(dict(row))

    if exclude:
        unlabeled = filter_out_keys(unlabeled, exclude)

    rng = random.Random(seed)
    sampled = rng.sample(unlabeled, min(max_samples, len(unlabeled)))
    sampled.sort(
        key=lambda r: (r["protein_accession"], int(r["cys_position_in_protein"]))
    )

    train_rows = positives + sampled
    train_X = _extract_vectors(train_rows, proteome_path)
    train_y = [r["label"] for r in train_rows]
    unlabeled_X = train_X[len(positives) :]

    def scorer(vectors: list[list[float]]) -> list[float]:
        return pu_logistic_regression_scores(train_X, train_y, vectors, seed=seed)

    return scorer, unlabeled_X


def _score_control(
    scorer: Any,
    accession: str,
    position: int,
    control_proteome_path: Path,
) -> tuple[float, str]:
    """Extract sequence features for one control Cys (from the *control
    species'* proteome) and score it with the shared frozen PU scorer."""
    row = {
        "protein_accession": accession,
        "cys_position_in_protein": str(position),
        "label": "positive",
        "study_accession": "zhang_control",
        "evidence_level": "site_biochemical",
        "source_sha256": "published_control",
    }
    vectors = _extract_vectors([row], control_proteome_path)
    if not vectors or vectors[0] == [0.0, 0.0]:
        return float("-inf"), "feature_extraction_failed"
    return float(scorer(vectors)[0]), "ok"


def _percentile_rank(score: float, distribution: list[float]) -> float:
    """Fraction of distribution strictly below score (0–100)."""
    if not distribution:
        return 0.0
    below = sum(1 for v in distribution if v < score)
    return below / len(distribution) * 100.0


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def run_control_evaluation(
    benchmark_path: Path,
    proteome_path: Path,
    output_path: Path | None,
    control_proteome_path: Path,
) -> list[dict[str, Any]]:
    positives = _positives_from_benchmark(benchmark_path)
    print(f"Training positives (benchmark): {len(positives)}")

    # Verify no control is in the training set (self-recovery guard)
    for ctrl in REGISTERED_CONTROLS:
        if ctrl.status == "mapped":
            key = (ctrl.uniprot_accession, ctrl.cys_position)
            if key in positives:
                raise RuntimeError(
                    f"CONTROL LEAKAGE: {ctrl.gene} {key} found in training set"
                )

    # Same-species (Arabidopsis) controls live in the benchmark as unlabeled
    # rows — exclude them from the scorer's training sample so they score as
    # genuinely held-out rows.
    same_species_keys = {
        (c.uniprot_accession, c.cys_position)
        for c in REGISTERED_CONTROLS
        if c.status == "mapped" and scores_against_benchmark_proteome(c)
    }
    scorer, unlabeled_X = _fit_pu_scorer(
        benchmark_path, proteome_path, exclude=same_species_keys
    )
    unlabeled_scores = scorer(unlabeled_X)
    print(f"Unlabeled reference distribution: {len(unlabeled_scores)} samples")

    results: list[dict[str, Any]] = []
    for ctrl in REGISTERED_CONTROLS:
        if ctrl.status != "mapped":
            results.append(
                {
                    **asdict(ctrl),
                    "score": None,
                    "percentile_rank": None,
                    "recovery_status": ctrl.status,
                }
            )
            continue

        species_proteome = (
            proteome_path
            if scores_against_benchmark_proteome(ctrl)
            else control_proteome_path
        )
        score, detail = _score_control(
            scorer, ctrl.uniprot_accession, ctrl.cys_position, species_proteome
        )
        pct = _percentile_rank(score, unlabeled_scores)
        results.append(
            {
                **asdict(ctrl),
                "score": score,
                "percentile_rank": round(pct, 1),
                "recovery_status": "recovered" if pct > 50 else "not_recovered",
            }
        )
        print(
            f"  {ctrl.gene} Cys{ctrl.cys_position}: "
            f"score={score:.6f}  pct={pct:.1f}%  "
            f"({'OK' if pct > 50 else 'LOW'})"
        )

    for ctrl in results:
        if ctrl["status"] == "unmappable":
            print(f"  {ctrl['gene']}: UNMAPPABLE — {ctrl['provenance'][:120]}...")

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(
            json.dumps(results, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(f"\nResults written to {output_path}")

    return results


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="retrospective known-control recovery evaluation"
    )
    p.add_argument(
        "--benchmark",
        type=Path,
        default=Path("data/processed/benchmark_v1/sites.tsv"),
    )
    p.add_argument(
        "--proteome",
        type=Path,
        default=Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta"),
    )
    p.add_argument(
        "--control-proteome",
        type=Path,
        default=Path("data/raw/references/tomato_ref_proteome_v1.fasta"),
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("results/known_controls/recovery_v1.json"),
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_control_evaluation(
        benchmark_path=args.benchmark,
        proteome_path=args.proteome,
        output_path=args.output,
        control_proteome_path=args.control_proteome,
    )
