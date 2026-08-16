#!/usr/bin/env python
"""Zhang-lab handoff: within-protein site-localization lists.

The task is site LOCALIZATION within known proteins — "which Cys of THIS
protein" — not a proteome-wide Top-K. This builder takes a scored candidate
table (every tomato Cys site, from ``top_k_candidates.tsv``), groups sites
by protein, and emits per-protein ranked Cys lists:

- the model's predicted site = the protein's rank-1 Cys;
- the within-protein controls = the protein's other Cys (rank-2..N), the
  natural matched comparator for a Cys→Ala functional readout;
- a confidence flag (score gap between rank-1 and rank-2) so the lab can
  prioritise proteins where the model's localization call is unambiguous.

Target proteins = those containing a site in the global top-K (default 200);
each such protein's FULL ranking is reported so the localization prediction
is visible, not just the single top site.

Usage::

    python scripts/build_zhang_within_protein_list.py \\
        --candidates results/candidates/multispecies_v2_candidate_release_v1/ \\
            top_k_candidates.tsv \\
        --output results/candidates/multispecies_v2_candidate_release_v1/ \\
            zhang_within_protein_list_v1.tsv \\
        --k 200
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import yaml  # noqa: E402

from plantpersulf.features.sequence import _load_proteome  # noqa: E402

DEFAULT_K = 200
CONTROLS_PER_PROTEIN = 3


def _window(sequence: str, position: int, radius: int = 10) -> str:
    start = max(0, position - 1 - radius)
    end = min(len(sequence), position - 1 + radius + 1)
    return sequence[start:end]


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=DEFAULT_K)
    parser.add_argument(
        "--controls-per-protein", type=int, default=CONTROLS_PER_PROTEIN
    )
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_argument_parser().parse_args(argv)

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    species_ref = next(
        item for item in cfg["reference_proteomes"] if item["species"] == "tomato"
    )
    proteome = _load_proteome(Path(str(species_ref["path"])))

    # site_key, position, score, uncertainty
    scored: list[tuple[str, int, float, float]] = []
    with args.candidates.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            scored.append(
                (
                    row["site_key"],
                    int(row["cys_position"]),
                    float(row["score"]),
                    float(row["uncertainty"]),
                )
            )
    ranked = sorted(scored, key=lambda item: item[2], reverse=True)
    top_k_keys = {(key, pos) for key, pos, _score, _unc in ranked[: args.k]}

    by_protein: dict[str, list[tuple[int, float, float]]] = defaultdict(list)
    for key, pos, score, unc in scored:
        by_protein[key].append((pos, score, unc))
    for positions in by_protein.values():
        positions.sort(key=lambda item: item[1], reverse=True)

    rows: list[dict[str, object]] = []
    for protein, positions in sorted(by_protein.items()):
        has_topk = any(
            (protein, pos) in top_k_keys for pos, _score, _unc in positions
        )
        if not has_topk:
            continue
        accession = protein.split("|", 1)[1]
        sequence = proteome[accession]
        predicted_pos, predicted_score, _ = positions[0]
        gap = positions[0][1] - positions[1][1] if len(positions) > 1 else None
        controls = positions[1 : 1 + args.controls_per_protein]
        rows.append(
            {
                "protein_accession": accession,
                "n_cys": len(positions),
                "predicted_site": predicted_pos,
                "predicted_window": _window(sequence, predicted_pos),
                "predicted_score": f"{predicted_score:.6g}",
                "rank1_minus_rank2_score_gap": (
                    f"{gap:.6g}" if gap is not None else ""
                ),
                "confidence_flag": (
                    "confident" if gap is not None and gap > 0.001 else "borderline"
                ),
                "within_protein_controls": ";".join(
                    f"{pos}:{_window(sequence, pos)}:{score:.4g}"
                    for pos, score, _unc in controls
                ),
                "top_site_in_global_topk": (
                    "yes" if (protein, predicted_pos) in top_k_keys else "no"
                ),
                "all_cys_ranked": ";".join(
                    f"{pos}:{score:.4g}" for pos, score, _unc in positions
                ),
            }
        )

    rows.sort(key=lambda r: -float(r["predicted_score"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "protein_accession",
        "n_cys",
        "predicted_site",
        "predicted_window",
        "predicted_score",
        "rank1_minus_rank2_score_gap",
        "confidence_flag",
        "within_protein_controls",
        "top_site_in_global_topk",
        "all_cys_ranked",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_confident = sum(1 for r in rows if r["confidence_flag"] == "confident")
    n_control_proteins = sum(1 for r in rows if r["within_protein_controls"])
    print(
        f"within-protein handoff: {len(rows)} proteins (from top-{args.k} "
        f"candidates), {n_confident} confident localization calls, "
        f"{n_control_proteins} with >=1 within-protein control"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
