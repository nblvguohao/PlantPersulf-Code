#!/usr/bin/env python
"""Select assayable target proteins for within-protein mutation validation.

The within-protein task needs target PROTEINS the wet lab can actually
assays. The global top-K candidate head is dominated by cysteine-rich small
proteins (metallothioneins, defensins/snakins, fragments, uncharacterised
entries) that have no clean functional readout — the model's high-cys-density
bias, not a usable validation set. This selector filters the tomato proteome
to proteins with a real functional annotation and a meaningful within-protein
contrast, then ranks the survivors by the model's within-protein rank-1 score.

Fail-closed exclusion rules (documented, applied BEFORE any outcome is known):
- no functional annotation (header contains "Uncharacterized" / "uncharacterised");
- "Fragment" entries (incomplete sequences — unreliable assay targets);
- metallothionein / MT2 (cysteine-rich metal-binding, no clean persulfidation
  functional readout);
- Snakin / defensin-like antimicrobial peptides (cysteine-rich small peptides);
- protein length < 150 (small peptides, dominant high-Cys-density class);
- < 3 Cys (within-protein contrast needs at least 3 sites).

Selection = remaining proteins ranked by within-protein rank-1 score
(descending), top-N reported with the full within-protein ranking.
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

DEFAULT_N = 25
MIN_CYS = 3
MIN_LENGTH = 150

_EXCLUDE_PATTERNS = (
    "uncharacter",
    "fragment",
    "metallothionein",
    "mt2",
    "snakin",
    "defensin",
    "cysteine-rich",
)


def _excluded(header: str) -> str | None:
    lowered = header.lower()
    for pattern in _EXCLUDE_PATTERNS:
        if pattern in lowered:
            return pattern
    return None


def _window(sequence: str, position: int, radius: int = 10) -> str:
    start = max(0, position - 1 - radius)
    end = min(len(sequence), position - 1 + radius + 1)
    return sequence[start:end]


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n", type=int, default=DEFAULT_N)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_argument_parser().parse_args(argv)

    cfg = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    species_ref = next(
        item for item in cfg["reference_proteomes"] if item["species"] == "tomato"
    )
    proteome_path = Path(str(species_ref["path"]))
    proteome = _load_proteome(proteome_path)

    # accession -> header description (first UniProt name field)
    descriptions: dict[str, str] = {}
    for line in proteome_path.open(encoding="utf-8"):
        if line.startswith(">"):
            header = line.strip()
            accession = header.split()[0].split("|")[1]
            name_field = header.split("|")[-1].split(" OS=")[0].strip()
            descriptions[accession] = name_field

    by_protein: dict[str, list[tuple[int, float]]] = defaultdict(list)
    with args.candidates.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accession = row["site_key"].split("|", 1)[1]
            by_protein[accession].append(
                (int(row["cys_position"]), float(row["score"]))
            )
    for positions in by_protein.values():
        positions.sort(key=lambda item: item[1], reverse=True)

    kept: list[dict[str, object]] = []
    excluded_counts: dict[str, int] = defaultdict(int)
    for accession, positions in by_protein.items():
        desc = descriptions.get(accession, "")
        reason = _excluded(desc)
        sequence = proteome[accession]
        if len(positions) < MIN_CYS:
            excluded_counts["<3 cys"] += 1
            continue
        if len(sequence) < MIN_LENGTH:
            excluded_counts["short (<150aa)"] += 1
            continue
        if reason is not None:
            excluded_counts[f"excluded: {reason}"] += 1
            continue
        pred_pos, pred_score = positions[0]
        gap = positions[0][1] - positions[1][1] if len(positions) > 1 else None
        kept.append(
            {
                "protein_accession": accession,
                "description": desc[:80],
                "n_cys": len(positions),
                "protein_length": len(sequence),
                "predicted_site": pred_pos,
                "predicted_window": _window(sequence, pred_pos),
                "predicted_score": f"{pred_score:.6g}",
                "rank1_minus_rank2_gap": f"{gap:.6g}" if gap is not None else "",
                "all_cys_ranked": ";".join(
                    f"{pos}:{score:.4g}" for pos, score in positions
                ),
            }
        )
    kept.sort(key=lambda r: -float(r["predicted_score"]))
    kept = kept[: args.n]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "protein_accession",
        "description",
        "protein_length",
        "n_cys",
        "predicted_site",
        "predicted_window",
        "predicted_score",
        "rank1_minus_rank2_gap",
        "all_cys_ranked",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(kept)

    print(f"excluded: {dict(excluded_counts)}")
    print(f"assayable proteins kept: {len(kept)} (top-{args.n})")
    for r in kept:
        print(
            f"  {r['protein_accession']:<14} {r['description'][:42]:<44} "
            f"n_cys={r['n_cys']:>2} pred={r['predicted_site']:>5} "
            f"gap={r['rank1_minus_rank2_gap']}"
        )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
