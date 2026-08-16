#!/usr/bin/env python
"""Build the Zhang-lab mutation-validation handoff list from a release.

Consumes the release's ``top_k_candidates.tsv`` (ranked site table) and
``matched_controls.tsv`` (candidate/control pairs) and produces ONE compact
handoff TSV for the wet-lab partner: every site with its protein context,
±10-residue sequence window around the Cys, structure availability, and the
matched controls that pair with each candidate.

The validation mode (PI decision 2026-08-16) is Cys→Ala substitution with a
functional readout (NOT MS detection), so the handoff includes the residue
context the lab needs to design the mutant and, for candidates, their
matched controls so the enrichment analysis has a defined comparator set.
K is fixed by the caller (default 200); the list is locked before the lab
receives any results.

Usage::

    python scripts/build_zhang_mutation_validation_list.py \\
        --release-dir results/candidates/multispecies_v2_candidate_release_v1 \\
        --output results/candidates/multispecies_v2_candidate_release_v1/ \\
            zhang_validation_list_v1.tsv \\
        --k 200
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import yaml  # noqa: E402

from plantpersulf.features.sequence import _load_proteome  # noqa: E402

DEFAULT_K = 200


def build_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--release-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=DEFAULT_K)
    return parser


def _window(sequence: str, position: int, radius: int = 10) -> str:
    """±``radius`` residue window around a 1-based position (clipped)."""
    start = max(0, position - 1 - radius)
    end = min(len(sequence), position - 1 + radius + 1)
    return sequence[start:end]


def main(argv: list[str] | None = None) -> None:
    args = build_argument_parser().parse_args(argv)

    candidates_path = args.release_dir / "top_k_candidates.tsv"
    controls_path = args.release_dir / "matched_controls.tsv"
    manifest_path = args.release_dir / "manifest.json"
    if not candidates_path.is_file() or not controls_path.is_file():
        raise RuntimeError(
            f"release artifacts missing in {args.release_dir}: need "
            "top_k_candidates.tsv and matched_controls.tsv"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    config_path = Path(str(manifest["source_config_path"]))
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    species_ref = next(
        item for item in cfg["reference_proteomes"] if item["species"] == "tomato"
    )
    proteome = _load_proteome(Path(str(species_ref["path"])))

    candidates: list[tuple[str, int, float, float]] = []
    with candidates_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            candidates.append(
                (
                    row["site_key"],
                    int(row["cys_position"]),
                    float(row["score"]),
                    float(row["uncertainty"]),
                )
            )
    candidates = candidates[: args.k]

    controls_by_candidate: dict[tuple[str, int], list[tuple[str, int]]] = {}
    with controls_path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            cand_key = (row["candidate_site_key"], int(row["candidate_cys_position"]))
            controls_by_candidate.setdefault(cand_key, []).append(
                (row["control_site_key"], int(row["control_cys_position"]))
            )

    def site_context(site_key: str, position: int) -> dict:
        _, accession = site_key.split("|", 1)
        sequence = proteome[accession]
        return {
            "protein_accession": accession,
            "protein_length": len(sequence),
            "sequence_window": _window(sequence, position),
            "residue": sequence[position - 1],
        }

    rows: list[dict[str, object]] = []
    for rank, (site_key, position, score, uncertainty) in enumerate(
        candidates, start=1
    ):
        context = site_context(site_key, position)
        controls = controls_by_candidate.get((site_key, position), [])
        rows.append(
            {
                "rank": rank,
                "site_key": site_key,
                "cys_position": position,
                "score": f"{score:.10g}",
                "uncertainty": f"{uncertainty:.6g}",
                "n_matched_controls": len(controls),
                "matched_controls": ";".join(
                    f"{key}|{pos}" for key, pos in controls
                ),
                **context,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "rank",
        "site_key",
        "cys_position",
        "protein_accession",
        "protein_length",
        "sequence_window",
        "residue",
        "score",
        "uncertainty",
        "n_matched_controls",
        "matched_controls",
    ]
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    n_with_controls = sum(1 for row in rows if row["n_matched_controls"] > 0)
    n_controls_total = sum(row["n_matched_controls"] for row in rows)
    print(f"handoff list: {len(rows)} candidates (K={args.k}), "
          f"{n_with_controls} with matched controls ({n_controls_total} pairs)")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
