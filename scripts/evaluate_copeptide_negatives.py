#!/usr/bin/env python
"""Diagnostic-only: does the frozen bundle separate co-peptide positives
from the registered co-peptide negatives?

The registered co-peptide evidence (``data/registry/copeptide_negatives_v1.tsv``,
see ``plantpersulf.evidence.copeptide_negatives``) contains the first explicit
negatives of the project: within ONE detected peptide — one spectrum, one
digestion, one enrichment — some Cys carry the persulfidation modification
and others do not, so detectability is matched by construction. For BRG3 the
three Cys are pairwise 3 residues apart (C206/C209/C212 on peptide
SSCMICLPCR), so a +/-10 sequence-window model sees nearly identical inputs
for all three and can only separate them if its features carry real signal;
for RNF144b the pair (C122 modified / C127 unmodified) shares the peptide
FYCPYKDCSAMLVNDSDEIVR.

This diagnostic scores EVERY Cys of the two proteins with the FROZEN
candidate-release bundle (structure branch masked, exactly the release's
tomato scoring semantics) and reports the within-protein ranks. It changes
no model, feature, Top-K, threshold or SAP; claim class is diagnostic_only.

The binary conclusion per peptide: are all positive Cys ranked above all
negative Cys of the same peptide?
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plantpersulf.evaluation.release_scoring import (  # noqa: E402
    cys_positions,
    score_sites,
)
from plantpersulf.evaluation.within_protein_ranking import rank_sites_desc  # noqa: E402
from plantpersulf.evidence.copeptide_negatives import (  # noqa: E402
    load_copeptide_registry,
)
from plantpersulf.features.sequence import _load_proteome  # noqa: E402
from plantpersulf.models.structure_ranker import StructureRankerBundle  # noqa: E402

BUNDLE_PATH = (
    _REPO_ROOT
    / "results"
    / "candidates"
    / "multispecies_v2_candidate_release_v1"
    / "model_weights"
    / "structure_ranker_bundle.pt"
)
PROTEOMES = {
    "tomato": (
        _REPO_ROOT
        / "data"
        / "raw"
        / "references"
        / "tomato_ref_proteome_v1.fasta"
    ),
    "arabidopsis": (
        _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
    ),
    "rice": (
        _REPO_ROOT
        / "data"
        / "raw"
        / "references"
        / "rice_proteome_v1"
        / "uniprot_rice_v1.fasta"
    ),
}
REGISTRY = _REPO_ROOT / "data" / "registry" / "copeptide_negatives_v1.tsv"
DEFAULT_OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "copeptide_negatives_v1.json"


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)

    rows = load_copeptide_registry(REGISTRY)
    proteins = sorted({(r["species"], r["protein_accession"]) for r in rows})
    print(f"copeptide registry: {len(rows)} rows, {len(proteins)} proteins")

    species_present = {r["species"] for r in rows}
    proteomes = {
        name: _load_proteome(path)
        for name, path in PROTEOMES.items()
        if name in species_present
    }
    bundle = StructureRankerBundle.load(BUNDLE_PATH)

    per_protein: dict[tuple[str, str], dict[int, float]] = {}

    for species, accession in proteins:
        sequence = proteomes[species][accession]
        sites = [
            (species, accession, position) for position in cys_positions(sequence)
        ]
        scored = score_sites(sites, bundle=bundle, proteomes=proteomes)
        per_protein[(species, accession)] = {
            position: scored[(species, accession, position)].score
            for position in cys_positions(sequence)
        }
        print(f"  scored {accession}: {len(scored)} Cys")

    # group registry rows by peptide
    by_peptide: dict[tuple[str, str, str, int, int], list[dict[str, str]]] = {}
    for record in rows:
        key = (
            record["species"],
            record["protein_accession"],
            record["peptide_sequence"],
            int(record["peptide_start"]),
            int(record["peptide_end"]),
        )
        by_peptide.setdefault(key, []).append(record)

    conclusions: list[dict[str, object]] = []
    for (species, accession, peptide, start, end), members in sorted(
        by_peptide.items()
    ):
        protein_scores = per_protein[(species, accession)]
        ranking = rank_sites_desc(protein_scores)
        rank_of = {position: rank for rank, position in enumerate(ranking, start=1)}
        positive = {int(r["cys_position"]) for r in members if r["state"] == "positive"}
        negative = {int(r["cys_position"]) for r in members if r["state"] == "negative"}
        undetermined = {
            int(r["cys_position"]) for r in members if r["state"] == "undetermined"
        }
        best_negative_rank = min((rank_of[p] for p in negative), default=None)
        worst_positive_rank = max((rank_of[p] for p in positive), default=None)
        separated = (
            best_negative_rank is not None
            and worst_positive_rank is not None
            and worst_positive_rank < best_negative_rank
        )
        member_rows = [
            {
                "site_id": r["site_id"],
                "cys_position": int(r["cys_position"]),
                "state": r["state"],
                "score": round(protein_scores[int(r["cys_position"])], 6),
                "within_protein_rank": rank_of[int(r["cys_position"])],
            }
            for r in sorted(members, key=lambda r: int(r["cys_position"]))
        ]
        conclusion = {
            "species": species,
            "protein_accession": accession,
            "peptide_sequence": peptide,
            "peptide_span": f"{start}-{end}",
            "n_protein_cys": len(protein_scores),
            "positive_sites": sorted(positive),
            "negative_sites": sorted(negative),
            "undetermined_sites": sorted(undetermined),
            "positive_ranks": {p: rank_of[p] for p in sorted(positive)},
            "negative_ranks": {p: rank_of[p] for p in sorted(negative)},
            "all_positives_above_all_negatives": separated,
            "site_rows": member_rows,
        }
        conclusions.append(conclusion)
        print(
            f"  {accession} {peptide}: positives {sorted(positive)} ranks "
            f"{[rank_of[p] for p in sorted(positive)]} vs negative {sorted(negative)} "
            f"ranks {[rank_of[p] for p in sorted(negative)]} "
            f"-> separated={separated}"
        )

    n_separated = sum(1 for c in conclusions if c["all_positives_above_all_negatives"])

    # --- permutation null ---------------------------------------------------
    # Within each peptide, shuffle the POS/NEG labels among the co-peptide Cys
    # (preserving each peptide's k_pos / k_neg, exchange only among that
    # peptide's own Cys) and recount "separated". The frozen sequence model
    # sees near-identical +/-10 windows for co-peptide Cys, so their scores
    # are exchangeable and this is the correct null.
    n_perm = 999
    rng = np.random.RandomState(20260815)
    null_counts: list[int] = []
    for _ in range(n_perm):
        count = 0
        for (species, accession, _peptide, _start, _end), members in (
            by_peptide.items()
        ):
            scores = per_protein[(species, accession)]
            positions = [
                int(r["cys_position"]) for r in members
            ]
            k_pos = sum(1 for r in members if r["state"] == "positive")
            rng.shuffle(positions)
            pos_set = set(positions[:k_pos])
            positive = {p for p in positions if p in pos_set}
            negative = {p for p in positions if p not in pos_set}
            if positive and negative:
                separated = (
                    min(scores[p] for p in positive) > max(scores[n] for n in negative)
                )
                count += separated
        null_counts.append(count)
    null_counts.sort()
    p_value = (sum(1 for value in null_counts if value >= n_separated) + 1) / (
        n_perm + 1
    )
    p_left = (sum(1 for value in null_counts if value <= n_separated) + 1) / (
        n_perm + 1
    )

    summary = {
        "track": "copeptide_negatives_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Frozen-bundle within-protein ranking of the registered co-peptide "
            "evidence. Explicit negatives from AGENTS.md-registered evidence; "
            "structure branch masked (release tomato semantics). Does not "
            "modify the model, features, Top-K, thresholds or SAP. Permutation "
            "null (B=999, composition-preserving within-peptide label shuffle) "
            "counts how many peptides separate beyond chance."
        ),
        "bundle": {"path": str(BUNDLE_PATH), "seed": bundle.seed},
        "registry": {
            "path": str(REGISTRY),
            "n_rows": len(rows),
        },
        "peptides": conclusions,
        "permutation": {
            "n_perm": n_perm,
            "seed": 20260815,
            "observed_separated": n_separated,
            "null_median": int(np.median(null_counts)),
            "null_p5": null_counts[len(null_counts) // 20],
            "null_p95": null_counts[
                min(len(null_counts) - 1, 19 * len(null_counts) // 20)
            ],
            "p_separate_beyond_chance": round(p_value, 4),
            "p_at_or_below_chance": round(p_left, 4),
        },
        "summary": {
            "n_peptides": len(conclusions),
            "n_with_all_positives_above_all_negatives": n_separated,
            "random_expectation": round(
                sum(
                    1
                    / math.comb(
                        len(c["positive_sites"]) + len(c["negative_sites"]),
                        len(c["positive_sites"]),
                    )
                    for c in conclusions
                ),
                3,
            ),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"\nwrote {args.output}")


if __name__ == "__main__":
    main()
