#!/usr/bin/env python
"""Leave-one-out + permutation validation of the regime-routed composite.

``evaluate_structure_regime_separation.py`` reported the signed composite as a
LABELED CEILING: signs were learned on all 12 proteins, so the burden (37
protein-wide, 32 regime-local) overstates what a real direction learner could
deliver out-of-sample. This diagnostic replaces the ceiling with an estimate:

- **leave-one-protein-out sign learning** — the signs for each held-out
  protein come ONLY from the other 11, then the held-out true site is ranked;
- two ranking scopes: ``protein`` (all Cys) and ``regime`` (Cys in the true
  site's own pLDDT bucket);
- two sign-training rules: ``global`` (all 11 others) and ``regime_grouped``
  (only same-regime others — the routing hypothesis, proposition 3);
- a **permutation null** (B=999) that reassigns the true site uniformly
  within each protein, re-learns signs, and records the LOO burden — the
  exchangeability reference for the observed value.

If the observed LOO burden sits in the left tail of the null, the structure
features carry real, out-of-sample rank signal for true-site identity.

Claim class ``diagnostic_only``; touches no frozen artifact; supervised
learning still not estimable at n=12.
"""

from __future__ import annotations

import json
import sys
import zlib
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plantpersulf.evaluation.known_controls import REGISTERED_CONTROLS  # noqa: E402
from plantpersulf.evaluation.structure_features import (  # noqa: E402
    STRUCTURE_FEATURE_NAMES,
    cys_structure_features,
    parse_pdb,
)
from plantpersulf.evaluation.structure_regime import (  # noqa: E402
    _in_regime_bucket,
    loo_composite_burden,
    permutation_null,
    regime_bucket,
    site_plddt_regime,
    subset_positions,
)
from plantpersulf.features.sequence import _load_proteome  # noqa: E402

ALPHAFOLD_DIR = _REPO_ROOT / "data" / "raw" / "alphafold"
TOMATO_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
)
ARABIDOPSIS_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
)
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "structure_regime_loo_v1.json"
N_PERM = 999
SEED = 20260815

SPECIES_PROTEOMES = {
    "Solanum lycopersicum": ("tomato", TOMATO_PROTEOME),
    "Arabidopsis thaliana": ("arabidopsis", ARABIDOPSIS_PROTEOME),
}


def _single_feature_burden(
    records: list[tuple[dict[int, dict[str, float]], int]],
    feature: str,
    scope_mode: str,
) -> int:
    """Total first-hit burden of ranking by one feature (desc), optionally
    within each true site's own pLDDT bucket. No learned signs — the pure
    unsupervised analogue that underlies the level-3 regime-local result."""
    total = 0
    for rows, true_position in records:
        values = {position: float(row[feature]) for position, row in rows.items()}
        if scope_mode == "protein":
            ranking = sorted(values, key=lambda p: values[p], reverse=True)
        else:
            bucket = regime_bucket(float(rows[true_position]["plddt"]))
            scope = subset_positions(rows, _in_regime_bucket(bucket))
            ranking = sorted(
                scope, key=lambda p: values[p], reverse=True
            )
        total += ranking.index(true_position) + 1
    return total


def main() -> None:
    controls = [c for c in REGISTERED_CONTROLS if c.status == "mapped"]
    proteomes: dict[str, dict[str, str]] = {}
    for _species, (name, path) in SPECIES_PROTEOMES.items():
        proteomes[name] = _load_proteome(path)

    records: list[tuple[dict[int, dict[str, float]], int, str]] = []
    for control in controls:
        accession = control.uniprot_accession
        species, _ = SPECIES_PROTEOMES[control.control_species]
        sequence = proteomes[species][accession]
        pdb_path = ALPHAFOLD_DIR / f"AF-{accession}-F1-model.pdb"
        residues = parse_pdb(pdb_path.read_text(encoding="utf-8"))
        if len(residues) != len(sequence):
            raise RuntimeError(
                f"{accession}: PDB {len(residues)} != proteome {len(sequence)}"
            )
        mismatches = [
            position
            for position in range(1, len(sequence) + 1)
            if residues[position - 1].type != sequence[position - 1]
        ]
        if mismatches:
            raise RuntimeError(
                f"{accession}: PDB/proteome mismatch at {mismatches[:10]}"
            )
        positions = [i + 1 for i, residue in enumerate(sequence) if residue == "C"]
        feature_rows = cys_structure_features(residues, positions)
        records.append((feature_rows, control.cys_position, accession))

    base = [(rows, true) for rows, true, _ in records]
    accessions = [accession for _, _, accession in records]
    regimes = [site_plddt_regime(rows, true) for rows, true, _ in records]

    # --- LOO variants -------------------------------------------------------
    variants = {
        "global_protein": loo_composite_burden(
            base, STRUCTURE_FEATURE_NAMES, scope_mode="protein"
        ),
        "global_regime": loo_composite_burden(
            base, STRUCTURE_FEATURE_NAMES, scope_mode="regime"
        ),
        "regime_grouped_regime": loo_composite_burden(
            base,
            STRUCTURE_FEATURE_NAMES,
            scope_mode="regime",
            group_key=lambda i: regimes[i],
        ),
    }

    # --- permutation nulls --------------------------------------------------
    rng = np.random.RandomState(SEED)
    nulls = {
        name: permutation_null(
            base,
            STRUCTURE_FEATURE_NAMES,
            scope_mode,
            N_PERM,
            rng,
        )
        for name, scope_mode in (
            ("global_protein", "protein"),
            ("global_regime", "regime"),
            ("regime_grouped_regime", "regime"),
        )
    }

    # --- single-feature observed + permutation nulls (no learned signs) -----
    # The level-3 finding (regime-local ranking helps every feature) is
    # unsupervised, so its exchangeability reference is a plain permutation
    # of true sites WITHOUT sign learning.
    single_feature: dict[str, dict[str, object]] = {}
    for feature in STRUCTURE_FEATURE_NAMES:
        for scope_mode in ("protein", "regime"):
            observed = _single_feature_burden(base, feature, scope_mode)
            key = f"{feature}:{scope_mode}"
            # deterministic per-key seed (hash() is randomized per process)
            null_rng = np.random.RandomState(SEED + zlib.crc32(key.encode()))
            null = []
            for _ in range(N_PERM):
                perm_records = [
                    (rows, int(null_rng.choice(list(rows)))) for rows, _ in base
                ]
                null.append(
                    _single_feature_burden(perm_records, feature, scope_mode)
                )
            null.sort()
            p_value = (
                sum(1 for value in null if value <= observed) + 1
            ) / (N_PERM + 1)
            single_feature[key] = {
                "observed_burden": observed,
                "random_baseline": round(
                    sum((len(rows) + 1) / 2.0 for rows, _ in base), 3
                ),
                "permutation_p": round(p_value, 4),
                "null_median": int(null[len(null) // 2]),
                "null_p5": null[len(null) // 20],
            }

    # --- p-values and left-tail position -----------------------------------
    def _pvalue(null: list[int], observed: int) -> float:
        # +1 on both sides: the observed value is itself one draw
        return (sum(1 for value in null if value <= observed) + 1) / (len(null) + 1)

    stats: dict[str, dict[str, object]] = {}
    for name, res in variants.items():
        observed = res["total_first_hit_burden"]
        null = nulls[name]
        stats[name] = {
            "observed_burden": observed,
            "random_baseline": res["total_random_burden"],
            "top1": res["top1"],
            "permutation_p": round(_pvalue(null, observed), 4),
            "null": {
                "median": int(np.median(null)),
                "p5": null[max(0, len(null) // 20)],
                "p95": null[min(len(null) - 1, 19 * len(null) // 20)],
                "min": null[0],
            },
        }

    # --- per-protein LOO ranks (global_regime = headline estimate) ----------
    per_protein_ranks = variants["global_regime"]["per_protein"]

    document = {
        "track": "structure_regime_loo_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "LOO + permutation validation of the regime-routed signed "
            "composite. Signs for each held-out protein are learned only from "
            "the other 11 (global) or from same-regime others only "
            "(regime_grouped, proposition-3 routing). scope_mode 'regime' "
            "ranks within the true site's own pLDDT bucket. The permutation "
            "null reassigns true sites uniformly within each protein and "
            "re-learns signs, so a left-tail observed burden is evidence of "
            "out-of-sample structure rank signal, not overfitting. The "
            "in-sample ceiling (37 protein-wide / 32 regime-local) is "
            "reported in structure_regime_separation_v1.json."
        ),
        "n_controls": len(records),
        "n_perm": N_PERM,
        "seed": SEED,
        "accessions": accessions,
        "site_regimes": regimes,
        "per_protein_loo_ranks_global_regime": {
            acc: variants["global_regime"]["per_protein"][index]
            for index, acc in enumerate(accessions)
        },
        "variants": variants,
        "permutation_nulls": nulls,
        "stats": stats,
        "single_feature": single_feature,
        "summary": {
            "headline": (
                "global_regime LOO burden (the unbiased estimate of the "
                "regime-local composite)"
            ),
            "observed_vs_random": {
                "observed": stats["global_regime"]["observed_burden"],
                "random": stats["global_regime"]["random_baseline"],
            },
            "permutation_p_global_regime": stats["global_regime"][
                "permutation_p"
            ],
            "permutation_p_global_protein": stats["global_protein"][
                "permutation_p"
            ],
            "permutation_p_regime_grouped": stats["regime_grouped_regime"][
                "permutation_p"
            ],
            "regime_routing_adds_value": (
                stats["regime_grouped_regime"]["observed_burden"]
                <= stats["global_regime"]["observed_burden"]
            ),
            "best_single_feature": min(
                single_feature,
                key=lambda k: single_feature[k]["permutation_p"],
            ),
            "single_feature_significant": [
                key
                for key, s in single_feature.items()
                if s["permutation_p"] <= 0.05
            ],
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )

    # --- console digest ------------------------------------------------------
    print(f"LOO + permutation (n=12, B={N_PERM})")
    print(f"{'variant':<24}{'obs':>5}{'rand':>6}{'top1':>6}{'nullMed':>9}{'p':>8}")
    for name, s in stats.items():
        print(
            f"{name:<24}{s['observed_burden']:>5}{s['random_baseline']:>6}"
            f"{s['top1']:>6}{s['null']['median']:>9}{s['permutation_p']:>8}"
        )
    print("\nper-protein LOO ranks (global signs, regime scope):")
    for index, acc in enumerate(accessions):
        rank = per_protein_ranks[index]["rank"]
        n_in = per_protein_ranks[index]["n_in_scope"]
        gene = next(c.gene for c in controls if c.uniprot_accession == acc)
        print(f"  {gene:<12} rank {rank}/{n_in}  ({regimes[index]})")
    print("\nsingle-feature permutation nulls (no learned signs):")
    print(f"{'feature:scope':<32}{'obs':>5}{'rand':>7}{'nullMed':>9}{'p':>8}")
    for key, s in sorted(single_feature.items()):
        marker = "  ***" if s["permutation_p"] <= 0.05 else ""
        print(
            f"{key:<32}{s['observed_burden']:>5}{s['random_baseline']:>7}"
            f"{s['null_median']:>9}{s['permutation_p']:>8}{marker}"
        )
    print(f"\nwrote {OUTPUT}")


if __name__ == "__main__":
    main()
