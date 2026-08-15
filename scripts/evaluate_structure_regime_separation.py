#!/usr/bin/env python
"""Systematic in-regime structure separation across all 12 mapped controls.

P2 (``evaluate_p2_structure_separation.py``) showed on two proteins that the
structure feature set separates co-peptide Cys *within a regime*: BRG3's RING
domain signed composite ranks true C206 first and gold-standard negative C209
last, while the same signs applied protein-wide collapse (N-terminal
disordered Cys share the "exposed" shape). This diagnostic extends that to
EVERY registered mapped control — 5 tomato + 7 Arabidopsis — using the
registered AlphaFold DB models (the 7 Arabidopsis structures were downloaded
and registered for this analysis; all PDB/proteome residue counts agree).

Evidence levels, in increasing degrees of learner assumption:

1. ``per_feature_aggregate`` — for each of the 7 structure features, the
   within-protein z of the true (registered representative) site, averaged
   across proteins (unbiased: computed per protein, then aggregated), plus
   Top-1 / Hit@2 / total first-hit burden vs the sum of random baselines.
2. ``regime_split`` / ``cluster_split`` — the same aggregation restricted to
   sites whose own structure is folded/linker/disordered (AF pLDDT prior) or
   metal-cluster-like (>=2 S-gamma within 8 A). Tests whether the exposure
   signal is regime-specific (proposition 3 "route, don't average").
3. ``regime_local_ranking`` — per feature, rank the true site among ONLY the
   Cys in its own pLDDT regime, vs ranking among all protein Cys. A routing
   win (regime-local burden < protein-wide burden) needs no learned signs.
4. ``composite_ceiling`` — signed z-composite with signs taken from the
   aggregate directions (n=12). An UPPER BOUND on what perfect direction
   learning at this sample could deliver — not unbiased evidence.

Group-variant sites (BRG3 206/212, DES1 44/205, RBOHD 825/890, SnRK2.6
131/137) are reported per protein; the aggregate uses the registered k=1
representative (release single-position-per-lineage convention).

Diagnostic-only: no fitting of frozen artifacts, no SAP changes, no new
pre-registration. Claim class ``diagnostic_only``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

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
    aggregate_site_z,
    composite_scores,
    ranking_stats,
    regime_bucket,
    site_plddt_regime,
    subset_positions,
)
from plantpersulf.evaluation.within_protein_ranking import (  # noqa: E402
    first_hit_rank,
    rank_sites_desc,
)
from plantpersulf.features.sequence import _load_proteome  # noqa: E402

ALPHAFOLD_DIR = _REPO_ROOT / "data" / "raw" / "alphafold"
TOMATO_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
)
ARABIDOPSIS_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
)
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "structure_regime_separation_v1.json"

SPECIES_PROTEOMES = {
    "Solanum lycopersicum": ("tomato", TOMATO_PROTEOME),
    "Arabidopsis thaliana": ("arabidopsis", ARABIDOPSIS_PROTEOME),
}

GROUP_SITE_VARIANTS: dict[str, tuple[int, ...]] = {
    "BRG3_H2S_UBIQUITINATION": (206, 212),
    "DES1_H2S_SELF_RBOHD_ABA": (44, 205),
    "RBOHD_H2S_ROS_ABA": (825, 890),
    "SNRK26_H2S_PHOSPHORYLATION": (131, 137),
}

# Methodology-design domain annotation carried over from P2 (not data-derived).
BRG3_RING_DOMAIN = (197, 231)


def _feature_values(
    feature_rows: dict[int, dict[str, float]], feature: str
) -> dict[int, float]:
    return {position: float(row[feature]) for position, row in feature_rows.items()}


def _aggregate_burden(
    records: list[tuple[dict[int, dict[str, float]], int]],
    feature: str,
    scope_positions: dict[int, list[int]] | None = None,
) -> dict[str, object]:
    """Burden aggregation over k=1 representative sites.

    ``scope_positions`` maps protein index -> Cys positions to rank within
    (regime-local). Falls back to all Cys when absent.
    """
    total = 0
    total_random = 0.0
    top1 = 0
    for index, (rows, true_position) in enumerate(records):
        scope = (
            scope_positions.get(index) if scope_positions is not None else list(rows)
        )
        values = _feature_values(rows, feature)
        ranking = sorted(scope, key=lambda position: values[position], reverse=True)
        rank = ranking.index(true_position) + 1
        total += rank
        total_random += (len(scope) + 1) / 2.0
        top1 += rank == 1
    return {
        "n": len(records),
        "top1": top1,
        "total_first_hit_burden": total,
        "total_random_burden": round(total_random, 3),
    }


def main() -> None:
    controls = [c for c in REGISTERED_CONTROLS if c.status == "mapped"]
    proteomes: dict[str, dict[str, str]] = {}
    for _species, (name, path) in SPECIES_PROTEOMES.items():
        proteomes[name] = _load_proteome(path)

    per_protein: dict[str, dict[str, object]] = {}
    records_k1: list[
        tuple[dict[int, dict[str, float]], int, str]
    ] = []  # (feature_rows, true_position, accession)
    regime_records: dict[str, list] = {"folded": [], "linker": [], "disordered": []}
    cluster_records: dict[str, list] = {"cluster_like": [], "isolated": []}

    for control in controls:
        accession = control.uniprot_accession
        species, _ = SPECIES_PROTEOMES[control.control_species]
        sequence = proteomes[species][accession]
        pdb_path = ALPHAFOLD_DIR / f"AF-{accession}-F1-model.pdb"
        if not pdb_path.exists():
            per_protein[accession] = {
                "gene": control.gene,
                "skipped": f"missing {pdb_path.name}",
            }
            continue
        residues = parse_pdb(pdb_path.read_text(encoding="utf-8"))
        if len(residues) != len(sequence):
            raise RuntimeError(
                f"{accession}: PDB {len(residues)} residues != proteome {len(sequence)}"
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

        true_site = control.cys_position
        site_regime = site_plddt_regime(feature_rows, true_site)
        regime_histogram = {
            bucket: len(
                subset_positions(
                    feature_rows,
                    lambda p, row, b=bucket: regime_bucket(float(row["plddt"])) == b,
                )
            )
            for bucket in ("folded", "linker", "disordered")
        }
        cluster_cys = subset_positions(
            feature_rows, lambda p, row: float(row["cys_count_8a"]) >= 2.0
        )

        # per-feature rank of the k=1 representative site
        feature_ranks: dict[str, dict[str, object]] = {}
        for feature in STRUCTURE_FEATURE_NAMES:
            values = _feature_values(feature_rows, feature)
            ranking = rank_sites_desc(values)
            rank = ranking.index(true_site) + 1
            feature_ranks[feature] = {
                "rank": rank,
                "n_cys": len(positions),
                "value": values[true_site],
                "site_regime": site_regime,
            }

        group_variants = GROUP_SITE_VARIANTS.get(control.mechanism_lineage_id)

        per_protein[accession] = {
            "gene": control.gene,
            "species": species,
            "mechanism_lineage_id": control.mechanism_lineage_id,
            "n_cys": len(positions),
            "true_site": true_site,
            "group_variants": list(group_variants) if group_variants else None,
            "site_regime": site_regime,
            "site_plddt": round(float(feature_rows[true_site]["plddt"]), 2),
            "regime_histogram": regime_histogram,
            "cluster_cys": cluster_cys,
            "feature_ranks": feature_ranks,
            "features": feature_rows,
        }
        records_k1.append((feature_rows, true_site, accession))
        regime_records[site_regime].append((feature_rows, true_site))
        cluster_records[
            "cluster_like" if true_site in cluster_cys else "isolated"
        ].append((feature_rows, true_site))

    # --- level 1: per-feature aggregate (k=1) ------------------------------
    per_feature_aggregate: dict[str, dict[str, object]] = {}
    for feature in STRUCTURE_FEATURE_NAMES:
        zagg = aggregate_site_z(
            [(rows, true_position) for rows, true_position, _ in records_k1],
            feature,
        )
        rstats = ranking_stats(
            [(rows, true_position) for rows, true_position, _ in records_k1],
            feature,
        )
        per_feature_aggregate[feature] = {
            "mean_z": round(float(zagg["mean_z"]), 4),
            "above_median": zagg["above_median"],
            "n": zagg["n"],
            "top1": rstats["top1"],
            "hit_at_2": rstats["hit_at_2"],
            "total_first_hit_burden": rstats["total_first_hit_burden"],
            "total_random_burden": rstats["total_random_burden"],
        }

    # --- level 2: regime / cluster splits ----------------------------------
    regime_split: dict[str, dict[str, object]] = {}
    for bucket, bucket_records in regime_records.items():
        if not bucket_records:
            regime_split[bucket] = {"n": 0}
            continue
        regime_split[bucket] = {
            "n": len(bucket_records),
            "features": {
                feature: {
                    "mean_z": round(
                        float(aggregate_site_z(bucket_records, feature)["mean_z"]),
                        4,
                    ),
                    "hit_at_2": ranking_stats(bucket_records, feature)["hit_at_2"],
                }
                for feature in STRUCTURE_FEATURE_NAMES
            },
        }
    cluster_split: dict[str, dict[str, object]] = {}
    for bucket, bucket_records in cluster_records.items():
        if not bucket_records:
            cluster_split[bucket] = {"n": 0}
            continue
        cluster_split[bucket] = {
            "n": len(bucket_records),
            "features": {
                feature: {
                    "mean_z": round(
                        float(aggregate_site_z(bucket_records, feature)["mean_z"]),
                        4,
                    ),
                    "hit_at_2": ranking_stats(bucket_records, feature)["hit_at_2"],
                }
                for feature in STRUCTURE_FEATURE_NAMES
            },
        }

    # --- level 3: regime-local ranking vs protein-wide (no signs) ----------
    regime_local_ranking: dict[str, dict[str, object]] = {}
    for feature in STRUCTURE_FEATURE_NAMES:
        scopes: dict[int, list[int]] = {}
        for index, (rows, true_position, _) in enumerate(records_k1):
            bucket = site_plddt_regime(rows, true_position)
            scopes[index] = subset_positions(
                rows,
                lambda p, row, b=bucket: regime_bucket(float(row["plddt"])) == b,
            )
        protein_wide = _aggregate_burden(
            [(rows, t) for rows, t, _ in records_k1], feature
        )
        regime_local = _aggregate_burden(
            [(rows, t) for rows, t, _ in records_k1], feature, scopes
        )
        regime_local_ranking[feature] = {
            "protein_wide": protein_wide,
            "regime_local": regime_local,
            "routing_help": (
                regime_local["total_first_hit_burden"]
                < protein_wide["total_first_hit_burden"]
            ),
        }

    # --- level 4: signed composite ceiling ---------------------------------
    signs: dict[str, float] = {
        feature: (
            1.0
            if per_feature_aggregate[feature]["mean_z"] > 0
            else -1.0
            if per_feature_aggregate[feature]["mean_z"] < 0
            else 0.0
        )
        for feature in STRUCTURE_FEATURE_NAMES
    }
    active_signs = {f: s for f, s in signs.items() if s != 0.0}

    def _composite_burden(
        use_regime_scope: bool,
    ) -> tuple[dict[str, object], dict[str, object]]:
        total = 0
        total_random = 0.0
        top1 = 0
        per_protein_ranks: dict[str, dict[str, object]] = {}
        for rows, true_position, accession in records_k1:
            scope = None
            if use_regime_scope:
                bucket = site_plddt_regime(rows, true_position)
                scope = subset_positions(
                    rows,
                    lambda p, row, b=bucket: regime_bucket(float(row["plddt"])) == b,
                )
            scores = composite_scores(rows, active_signs, positions=scope)
            ranking = rank_sites_desc(scores)
            rank = ranking.index(true_position) + 1
            total += rank
            total_random += (len(scores) + 1) / 2.0
            top1 += rank == 1
            per_protein_ranks[accession] = {
                "rank": rank,
                "n_in_scope": len(scores),
                "ranking": ranking,
            }
        return (
            {
                "top1": top1,
                "total_first_hit_burden": total,
                "total_random_burden": round(total_random, 3),
            },
            per_protein_ranks,
        )

    pw_agg, pw_ranks = _composite_burden(use_regime_scope=False)
    rl_agg, rl_ranks = _composite_burden(use_regime_scope=True)
    composite_ceiling = {
        "signs": signs,
        "n_active_features": len(active_signs),
        "note": (
            "Signs from the n=12 aggregate within-protein directions — an "
            "UPPER BOUND on what perfect direction learning at this sample "
            "could deliver, not unbiased evidence. Regime-local scope ranks "
            "only within the true site's own pLDDT bucket."
        ),
        "protein_wide": pw_agg,
        "regime_local": rl_agg,
        "per_protein_protein_wide": pw_ranks,
        "per_protein_regime_local": rl_ranks,
    }

    # --- BRG3 RING domain continuity (P2 level-5 method, exposure signs) ----
    # Reproduce P2's RING-local signed composite faithfully: z-normalisation
    # and signs BOTH within the RING domain (C197-231), using the exposure
    # directions the P2 separation table found (RSA high, contact low,
    # Coulomb low = free/exposed Cys). This is the continuity check against
    # p2_structure_separation_v1.json, NOT a re-estimate from the aggregate.
    brg3 = per_protein["A0A3Q7EW23"]
    brg3_rows = brg3["features"]
    ring_positions = subset_positions(
        brg3_rows,
        lambda p, row: BRG3_RING_DOMAIN[0] <= p <= BRG3_RING_DOMAIN[1],
    )
    p2_ring_signs = {
        "rsa_relative": 1.0,
        "contact_number_10a": -1.0,
        "coulomb_potential_sg": -1.0,
    }
    ring_scores = composite_scores(brg3_rows, p2_ring_signs, positions=ring_positions)
    ring_ranking = rank_sites_desc(ring_scores)
    brg3_ring = {
        "domain": list(BRG3_RING_DOMAIN),
        "positions": ring_positions,
        "method": (
            "P2 level-5 reproduction: RING-local z and signs from the P2 "
            "separation table (rsa +1, contact -1, coulomb -1). Continuity "
            "check against p2_structure_separation_v1.json."
        ),
        "true_sites": [206, 212],
        "gold_standard_negative": 209,
        "ranking": ring_ranking,
        "true_ranks": {
            position: ring_ranking.index(position) + 1
            for position in (206, 212)
            if position in ring_positions
        },
        "negative_rank": ring_ranking.index(209) + 1,
        "first_hit_burden": first_hit_rank(ring_ranking, {206, 212}),
        "random_baseline_burden": (len(ring_positions) + 1) / 3.0,
    }
    # Contrast: the same RING positions under the GLOBAL aggregate signs —
    # makes the regime-sign tension explicit (C209 rises when burial is the
    # preferred direction).
    ring_global_scores = composite_scores(
        brg3_rows, active_signs, positions=ring_positions
    )
    ring_global_ranking = rank_sites_desc(ring_global_scores)
    brg3_ring["global_signs_ranking"] = ring_global_ranking
    brg3_ring["global_signs_negative_rank"] = (
        ring_global_ranking.index(209) + 1 if 209 in ring_global_ranking else None
    )

    summary = {
        "n_controls": len(controls),
        "n_with_structure": sum(1 for r in per_protein.values() if "features" in r),
        "n_idr_sites": regime_split.get("disordered", {}).get("n", 0)
        + regime_split.get("linker", {}).get("n", 0),
        "n_folded_sites": regime_split.get("folded", {}).get("n", 0),
        "n_cluster_like_sites": cluster_split.get("cluster_like", {}).get("n", 0),
        "best_single_feature": max(
            per_feature_aggregate,
            key=lambda f: per_feature_aggregate[f]["hit_at_2"],
        ),
        "best_feature_burden_ratio": min(
            per_feature_aggregate[f]["total_first_hit_burden"]
            / per_feature_aggregate[f]["total_random_burden"]
            for f in per_feature_aggregate
        ),
        "routing_help_features": [
            f for f, r in regime_local_ranking.items() if r["routing_help"]
        ],
        "composite_ceiling": {
            "protein_wide_burden": composite_ceiling["protein_wide"][
                "total_first_hit_burden"
            ],
            "regime_local_burden": composite_ceiling["regime_local"][
                "total_first_hit_burden"
            ],
        },
        "brg3_ring_true_rank": brg3_ring["true_ranks"].get(206),
        "brg3_ring_negative_rank": brg3_ring["negative_rank"],
        "brg3_ring_global_signs_negative_rank": brg3_ring["global_signs_negative_rank"],
    }

    document = {
        "track": "structure_regime_separation_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Systematic in-regime structure separation for all 12 registered "
            "mapped controls (5 tomato + 7 Arabidopsis), registered AFDB v6 "
            "models, PDB/proteome residue counts all verified. Features: "
            "pLDDT, RSA, S-gamma cluster geometry, positive-residue "
            "microenvironment, Coulomb potential, contact number. No fitting "
            "of frozen artifacts; metal coordination is a geometry proxy "
            "(AFDB models contain no metal ions); supervised learning is "
            "still not estimable at n=12 labelled sites — ranking/burden "
            "aggregates are the evidence, the signed composite is a labeled "
            "ceiling. HEADLINE: within-protein direction of the true site is "
            "REGIME-DEPENDENT — folded-regime sites are structurally buried/"
            "packed (contact_number mean_z +0.72, hit@2 6/9), while "
            "disordered-regime sites are exposed (RSA / nearest-Sg-distance "
            "hit@2 3/3). A single global-sign composite is folded-biased and "
            "cannot serve both regimes; regime-local ranking (no signs) "
            "improves every feature — the proposition-3 routing test at "
            "n=12. The BRG3 RING exposure-sign result from P2 is reproduced "
            "for continuity, and its fragility under global burial signs is "
            "reported explicitly (C209 rank 7 under global signs vs rank 9 "
            "under RING exposure signs)."
        ),
        "coverage": {
            "n_controls": len(controls),
            "n_with_structure": summary["n_with_structure"],
            "structures": {
                accession: "AFDB v6 (registered alphafold_structures.tsv)"
                for accession in per_protein
            },
        },
        "per_protein": per_protein,
        "per_feature_aggregate": per_feature_aggregate,
        "regime_split": regime_split,
        "cluster_split": cluster_split,
        "regime_local_ranking": regime_local_ranking,
        "composite_ceiling": composite_ceiling,
        "brg3_ring_domain": brg3_ring,
        "supervised_note": (
            "Supervised gradient boosting is still NOT estimable: 12 labelled "
            "sites across ~100 Cys. The n=12 aggregate direction (per_feature "
            "aggregate) and the regime routing comparison are the evidence; "
            "the signed composite is a ceiling, not an estimate."
        ),
        "summary": summary,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )

    # --- console digest ------------------------------------------------------
    print("coverage:", summary["n_with_structure"], "/", summary["n_controls"])
    print("\nper-feature aggregate (k=1, n=12):")
    print(
        f"  {'feature':<28}{'mean_z':>8}{'aboveMed':>9}{'top1':>6}"
        f"{'hit@2':>7}{'burden':>8}{'randB':>8}"
    )
    for feature, agg in per_feature_aggregate.items():
        print(
            f"  {feature:<28}{agg['mean_z']:>8}{agg['above_median']:>9}"
            f"{agg['top1']:>6}{agg['hit_at_2']:>7}"
            f"{agg['total_first_hit_burden']:>8}{agg['total_random_burden']:>8}"
        )
    print("\nregime split (site pLDDT):")
    for bucket, rec in regime_split.items():
        if not rec.get("n"):
            print(f"  {bucket}: n=0")
            continue
        hits = {f: v["hit_at_2"] for f, v in rec["features"].items()}
        print(f"  {bucket}: n={rec['n']} hit@2 {hits}")
    print("\nregime-local vs protein-wide single-feature ranking (no signs):")
    for feature, rec in regime_local_ranking.items():
        pw = rec["protein_wide"]["total_first_hit_burden"]
        rl = rec["regime_local"]["total_first_hit_burden"]
        flag = "  <-- routing helps" if rec["routing_help"] else ""
        print(f"  {feature:<28} pw={pw:>3} rl={rl:>3}{flag}")
    print("\ncomposite ceiling (signed, labels from aggregate):")
    print(
        f"  protein-wide: {composite_ceiling['protein_wide']} "
        f"regime-local: {composite_ceiling['regime_local']}"
    )
    print(f"  signs: { {f: signs[f] for f in signs if signs[f] != 0.0} }")
    print(
        f"\nBRG3 RING domain (C197-231): true C206/C212 ranks "
        f"{brg3_ring['true_ranks']}, negative C209 rank "
        f"{brg3_ring['negative_rank']}, burden "
        f"{brg3_ring['first_hit_burden']} (random "
        f"{round(brg3_ring['random_baseline_burden'], 2)})"
    )
    print(
        f"  (global burial signs within RING: C209 rank "
        f"{brg3_ring['global_signs_negative_rank']} — regime-sign tension "
        f"explicit)"
    )
    print("\nper-protein best structure-feature rank of the true site:")
    for _accession, record in per_protein.items():
        if "features" not in record:
            continue
        best = min(
            record["feature_ranks"],
            key=lambda f: record["feature_ranks"][f]["rank"],
        )
        best_rank = record["feature_ranks"][best]["rank"]
        print(
            f"  {record['gene']:<12} Cys{record['true_site']:<4} "
            f"regime={record['site_regime']:<10} "
            f"best_feature={best} (rank {best_rank}/{record['n_cys']})"
        )
    print(f"\nwrote {OUTPUT}")


if __name__ == "__main__":
    main()
