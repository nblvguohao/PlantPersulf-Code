#!/usr/bin/env python
"""Diagnostic A: does regime-routed MoE beat the single contact feature?

The v2 regime-routing design (``2026-08-15-regime-routing-moe-design.md``)
proposes a MoE structure branch whose gate is driven by pLDDT and whose
experts use fixed structural-prior directions. Before any of that is trained
(which would require a new release + preregistration), this diagnostic asks
the minimal question on the 12 registered known controls: **does routing a
site to its regime expert and scoring with that expert's prior beat simply
ranking by contact_number_10a?**

Two routing scopes are reported (both with B=999 exchangeability nulls, seed
20260815):

- ``protein``: routed scores across regimes compete on one scale — the
  stricter bar, since the folded and disordered experts use different
  features and different z-pools.
- ``regime``: each site is ranked only within its own pLDDT bucket — the MoE
  claim itself (regime-local routing).

Expert directions are fixed structural priors, so there is no LOO sign
learning and no ceiling: the permutation null is the honest reference.

Reference result (``structure_regime_loo_v1.json``): contact_number_10a
single feature is the ONLY LOO-significant feature — burden 36 protein-wide
(p=0.018) / 30 regime-local (p=0.032) vs random 58.5. The routed composite
must at least match those numbers to justify the architecture.

Claim class ``diagnostic_only``; touches no frozen artifact, trains nothing.
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
from plantpersulf.evaluation.regime_moe import (  # noqa: E402
    EXPERT_DIRECTIONS,
    routed_burden,
    routed_permutation_null,
)
from plantpersulf.evaluation.structure_features import (  # noqa: E402
    cys_structure_features,
    parse_pdb,
)
from plantpersulf.evaluation.structure_regime import (  # noqa: E402
    _in_regime_bucket,
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
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "regime_moe_v1.json"
N_PERM = 999
SEED = 20260815

SPECIES_PROTEOMES = {
    "Solanum lycopersicum": ("tomato", TOMATO_PROTEOME),
    "Arabidopsis thaliana": ("arabidopsis", ARABIDOPSIS_PROTEOME),
}

CONTACT = "contact_number_10a"


def _contact_burden(
    records: list[tuple[dict[int, dict[str, float]], int]],
    scope_mode: str,
) -> int:
    """Total first-hit burden of ranking by contact_number_10a alone (the
    reference single feature), optionally regime-local."""
    total = 0
    for rows, true_position in records:
        values = {position: float(row[CONTACT]) for position, row in rows.items()}
        if scope_mode == "protein":
            ranking = sorted(values, key=lambda p: values[p], reverse=True)
        else:
            bucket = regime_bucket(float(rows[true_position]["plddt"]))
            scope = subset_positions(rows, _in_regime_bucket(bucket))
            ranking = sorted(scope, key=lambda p: values[p], reverse=True)
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
    regime_counts = {bucket: regimes.count(bucket) for bucket in sorted(set(regimes))}

    # --- observed burdens ---------------------------------------------------
    routed: dict[str, object] = {}
    contact: dict[str, object] = {}
    for scope_mode in ("protein", "regime"):
        routed[scope_mode] = routed_burden(
            base, EXPERT_DIRECTIONS, scope_mode=scope_mode
        )
        observed_contact = _contact_burden(base, scope_mode)
        rng = np.random.RandomState(
            SEED + zlib.crc32(f"contact:{scope_mode}".encode())
        )
        contact_null = []
        for _ in range(N_PERM):
            perm_records = [
                (rows, int(rng.choice(list(rows)))) for rows, _ in base
            ]
            contact_null.append(_contact_burden(perm_records, scope_mode))
        contact_null.sort()
        contact[scope_mode] = {
            "observed_burden": observed_contact,
            "random_baseline": round(
                sum((len(rows) + 1) / 2.0 for rows, _ in base), 3
            ),
            "permutation_p": round(
                (
                    sum(1 for value in contact_null if value <= observed_contact) + 1
                )
                / (N_PERM + 1),
                4,
            ),
            "null_median": int(np.median(contact_null)),
        }

    # --- routed permutation nulls -------------------------------------------
    routed_nulls: dict[str, list[int]] = {}
    for scope_mode in ("protein", "regime"):
        rng = np.random.RandomState(
            SEED + zlib.crc32(f"routed:{scope_mode}".encode())
        )
        routed_nulls[scope_mode] = routed_permutation_null(
            base, EXPERT_DIRECTIONS, scope_mode, N_PERM, rng
        )

    def _pvalue(null: list[int], observed: int) -> float:
        return (sum(1 for value in null if value <= observed) + 1) / (len(null) + 1)

    routed_stats: dict[str, object] = {}
    for scope_mode in ("protein", "regime"):
        res = routed[scope_mode]
        observed = res["total_first_hit_burden"]
        null = routed_nulls[scope_mode]
        routed_stats[scope_mode] = {
            "observed_burden": observed,
            "random_baseline": res["total_random_burden"],
            "top1": res["top1"],
            "permutation_p": round(_pvalue(null, observed), 4),
            "null_median": int(np.median(null)),
            "null_p5": null[max(0, len(null) // 20)],
            "per_protein": res["per_protein"],
        }

    document = {
        "track": "regime_moe_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Diagnostic A of the v2 regime-routing MoE design: does routing "
            "each site to its pLDDT regime expert (fixed structural-prior "
            "directions) beat the single contact_number_10a feature? Expert "
            "directions are written down as structural hypotheses, never "
            "learned from labels, so there is no LOO sign-learning bias and "
            "no ceiling. scope 'regime' ranks within each site's own pLDDT "
            "bucket (the MoE claim); 'protein' makes routed scores across "
            "regimes compete on one scale (stricter). Reference: contact "
            "single feature burden 36 protein / 30 regime vs random 58.5 "
            "(structure_regime_loo_v1.json)."
        ),
        "n_controls": len(records),
        "n_perm": N_PERM,
        "seed": SEED,
        "accessions": accessions,
        "site_regimes": regimes,
        "regime_counts": regime_counts,
        "expert_directions": EXPERT_DIRECTIONS,
        "routed": routed,
        "routed_stats": routed_stats,
        "contact_single_feature": contact,
        "summary": {
            "headline": (
                "regime-local routed burden vs contact single feature: "
                "does the MoE routing add value over the only significant "
                "structure feature?"
            ),
            "routed_regime": {
                "observed": routed_stats["regime"]["observed_burden"],
                "random": routed_stats["regime"]["random_baseline"],
                "permutation_p": routed_stats["regime"]["permutation_p"],
            },
            "contact_regime": contact["regime"],
            "routed_protein": {
                "observed": routed_stats["protein"]["observed_burden"],
                "random": routed_stats["protein"]["random_baseline"],
                "permutation_p": routed_stats["protein"]["permutation_p"],
            },
            "contact_protein": contact["protein"],
            "routing_beats_contact_regime": (
                routed_stats["regime"]["observed_burden"]
                <= contact["regime"]["observed_burden"]
            ),
            "routing_beats_contact_protein": (
                routed_stats["protein"]["observed_burden"]
                <= contact["protein"]["observed_burden"]
            ),
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )

    # --- console digest -----------------------------------------------------
    print(f"Regime-routed MoE vs contact single feature (n={len(records)}, B={N_PERM})")
    print(f"regime counts: {regime_counts}")
    print(
        f"{'variant':<24}{'observed':>9}{'random':>9}{'nullMed':>9}{'p':>7}"
    )
    for label, _kind, mode in (
        ("contact protein", "contact", "protein"),
        ("contact regime", "contact", "regime"),
        ("routed protein", "routed", "protein"),
        ("routed regime", "routed", "regime"),
    ):
        if label.startswith("contact"):
            res = contact[mode]
            observed = res["observed_burden"]
            p = res["permutation_p"]
            null_med = res["null_median"]
        else:
            res = routed_stats[mode]
            observed = res["observed_burden"]
            p = res["permutation_p"]
            null_med = res["null_median"]
        print(
            f"{label:<24}{observed:>9}{res['random_baseline']:>9}{null_med:>9}{p:>7}"
        )
    print("\nper-protein routed regime-local ranks:")
    for index, acc in enumerate(accessions):
        rank = routed_stats["regime"]["per_protein"][index]
        print(
            f"  {acc:<12} regime={regimes[index]:<10} rank={rank['rank']}"
            f"/{rank['n_in_scope']}"
        )


if __name__ == "__main__":
    main()
