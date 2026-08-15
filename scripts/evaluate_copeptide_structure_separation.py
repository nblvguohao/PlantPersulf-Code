#!/usr/bin/env python
"""Co-peptide structure separation: POS vs NEG Cys within the same peptide.

The frozen sequence model cannot separate the co-peptide negatives (9
peptides, 8/9 unseparated — the +/-10 window sees nearly identical inputs).
This diagnostic asks the structure side of the same question, with the
strongest matched-detection design the project has: within one peptide, one
spectrum, one digestion, one enrichment, do the MS-confirmed modified Cys
separate from their unmodified siblings on any structure feature?

Design:
- **units** = the 9 co-peptide groups from ``copeptide_negatives_v1.tsv``
  (2 tomato gold-standard: BRG3 SSCMICLPCR, RNF144B FYCPYKDCSAMLVNDSDEIVR;
  7 Arabidopsis from PXD024061, peptide-localised with site-determining ion
  coverage).
- **per-peptide contrast** = for each feature, is the modified set CLEANLY
  above (``pos_higher``), CLEANLY below (``pos_lower``), or overlapping
  (``mixed``) the unmodified set in that same peptide?
- **aggregate** = across the 9 peptides, how many separate in each
  direction, per feature.
- **permutation null** = B=999 within-peptide shuffles of the POS/NEG labels
  (composition preserved: same k_pos / k_neg per group, exchange only among
  that peptide's own Cys), giving the null count for each direction.

Direction is deliberately left free. The known-control analysis (n=12, LOO)
validated burial (``contact_number_10a`` high) for functional sites ACROSS
proteins — but the BRG3 RING co-peptide case shows the modified Cys there is
the LESS buried of the cluster. The within-peptide contrast is the direct
probe of which direction holds when detection is matched; either clean
separation beyond permutation chance is a positive for the structure
hypothesis.

Claim class ``diagnostic_only``; touches no frozen artifact.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from plantpersulf.evaluation.copeptide_structure import (  # noqa: E402
    FEATURE_NAMES,
    aggregate_contrast,
    contrast_permutation_null,
    peptide_contrast,
)
from plantpersulf.evaluation.structure_features import (  # noqa: E402
    cys_structure_features,
    parse_pdb,
)
from plantpersulf.features.sequence import _load_proteome  # noqa: E402

ALPHAFOLD_DIR = _REPO_ROOT / "data" / "raw" / "alphafold"
TOMOTO_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
)
ARABIDOPSIS_PROTEOME = (
    _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
)
CO_PEPTIDE_REGISTRY = _REPO_ROOT / "data" / "registry" / "copeptide_negatives_v1.tsv"
OUTPUT = (
    _REPO_ROOT
    / "results"
    / "diagnostics"
    / "copeptide_structure_separation_v1.json"
)
N_PERM = 999
SEED = 20260815

PROTEOMES = {
    "tomato": TOMOTO_PROTEOME,
    "arabidopsis": ARABIDOPSIS_PROTEOME,
}


def _load_copeptide_groups() -> list[dict]:
    """Read the registry, grouping rows into peptides keyed by accession +
    peptide sequence. Each group: accession, species, peptide, positive and
    negative positions (the peptide's own Cys)."""
    rows: list[dict] = []
    with CO_PEPTIDE_REGISTRY.open(encoding="utf-8") as handle:
        header = handle.readline().strip().split("\t")
        for line in handle:
            if not line.strip():
                continue
            values = line.strip().split("\t")
            rows.append(dict(zip(header, values, strict=False)))

    groups_by_key: dict[tuple[str, str], dict] = {}
    for row in rows:
        key = (row["protein_accession"], row["peptide_sequence"])
        group = groups_by_key.setdefault(
            key,
            {
                "accession": row["protein_accession"],
                "species": row["species"],
                "peptide": row["peptide_sequence"],
                "positive_positions": [],
                "negative_positions": [],
            },
        )
        target = (
            group["positive_positions"]
            if row["state"] == "positive"
            else group["negative_positions"]
        )
        target.append(int(row["cys_position"]))
    groups = list(groups_by_key.values())
    for group in groups:
        group["positive_positions"].sort()
        group["negative_positions"].sort()
    return groups


def _load_feature_rows(groups: list[dict]) -> dict[str, dict[int, dict[str, float]]]:
    """Compute structure features for every Cys of each group's protein."""
    proteomes: dict[str, dict[str, str]] = {}
    for species, path in PROTEOMES.items():
        proteomes[species] = _load_proteome(path)

    result: dict[str, dict[int, dict[str, float]]] = {}
    seen: set[str] = set()
    for group in groups:
        accession = group["accession"]
        if accession in seen:
            continue
        seen.add(accession)
        sequence = proteomes[group["species"]][accession]
        pdb_path = ALPHAFOLD_DIR / f"AF-{accession}-F1-model.pdb"
        if not pdb_path.exists():
            raise RuntimeError(f"no AFDB structure for {accession}: {pdb_path}")
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
        result[accession] = cys_structure_features(residues, positions)
    return result


def _pvalue(null: list[int], observed: int) -> float:
    """Right-tail p for a count direction (observed is one draw)."""
    return (sum(1 for value in null if value >= observed) + 1) / (len(null) + 1)


def main() -> None:
    groups = _load_copeptide_groups()
    feature_rows = _load_feature_rows(groups)

    for group in groups:
        group["features"] = feature_rows[group["accession"]]

    aggregate = aggregate_contrast(groups)
    n_groups = len(groups)

    # --- permutation nulls (B=999 per feature, composition-preserving) ------
    per_feature: dict[str, dict] = {}
    for feature in FEATURE_NAMES:
        observed = aggregate[feature]
        null = contrast_permutation_null(
            groups, feature, N_PERM, np.random.RandomState(SEED)
        )
        per_feature[feature] = {
            "observed": observed,
            "null_pos_higher": null["pos_higher"],
            "null_pos_lower": null["pos_lower"],
            "p_pos_higher": round(
                _pvalue(null["pos_higher"], observed["pos_higher"]), 4
            ),
            "p_pos_lower": round(
                _pvalue(null["pos_lower"], observed["pos_lower"]), 4
            ),
            "null_pos_higher_median": int(np.median(null["pos_higher"])),
            "null_pos_lower_median": int(np.median(null["pos_lower"])),
        }

    # --- per-peptide detail for the report ----------------------------------
    per_peptide: list[dict] = []
    for group in groups:
        contrast = peptide_contrast(
            group["features"],
            group["positive_positions"],
            group["negative_positions"],
        )
        contact = {
            p: round(float(group["features"][p]["contact_number_10a"]), 3)
            for p in group["positive_positions"] + group["negative_positions"]
        }
        per_peptide.append(
            {
                "accession": group["accession"],
                "species": group["species"],
                "peptide": group["peptide"],
                "positive_positions": group["positive_positions"],
                "negative_positions": group["negative_positions"],
                "contrast": contrast,
                "contact_number_10a": contact,
            }
        )

    significant = {
        f: pf
        for f, pf in per_feature.items()
        if min(pf["p_pos_higher"], pf["p_pos_lower"]) <= 0.05
    }

    document = {
        "track": "copeptide_structure_separation_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Within-peptide structural contrast on the explicit co-peptide "
            "negative pool: for each of the 9 peptides (2 tomato "
            "gold-standard, 7 Arabidopsis PXD024061), is the modified set "
            "cleanly above (pos_higher), cleanly below (pos_lower), or "
            "overlapping (mixed) the unmodified set on each structure "
            "feature? Direction is left free — the known-control burial "
            "signal (contact_number) was established ACROSS proteins and may "
            "run the other way WITHIN a peptide (BRG3 RING: modified C206 is "
            "the cluster's exposed Cys). Permutation null preserves each "
            "peptide's composition (k_pos/k_neg, exchange only among that "
            "peptide's own Cys)."
        ),
        "n_peptides": n_groups,
        "n_perm": N_PERM,
        "seed": SEED,
        "per_peptide": per_peptide,
        "aggregate": aggregate,
        "per_feature": per_feature,
        "summary": {
            "significant_features": sorted(significant),
            "best_direction": {
                f: (
                    "pos_higher"
                    if per_feature[f]["p_pos_higher"]
                    <= per_feature[f]["p_pos_lower"]
                    else "pos_lower"
                )
                for f in significant
            },
            "contact_number_10a": per_feature["contact_number_10a"],
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )

    # --- console digest ------------------------------------------------------
    print(f"Co-peptide structure separation (n={n_groups} peptides, B={N_PERM})")
    print(
        f"{'feature':<22}{'posHi':>6}{'posLo':>6}{'mix':>5}"
        f"{'nullMedHi':>10}{'pHi':>7}{'pLo':>7}"
    )
    for feature in FEATURE_NAMES:
        pf = per_feature[feature]
        marker = "  **" if min(pf["p_pos_higher"], pf["p_pos_lower"]) <= 0.05 else ""
        print(
            f"{feature:<22}{pf['observed']['pos_higher']:>6}"
            f"{pf['observed']['pos_lower']:>6}{pf['observed']['mixed']:>5}"
            f"{pf['null_pos_higher_median']:>10}"
            f"{pf['p_pos_higher']:>7}{pf['p_pos_lower']:>7}{marker}"
        )
    print("\nper-peptide contrast + contact (all 7 features):")
    for pp in per_peptide:
        direction = pp["contrast"]["contact_number_10a"]
        print(
            f"  {pp['accession']:<12} POS={pp['positive_positions']} "
            f"NEG={pp['negative_positions']} contact={pp['contact_number_10a']} "
            f"=> contact {direction}"
        )
    print(f"\nwrote {OUTPUT}")


if __name__ == "__main__":
    main()
