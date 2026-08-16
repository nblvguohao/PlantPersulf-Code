#!/usr/bin/env python
"""BRG3 metal-coordination proxy — the S-gamma-layer "last shot".

BRG3 (A0A3Q7EW23) carries a C3HC4-type RING with 9 Cys within 35 aa; a
canonical RING needs only 7 Cys + 1 His as metal ligands, so at least two of
the cluster's Cys are non-ligands. The co-peptide group SSCMICLPCR labels
Cys206 and Cys212 persulfidated (POS) and Cys209 unmodified (NEG). The
metal-coordination proxy (``metal_coordination_sg_3a``) counts N/O/S atoms of
His/Cys/Asp/Glu residues within 3 A of each Cys S-gamma — a surrogate for Zn
coordination (AFDB models carry no metal ions). This diagnostic asks whether
the proxy is a binary discriminator invisible to every existing feature:
modified (non-ligand, exposed) Cys should read LOW, unmodified (ligand) Cys
HIGH. The same composition-preserving permutation null is applied across the
14-group / 3-species co-peptide pool to see whether any separation is
BRG3-local or general.

claim_class diagnostic_only; touches no frozen artifact.
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
CO_PEPTIDE_REGISTRY = _REPO_ROOT / "data" / "registry" / "copeptide_negatives_v1.tsv"
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "brg3_metal_coordination_v1.json"
METAL_FEATURE = "metal_coordination_sg_3a"
N_PERM = 999
SEED = 20260816
BRG3_ACCESSION = "A0A3Q7EW23"
BRG3_POSITIVES = (206, 212)
BRG3_NEGATIVES = (209,)

PROTEOMES = {
    "tomato": (
        _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
    ),
    "arabidopsis": (
        _REPO_ROOT
        / "data"
        / "raw"
        / "references"
        / "arabidopsis_ref_proteome_v2.fasta"
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


def _load_copeptide_groups() -> list[dict]:
    """Read the registry, grouping rows into peptides keyed by accession +
    peptide sequence (same dedup as the co-peptide separation diagnostic)."""
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
    # Dedup by peptide sequence, keeping the first protein instance (the rice
    # IIPTPNC peptide is carried by 4 paralogous proteins; the peptide is the
    # unit of the diagnostic, so each distinct peptide counts once).
    groups_by_peptide: dict[str, dict] = {}
    for group in groups_by_key.values():
        groups_by_peptide.setdefault(group["peptide"], group)
    groups = list(groups_by_peptide.values())
    for group in groups:
        group["positive_positions"].sort()
        group["negative_positions"].sort()
    return groups


def _load_feature_rows(
    groups: list[dict],
) -> tuple[dict[str, dict[int, dict[str, float]]], list[str]]:
    """Structure features for every Cys of each group's protein (incl. the
    metal-coordination proxy). Missing/mismatched structures are skipped and
    their reasons returned, never silently dropped."""
    proteomes: dict[str, dict[str, str]] = {}
    for species, path in PROTEOMES.items():
        proteomes[species] = _load_proteome(path)

    result: dict[str, dict[int, dict[str, float]]] = {}
    skipped: list[str] = []
    seen: set[str] = set()
    for group in groups:
        accession = group["accession"]
        if accession in seen:
            continue
        seen.add(accession)
        sequence = proteomes[group["species"]][accession]
        pdb_path = ALPHAFOLD_DIR / f"AF-{accession}-F1-model.pdb"
        if not pdb_path.exists():
            skipped.append(f"{accession}: no AFDB structure {pdb_path.name}")
            continue
        residues = parse_pdb(pdb_path.read_text(encoding="utf-8"))
        if len(residues) != len(sequence):
            skipped.append(
                f"{accession}: PDB {len(residues)} != proteome {len(sequence)}"
            )
            continue
        mismatches = [
            position
            for position in range(1, len(sequence) + 1)
            if residues[position - 1].type != sequence[position - 1]
        ]
        if mismatches:
            skipped.append(
                f"{accession}: PDB/proteome mismatch at {mismatches[:10]}"
            )
            continue
        positions = [i + 1 for i, residue in enumerate(sequence) if residue == "C"]
        result[accession] = cys_structure_features(residues, positions)
    return result, skipped


def _pvalue(null: list[int], observed: int) -> float:
    """Right-tail p for a count direction (observed is one draw)."""
    return (sum(1 for value in null if value >= observed) + 1) / (len(null) + 1)


def main() -> None:
    groups = _load_copeptide_groups()
    feature_rows, skipped = _load_feature_rows(groups)

    evaluated = [g for g in groups if g["accession"] in feature_rows]
    for group in evaluated:
        group["features"] = feature_rows[group["accession"]]
    for reason in skipped:
        print(f"  skipped: {reason}")

    # --- aggregate contrast + permutation null for the metal feature --------
    aggregate = aggregate_contrast(evaluated, (METAL_FEATURE,))
    observed = aggregate[METAL_FEATURE]
    null = contrast_permutation_null(
        evaluated, METAL_FEATURE, N_PERM, np.random.RandomState(SEED)
    )
    p_hi = _pvalue(null["pos_higher"], observed["pos_higher"])
    p_lo = _pvalue(null["pos_lower"], observed["pos_lower"])

    # --- BRG3 detail ---------------------------------------------------------
    brg3 = feature_rows[BRG3_ACCESSION]
    trio = {
        p: {k: v for k, v in brg3[p].items()}
        for p in (*BRG3_POSITIVES, *BRG3_NEGATIVES)
    }
    contrast_brg3 = peptide_contrast(
        brg3,
        list(BRG3_POSITIVES),
        list(BRG3_NEGATIVES),
        (METAL_FEATURE,),
    )[METAL_FEATURE]
    metal_all = {
        p: float(brg3[p][METAL_FEATURE]) for p in sorted(brg3)
    }

    document = {
        "track": "brg3_metal_coordination_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Metal-coordination proxy (metal_coordination_sg_3a: N/O/S atoms "
            "of His/Cys/Asp/Glu residues within 3 A of each Cys S-gamma) as a "
            "binary discriminator for the BRG3 C3HC4 RING co-peptide group "
            "(SSCMICLPCR: Cys206/Cys212 persulfidated POS, Cys209 unmodified "
            "NEG). AFDB models carry no metal ions, so Zn coordination is "
            "proxied by ligand-atom density around the sulphur. Modified "
            "(non-ligand, exposed) Cys should read LOW; unmodified (ligand) "
            "Cys HIGH. Same composition-preserving permutation null across the "
            "14-group / 3-species pool checks whether any separation is "
            "BRG3-local or general."
        ),
        "metal_feature": METAL_FEATURE,
        "radius_angstrom": 3.0,
        "n_peptides_pool": len(evaluated),
        "n_dropped_no_structure": len(skipped),
        "aggregate": observed,
        "p_pos_higher": round(p_hi, 4),
        "p_pos_lower": round(p_lo, 4),
        "null_pos_higher_median": int(np.median(null["pos_higher"])),
        "null_pos_lower_median": int(np.median(null["pos_lower"])),
        "brg3_accession": BRG3_ACCESSION,
        "brg3_contrast": contrast_brg3,
        "brg3_positive_positions": list(BRG3_POSITIVES),
        "brg3_negative_positions": list(BRG3_NEGATIVES),
        "brg3_trio_metal": {
            str(p): float(trio[p][METAL_FEATURE])
            for p in (*BRG3_POSITIVES, *BRG3_NEGATIVES)
        },
        "brg3_all_cys_metal": {str(p): v for p, v in metal_all.items()},
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    # --- console digest ------------------------------------------------------
    print(
        f"metal_coordination_sg_3a across co-peptide pool "
        f"(n={len(evaluated)} peptides, B={N_PERM})"
    )
    print(
        f"  posHi={observed['pos_higher']} posLo={observed['pos_lower']} "
        f"mixed={observed['mixed']}  pHi={p_hi:.3f} pLo={p_lo:.3f} "
        f"(nullMed {int(np.median(null['pos_higher']))}/"
        f"{int(np.median(null['pos_lower']))})"
    )
    print("\nBRG3 (A0A3Q7EW23) co-peptide trio SSCMICLPCR "
          f"POS={list(BRG3_POSITIVES)} NEG={list(BRG3_NEGATIVES)}:")
    print(f"  => metal contrast {contrast_brg3}")
    for p in (*BRG3_POSITIVES, *BRG3_NEGATIVES):
        label = "POS" if p in BRG3_POSITIVES else "NEG"
        print(
            f"    Cys{p:<4} ({label})  metal={trio[p][METAL_FEATURE]:.0f}  "
            f"cn_sg6a={trio[p]['contact_number_sg_6a']:.0f}  "
            f"rsa={trio[p]['rsa_relative']:.3f}"
        )
    print("\nAll BRG3 Cys, metal_coordination_sg_3a (ascending):")
    for p, v in sorted(metal_all.items(), key=lambda kv: kv[1]):
        mark = ""
        if p in BRG3_POSITIVES:
            mark = "  <== POS (modified)"
        elif p in BRG3_NEGATIVES:
            mark = "  <== NEG (unmodified)"
        print(f"    Cys{p:<4}  metal={v:.0f}{mark}")
    print(f"\nwrote {OUTPUT}")


if __name__ == "__main__":
    main()
