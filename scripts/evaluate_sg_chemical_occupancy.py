#!/usr/bin/env python
"""Chemical-occupancy last shot: is the modified Cys S-gamma *occupied*?

The packing-density class (contact number / SASA / depth) is structurally
blind to chemical occupancy: a Zn-coordinated Cys and a free Cys can have
identical packing density, because every ligand of a metal site sits in the
same crowded pocket. contact_number at any radius cannot see the difference —
not a resolution issue, a wrong-question issue. This diagnostic tests the
occupancy class directly:

- **BRG3 (pre-declared single case)**: cluster the RING-region S-gamma atoms
  (Cys 197-231) geometrically into the two Zn sites. A tetrahedral Zn site's
  ligand S-gamma atoms sit ~3.5-4 A apart (AF has no Zn ion, but a RING folds
  around its Zn and at high pLDDT the holo ligand geometry is reproduced), so
  S-gamma mutual-distance clustering reads out the sites. Direction is NOT
  fixed: coordinated Cys could be locked (no room to modify) OR could be the
  low-pKa thiolates a persulfidation "zinc-finger knockout" targets. Either
  direction — if coordination state separates C206/C212 (modified) from
  C209 (unmodified) — is a discriminator invisible to every tested feature.
- **Pool (per-group, no aggregate permutation)**: the occupancy feature
  (nearest S-gamma distance) has near-zero variance in groups with no metal
  site or disulfide; an aggregate permutation null is structurally incapable
  of finding a signal that lives in the few occupied groups. So groups are
  first filtered to those where the feature actually varies (a Cys with a
  close S-gamma partner), then reported one by one.

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

from plantpersulf.evaluation.structure_features import (  # noqa: E402
    cys_structure_features,
    parse_pdb,
)
from plantpersulf.features.sequence import _load_proteome  # noqa: E402

ALPHAFOLD_DIR = _REPO_ROOT / "data" / "raw" / "alphafold"
CO_PEPTIDE_REGISTRY = _REPO_ROOT / "data" / "registry" / "copeptide_negatives_v1.tsv"
OUTPUT = _REPO_ROOT / "results" / "diagnostics" / "sg_chemical_occupancy_v1.json"
N_PERM = 999
SEED = 20260816

BRG3_ACCESSION = "A0A3Q7EW23"
BRG3_RING_START = 197
BRG3_RING_END = 231
BRG3_POSITIVES = (206, 212)
BRG3_NEGATIVES = (209,)
# Zn-cluster ligand-ligand S-gamma range (tetrahedral Zn: S-S ~3.5-4.0 A).
ZINC_CLUSTER_CUTOFF = 4.5
# A Cys is "occupied" when its nearest other S-gamma is this close (disulfide
# ~2 A, Zn-cluster partner ~3.5-4 A).
OCCUPANCY_CUTOFF = 4.0

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
                "domain": row.get("source_detail", ""),
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
    """Structure features for every Cys of each group's protein."""
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


def _brg3_ring_sg(residues: list, ring_start: int, ring_end: int) -> dict[int, dict]:
    """RING-region Cys S-gamma coordinates, CA pLDDT and residue type."""
    out: dict[int, dict] = {}
    for residue in residues:
        if residue.type != "C":
            continue
        if not (ring_start <= residue.resseq <= ring_end):
            continue
        sg = residue.atom("SG")
        ca = residue.atom("CA")
        if sg is None or ca is None:
            continue
        out[residue.resseq] = {
            "sg_xyz": sg.xyz,
            "plddt": float(ca.plddt),
        }
    return out


def _sg_distance_matrix(
    ring: dict[int, dict],
) -> tuple[dict[int, dict[int, float]], np.ndarray]:
    """Pairwise S-gamma distances among RING Cys (and a matrix for clustering)."""
    positions = sorted(ring)
    index = {p: i for i, p in enumerate(positions)}
    matrix = np.full((len(positions), len(positions)), np.inf)
    for i, p in enumerate(positions):
        for j, q in enumerate(positions):
            if i == j:
                continue
            delta = ring[p]["sg_xyz"] - ring[q]["sg_xyz"]
            matrix[i, j] = float(np.sqrt(np.dot(delta, delta)))
    distances: dict[int, dict[int, float]] = {}
    for i, p in enumerate(positions):
        distances[p] = {
            q: round(float(matrix[i, index[q]]), 2)
            for q in positions
            if q != p
        }
    return distances, matrix


def _cluster_zinc_sites(
    positions: list[int],
    matrix: np.ndarray,
    cutoff: float,
) -> tuple[list[list[int]], list[int]]:
    """Connected components of the S-gamma distance graph at ``cutoff`` A.

    A tetrahedral Zn site's ligand S-gamma atoms are all ~3.5-4 A apart, so
    they form one dense connected component; non-ligand Cys (no partner within
    the cluster range) stay isolated. Returns component members and isolated
    Cys as **sequence positions** (not matrix indices).
    """
    n = len(positions)
    adjacency = [
        [j for j in range(n) if j != i and matrix[i, j] <= cutoff]
        for i in range(n)
    ]
    seen: set[int] = set()
    index_components: list[list[int]] = []
    for start in range(n):
        if start in seen:
            continue
        stack = [start]
        comp: set[int] = set()
        while stack:
            node = stack.pop()
            if node in comp:
                continue
            comp.add(node)
            for neighbor in adjacency[node]:
                if neighbor not in comp:
                    stack.append(neighbor)
        seen |= comp
        index_components.append(sorted(comp))
    index_components = [c for c in index_components if len(c) >= 2]
    component_positions = [
        [positions[i] for i in component] for component in index_components
    ]
    flat = set(sum(component_positions, []))
    isolated = [p for p in positions if p not in flat]
    return component_positions, isolated


def main() -> None:
    groups = _load_copeptide_groups()
    feature_rows, skipped = _load_feature_rows(groups)

    # --- BRG3 single-case: RING S-gamma clustering + pLDDT -------------------
    pdb_path = ALPHAFOLD_DIR / f"AF-{BRG3_ACCESSION}-F1-model.pdb"
    residues = parse_pdb(pdb_path.read_text(encoding="utf-8"))
    ring = _brg3_ring_sg(residues, BRG3_RING_START, BRG3_RING_END)
    ring_positions = sorted(ring)
    distances, matrix = _sg_distance_matrix(ring)
    components, isolated = _cluster_zinc_sites(
        ring_positions, matrix, ZINC_CLUSTER_CUTOFF
    )
    brg3_features = feature_rows[BRG3_ACCESSION]

    ring_report = {
        str(p): {
            "plddt": round(ring[p]["plddt"], 1),
            "nearest_sg_distance": brg3_features[p]["nearest_sg_distance"],
            "nearest_sg_to": min(distances[p], key=distances[p].get),
        }
        for p in ring_positions
    }
    brg3_trio = {
        p: {
            "label": "POS" if p in BRG3_POSITIVES else "NEG",
            "nearest_sg_distance": brg3_features[p]["nearest_sg_distance"],
            "occupied": (
                brg3_features[p]["nearest_sg_distance"] <= OCCUPANCY_CUTOFF
            ),
            "in_zinc_cluster": any(p in comp for comp in components),
        }
        for p in (*BRG3_POSITIVES, *BRG3_NEGATIVES)
    }

    # --- Pool: per-group occupancy, variance-filtered, no merging -----------
    per_group: list[dict] = []
    for group in groups:
        if group["accession"] not in feature_rows:
            continue
        features = feature_rows[group["accession"]]
        positions = group["positive_positions"] + group["negative_positions"]
        nearest = {
            p: round(float(features[p]["nearest_sg_distance"]), 2) for p in positions
        }
        occupied = {p: nearest[p] <= OCCUPANCY_CUTOFF for p in positions}
        has_variance = any(occupied.values())
        if not has_variance:
            continue
        pos_vals = [nearest[p] for p in group["positive_positions"]]
        neg_vals = [nearest[p] for p in group["negative_positions"]]
        if min(pos_vals) > max(neg_vals):
            direction = "pos_farther"  # modified Cys further from any S-gamma partner
        elif max(pos_vals) < min(neg_vals):
            direction = "pos_closer"  # modified Cys closer to an S-gamma partner
        else:
            direction = "mixed"
        per_group.append(
            {
                "accession": group["accession"],
                "species": group["species"],
                "peptide": group["peptide"],
                "domain": group["domain"],
                "positive_positions": group["positive_positions"],
                "negative_positions": group["negative_positions"],
                "nearest_sg_distance": nearest,
                "occupied": {str(p): v for p, v in occupied.items()},
                "direction": direction,
            }
        )

    document = {
        "track": "sg_chemical_occupancy_v1",
        "claim_class": "diagnostic_only",
        "note": (
            "Chemical-occupancy class (is the S-gamma already engaged in a "
            "disulfide or Zn cluster), tested per-group WITHOUT aggregate "
            "permutation: the feature has near-zero variance in the many "
            "peptides with no metal site or disulfide, so groups are first "
            f"filtered to those where a Cys has a close S-gamma partner "
            f"(<= {OCCUPANCY_CUTOFF} A), then reported one by one. BRG3 is a "
            "pre-declared single-case prediction test; direction is left free "
            "(coordinated Cys locked vs low-pKa thiolate target), either "
            "clean separation of C206/C212 (modified) from C209 (unmodified) "
            "counts as a discriminator invisible to the packing-density "
            "features."
        ),
        "zinc_cluster_cutoff_angstrom": ZINC_CLUSTER_CUTOFF,
        "occupancy_cutoff_angstrom": OCCUPANCY_CUTOFF,
        "brg3_accession": BRG3_ACCESSION,
        "brg3_ring_region": [BRG3_RING_START, BRG3_RING_END],
        "brg3_ring_cys": [str(p) for p in ring_positions],
        "brg3_zinc_sites_clusters": [
            [str(p) for p in comp] for comp in components
        ],
        "brg3_ring_isolated_cys": [str(p) for p in isolated],
        "brg3_ring_report": ring_report,
        "brg3_trio": {str(p): v for p, v in brg3_trio.items()},
        "pool": {
            "n_groups_total": len(groups),
            "n_groups_with_occupancy_variance": len(per_group),
            "groups": per_group,
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(document, indent=2, ensure_ascii=False, default=str) + "\n",
        encoding="utf-8",
    )

    # --- console -------------------------------------------------------------
    print(f"BRG3 RING ({BRG3_RING_START}-{BRG3_RING_END}) S-gamma Zn clustering "
          f"(cutoff {ZINC_CLUSTER_CUTOFF} A)")
    print(f"  RING Cys: {ring_positions}")
    print(f"  Zn-site clusters: {components}")
    print(f"  isolated (no cluster partner): {isolated}")
    print("\n  per-Cys  pLDDT  nearestSG  partner")
    for p in ring_positions:
        rp = ring_report[str(p)]
        if p in BRG3_POSITIVES:
            mark = "  <== POS"
        elif p in BRG3_NEGATIVES:
            mark = "  <== NEG"
        else:
            mark = ""
        print(
            f"    Cys{p:<4} {rp['plddt']:>5}  {rp['nearest_sg_distance']:>8.2f}  "
            f"Cys{rp['nearest_sg_to']:<4}{mark}"
        )
    print("\n  BRG3 co-peptide trio (POS=C206,C212 / NEG=C209):")
    for p in (*BRG3_POSITIVES, *BRG3_NEGATIVES):
        t = brg3_trio[p]
        print(
            f"    Cys{p} ({t['label']}): nearestSG={t['nearest_sg_distance']:.2f} "
            f"occupied={t['occupied']} in_zinc_cluster={t['in_zinc_cluster']}"
        )

    print(
        f"\nPool: {len(per_group)}/{len(groups)} groups have occupancy variance "
        f"(a Cys with nearest S-gamma <= {OCCUPANCY_CUTOFF} A)"
    )
    for g in per_group:
        print(
            f"  {g['species']:<11} {g['accession']:<11} {g['peptide'][:28]:<28} "
            f"POS={g['positive_positions']} NEG={g['negative_positions']} "
            f"nearestSG={g['nearest_sg_distance']} => {g['direction']}"
        )
    print(f"\nwrote {OUTPUT}")


if __name__ == "__main__":
    main()
