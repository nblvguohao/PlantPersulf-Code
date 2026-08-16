"""Diagnostic structure-feature pipeline for the P2 separation gate.

The methodology design's P2 early gate asks whether structure-derived
features (solvent accessibility, metal-cluster geometry, electrostatic
microenvironment) can separate co-peptide Cys that a sequence-window model
cannot (BRG3 C206/C212 modified vs C209 unmodified). This module computes
those features from AlphaFold DB PDB files with pure numpy — no DSSP /
FreeSASA / APBS dependency — at the precision needed for a gate test:

- ``plddt`` (AF B-factor column) — disorder proxy for the IDR regime;
- ``rsa_relative`` — Shrake-Rupley side-chain SASA / reference A^2 (Tien
  et al. 2013 values);
- ``cys_count_8a`` / ``nearest_sg_distance`` — metal-cluster geometry
  PROXY: AFDB models contain no metal ions, so the RING/Zn-finger cluster
  is represented by S-gamma neighbourhood density (limitation recorded;
  Metal3D / AF3-with-Zn is a follow-up);
- ``positive_residue_count_6a`` — Arg/Lys/His side chains within 6 A of
  the S-gamma (electrostatic microenvironment, methodology L1 table);
- ``coulomb_potential_sg`` — Coulomb potential at the S-gamma from
  charged side-chain partial charges (dielectric folded into a constant;
  APBS Poisson-Boltzmann is the eventual rigorous feature);
- ``contact_number_10a`` — C-alpha neighbours within 10 A (packing
  density);
- ``contact_number_sg_6a`` — heavy atoms (non-hydrogen) of other residues
  within 6 A of the Cys S-gamma (atom-level local packing around the
  sulphur; the S-gamma-layer replacement for the C-alpha contact count);
- ``metal_coordination_sg_3a`` — N/O/S atoms of His/Cys/Asp/Glu residues
  within 3 A of the Cys S-gamma (metal-coordination proxy: AFDB models carry
  no metal ions, so Zn coordination is represented by ligand-atom density
  around the sulphur).

Diagnostic-only: no fitting, no mutation of frozen artifacts.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

import numpy as np

# --- atom geometry ----------------------------------------------------------

VAN_DER_WAALS_RADIUS = {
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "S": 1.8,
    "H": 1.2,
}

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# side-chain atoms per residue (excludes main-chain N/CA/C/O/OXT)
_SIDE_CHAIN = {
    "A": ("CB",), "R": ("CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"),
    "N": ("CB", "CG", "OD1", "ND2"), "D": ("CB", "CG", "OD1", "OD2"),
    "C": ("CB", "SG"), "E": ("CB", "CG", "CD", "OE1", "OE2"),
    "Q": ("CB", "CG", "CD", "OE1", "NE2"), "G": (),
    "H": ("CB", "CG", "ND1", "CD2", "CE1", "NE2"), "I": ("CB", "CG1", "CG2", "CD1"),
    "L": ("CB", "CG", "CD1", "CD2"), "K": ("CB", "CG", "CD", "CE", "NZ"),
    "M": ("CB", "CG", "SD", "CE"), "F": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "P": ("CB", "CG", "CD"), "S": ("CB", "OG"), "T": ("CB", "OG1", "CG2"),
    "W": ("CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "Y": ("CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"),
    "V": ("CB", "CG1", "CG2"),
}

# side-chain relative solvent accessibility reference (A^2), Tien et al. 2013
RSA_REFERENCE = {"C": 104.6, "M": 194.1, "F": 199.5, "I": 182.2, "L": 183.1,
                 "V": 144.1, "W": 249.4, "Y": 212.2, "T": 116.2, "S": 92.9,
                 "R": 225.6, "K": 202.5, "H": 180.7, "D": 156.6, "E": 183.2,
                 "N": 134.2, "Q": 172.4, "P": 132.7, "A": 71.6, "G": 52.5}

# simplified AMBER-style partial charges (e) on side-chain atoms
# (keys are one-letter residue codes, matching ``Residue.type``)
_PARTIAL_CHARGES = {
    "R": {"NE": 1 / 3, "NH1": 1 / 3, "NH2": 1 / 3},
    "K": {"NZ": 1.0},
    "D": {"OD1": -0.5, "OD2": -0.5},
    "E": {"OE1": -0.5, "OE2": -0.5},
}

POSITIVE_RESIDUES = ("R", "K", "H")

S_GAMMA_HEAVY_ATOM_RADIUS_ANGSTROM = 6.0
METAL_COORDINATION_RADIUS_ANGSTROM = 3.0
METAL_LIGAND_RESIDUES = ("H", "C", "D", "E")

STRUCTURE_FEATURE_NAMES = (
    "plddt",
    "rsa_relative",
    "cys_count_8a",
    "nearest_sg_distance",
    "positive_residue_count_6a",
    "coulomb_potential_sg",
    "contact_number_10a",
    "contact_number_sg_6a",
    "metal_coordination_sg_3a",
)


@dataclass(frozen=True)
class Atom:
    """One PDB atom (backbone or side chain)."""

    name: str
    xyz: np.ndarray[Any, Any]
    plddt: float
    element: str

    @property
    def radius(self) -> float:
        return VAN_DER_WAALS_RADIUS.get(self.element, 1.7)


@dataclass(frozen=True)
class Residue:
    """One residue of a PDB chain, with its atoms in file order."""

    resseq: int
    type: str
    atoms: tuple[Atom, ...]

    def atom(self, name: str) -> Atom | None:
        for atom in self.atoms:
            if atom.name == name:
                return atom
        return None


def _element_from_name(name: str) -> str:
    stripped = name.strip()
    if not stripped:
        return "C"
    if len(stripped) >= 2 and stripped[0].isalpha() and stripped[1].islower():
        return stripped[:2].capitalize()
    return stripped[0]


def parse_pdb(text: str) -> list[Residue]:
    """Parse a PDB model into ordered residues (single chain, first altloc).

    Ignores HETATM, water and alternate conformations other than the first.
    Residues are ordered by first appearance (chain-major, resseq-ordered).
    """
    current: dict[tuple[str, int], list[Atom]] = {}
    order: list[tuple[str, int]] = []
    resname_by_key: dict[tuple[str, int], str] = {}
    for line in text.splitlines():
        if not line.startswith("ATOM"):
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        chain = line[21]
        resseq = int(line[22:26])
        key = (chain, resseq)
        if key not in current:
            current[key] = []
            order.append(key)
            resname_by_key[key] = line[17:20].strip()
        current[key].append(
            Atom(
                name=line[12:16].strip(),
                xyz=np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                ),
                plddt=float(line[60:66]),
                element=_element_from_name(line[12:16]),
            )
        )
    residues: list[Residue] = []
    for key in order:
        residues.append(
            Residue(
                resseq=key[1],
                type=AA3_TO_AA1.get(resname_by_key[key], "X"),
                atoms=tuple(current[key]),
            )
        )
    return residues


# --- Shrake-Rupley SASA -----------------------------------------------------

def _fibonacci_sphere(n_points: int) -> np.ndarray[Any, Any]:
    """Approximately uniform points on the unit sphere."""
    indices = np.arange(n_points, dtype=np.float64)
    phi = math.pi * (3.0 - math.sqrt(5.0))
    y = 1.0 - 2.0 * (indices + 0.5) / n_points
    radius = np.sqrt(np.maximum(0.0, 1.0 - y * y))
    theta = phi * indices
    return np.stack([radius * np.cos(theta), y, radius * np.sin(theta)], axis=1)


def sasa_shrake_rupley(
    coords: np.ndarray[Any, Any],
    radii: np.ndarray[Any, Any],
    probe: float = 1.4,
    n_points: int = 512,
) -> np.ndarray[Any, Any]:
    """Solvent-accessible surface area per atom (Shrake-Rupley, numpy).

    ``coords`` (n, 3) and ``radii`` (n,) in Angstrom. Points of the
    probe-expanded sphere of each atom are occluded if they fall inside the
    probe-expanded sphere of any other atom.
    """
    coords = np.asarray(coords, dtype=np.float64)
    radii = np.asarray(radii, dtype=np.float64)
    n_atoms = len(coords)
    unit_points = _fibonacci_sphere(n_points)
    sasa = np.zeros(n_atoms)
    for i in range(n_atoms):
        center = coords[i]
        radius = radii[i] + probe
        points = center + radius * unit_points  # (n_points, 3)
        # neighbours within reach of the probe-expanded surface
        deltas = coords - center
        distances = np.sqrt(np.einsum("ij,ij->i", deltas, deltas))
        reach = 2 * probe + radii[i] + radii
        neighbours = np.where((distances < reach) & (np.arange(n_atoms) != i))[0]
        if neighbours.size == 0:
            sasa[i] = 4 * math.pi * radius**2
            continue
        visible = np.ones(n_points, dtype=bool)
        for j in neighbours:
            occlusion = radii[j] + probe
            point_deltas = points - coords[j]
            point_dist = np.sqrt(
                np.einsum("ij,ij->i", point_deltas, point_deltas)
            )
            visible &= point_dist >= occlusion - 1e-6
        sasa[i] = visible.sum() / n_points * 4 * math.pi * radius**2
    return sasa


# --- Coulomb potential ------------------------------------------------------

class CoulombPotentialEstimator:
    """Coulomb potential at a probe point from charged side-chain atoms.

    Dielectric constant and the 1/(4 pi epsilon_0) factor are folded into a
    unit scale: only RELATIVE values are used. ``cutoff`` A distance
    truncation.
    """

    def __init__(self, cutoff: float = 20.0) -> None:
        self.cutoff = cutoff

    def potential_at(
        self,
        probe: tuple[float, float, float],
        charges: Iterable[tuple[float, tuple[float, float, float]]],
    ) -> float:
        probe_array = np.asarray(probe, dtype=np.float64)
        potential = 0.0
        for charge, position in charges:
            delta = np.asarray(position, dtype=np.float64) - probe_array
            distance = float(np.sqrt(np.dot(delta, delta)))
            if distance > self.cutoff:
                continue
            potential += charge / max(distance, 1e-6)
        return potential


# --- per-Cys feature assembly ----------------------------------------------

def _side_chain_atoms(residue: Residue) -> list[Atom]:
    return [
        atom
        for atom in residue.atoms
        if atom.name in _SIDE_CHAIN.get(residue.type, ())
    ]


def _heavy_atom_contacts_within(
    residues: list[Residue],
    probe_xyz: np.ndarray[Any, Any],
    self_resseq: int,
    radius_angstrom: float,
) -> int:
    """Count heavy atoms (non-hydrogen) of OTHER residues within a radius of
    a probe point (the Cys S-gamma). The same residue's own atoms are excluded
    so the count measures local packing *around* the sulphur rather than the
    residue's constant internal geometry (its own C-beta is ~1.8 A away)."""
    count = 0
    for residue in residues:
        if residue.resseq == self_resseq:
            continue
        for atom in residue.atoms:
            if atom.element == "H":
                continue
            delta = atom.xyz - probe_xyz
            if float(np.sqrt(np.dot(delta, delta))) <= radius_angstrom:
                count += 1
    return count


def _metal_coordination_count(
    residues: list[Residue],
    probe_xyz: np.ndarray[Any, Any],
    self_resseq: int,
    radius_angstrom: float = METAL_COORDINATION_RADIUS_ANGSTROM,
) -> int:
    """Count N/O/S atoms of His/Cys/Asp/Glu residues within a radius of the
    Cys S-gamma (metal-coordination proxy). AFDB models carry no metal ions,
    so a Zn-coordinating Cys is recognised by the density of nearby ligand
    atoms (His ND1/NE2, Cys SG, Asp OD1/OD2, Glu OE1/OE2) around its sulphur.
    The same residue's own atoms are excluded."""
    count = 0
    for residue in residues:
        if residue.type not in METAL_LIGAND_RESIDUES:
            continue
        if residue.resseq == self_resseq:
            continue
        for atom in residue.atoms:
            if atom.element not in ("N", "O", "S"):
                continue
            delta = atom.xyz - probe_xyz
            if float(np.sqrt(np.dot(delta, delta))) <= radius_angstrom:
                count += 1
    return count


def cys_structure_features(
    residues: list[Residue], positions: Iterable[int]
) -> dict[int, dict[str, float]]:
    """Structure features for every requested Cys position (1-based).

    Raises ``ValueError`` when a requested position is not a Cys in the
    structure. Features are relative/geometric, so no reference proteome is
    needed beyond the structure itself.
    """
    by_resseq = {residue.resseq: residue for residue in residues}
    missing = [p for p in positions if p not in by_resseq]
    if missing:
        raise ValueError(f"positions missing from structure: {missing}")

    all_atoms = [atom for residue in residues for atom in residue.atoms]
    coords = np.array([atom.xyz for atom in all_atoms])
    radii = np.array([atom.radius for atom in all_atoms])
    sasa = sasa_shrake_rupley(coords, radii)

    # atom index lookup by residue
    atom_index: dict[tuple[int, str], Atom] = {}
    offset = 0
    for residue in residues:
        for atom in residue.atoms:
            atom_index[(residue.resseq, atom.name)] = atom
        offset += len(residue.atoms)

    # per-residue side-chain SASA (sum of its side-chain atom SASAs)
    side_sasa: dict[int, float] = {}
    offset = 0
    for residue in residues:
        side_sasa[residue.resseq] = sum(
            sasa[offset + index]
            for index, atom in enumerate(residue.atoms)
            if atom.name in _SIDE_CHAIN.get(residue.type, ())
        )
        offset += len(residue.atoms)

    # Cys S-gamma coordinates of every Cys
    sg_by_resseq: dict[int, np.ndarray[Any, Any]] = {}
    for residue in residues:
        if residue.type != "C":
            continue
        first_sg = residue.atom("SG")
        if first_sg is not None:
            sg_by_resseq[residue.resseq] = first_sg.xyz

    ca_by_resseq: dict[int, np.ndarray[Any, Any]] = {}
    for residue in residues:
        ca = residue.atom("CA")
        if ca is not None:
            ca_by_resseq[residue.resseq] = ca.xyz

    features: dict[int, dict[str, float]] = {}
    for position in positions:
        residue = by_resseq[position]
        if residue.type != "C":
            raise ValueError(f"position {position} is {residue.type}, not Cys")
        sg = residue.atom("SG")
        if sg is None:
            raise ValueError(f"position {position} Cys has no SG atom")
        ca = residue.atom("CA")
        plddt = ca.plddt if ca is not None else sg.plddt

        sg_deltas = {
            other: np.sqrt(np.dot(sg.xyz - other_sg, sg.xyz - other_sg))
            for other, other_sg in sg_by_resseq.items()
            if other != position
        }
        cys_count_8a = sum(1 for d in sg_deltas.values() if d <= 8.0)
        nearest_sg_distance = (
            min(sg_deltas.values()) if sg_deltas else 99.0
        )

        # positive residues within 6 A of SG (side-chain atoms)
        positive_count = 0
        for neighbor in residues:
            if neighbor.type not in POSITIVE_RESIDUES:
                continue
            close = any(
                np.sqrt(np.dot(atom.xyz - sg.xyz, atom.xyz - sg.xyz)) <= 6.0
                for atom in neighbor.atoms
            )
            if close:
                positive_count += 1

        # Coulomb potential at SG
        charges: list[tuple[float, tuple[float, float, float]]] = []
        for neighbor in residues:
            for name, charge in _PARTIAL_CHARGES.get(neighbor.type, {}).items():
                charge_atom = neighbor.atom(name)
                if charge_atom is not None:
                    charges.append((charge, tuple(charge_atom.xyz)))
        estimator = CoulombPotentialEstimator()
        coulomb = estimator.potential_at(tuple(sg.xyz), charges)

        # contact number: CA within 10 A (excluding self)
        contact_number = sum(
            1
            for other, other_ca in ca_by_resseq.items()
            if other != position
            and np.sqrt(np.dot(ca.xyz - other_ca, ca.xyz - other_ca)) <= 10.0
        ) if ca is not None else 0

        reference = RSA_REFERENCE.get(residue.type, 100.0)
        contact_number_sg_6a = float(
            _heavy_atom_contacts_within(
                residues,
                sg.xyz,
                position,
                S_GAMMA_HEAVY_ATOM_RADIUS_ANGSTROM,
            )
        )
        metal_coordination_sg_3a = float(
            _metal_coordination_count(residues, sg.xyz, position)
        )
        features[position] = {
            "plddt": float(plddt),
            "rsa_relative": float(side_sasa[position]) / reference,
            "cys_count_8a": float(cys_count_8a),
            "nearest_sg_distance": float(nearest_sg_distance),
            "positive_residue_count_6a": float(positive_count),
            "coulomb_potential_sg": float(coulomb),
            "contact_number_10a": float(contact_number),
            "contact_number_sg_6a": contact_number_sg_6a,
            "metal_coordination_sg_3a": metal_coordination_sg_3a,
        }
    return features


def extract_cys_structure_features(
    residues: list[Residue], positions: Iterable[int]
) -> dict[int, dict[str, float]]:
    """Cys structure features (kept as a thin named alias for callers)."""
    return cys_structure_features(residues, positions)
