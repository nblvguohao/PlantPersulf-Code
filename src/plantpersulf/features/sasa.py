r"""Real Shrake-Rupley solvent-accessible surface area (SASA).

Complements ``features/structure.py``'s ``contact_number_proxy`` (explicitly
documented there as NOT a rigorous SASA calculation — a contact-count
correlate of burial) with an actual Shrake-Rupley rolling-probe calculation
on the full (non-CA-only) heavy-atom set. Added 2026-07-23 in response to a
self-review finding: the structural-context conservation claim
(§5.5.2, ``docs/phase_z_evidence_audit.md``) relied only on the proxy.

**Method**: for each target atom, ``n_sphere_points`` are placed on a sphere
of radius ``vdW(atom) + probe_radius`` via a deterministic Fibonacci-lattice
distribution (no randomness — SASA is exactly reproducible given the same
point count). A point is "buried" if it falls within any neighboring atom's
own probe-inflated van der Waals sphere; the exposed fraction times the
sphere's surface area is that atom's SASA contribution. A residue's SASA is
the sum over its atoms.

**Scope note**: this project's persulfidation-site comparisons are always
Cys-vs-Cys (persulfidated vs. other cysteines in the *same* protein — see
``evaluation/cross_species_structural_context.py``), so no cross-residue-type
relative-SASA normalization (e.g. Tien et al. 2013 reference maxima) is
needed or implemented; raw absolute SASA in Å² is already directly
comparable within this same-residue-type comparison.

**Performance**: SASA for only the requested atoms (typically one cysteine
residue, ~6 atoms) is computed against a spatial grid of the full
structure's atoms, not against every atom in the protein — this keeps the
per-residue query tractable across thousands of residues without an
external spatial-indexing dependency.

**No dependency added**: pure Python (``math`` only), consistent with this
project's existing from-scratch numerical code (``permutation.py``,
``cross_species_conservation.py``'s hypergeometric tail, and
``structure.py`` itself).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# Bondi (1964) van der Waals radii, Angstrom — the standard minimal set for
# the heavy atoms present in AlphaFold DB's hydrogen-free PDB output.
VDW_RADII_ANGSTROM: dict[str, float] = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
}
PROBE_RADIUS_ANGSTROM = 1.4  # water probe, standard Shrake-Rupley choice
DEFAULT_N_SPHERE_POINTS = 92


@dataclass(frozen=True)
class Atom:
    chain_id: str
    res_seq: int
    res_name: str
    atom_name: str
    element: str
    x: float
    y: float
    z: float


def parse_all_atoms(pdb_text: str) -> tuple[Atom, ...]:
    """Parse every ``ATOM`` record (not just C-alpha) from a standard PDB
    file, using the element-symbol field (wwPDB columns 77-78) rather than
    inferring element from the atom name (which is ambiguous — e.g. ``CA``
    is the alpha-carbon atom name, not the element calcium).

    Atoms whose element is not in ``VDW_RADII_ANGSTROM`` (i.e. not C/N/O/S
    — AlphaFold DB output has no hydrogens and no heteroatoms in the
    standard amino-acid chain) are skipped rather than guessed.
    """
    atoms: list[Atom] = []
    for line in pdb_text.splitlines():
        if len(line) < 78 or line[0:6] != "ATOM  ":
            continue
        try:
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21]
            res_seq = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        element = line[76:78].strip().upper()
        if element not in VDW_RADII_ANGSTROM:
            continue
        atoms.append(
            Atom(
                chain_id=chain_id,
                res_seq=res_seq,
                res_name=res_name,
                atom_name=atom_name,
                element=element,
                x=x,
                y=y,
                z=z,
            )
        )
    return tuple(atoms)


def _fibonacci_sphere_points(
    n_points: int,
) -> tuple[tuple[float, float, float], ...]:
    """Deterministic, near-evenly-distributed unit-sphere points via the
    Fibonacci lattice — no RNG, so SASA is exactly reproducible for a fixed
    ``n_points``."""
    if n_points < 2:
        raise ValueError("n_points must be >= 2")
    points: list[tuple[float, float, float]] = []
    golden_angle = math.pi * (3.0 - math.sqrt(5.0))
    for i in range(n_points):
        y = 1.0 - (i / (n_points - 1)) * 2.0
        radius_at_y = math.sqrt(max(0.0, 1.0 - y * y))
        theta = golden_angle * i
        x = math.cos(theta) * radius_at_y
        z = math.sin(theta) * radius_at_y
        points.append((x, y, z))
    return tuple(points)


class _SpatialGrid:
    """Uniform-cell spatial hash for atom neighbor queries.

    Cell size is fixed to ``cutoff`` so that any atom within ``cutoff`` of a
    query atom is guaranteed to lie in one of the 27 cells (3x3x3) centered
    on the query atom's own cell — no true neighbor is ever missed.
    """

    def __init__(self, atoms: tuple[Atom, ...], cutoff: float) -> None:
        self._cutoff = cutoff
        self._cells: dict[tuple[int, int, int], list[int]] = {}
        for i, atom in enumerate(atoms):
            key = self._cell_key(atom.x, atom.y, atom.z)
            self._cells.setdefault(key, []).append(i)

    def _cell_key(self, x: float, y: float, z: float) -> tuple[int, int, int]:
        c = self._cutoff
        return (int(math.floor(x / c)), int(math.floor(y / c)), int(math.floor(z / c)))

    def neighbor_indices(self, x: float, y: float, z: float) -> list[int]:
        cx, cy, cz = self._cell_key(x, y, z)
        result: list[int] = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    result.extend(self._cells.get((cx + dx, cy + dy, cz + dz), ()))
        return result


def compute_atom_sasa(
    atoms: tuple[Atom, ...],
    target_indices: list[int],
    n_sphere_points: int = DEFAULT_N_SPHERE_POINTS,
    probe_radius: float = PROBE_RADIUS_ANGSTROM,
) -> dict[int, float]:
    """Shrake-Rupley SASA (Angstrom^2) for the atoms at ``target_indices``.

    Burial is checked against every atom in ``atoms`` that falls within
    reach (via a spatial grid, not a full O(n^2) scan) — including atoms
    outside the target set, so residues from other parts of the folded
    structure correctly bury a target atom's surface.
    """
    max_inflated_radius = max(VDW_RADII_ANGSTROM.values()) + probe_radius
    cutoff = 2.0 * max_inflated_radius
    grid = _SpatialGrid(atoms, cutoff)
    sphere_points = _fibonacci_sphere_points(n_sphere_points)
    cutoff_sq = cutoff * cutoff

    results: dict[int, float] = {}
    for idx in target_indices:
        atom = atoms[idx]
        r_i = VDW_RADII_ANGSTROM[atom.element] + probe_radius

        neighbor_data: list[tuple[float, float, float, float]] = []
        for j in grid.neighbor_indices(atom.x, atom.y, atom.z):
            if j == idx:
                continue
            other = atoms[j]
            dx = other.x - atom.x
            dy = other.y - atom.y
            dz = other.z - atom.z
            if dx * dx + dy * dy + dz * dz <= cutoff_sq:
                r_j = VDW_RADII_ANGSTROM[other.element] + probe_radius
                neighbor_data.append((other.x, other.y, other.z, r_j * r_j))

        exposed = 0
        for px, py, pz in sphere_points:
            wx = atom.x + px * r_i
            wy = atom.y + py * r_i
            wz = atom.z + pz * r_i
            buried = False
            for nx, ny, nz, r_sq in neighbor_data:
                ddx = wx - nx
                ddy = wy - ny
                ddz = wz - nz
                if ddx * ddx + ddy * ddy + ddz * ddz < r_sq:
                    buried = True
                    break
            if not buried:
                exposed += 1

        fraction_exposed = exposed / n_sphere_points
        results[idx] = fraction_exposed * 4.0 * math.pi * r_i * r_i
    return results


@dataclass(frozen=True)
class CysSasaResult:
    residue_sasa: float
    sg_sasa: float | None  # None if the SG atom is missing from the record


def cys_residue_sasa(
    atoms: tuple[Atom, ...],
    chain_id: str,
    res_seq: int,
    n_sphere_points: int = DEFAULT_N_SPHERE_POINTS,
) -> CysSasaResult | None:
    """Whole-residue and side-chain-sulfur (SG) SASA for one cysteine.

    Returns ``None`` if no CYS residue is found at ``(chain_id, res_seq)``
    in ``atoms`` (coordinate-mapping failure — never imputed).
    """
    residue_indices = [
        i
        for i, a in enumerate(atoms)
        if a.chain_id == chain_id and a.res_seq == res_seq and a.res_name == "CYS"
    ]
    if not residue_indices:
        return None

    per_atom = compute_atom_sasa(
        atoms, residue_indices, n_sphere_points=n_sphere_points
    )
    residue_total = sum(per_atom.values())

    sg_sasa: float | None = None
    for i in residue_indices:
        if atoms[i].atom_name == "SG":
            sg_sasa = per_atom[i]
            break

    return CysSasaResult(residue_sasa=residue_total, sg_sasa=sg_sasa)
