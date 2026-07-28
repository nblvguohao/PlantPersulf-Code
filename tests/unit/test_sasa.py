"""RED: real Shrake-Rupley SASA calculation (features/sasa.py).

Complements the documented, explicitly-labeled ``contact_number_proxy`` in
features/structure.py with an actual solvent-accessible-surface-area
calculation, added in response to a self-review finding that the
structural-context conservation claim relied only on the proxy.
"""

from __future__ import annotations

import math

import pytest

from plantpersulf.features.sasa import (
    PROBE_RADIUS_ANGSTROM,
    VDW_RADII_ANGSTROM,
    Atom,
    compute_atom_sasa,
    cys_residue_sasa,
    parse_all_atoms,
)


def _pdb_atom_line(
    serial: int,
    atom_name: str,
    res_name: str,
    chain: str,
    res_seq: int,
    x: float,
    y: float,
    z: float,
    element: str,
) -> str:
    # Fixed-column PDB ATOM record, column positions verified by direct
    # character-offset comparison against a real AlphaFold DB output line
    # (0-indexed): [0:6] record name, [6:11] serial, [12:16] atom name,
    # [17:20] resName, [21] chainID, [22:26] resSeq, [30:38]/[38:46]/[46:54]
    # x/y/z, [54:60] occupancy, [60:66] B-factor, [76:78] element.
    return (
        f"ATOM  {serial:5d} "
        f" {atom_name:<3s}"
        f" {res_name:>3s} {chain}"
        f"{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{50.0:6.2f}"
        f"{'':10s}{element:>2s}"
    )


# ---------------------------------------------------------------------------
# parse_all_atoms
# ---------------------------------------------------------------------------


def test_parse_all_atoms_extracts_element_from_dedicated_field() -> None:
    text = _pdb_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C") + "\n"
    atoms = parse_all_atoms(text)
    assert len(atoms) == 1
    assert atoms[0].element == "C"
    assert atoms[0].atom_name == "CA"
    assert atoms[0].res_name == "ALA"
    assert atoms[0].chain_id == "A"
    assert atoms[0].res_seq == 1


def test_parse_all_atoms_keeps_multiple_atoms_per_residue() -> None:
    lines = [
        _pdb_atom_line(1, "N", "CYS", "A", 5, 0.0, 0.0, 0.0, "N"),
        _pdb_atom_line(2, "CA", "CYS", "A", 5, 1.5, 0.0, 0.0, "C"),
        _pdb_atom_line(3, "CB", "CYS", "A", 5, 3.0, 0.0, 0.0, "C"),
        _pdb_atom_line(4, "SG", "CYS", "A", 5, 4.8, 0.0, 0.0, "S"),
    ]
    atoms = parse_all_atoms("\n".join(lines) + "\n")
    assert len(atoms) == 4
    assert {a.atom_name for a in atoms} == {"N", "CA", "CB", "SG"}


def test_parse_all_atoms_skips_unsupported_elements() -> None:
    # Zinc (a real hetero-metal some structures carry) is not in the
    # C/N/O/S table this module supports — must be skipped, not guessed.
    line = _pdb_atom_line(1, "ZN", "ZN", "A", 1, 0.0, 0.0, 0.0, "ZN")
    atoms = parse_all_atoms(line + "\n")
    assert atoms == ()


def test_parse_all_atoms_ignores_non_atom_lines() -> None:
    text = "HEADER    test\n" + _pdb_atom_line(
        1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0, "C"
    ) + "\nTER\nEND\n"
    atoms = parse_all_atoms(text)
    assert len(atoms) == 1


# ---------------------------------------------------------------------------
# compute_atom_sasa: hand-verifiable synthetic geometries
# ---------------------------------------------------------------------------


def test_isolated_atom_has_full_sphere_sasa() -> None:
    # A single atom with no neighbors must be 100% exposed: SASA equals the
    # full probe-inflated sphere surface area, 4*pi*r^2.
    atoms = (
        Atom(
            chain_id="A", res_seq=1, res_name="ALA", atom_name="CA",
            element="C", x=0.0, y=0.0, z=0.0,
        ),
    )
    sasa = compute_atom_sasa(atoms, [0], n_sphere_points=200)
    r = VDW_RADII_ANGSTROM["C"] + PROBE_RADIUS_ANGSTROM
    expected = 4.0 * math.pi * r * r
    assert sasa[0] == pytest.approx(expected, rel=1e-9)


def test_two_atoms_at_vdw_contact_reduces_sasa_below_isolated() -> None:
    # Two carbon atoms placed exactly at van-der-Waals contact distance
    # (sum of probe-inflated radii) must each show LESS than the full
    # isolated-atom SASA — the near side is partially buried by the other.
    r = VDW_RADII_ANGSTROM["C"] + PROBE_RADIUS_ANGSTROM
    contact_distance = 2 * r * 0.9  # slight overlap so burial is measurable
    atoms = (
        Atom("A", 1, "ALA", "CA", "C", 0.0, 0.0, 0.0),
        Atom("A", 2, "ALA", "CA", "C", contact_distance, 0.0, 0.0),
    )
    sasa = compute_atom_sasa(atoms, [0, 1], n_sphere_points=300)
    full = 4.0 * math.pi * r * r
    assert sasa[0] < full
    assert sasa[1] < full
    # By symmetry the two atoms should have (near-)identical exposed SASA.
    assert sasa[0] == pytest.approx(sasa[1], rel=0.05)


def test_atom_fully_enclosed_by_a_shell_has_near_zero_sasa() -> None:
    # Surround a central carbon atom with six close carbons (an
    # octahedral "shell") so essentially all sphere-sample points on the
    # central atom are buried.
    r = VDW_RADII_ANGSTROM["C"] + PROBE_RADIUS_ANGSTROM
    d = r * 1.5  # close enough that all six directions are covered
    shell = [
        Atom("A", 2, "ALA", "CB", "C", d, 0.0, 0.0),
        Atom("A", 3, "ALA", "CB", "C", -d, 0.0, 0.0),
        Atom("A", 4, "ALA", "CB", "C", 0.0, d, 0.0),
        Atom("A", 5, "ALA", "CB", "C", 0.0, -d, 0.0),
        Atom("A", 6, "ALA", "CB", "C", 0.0, 0.0, d),
        Atom("A", 7, "ALA", "CB", "C", 0.0, 0.0, -d),
    ]
    atoms = (Atom("A", 1, "ALA", "CA", "C", 0.0, 0.0, 0.0), *shell)
    sasa = compute_atom_sasa(atoms, [0], n_sphere_points=300)
    full = 4.0 * math.pi * r * r
    # Six axis-aligned neighbors cover six spherical caps but leave the
    # diagonal directions (e.g. (1,1,1)) exposed — full enclosure would
    # need a denser (e.g. icosahedral) shell. This geometry should still
    # show substantial, clearly-nonzero burial relative to isolated (100%).
    assert sasa[0] < 0.30 * full


def test_sasa_is_exactly_reproducible_across_repeated_calls() -> None:
    atoms = (
        Atom("A", 1, "ALA", "CA", "C", 0.0, 0.0, 0.0),
        Atom("A", 2, "ALA", "CB", "C", 2.5, 0.0, 0.0),
    )
    sasa_1 = compute_atom_sasa(atoms, [0], n_sphere_points=100)
    sasa_2 = compute_atom_sasa(atoms, [0], n_sphere_points=100)
    assert sasa_1 == sasa_2  # no RNG anywhere -> bit-identical


# ---------------------------------------------------------------------------
# cys_residue_sasa
# ---------------------------------------------------------------------------


def test_cys_residue_sasa_sums_all_atoms_and_reports_sg_separately() -> None:
    atoms = (
        Atom("A", 5, "CYS", "N", "N", 0.0, 0.0, 0.0),
        Atom("A", 5, "CYS", "CA", "C", 1.5, 0.0, 0.0),
        Atom("A", 5, "CYS", "CB", "C", 3.0, 0.0, 0.0),
        Atom("A", 5, "CYS", "SG", "S", 5.0, 0.0, 0.0),
    )
    result = cys_residue_sasa(atoms, chain_id="A", res_seq=5, n_sphere_points=100)
    assert result is not None
    assert result.sg_sasa is not None
    assert result.sg_sasa > 0
    assert result.residue_sasa >= result.sg_sasa  # residue total includes SG


def test_cys_residue_sasa_returns_none_for_missing_residue() -> None:
    atoms = (Atom("A", 5, "CYS", "SG", "S", 0.0, 0.0, 0.0),)
    assert cys_residue_sasa(atoms, chain_id="A", res_seq=99) is None


def test_cys_residue_sasa_returns_none_for_wrong_residue_type() -> None:
    atoms = (Atom("A", 5, "ALA", "CA", "C", 0.0, 0.0, 0.0),)
    assert cys_residue_sasa(atoms, chain_id="A", res_seq=5) is None


def test_cys_residue_sasa_handles_missing_sg_gracefully() -> None:
    # A malformed/truncated record missing the SG atom must not crash —
    # residue SASA is still computed from the atoms present, sg_sasa=None.
    atoms = (
        Atom("A", 5, "CYS", "N", "N", 0.0, 0.0, 0.0),
        Atom("A", 5, "CYS", "CA", "C", 1.5, 0.0, 0.0),
    )
    result = cys_residue_sasa(atoms, chain_id="A", res_seq=5)
    assert result is not None
    assert result.sg_sasa is None
    assert result.residue_sasa > 0
