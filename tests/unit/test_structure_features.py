"""Unit tests for the diagnostic structure-feature pipeline (P2 gate).

Covers PDB parsing (residue order, coordinates, pLDDT), the numpy
Shrake-Rupley SASA implementation, S-gamma geometry, positive-residue
counts, the Coulomb potential approximation and per-Cys feature assembly.
The pipeline is diagnostic-only: it feeds the P2 structure-separation
test, not the frozen bundle.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from plantpersulf.evaluation.structure_features import (
    CoulombPotentialEstimator,
    cys_structure_features,
    extract_cys_structure_features,
    parse_pdb,
    sasa_shrake_rupley,
)

PDB_SNIPPET = (
    """ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 85.30           N
ATOM      2  CA  MET A   1       1.452   0.000   0.000  1.00 84.90           C
ATOM      3  C   MET A   1       2.043   1.411   0.000  1.00 84.50           C
ATOM      4  O   MET A   1       3.245   1.598   0.000  1.00 84.10           O
ATOM      5  CB  MET A   1       1.980  -0.736   1.230  1.00 84.00           C
ATOM      6  SG  MET A   1       3.462  -0.236   1.662  1.00 83.60           S
ATOM      7  CE  MET A   1       4.195  -1.671   2.481  1.00 83.20           C
ATOM      8  N   CYS A   2       1.352   2.411   0.000  1.00 70.00           N
ATOM      9  CA  CYS A   2       1.822   3.786   0.000  1.00 69.50           C
ATOM     10  C   CYS A   2       0.775   4.783   0.000  1.00 69.00           C
ATOM     11  O   CYS A   2       0.969   5.980   0.000  1.00 68.60           O
ATOM     12  CB  CYS A   2       2.678   4.057   1.241  1.00 68.20           C
ATOM     13  SG  CYS A   2       3.315   5.748   1.241  1.00 67.80           S
"""
)


def test_parse_pdb_residue_order_and_coordinates() -> None:
    residues = parse_pdb(PDB_SNIPPET)
    assert [r.type for r in residues] == ["M", "C"]
    assert residues[1].resseq == 2
    sg = [a for a in residues[1].atoms if a.name == "SG"][0]
    assert np.allclose(sg.xyz, (3.315, 5.748, 1.241))
    assert sg.plddt == pytest.approx(67.8)


def test_parse_pdb_raises_on_non_cys_position() -> None:
    residues = parse_pdb(PDB_SNIPPET)
    with pytest.raises(ValueError, match="not Cys"):
        extract_cys_structure_features(residues, positions=(1,))  # position 1 is M


def test_sasa_of_isolated_atom_is_sphere_area() -> None:
    coords = np.array([[0.0, 0.0, 0.0]])
    radii = np.array([1.7])
    sasa = sasa_shrake_rupley(coords, radii, probe=0.0, n_points=4096)
    expected = 4 * math.pi * 1.7**2
    assert sasa[0] == pytest.approx(expected, rel=0.02)


def test_sasa_probe_expands_radius() -> None:
    coords = np.array([[0.0, 0.0, 0.0]])
    radii = np.array([1.7])
    sasa = sasa_shrake_rupley(coords, radii, probe=1.4, n_points=4096)
    expected = 4 * math.pi * (1.7 + 1.4) ** 2
    assert sasa[0] == pytest.approx(expected, rel=0.02)


def test_sasa_two_atoms_occlude_each_other() -> None:
    coords = np.array([[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]])
    radii = np.array([1.7, 1.7])
    sasa = sasa_shrake_rupley(coords, radii, probe=0.0, n_points=4096)
    # distance 2.5 < 2*1.7: partial occlusion
    assert 0.0 < sasa[0] < 4 * math.pi * 1.7**2


def test_coulomb_potential_basic_physics() -> None:
    estimator = CoulombPotentialEstimator()
    # positive charge at distance 5 from the probe point
    potential = estimator.potential_at(
        probe=(0.0, 0.0, 0.0), charges=((1.0, (5.0, 0.0, 0.0)),)
    )
    assert potential > 0
    # negative charge flips sign
    potential_neg = estimator.potential_at(
        probe=(0.0, 0.0, 0.0), charges=((-1.0, (5.0, 0.0, 0.0)),)
    )
    assert potential_neg < 0
    # closer charge dominates
    potential_close = estimator.potential_at(
        probe=(0.0, 0.0, 0.0), charges=((1.0, (2.0, 0.0, 0.0)),)
    )
    assert potential_close > potential


def test_extract_cys_features_shape_and_keys() -> None:
    residues = parse_pdb(PDB_SNIPPET)
    features = cys_structure_features(residues, positions=(2,))
    required = {
        "plddt",
        "rsa_relative",
        "cys_count_8a",
        "nearest_sg_distance",
        "positive_residue_count_6a",
        "coulomb_potential_sg",
        "contact_number_10a",
    }
    assert required <= set(features[2])
    assert 0.0 <= features[2]["rsa_relative"] <= 2.0
    # per-residue pLDDT is read from the CA atom
    assert features[2]["plddt"] == pytest.approx(69.5)
