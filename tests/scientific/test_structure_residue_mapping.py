"""RED (Task 7): AlphaFold PDB residue mapping + explicit missingness mask.

Parses the well-documented, standard fixed-column PDB ``ATOM`` record format
(wwPDB format v3.3, sect. 9: columns 7-11 serial, 13-16 atom name, 18-20
resName, 22 chainID, 23-26 resSeq, 31-38/39-46/47-54 x/y/z, 61-66
tempFactor). AlphaFold stores per-residue pLDDT in the tempFactor
(B-factor) column, replicated across every atom of a residue.

The fixture below is a hand-built, *clearly synthetic* five-residue chain
(reusing the toy "MVCGK" sequence already used in
``tests/scientific/test_esm2_features.py`` for the same reason: it exists
only to exercise code correctness, not to claim anything about a real
protein) with Cα atoms spaced at the textbook ~3.8 A trans-peptide bond
distance. No network access or real AlphaFold download is used or required
by this test.

Expected RED: ``plantpersulf.features.structure`` does not exist yet.
"""

from __future__ import annotations

from plantpersulf.features.structure import (  # RED: module missing
    LOW_PLDDT_THRESHOLD,
    STRUCTURE_FEATURE_FIELDS,
    CysStructureFeature,
    extract_cys_structure_features,
    missing_structure_feature,
    parse_alphafold_pdb,
    read_structure_features_tsv,
    write_structure_features_tsv,
)


def _atom_line(
    serial: int,
    atom_name: str,
    res_name: str,
    chain: str,
    res_seq: int,
    x: float,
    y: float,
    z: float,
    b_factor: float,
    element: str,
) -> str:
    """Build one fixed-column PDB ATOM record per the wwPDB v3.3 spec."""
    line = [" "] * 80
    line[0:6] = list("ATOM  ")
    line[6:11] = list(f"{serial:>5}")
    line[12:16] = list(f"{atom_name:<4}")
    line[17:20] = list(f"{res_name:>3}")
    line[21] = chain
    line[22:26] = list(f"{res_seq:>4}")
    line[30:38] = list(f"{x:8.3f}")
    line[38:46] = list(f"{y:8.3f}")
    line[46:54] = list(f"{z:8.3f}")
    line[54:60] = list(f"{1.0:6.2f}")
    line[60:66] = list(f"{b_factor:6.2f}")
    line[76:78] = list(f"{element:>2}")
    return "".join(line)


# Synthetic "MVCGK" chain: one CA per residue, spaced 3.8 A apart on the
# x-axis (the textbook Ca-Ca trans-peptide-bond spacing), each carrying a
# distinct per-residue pLDDT in the B-factor column. This is a
# code-correctness fixture, not a claim about a real structure.
_RESIDUES = [
    (1, "MET", 0.0, 90.0),
    (2, "VAL", 3.8, 85.0),
    (3, "CYS", 7.6, 40.0),  # deliberately low-confidence (<50)
    (4, "GLY", 11.4, 88.0),
    (5, "LYS", 15.2, 91.0),
]
_PDB_TEXT = "\n".join(
    _atom_line(idx, "CA", res_name, "A", seq, x, 0.0, 0.0, plddt, "C")
    for idx, (seq, res_name, x, plddt) in enumerate(_RESIDUES, start=1)
)


def test_parses_per_residue_plddt_from_bfactor_column() -> None:
    residues = parse_alphafold_pdb(_PDB_TEXT)

    by_seq = {r.res_seq: r for r in residues}
    assert by_seq[1].plddt == 90.0
    assert by_seq[3].plddt == 40.0
    assert by_seq[3].res_name == "CYS"


def test_low_plddt_threshold_matches_alphafold_convention() -> None:
    assert LOW_PLDDT_THRESHOLD == 50.0


def test_extract_features_flags_low_confidence_residue() -> None:
    features = extract_cys_structure_features(_PDB_TEXT, "TESTP", [3])

    assert len(features) == 1
    feature = features[0]
    assert feature.has_structure is True
    assert feature.plddt == 40.0
    assert feature.low_plddt is True


def test_contact_number_proxy_counts_nearby_residues_not_self() -> None:
    features = extract_cys_structure_features(
        _PDB_TEXT, "TESTP", [3], contact_radius_angstrom=10.0
    )

    feature = features[0]
    # residues 1,2,4,5 are all within 10 A of residue 3's Ca (7.6,3.8,3.8,7.6)
    assert feature.contact_number_proxy == 4.0


def test_contact_number_proxy_shrinks_with_smaller_radius() -> None:
    features = extract_cys_structure_features(
        _PDB_TEXT, "TESTP", [3], contact_radius_angstrom=5.0
    )

    # only residues 2 and 4 (3.8 A away) fall within a 5 A shell
    assert features[0].contact_number_proxy == 2.0


def test_cys_position_not_matching_cys_residue_is_missing_not_a_crash() -> None:
    # position 1 is MET, not CYS: a coordinate-mapping failure must yield an
    # explicit missing record, never raise and never fabricate a value.
    features = extract_cys_structure_features(_PDB_TEXT, "TESTP", [1, 3])

    by_pos = {f.cys_position: f for f in features}
    assert by_pos[1].has_structure is False
    assert by_pos[1].plddt is None
    assert by_pos[1].contact_number_proxy is None
    assert by_pos[3].has_structure is True


def test_missing_structure_feature_uses_explicit_mask_not_average() -> None:
    feature = missing_structure_feature("NOSTRUCT", 7)

    assert feature.has_structure is False
    assert feature.plddt is None
    assert feature.low_plddt is None
    assert feature.contact_number_proxy is None
    assert isinstance(feature, CysStructureFeature)


def test_structure_feature_tsv_round_trip_preserves_missing_mask(tmp_path) -> None:  # type: ignore[no-untyped-def]
    rows = extract_cys_structure_features(_PDB_TEXT, "TESTP", [3]) + [
        missing_structure_feature("TESTP", 999)
    ]
    path = tmp_path / "structure_features.tsv"
    write_structure_features_tsv(rows, path)

    header = path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert tuple(header) == STRUCTURE_FEATURE_FIELDS

    round_tripped = read_structure_features_tsv(path)
    by_pos = {r.cys_position: r for r in round_tripped}
    assert by_pos[3].has_structure is True
    assert by_pos[3].plddt == 40.0
    assert by_pos[999].has_structure is False
    assert by_pos[999].plddt is None
    assert by_pos[999].contact_number_proxy is None
