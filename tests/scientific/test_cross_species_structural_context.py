"""RED: structural context of persulfidation targeting within proteins,
compared across species (mechanistic companion to the family-level
conservation test in test_cross_species_conservation.py).

PU semantics: non-persulfidated cysteines in a persulfidated protein are
labelled "unlabeled" (not detected), never "negative" (not confirmed
absent) — consistent with the project-wide convention in
evaluation/metrics.py.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.evaluation.cross_species_structural_context import (
    StructuralContextRow,
    StructuralContextRowSasa,
    build_structural_context_rows,
    build_structural_context_rows_sasa,
    filter_rows_by_plddt,
    filter_sasa_rows_by_plddt,
    mean_difference_metric,
    plddt_by_key,
    plddt_context_test,
    sasa_structural_context_test,
    structural_context_test,
)
from plantpersulf.features.structure import CysStructureFeature

# ---------------------------------------------------------------------------
# mean_difference_metric
# ---------------------------------------------------------------------------


def test_mean_difference_metric_basic() -> None:
    scored = [
        (1.0, "positive"),
        (3.0, "positive"),
        (5.0, "unlabeled"),
        (7.0, "unlabeled"),
    ]
    # mean(positive)=2.0, mean(unlabeled)=6.0 -> diff = -4.0
    assert mean_difference_metric(scored) == pytest.approx(-4.0)


def test_mean_difference_metric_none_when_group_empty() -> None:
    assert mean_difference_metric([(1.0, "positive")]) is None
    assert mean_difference_metric([(1.0, "unlabeled")]) is None
    assert mean_difference_metric([]) is None


# ---------------------------------------------------------------------------
# build_structural_context_rows: synthetic PDB text
# ---------------------------------------------------------------------------

_PDB_HEADER = ""


def _pdb_ca_line(res_seq: int, res_name: str, x: float, plddt: float) -> str:
    # Fixed-column PDB ATOM record, CA atom, chain A.
    return (
        f"ATOM  {res_seq:5d}  CA  {res_name:>3s} A{res_seq:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}{1.0:6.2f}{plddt:6.2f}"
    )


def _make_synthetic_pdb(sequence: str, plddt_by_pos: dict[int, float]) -> str:
    three_letter = {"C": "CYS", "M": "MET", "A": "ALA"}
    lines = []
    for i, aa in enumerate(sequence, start=1):
        lines.append(
            _pdb_ca_line(
                i, three_letter.get(aa, "ALA"), float(i), plddt_by_pos.get(i, 90.0)
            )
        )
    return "\n".join(lines) + "\n"


def test_build_rows_labels_positive_and_unlabeled_cys(tmp_path: Path) -> None:
    # Protein "P1": MACMAC -> Cys at 3 and 6. Position 3 is persulfidated.
    seq = "MACMAC"
    pdb_text = _make_synthetic_pdb(seq, {})
    pdb_path = tmp_path / "P1.pdb"
    pdb_path.write_text(pdb_text, encoding="utf-8")

    rows = build_structural_context_rows(
        persulfidated_keys={("P1", 3)},
        proteome={"P1": seq},
        structure_dir=tmp_path,
        accession_to_structure_file={"P1": pdb_path},
    )
    by_pos = {r.cys_position: r.label for r in rows}
    assert by_pos == {3: "positive", 6: "unlabeled"}


def test_build_rows_skips_proteins_without_structure_file(tmp_path: Path) -> None:
    rows = build_structural_context_rows(
        persulfidated_keys={("P1", 3)},
        proteome={"P1": "MACMAC"},
        structure_dir=tmp_path,
        accession_to_structure_file={},  # no structure registered for P1
    )
    assert rows == []


def test_build_rows_only_covers_proteins_with_a_positive(tmp_path: Path) -> None:
    seq = "MACMAC"
    pdb_path = tmp_path / "P2.pdb"
    pdb_path.write_text(_make_synthetic_pdb(seq, {}), encoding="utf-8")

    # P2 has a structure but no persulfidated Cys at all -> not in
    # persulfidated_keys for ANY position -> contributes no rows.
    rows = build_structural_context_rows(
        persulfidated_keys={("P1", 3)},  # P1, not P2
        proteome={"P2": seq},
        structure_dir=tmp_path,
        accession_to_structure_file={"P2": pdb_path},
    )
    assert rows == []


def test_build_rows_are_sorted_by_accession_regardless_of_set_hash_order(
    tmp_path: Path,
) -> None:
    # Row order must depend only on sorted(accession), never on Python's
    # per-process-randomized set-iteration order — otherwise the
    # downstream permutation test's exact p-value would not be
    # reproducible across process restarts at a fixed seed (2026-07-23
    # self-review finding).
    seq = "MACMAC"
    paths = {}
    for name in ("Zebra", "Apple", "Mango"):
        p = tmp_path / f"{name}.pdb"
        p.write_text(_make_synthetic_pdb(seq, {}), encoding="utf-8")
        paths[name] = p

    rows = build_structural_context_rows(
        persulfidated_keys={("Zebra", 3), ("Apple", 3), ("Mango", 3)},
        proteome={"Zebra": seq, "Apple": seq, "Mango": seq},
        structure_dir=tmp_path,
        accession_to_structure_file=paths,
    )
    accessions_in_order = [r.protein_accession for r in rows]
    assert accessions_in_order == sorted(accessions_in_order)


# ---------------------------------------------------------------------------
# structural_context_test
# ---------------------------------------------------------------------------


def _row(accession: str, position: int, label: str, contact: float):
    from plantpersulf.evaluation.cross_species_structural_context import (
        StructuralContextRow,
    )

    feature = CysStructureFeature(
        protein_accession=accession,
        cys_position=position,
        has_structure=True,
        plddt=90.0,
        low_plddt=False,
        contact_number_proxy=contact,
    )
    return StructuralContextRow(
        protein_accession=accession, cys_position=position, label=label, feature=feature
    )


def test_structural_context_test_reports_signed_direction() -> None:
    # Positive Cys have systematically LOWER contact number (more exposed).
    rows = [
        _row("P1", 1, "positive", 2.0),
        _row("P1", 2, "unlabeled", 10.0),
        _row("P2", 1, "positive", 3.0),
        _row("P2", 2, "unlabeled", 11.0),
    ]
    result = structural_context_test("TestSpecies", rows, n_perm=200, seed=0)
    assert result.n_positive_cys == 2
    assert result.n_unlabeled_cys == 2
    assert result.mean_contact_positive == pytest.approx(2.5)
    assert result.mean_contact_unlabeled == pytest.approx(10.5)
    assert result.mean_diff == pytest.approx(-8.0)
    assert result.permutation.p_value < 0.5


def test_structural_context_test_no_signal_gives_weak_p_value() -> None:
    rows = [
        _row("P1", 1, "positive", 5.0),
        _row("P1", 2, "unlabeled", 5.1),
        _row("P2", 1, "positive", 4.9),
        _row("P2", 2, "unlabeled", 5.0),
    ]
    result = structural_context_test("TestSpecies", rows, n_perm=200, seed=0)
    assert result.permutation.p_value > 0.1


# ---------------------------------------------------------------------------
# Real SASA pathway (2026-07-23 self-review cross-check): independent of
# the contact_number_proxy pathway above, using features/sasa.py's genuine
# Shrake-Rupley calculation on full heavy-atom geometry.
# ---------------------------------------------------------------------------


def _full_atom_pdb_line(
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
    # Same verified column layout as tests/unit/test_sasa.py.
    return (
        f"ATOM  {serial:5d} "
        f" {atom_name:<3s}"
        f" {res_name:>3s} {chain}"
        f"{res_seq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{50.0:6.2f}"
        f"{'':10s}{element:>2s}"
    )


def _make_full_atom_cys_pdb(residues: dict[int, str]) -> str:
    """Build a minimal multi-residue PDB with full backbone + side-chain
    atoms for CYS/ALA, spaced 3.8 Å apart along x (typical CA-CA spacing),
    each residue's own atoms offset slightly so every residue gets a
    resolvable, non-degenerate local geometry."""
    three_letter = {"C": "CYS", "A": "ALA"}
    lines = []
    serial = 1
    for res_seq, aa in sorted(residues.items()):
        base_x = res_seq * 3.8
        res_name = three_letter[aa]
        atoms = [
            ("N", -0.5, 0.0, "N"),
            ("CA", 0.0, 0.0, "C"),
            ("C", 0.5, 0.3, "C"),
            ("O", 0.5, 1.3, "O"),
            ("CB", 0.0, -1.5, "C"),
        ]
        if aa == "C":
            atoms.append(("SG", 0.0, -3.3, "S"))
        for atom_name, dy, dz, element in atoms:
            lines.append(
                _full_atom_pdb_line(
                    serial,
                    atom_name,
                    res_name,
                    "A",
                    res_seq,
                    base_x,
                    dy,
                    dz,
                    element,
                )
            )
            serial += 1
    return "\n".join(lines) + "\n"


def test_build_rows_sasa_labels_positive_and_unlabeled(tmp_path: Path) -> None:
    # Two Cys (positions 3, 6), one Ala spacer between them.
    pdb_text = _make_full_atom_cys_pdb({3: "C", 4: "A", 6: "C"})
    pdb_path = tmp_path / "P1.pdb"
    pdb_path.write_text(pdb_text, encoding="utf-8")
    seq = "MMCAMC"  # Cys at 1-indexed positions 3 and 6

    rows = build_structural_context_rows_sasa(
        persulfidated_keys={("P1", 3)},
        proteome={"P1": seq},
        accession_to_structure_file={"P1": pdb_path},
        n_sphere_points=50,
    )
    by_pos = {r.cys_position: r.label for r in rows}
    assert by_pos == {3: "positive", 6: "unlabeled"}
    for row in rows:
        assert row.residue_sasa > 0
        assert row.sg_sasa is not None
        assert row.sg_sasa > 0


def test_build_rows_sasa_skips_proteins_without_structure(tmp_path: Path) -> None:
    rows = build_structural_context_rows_sasa(
        persulfidated_keys={("P1", 3)},
        proteome={"P1": "MMCAMC"},
        accession_to_structure_file={},
    )
    assert rows == []


def test_build_rows_sasa_are_sorted_by_accession(tmp_path: Path) -> None:
    seq = "MMCAMC"
    paths = {}
    for name in ("Zebra", "Apple", "Mango"):
        p = tmp_path / f"{name}.pdb"
        p.write_text(_make_full_atom_cys_pdb({3: "C", 6: "C"}), encoding="utf-8")
        paths[name] = p

    rows = build_structural_context_rows_sasa(
        persulfidated_keys={("Zebra", 3), ("Apple", 3), ("Mango", 3)},
        proteome={"Zebra": seq, "Apple": seq, "Mango": seq},
        accession_to_structure_file=paths,
        n_sphere_points=50,
    )
    accessions_in_order = [r.protein_accession for r in rows]
    assert accessions_in_order == sorted(accessions_in_order)


def test_sasa_structural_context_test_residue_metric() -> None:
    from plantpersulf.evaluation.cross_species_structural_context import (
        StructuralContextRowSasa,
    )

    rows = [
        StructuralContextRowSasa("P1", 1, "positive", residue_sasa=10.0, sg_sasa=2.0),
        StructuralContextRowSasa("P1", 2, "unlabeled", residue_sasa=50.0, sg_sasa=20.0),
        StructuralContextRowSasa("P2", 1, "positive", residue_sasa=12.0, sg_sasa=3.0),
        StructuralContextRowSasa("P2", 2, "unlabeled", residue_sasa=48.0, sg_sasa=19.0),
    ]
    result = sasa_structural_context_test(
        "TestSpecies", rows, metric="residue_sasa", n_perm=200, seed=0
    )
    assert result.metric_name == "residue_sasa"
    assert result.mean_sasa_positive == pytest.approx(11.0)
    assert result.mean_sasa_unlabeled == pytest.approx(49.0)
    assert result.mean_diff < 0  # positives have LOWER SASA -> more buried
    assert result.permutation.p_value < 0.5


def test_sasa_structural_context_test_sg_metric_excludes_missing_sg() -> None:
    from plantpersulf.evaluation.cross_species_structural_context import (
        StructuralContextRowSasa,
    )

    rows = [
        StructuralContextRowSasa("P1", 1, "positive", residue_sasa=10.0, sg_sasa=None),
        StructuralContextRowSasa("P1", 2, "unlabeled", residue_sasa=50.0, sg_sasa=20.0),
    ]
    result = sasa_structural_context_test(
        "TestSpecies", rows, metric="sg_sasa", n_perm=50, seed=0
    )
    assert result.n_positive_cys == 0  # the only positive row had sg_sasa=None
    assert result.n_unlabeled_cys == 1


def test_sasa_structural_context_test_rejects_unknown_metric() -> None:
    with pytest.raises(ValueError, match="unknown metric"):
        sasa_structural_context_test("TestSpecies", [], metric="bogus")


# ---------------------------------------------------------------------------
# Real-data smoke test (skipped unless the bulk AlphaFold download + all
# three species' evidence tables are present)
# ---------------------------------------------------------------------------

ALPHAFOLD_REGISTRY = Path("data/registry/alphafold_structures.tsv")


@pytest.mark.skipif(
    not ALPHAFOLD_REGISTRY.is_file(),
    reason="bulk AlphaFold structure registry not present",
)
def test_real_alphafold_registry_has_bulk_download_scale() -> None:
    import csv

    with ALPHAFOLD_REGISTRY.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    # The bulk download targets ~2048 persulfidated proteins across three
    # species; this is a loose scale sanity check, not an exact count (the
    # download may still be in progress or include not_found/isoform gaps).
    assert len(rows) > 100


# ---------------------------------------------------------------------------
# Model-confidence control (2026-08-14)
#
# The accessibility comparison is only meaningful where AlphaFold is
# confident: a very-low-pLDDT region is predicted as an extended chain, which
# reads as "exposed" no matter what the real structure does. Mass spectrometry
# also favours flexible, protease-accessible regions, so "positives sit in
# low-confidence regions" is a live alternative explanation for any exposure
# result and has to be measured, not assumed away.
# ---------------------------------------------------------------------------


def _plddt_row(accession: str, position: int, label: str, contact: float, plddt: float):
    return StructuralContextRow(
        protein_accession=accession,
        cys_position=position,
        label=label,
        feature=CysStructureFeature(
            protein_accession=accession,
            cys_position=position,
            has_structure=True,
            plddt=plddt,
            low_plddt=plddt < 70.0,
            contact_number_proxy=contact,
        ),
    )


def test_plddt_context_test_detects_lower_confidence_at_positives() -> None:
    rows = [_plddt_row("P1", i, "positive", 10.0, 30.0) for i in range(1, 11)]
    rows += [_plddt_row("P1", i, "unlabeled", 10.0, 90.0) for i in range(11, 31)]
    result = plddt_context_test("Species", rows, n_perm=200, seed=1)
    assert result.n_positive_cys == 10
    assert result.n_unlabeled_cys == 20
    assert result.mean_plddt_positive == pytest.approx(30.0)
    assert result.mean_plddt_unlabeled == pytest.approx(90.0)
    assert result.mean_diff == pytest.approx(-60.0)
    assert result.permutation.p_value < 0.05


def test_plddt_context_test_is_null_when_confidence_matches() -> None:
    rows = [_plddt_row("P1", i, "positive", 10.0, 80.0) for i in range(1, 11)]
    rows += [_plddt_row("P1", i, "unlabeled", 10.0, 80.0) for i in range(11, 31)]
    result = plddt_context_test("Species", rows, n_perm=200, seed=1)
    assert result.mean_diff == pytest.approx(0.0)
    assert result.permutation.p_value == pytest.approx(1.0)


def test_filter_rows_by_plddt_keeps_only_confident_residues() -> None:
    rows = [
        _plddt_row("P1", 1, "positive", 10.0, 95.0),
        _plddt_row("P1", 2, "unlabeled", 12.0, 69.9),
        _plddt_row("P1", 3, "unlabeled", 14.0, 70.0),
    ]
    kept = filter_rows_by_plddt(rows, min_plddt=70.0)
    assert [r.cys_position for r in kept] == [1, 3]


def test_filter_rows_by_plddt_drops_rows_without_a_confidence_value() -> None:
    missing = StructuralContextRow(
        protein_accession="P1",
        cys_position=9,
        label="positive",
        feature=CysStructureFeature(
            protein_accession="P1",
            cys_position=9,
            has_structure=True,
            plddt=None,
            low_plddt=None,
            contact_number_proxy=10.0,
        ),
    )
    assert filter_rows_by_plddt([missing], min_plddt=70.0) == []


def test_plddt_by_key_maps_every_scored_cysteine() -> None:
    rows = [
        _plddt_row("P1", 1, "positive", 10.0, 95.0),
        _plddt_row("P2", 4, "unlabeled", 8.0, 40.0),
    ]
    mapping = plddt_by_key(rows)
    assert mapping == {("P1", 1): 95.0, ("P2", 4): 40.0}


def test_filter_sasa_rows_by_plddt_uses_the_shared_key_map() -> None:
    sasa_rows = [
        StructuralContextRowSasa(
            protein_accession="P1",
            cys_position=1,
            label="positive",
            residue_sasa=40.0,
            sg_sasa=20.0,
        ),
        StructuralContextRowSasa(
            protein_accession="P2",
            cys_position=4,
            label="unlabeled",
            residue_sasa=10.0,
            sg_sasa=5.0,
        ),
        # No pLDDT entry at all -> dropped, never assumed confident.
        StructuralContextRowSasa(
            protein_accession="P3",
            cys_position=7,
            label="unlabeled",
            residue_sasa=1.0,
            sg_sasa=0.5,
        ),
    ]
    mapping = {("P1", 1): 95.0, ("P2", 4): 40.0}
    kept = filter_sasa_rows_by_plddt(sasa_rows, mapping, min_plddt=70.0)
    assert [(r.protein_accession, r.cys_position) for r in kept] == [("P1", 1)]
