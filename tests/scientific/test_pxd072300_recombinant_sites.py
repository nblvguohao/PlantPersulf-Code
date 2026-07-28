"""RED: PXD072300 (recombinant rice protein persulfidation) control parser.

The parser reads PEAKS protein-summary CSVs deposited on PRIDE and extracts
coordinate-verified persulfidation control sites for purified recombinant
proteins. These are SAME-PAPER, SAME-LAB orthogonal in-vitro validations —
not an independent unit for Gate 2 conditions 1 or 3 — used only as an
internal-consistency probe for the rice cross-species transfer track.

Reverse-engineered format quirks under test:
* protein blocks are located by AC match, not row order (contaminants can
  outrank the real target protein);
* the ``Positions`` field gives the position of the residue BEFORE the
  peptide (peptide start = field value + 1);
* only rows with a ``Carbamidomethyl[C]`` entry in ``Modification`` mark a
  persulfidation-derived site; duplicate unmodified rows are ignored.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.proteomics.pxd072300_recombinant_sites import (
    RECOMBINANT_PROTEINS,
    parse_recombinant_control_sites,
)

PROTEIN_HEADER = [
    "ID", "AC", "Score", "Q-Value", "Coverage", "No.Peptide",
    "No.Sameset", "No.Subset", "Have_Distinct_Pep", "Description",
]
PEPTIDE_HEADER = [
    "", "", "ID", "Sequence", "Calc.MH+", "Mass_Shift(Exp.-Calc.)",
    "Raw_Score", "Final_Score", "Modification", "Specificity", "Proteins",
    "Positions", "Label", "Target/Decoy", "Miss.Clv.Sites",
    "Avg.Frag.Mass.Shift", "File_Name", "Charge", "Spec_Num",
]


def _protein_row(protein_id: str, ac: str, desc: str) -> list[str]:
    return [protein_id, ac, "1.5", "0", "0.9", "10", "0", "0", "1", desc]


def _peptide_row(
    pep_id: str,
    seq: str,
    pos_before: int,
    before_res: str,
    after_res: str,
    modification: str = "",
) -> list[str]:
    return [
        "", "", pep_id, seq, "1000.0", "0.001", "20.0", "1e-5",
        modification, "3", "FAKE/", f"{pos_before},{before_res},{after_res}/",
        "1|", "target", "0", "0.05", "fake.dta", "2", "1",
    ]


def _write_csv(path: Path, rows: list[list[str]]) -> None:
    import csv

    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerows(rows)


MINI_PROTEOME = {
    # 20 aa, Cys at positions 7 and 15 (1-indexed)
    "Q6ZJ08": "MPGFHTCVGSGGERCLLPQR",
}

MDHAR4_ONLY = {"MDHAR4": {"internal_ac": "MDAR4", "uniprot_accession": "Q6ZJ08"}}


def _mdhar4_csv(tmp_path: Path, peptide_rows: list[list[str]]) -> Path:
    path = tmp_path / "Persulfidation_MDHAR4.csv"
    rows = [
        PROTEIN_HEADER,
        PEPTIDE_HEADER,
        _protein_row("1", "MDAR4", ">MDAR4"),
        *peptide_rows,
        _protein_row("2", "CON_KRT1", "contaminant"),
        _peptide_row("99", "AAAAAAA", 1, "M", "A"),
    ]
    _write_csv(path, rows)
    return path


def test_registry_lists_all_ten_recombinant_proteins() -> None:
    assert len(RECOMBINANT_PROTEINS) == 10
    assert "MDHAR4" in RECOMBINANT_PROTEINS
    assert "TKT" in RECOMBINANT_PROTEINS


def test_verified_carbamidomethyl_site_is_extracted(tmp_path: Path) -> None:
    # Peptide "PGFHTCVGSGGER" starts right after position 1 (M), i.e. at
    # position 2; sequence[1:14] == "PGFHTCVGSGGER". Positions field must
    # give the PRECEDING residue's position: 1 (M).
    seq = MINI_PROTEOME["Q6ZJ08"]
    pep = "PGFHTCVGSGGER"
    assert seq[1 : 1 + len(pep)] == pep
    _mdhar4_csv(
        tmp_path,
        [_peptide_row("1", pep, 1, "M", "L", modification="6,Carbamidomethyl[C]S;")],
    )
    sites = parse_recombinant_control_sites(
        tmp_path, MINI_PROTEOME, proteins=MDHAR4_ONLY
    )
    assert len(sites) == 1
    assert sites[0].gene == "MDHAR4"
    assert sites[0].uniprot_accession == "Q6ZJ08"
    assert sites[0].cys_position == 7  # position 2 (start) + index 6 - 1
    assert sites[0].status == "mapped"


def test_unmodified_duplicate_row_is_ignored(tmp_path: Path) -> None:
    pep = "PGFHTCVGSGGER"
    _mdhar4_csv(
        tmp_path,
        [
            _peptide_row("1", pep, 1, "M", "L", modification=""),
            _peptide_row("2", pep, 1, "M", "L", modification=""),
        ],
    )
    sites = parse_recombinant_control_sites(
        tmp_path, MINI_PROTEOME, proteins=MDHAR4_ONLY
    )
    assert len(sites) == 1
    assert sites[0].status == "no_sites_detected"
    assert sites[0].cys_position == 0


def test_coordinate_mismatch_row_is_ignored(tmp_path: Path) -> None:
    # Claimed peptide does not actually occur at the given position.
    _mdhar4_csv(
        tmp_path,
        [
            _peptide_row(
                "1", "ZZZZZZZZZZZZZ", 1, "M", "L",
                modification="6,Carbamidomethyl[C]S;",
            ),
        ],
    )
    sites = parse_recombinant_control_sites(
        tmp_path, MINI_PROTEOME, proteins=MDHAR4_ONLY
    )
    assert len(sites) == 1
    assert sites[0].status == "no_sites_detected"


def test_unresolved_accession_is_unmappable(tmp_path: Path) -> None:
    proteins = {"TKT": {"internal_ac": "transketolase", "uniprot_accession": ""}}
    sites = parse_recombinant_control_sites(tmp_path, {}, proteins=proteins)
    assert len(sites) == 1
    assert sites[0].status == "unmappable"
    assert sites[0].gene == "TKT"


def test_missing_csv_file_is_unmappable(tmp_path: Path) -> None:
    sites = parse_recombinant_control_sites(
        tmp_path, MINI_PROTEOME, proteins=MDHAR4_ONLY
    )
    assert len(sites) == 1
    assert sites[0].status == "unmappable"


def test_protein_block_located_by_ac_not_row_order(tmp_path: Path) -> None:
    """The target protein is the SECOND block (a contaminant ranks first)."""
    pep = "PGFHTCVGSGGER"
    path = tmp_path / "Persulfidation_MDHAR4.csv"
    rows = [
        PROTEIN_HEADER,
        PEPTIDE_HEADER,
        _protein_row("1", "CON_KRT1", "contaminant, ranks first"),
        _peptide_row("50", "AAAAAAA", 1, "M", "A"),
        _protein_row("2", "MDAR4", ">MDAR4"),
        _peptide_row("1", pep, 1, "M", "L", modification="6,Carbamidomethyl[C]S;"),
    ]
    _write_csv(path, rows)
    sites = parse_recombinant_control_sites(
        tmp_path, MINI_PROTEOME, proteins=MDHAR4_ONLY
    )
    assert len(sites) == 1
    assert sites[0].status == "mapped"
    assert sites[0].cys_position == 7


# --- real-data test: runs against the registered PXD072300 downloads ------

REAL_CSV_DIR = Path("data/raw/supplements/PXD072300")
REAL_PROTEOME = Path("data/raw/references/rice_proteome_v1/uniprot_rice_v1.fasta")


@pytest.mark.skipif(
    not REAL_CSV_DIR.is_dir() or not REAL_PROTEOME.is_file(),
    reason="registered PXD072300 downloads or rice proteome not present",
)
def test_real_pxd072300_sites_verify_against_rice_proteome() -> None:
    from plantpersulf.features.sequence import _load_proteome

    proteome = _load_proteome(REAL_PROTEOME)
    sites = parse_recombinant_control_sites(REAL_CSV_DIR, proteome)

    mapped = [s for s in sites if s.status == "mapped"]
    unmappable = [s for s in sites if s.status == "unmappable"]

    # 23 verified sites across 8 mapped proteins (per manual audit); TKT is
    # unmappable (best candidate covers only 60/84 unique peptides).
    assert len(mapped) == 23
    assert any(s.gene == "TKT" for s in unmappable)

    for site in mapped:
        seq = proteome[site.uniprot_accession]
        assert seq[site.cys_position - 1] == "C"
