"""RED: PXD063170 (Magnaporthe oryzae S-sulfhydration) site-table parser.

The parser consumes the TSV export of the paper's Supplementary Data 1 and
applies the same integrity rules as the Arabidopsis benchmark: every reported
site must be a confidently localized cysteine whose coordinate matches the
reference proteome sequence; anything else is dropped and counted, never
silently repaired.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.proteomics.pxd063170_sites import (
    load_ensembl_fungi_proteome,
    parse_pxd063170_sites,
)

HEADER = (
    "Protein accession\tPosition\tAmino acid\tCSE_OE/WT Ratio\t"
    "Localization probability\tModified sequence\n"
)


def _write_tsv(path: Path, rows: list[list[str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(HEADER)
        for row in rows:
            handle.write("\t".join(row) + "\n")


PROTEOME = {
    "MGG_00001T0": "MACDEFGHIKLMNCQRSTVWY",  # Cys at 3 and 14
    "MGG_00002T0": "MADEFGCCLM",  # Cys at 7, 8
}


def test_valid_localized_site_is_accepted(tmp_path: Path) -> None:
    tsv = tmp_path / "sites.tsv"
    _write_tsv(tsv, [["MGG_00001T0", "3", "C", "1.5", "0.99", "AAC(1)DEF"]])
    table = parse_pxd063170_sites(tsv, PROTEOME, min_localization=0.75)
    assert [(s.protein_accession, s.cys_position) for s in table.sites] == [
        ("MGG_00001T0", 3)
    ]
    assert table.dropped_low_localization == 0
    assert table.dropped_coordinate_mismatch == 0
    assert table.dropped_missing_accession == 0


def test_low_localization_site_is_dropped(tmp_path: Path) -> None:
    tsv = tmp_path / "sites.tsv"
    _write_tsv(
        tsv,
        [
            ["MGG_00001T0", "3", "C", "1.5", "0.30", "AAC(1)DEF"],
            ["MGG_00001T0", "14", "C", "1.5", "0.80", "NC(1)QRST"],
        ],
    )
    table = parse_pxd063170_sites(tsv, PROTEOME, min_localization=0.75)
    assert [(s.protein_accession, s.cys_position) for s in table.sites] == [
        ("MGG_00001T0", 14)
    ]
    assert table.dropped_low_localization == 1


def test_coordinate_mismatch_is_dropped_not_repaired(tmp_path: Path) -> None:
    tsv = tmp_path / "sites.tsv"
    _write_tsv(tsv, [["MGG_00001T0", "2", "C", "1.5", "0.99", "ACD(1)EFG"]])
    table = parse_pxd063170_sites(tsv, PROTEOME)
    assert table.sites == ()
    assert table.dropped_coordinate_mismatch == 1


def test_missing_accession_is_counted(tmp_path: Path) -> None:
    tsv = tmp_path / "sites.tsv"
    _write_tsv(tsv, [["MGG_99999T0", "5", "C", "1.5", "0.99", "AAC(1)DEF"]])
    table = parse_pxd063170_sites(tsv, PROTEOME)
    assert table.sites == ()
    assert table.dropped_missing_accession == 1


def test_non_cysteine_row_is_rejected(tmp_path: Path) -> None:
    tsv = tmp_path / "sites.tsv"
    _write_tsv(tsv, [["MGG_00001T0", "4", "D", "1.5", "0.99", "AACD(1)EFG"]])
    with pytest.raises(RuntimeError, match="non-cysteine"):
        parse_pxd063170_sites(tsv, PROTEOME)


def test_invalid_columns_are_rejected(tmp_path: Path) -> None:
    tsv = tmp_path / "sites.tsv"
    tsv.write_text("a\tb\n1\t2\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="invalid columns"):
        parse_pxd063170_sites(tsv, PROTEOME)


def test_duplicate_sites_collapse_to_best_localization(tmp_path: Path) -> None:
    tsv = tmp_path / "sites.tsv"
    _write_tsv(
        tsv,
        [
            ["MGG_00002T0", "7", "C", "1.5", "0.80", "FGC(1)CLM"],
            ["MGG_00002T0", "7", "C", "1.5", "0.95", "FGC(1)CLM"],
        ],
    )
    table = parse_pxd063170_sites(tsv, PROTEOME)
    assert len(table.sites) == 1
    assert table.sites[0].localization_probability == pytest.approx(0.95)


def test_load_ensembl_fungi_proteome(tmp_path: Path) -> None:
    fasta = tmp_path / "pep.fa"
    fasta.write_text(
        ">MGG_00001T0 pep chromosome:MG8:1:100:200:1 gene:MGG_00001\n"
        "MACDEFGHIK\nLMNCQRSTVWY\n"
        ">MGG_00002T0 pep chromosome:MG8:1:300:400:-1 gene:MGG_00002\n"
        "MADEFGCCLM\n",
        encoding="utf-8",
    )
    proteome = load_ensembl_fungi_proteome(fasta)
    assert proteome == PROTEOME


# --- real-data test: runs against the registered downloads -----------------

REAL_TSV = Path("data/raw/PXD063170/PXD063170_sites_moesm3.tsv")
REAL_PROTEOME = Path("data/raw/PXD063170/Magnaporthe_oryzae.MG8.pep.all.fa")


@pytest.mark.skipif(
    not REAL_TSV.is_file() or not REAL_PROTEOME.is_file(),
    reason="registered PXD063170 downloads not present",
)
def test_real_pxd063170_sites_verify_against_proteome() -> None:
    proteome = load_ensembl_fungi_proteome(REAL_PROTEOME)
    table = parse_pxd063170_sites(REAL_TSV, proteome, min_localization=0.75)
    # The table ships 1482 unique site rows; confidently localized ones are
    # the large majority, and coordinate mismatches must be rare and counted.
    assert len(table.sites) > 1000
    assert table.dropped_missing_accession == 0
    assert table.dropped_coordinate_mismatch <= 20
