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
    PXD063170SiteTable,
    build_cross_species_eval_rows,
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

REAL_TSV = Path("data/raw/supplements/PXD063170/PXD063170_sites_moesm3.tsv")
REAL_PROTEOME = Path("data/raw/supplements/PXD063170/Magnaporthe_oryzae.MG8.pep.all.fa")


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


# --- cross-species evaluation-row construction (PU semantics) ---------------


def _parsed_table(tmp_path: Path) -> PXD063170SiteTable:
    tsv = tmp_path / "sites.tsv"
    _write_tsv(
        tsv,
        [
            ["MGG_00001T0", "3", "C", "1.5", "0.99", "AAC(1)DEF"],
            ["MGG_00002T0", "7", "C", "1.5", "0.90", "FGC(1)CLM"],
        ],
    )
    return parse_pxd063170_sites(tsv, PROTEOME)


def test_eval_rows_mark_sites_positive_and_other_cysteines_unlabeled(
    tmp_path: Path,
) -> None:
    table = _parsed_table(tmp_path)
    rows = build_cross_species_eval_rows(
        table, PROTEOME, study_accession="PXD063170", source_sha256="deadbeef"
    )
    by_key = {(r["protein_accession"], r["cys_position_in_protein"]): r for r in rows}
    pos = by_key[("MGG_00001T0", "3")]
    assert pos["label"] == "positive"
    assert pos["study_accession"] == "PXD063170"
    assert pos["evidence_level"] == "site_ms"
    assert pos["source_sha256"] == "deadbeef"
    # MGG_00001T0 Cys 14 and MGG_00002T0 Cys 8 are cysteines but not sites.
    for key in (("MGG_00001T0", "14"), ("MGG_00002T0", "8")):
        row = by_key[key]
        assert row["label"] == "unlabeled"
        assert row["study_accession"] == ""
        assert row["evidence_level"] == ""
        assert row["source_sha256"] == ""


def test_eval_rows_cover_every_proteome_cysteine_exactly_once(
    tmp_path: Path,
) -> None:
    table = _parsed_table(tmp_path)
    rows = build_cross_species_eval_rows(table, PROTEOME)
    n_cys = sum(seq.count("C") for seq in PROTEOME.values())
    assert len(rows) == n_cys
    keys = [(r["protein_accession"], r["cys_position_in_protein"]) for r in rows]
    assert len(set(keys)) == len(keys)
    # deterministic order: sorted by (accession, numeric position)
    assert keys == sorted(keys, key=lambda k: (k[0], int(k[1])))


def test_eval_rows_never_emit_non_cysteine_positions(tmp_path: Path) -> None:
    table = _parsed_table(tmp_path)
    rows = build_cross_species_eval_rows(table, PROTEOME)
    for row in rows:
        seq = PROTEOME[row["protein_accession"]]
        assert seq[int(row["cys_position_in_protein"]) - 1] == "C"


@pytest.mark.skipif(
    not REAL_TSV.is_file() or not REAL_PROTEOME.is_file(),
    reason="registered PXD063170 downloads not present",
)
def test_real_eval_rows_background_matches_proteome_cysteine_count() -> None:
    proteome = load_ensembl_fungi_proteome(REAL_PROTEOME)
    table = parse_pxd063170_sites(REAL_TSV, proteome, min_localization=0.75)
    rows = build_cross_species_eval_rows(table, proteome)
    n_cys = sum(seq.count("C") for seq in proteome.values())
    assert len(rows) == n_cys
    n_pos = sum(1 for r in rows if r["label"] == "positive")
    assert n_pos == len(table.sites)
