"""kiae271 (SlLCD1-OE vs WT tomato leaf) differential persulfidation parser.

Supplementary Dataset S1 is mistitled ("SlWRKY6 site analysis") but actually
holds the paper's full 119-row differential persulfidation screen. Every
reported site must be a coordinate-verified cysteine in the tomato reference
proteome, clear a localization-probability bar, and have nonzero intensity
in at least one condition; anything else is dropped and counted.
"""

from __future__ import annotations

from pathlib import Path

from plantpersulf.proteomics.kiae271_sites import (
    REGULATION_BOTH,
    REGULATION_LCD_GAIN,
    REGULATION_WT_ONLY,
    parse_kiae271_sites,
)

MINI_PROTEOME = {
    "A0A3Q7F586": "M" * 395 + "C" + "SSNMATISASAPFPTVTLDLTAQNPNAALPNYHQRINQANP",
    "Q6H3X6": "M" * 132 + "C" + "YYFGKDGKFVCEGESDEPKANMYPVM",
    "A0A3Q7XXXX": "MSEQNCE",  # Cys at position 6
}

HEADER = (
    "Proteins",
    "Positions within proteins",
    "Leading proteins",
    "Protein",
    "Fasta headers",
    "Localization prob",
    "PEP",
    "Score",
    "Score for localization",
    "Number of S(C)",
    "Amino acid",
    "Sequence window",
    "S(C) Probabilities",
    "Position in peptide",
    "Charge",
    "Intensity LCD_1",
    "Intensity LCD_2",
    "Intensity WT_1",
    "Intensity WT_2",
    "Intensity LCD",
    "Intensity WT",
)


def _write_xlsx(path: Path, rows: list[tuple]) -> None:
    import openpyxl

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "Supplementary Dataset S1"
    ws.append(("Title row (ignored)",))
    ws.append(HEADER)
    for row in rows:
        ws.append(row)
    wb.save(path)


def _row(
    proteins: str,
    positions: str,
    leading: str,
    amino: str = "C",
    loc_prob: float = 1.0,
    intensity_lcd: float = 0.0,
    intensity_wt: float = 0.0,
) -> tuple:
    return (
        proteins,
        positions,
        leading,
        leading,
        f"{leading} description",
        loc_prob,
        0.001,
        100,
        100,
        1,
        amino,
        "SEQWINDOW",
        "1",
        "PEP",
        2,
        intensity_lcd,
        0,
        intensity_wt,
        0,
        intensity_lcd,
        intensity_wt,
    )


def test_lcd_only_site_is_gain_regulation(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "396",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                intensity_lcd=221591.0,
                intensity_wt=0.0,
            )
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    assert table.total_verified_sites == 1
    site = table.sites[0]
    assert (site.protein_accession, site.cys_position) == ("A0A3Q7F586", 396)
    assert site.regulation == REGULATION_LCD_GAIN


def test_wt_only_and_both_regulation_classified(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7XXXX|X_SOLLC",
                "6",
                "tr|A0A3Q7XXXX|X_SOLLC",
                intensity_lcd=0.0,
                intensity_wt=500.0,
            ),
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "396",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                intensity_lcd=100.0,
                intensity_wt=200.0,
            ),
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    by_key = {(s.protein_accession, s.cys_position): s for s in table.sites}
    assert by_key[("A0A3Q7XXXX", 6)].regulation == REGULATION_WT_ONLY
    assert by_key[("A0A3Q7F586", 396)].regulation == REGULATION_BOTH


def test_multi_protein_position_alignment_uses_leading_protein_index(
    tmp_path: Path,
) -> None:
    """Positions within proteins is parallel to Proteins, not to Leading
    proteins' assumed index-0 — the leading protein here is 2nd in the list,
    with the 2nd position; a naive index-0 read would wrongly grab pos 133
    for Q6H3X6 (whose real, verified position is 133 in this fixture)."""
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC;tr|Q6H3X6|Q6H3X6_SOLLC",
                "396;133",
                "tr|Q6H3X6|Q6H3X6_SOLLC",
                intensity_lcd=1000.0,
            )
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    assert table.total_verified_sites == 1
    assert (table.sites[0].protein_accession, table.sites[0].cys_position) == (
        "Q6H3X6",
        133,
    )


def test_low_localization_probability_is_dropped(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "396",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                loc_prob=0.5,
                intensity_lcd=1000.0,
            )
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME, localization_min=0.75)
    assert table.total_verified_sites == 0
    assert table.dropped_low_localization == 1


def test_no_quantified_intensity_is_dropped(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "396",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                intensity_lcd=0.0,
                intensity_wt=0.0,
            )
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.dropped_no_quantified_intensity == 1


def test_non_cysteine_amino_acid_is_dropped(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "396",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                amino="S",
                intensity_lcd=1000.0,
            )
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.dropped_not_cysteine == 1


def test_coordinate_mismatch_is_dropped_not_repaired(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "5",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                intensity_lcd=1000.0,
            )
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.dropped_coordinate_mismatch == 1


def test_missing_accession_is_dropped(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|Q9NOTREAL|NOTREAL_SOLLC",
                "10",
                "tr|Q9NOTREAL|NOTREAL_SOLLC",
                intensity_lcd=1000.0,
            )
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.dropped_missing_accession == 1


def test_duplicate_accession_position_counted_once(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "396",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                intensity_lcd=1000.0,
            ),
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "396",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                intensity_lcd=500.0,
            ),
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    assert table.total_verified_sites == 1
    assert table.dropped_duplicate == 1


def test_rows_total_counts_every_data_row(tmp_path: Path) -> None:
    xlsx = tmp_path / "DSs.xlsx"
    _write_xlsx(
        xlsx,
        [
            _row(
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                "396",
                "tr|A0A3Q7F586|A0A3Q7F586_SOLLC",
                intensity_lcd=1000.0,
            ),
            _row(
                "tr|Q9NOTREAL|X_SOLLC", "10", "tr|Q9NOTREAL|X_SOLLC", intensity_lcd=1.0
            ),
        ],
    )
    table = parse_kiae271_sites(xlsx, MINI_PROTEOME)
    assert table.rows_total == 2
    assert table.total_verified_sites == 1
    assert table.dropped_missing_accession == 1
