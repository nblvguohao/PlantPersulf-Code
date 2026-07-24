"""RED: PXD072089 (Oryza sativa persulfidome) site-table parser.

The parser consumes SD01 (peptide-level) + SD04 (site-level) from the
PNAS supplementary data, optionally joined by SS-all-peptides.tsv (the
submitter-deposited MaxQuant peptide table on PRIDE), and derives a
coordinate-verified site union. Every reported site must be a cysteine
whose coordinate matches the rice reference proteome sequence; anything
else is dropped and counted, never silently repaired.

PU semantics: positives = the verified persulfidation sites (897 from
SD01+SD04 alone, 929 once SS-all-peptides is joined in); unlabeled =
every other cysteine in the rice proteome.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.proteomics.pxd072089_sites import (
    PXD072089SiteTable,
    build_cross_species_eval_rows,
    parse_pxd072089_sites,
)

# Minimal test proteome with UniProt-style headers
MINI_PROTEOME = {
    "Q6YZ09": "MNYFAFTPLLPSVTCGTVEPR",  # Cys at 15
    "Q8S5T1": "MAAVCEMPFATVASDDLGGVGGTCVLR",  # Cys at 5, 24
    "A3BE48": "MLVGSVPLVHADDVCDALVFCMDQPSLAGR",  # Cys at 15, 22
}


def _write_xlsx(path: Path, header: tuple[str, ...], rows: list[tuple]) -> None:
    """Write a minimal xlsx with a title row + header + data rows."""
    import openpyxl

    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(("Title row (ignored)",))
    ws.append(header)
    for row in rows:
        ws.append(row)
    wb.save(path)


# ---------------------------------------------------------------------------
# SD01 + SD04 synthetic data
# ---------------------------------------------------------------------------

SD01_HEADER = (
    "Sequence",
    "Amino acid",
    "Position",
    "Cysteine quantity",
    "Unique (Groups)",
    "Charges",
    "Start position",
    "End position",
    "Leading razor protein",
    "Protein description",
)

SD04_HEADER = (
    "Peptide sequence",
    "Protein_ID",
    "leading razor protein",
    "-SSH Cys position",
    "Protein description",
    "Domain",
)


def _make_sd01(path: Path, rows: list[tuple]) -> None:
    _write_xlsx(path, SD01_HEADER, rows)


def _make_sd04(path: Path, rows: list[tuple]) -> None:
    _write_xlsx(path, SD04_HEADER, rows)


# ---------------------------------------------------------------------------
# Basic parsing tests
# ---------------------------------------------------------------------------


def test_sd01_single_cys_site_is_accepted(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(
        sd01,
        [
            (
                "VTCGTVEPR",
                "C",
                "15",
                "1",
                "yes",
                "2",
                "12",
                "20",
                "Q6YZ09",
                "NADH protein",
            ),
        ],
    )
    _make_sd04(sd04, [])  # empty
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 1
    assert table.total_verified_proteins == 1
    assert table.sites[0].protein_accession == "Q6YZ09"
    assert table.sites[0].cys_position == 15
    assert table.sites[0].source == "sd01"
    assert table.sd01_single_cys_parsed == 1


def test_sd04_site_level_row_is_accepted(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(sd01, [])
    _make_sd04(
        sd04,
        [
            (
                "AAVCEMPFATVASDDLGGVGGTCVLR",
                "Q8S5T1;Q7G6Z6",
                "Q8S5T1",
                5,
                "Glutathione reductase",
                "FAD/NAD(P)-binding",
            ),
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 1
    assert table.sites[0].protein_accession == "Q8S5T1"
    assert table.sites[0].cys_position == 5
    assert table.sites[0].source == "sd04"


def test_sd01_and_sd04_overlap_is_marked_combined(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(
        sd01,
        [
            (
                "MAAVCEMPFATVASDDLGGVGGTCVLR",
                "C",
                "5",
                "2",
                "yes",
                "3",
                "1",
                "27",
                "Q8S5T1",
                "GR protein",
            ),
        ],
    )
    _make_sd04(
        sd04,
        [
            (
                "AAVCEMPFATVASDDLGGVGGTCVLR",
                "Q8S5T1;Q7G6Z6",
                "Q8S5T1",
                5,
                "Glutathione reductase",
                "FAD/NAD(P)-binding",
            ),
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 1
    assert table.sites[0].source == "sd01+sd04"


# ---------------------------------------------------------------------------
# Multi-Cys peptides
# ---------------------------------------------------------------------------


def test_sd01_multi_cys_peptide_is_skipped(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(
        sd01,
        [
            (
                "AAVCEMPFATVASDDLGGVGGTCVLR",
                "C",
                "5,25",
                "2",
                "yes",
                "3",
                "1",
                "27",
                "Q8S5T1",
                "GR protein",
            ),
        ],
    )
    _make_sd04(sd04, [])
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.sd01_multi_cys_skipped == 1


def test_sd04_resolves_multi_cys_peptide_individual_sites(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(sd01, [])
    _make_sd04(
        sd04,
        [
            (
                "AAVCEMPFATVASDDLGGVGGTCVLR",
                "Q8S5T1;Q7G6Z6",
                "Q8S5T1",
                5,
                "GR protein",
                "FAD/NAD(P)-binding",
            ),
            (
                "AAVCEMPFATVASDDLGGVGGTCVLR",
                "Q8S5T1;Q7G6Z6",
                "Q8S5T1",
                24,
                "GR protein",
                "FAD/NAD(P)-binding",
            ),
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 2
    positions = {s.cys_position for s in table.sites}
    assert positions == {5, 24}


# ---------------------------------------------------------------------------
# SS-all-peptides.tsv (optional third source)
# ---------------------------------------------------------------------------

SS_FIELDS = ("Sequence", "C Count", "Leading razor protein", "Start position")


def _write_ss_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(SS_FIELDS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(f, "")) for f in SS_FIELDS) + "\n")


def test_ss_single_cys_peptide_is_accepted(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    ss = tmp_path / "ss.tsv"
    _make_sd01(sd01, [])
    _make_sd04(sd04, [])
    _write_ss_tsv(
        ss,
        [
            {
                "Sequence": "VTCGTVEPR",
                "C Count": "1",
                "Leading razor protein": "Q6YZ09",
                "Start position": "13",  # Cys at rel-idx 2 -> abs 15
            },
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME, ss_all_peptides_path=ss)
    assert table.total_verified_sites == 1
    assert table.sites[0].protein_accession == "Q6YZ09"
    assert table.sites[0].cys_position == 15
    assert table.sites[0].source == "ss"
    assert table.ss_single_cys_parsed == 1


def test_ss_multi_cys_peptide_is_skipped(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    ss = tmp_path / "ss.tsv"
    _make_sd01(sd01, [])
    _make_sd04(sd04, [])
    _write_ss_tsv(
        ss,
        [
            {
                "Sequence": "VCDALVFCMDQ",
                "C Count": "2",
                "Leading razor protein": "A3BE48",
                "Start position": "14",
            },
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME, ss_all_peptides_path=ss)
    assert table.total_verified_sites == 0
    assert table.ss_multi_cys_skipped == 1


def test_ss_missing_accession_is_counted(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    ss = tmp_path / "ss.tsv"
    _make_sd01(sd01, [])
    _make_sd04(sd04, [])
    _write_ss_tsv(
        ss,
        [
            {
                "Sequence": "ACGTACGT",
                "C Count": "1",
                "Leading razor protein": "P99999",
                "Start position": "1",
            },
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME, ss_all_peptides_path=ss)
    assert table.total_verified_sites == 0
    assert table.ss_dropped_missing_accession == 1


def test_ss_overlapping_with_sd01_and_sd04_is_combined(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    ss = tmp_path / "ss.tsv"
    _make_sd01(
        sd01,
        [
            (
                "VTCGTVEPR",
                "C",
                "15",
                "1",
                "yes",
                "2",
                "12",
                "20",
                "Q6YZ09",
                "NADH protein",
            ),
        ],
    )
    _make_sd04(sd04, [])
    _write_ss_tsv(
        ss,
        [
            {
                "Sequence": "VTCGTVEPR",
                "C Count": "1",
                "Leading razor protein": "Q6YZ09",
                "Start position": "13",
            },
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME, ss_all_peptides_path=ss)
    assert table.total_verified_sites == 1
    assert table.sites[0].source == "sd01+ss"


def test_ss_default_absent_leaves_ss_counters_at_zero(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(sd01, [])
    _make_sd04(sd04, [])
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.ss_rows_total == 0
    assert table.ss_single_cys_parsed == 0
    assert table.ss_multi_cys_skipped == 0
    assert table.ss_dropped_missing_accession == 0
    assert table.ss_dropped_coordinate_mismatch == 0


# ---------------------------------------------------------------------------
# Coordinate verification
# ---------------------------------------------------------------------------


def test_missing_accession_is_counted_in_sd01(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(
        sd01,
        [
            (
                "ACGTACGT",
                "C",
                "5",
                "1",
                "yes",
                "2",
                "1",
                "8",
                "P99999",
                "Deleted protein",
            ),
        ],
    )
    _make_sd04(sd04, [])
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.sd01_dropped_missing_accession == 1


def test_missing_accession_is_counted_in_sd04(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(sd01, [])
    _make_sd04(
        sd04,
        [
            (
                "ACGTACGT",
                "P99999",
                "P99999",
                5,
                "Deleted protein",
                "Domain",
            ),
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.sd04_dropped_missing_accession == 1


def test_coordinate_mismatch_is_dropped_in_sd01(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(
        sd01,
        [
            (
                "MTCGTVEPR",  # Cys at position 3
                "C",
                "1",  # position 1 is M, not C
                "1",
                "yes",
                "2",
                "1",
                "9",
                "Q6YZ09",
                "mismatched",
            ),
        ],
    )
    _make_sd04(sd04, [])
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.sd01_dropped_coordinate_mismatch == 1


def test_coordinate_mismatch_is_dropped_in_sd04(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(sd01, [])
    _make_sd04(
        sd04,
        [
            (
                "AAVCEMPFATVASDDLGGVGGTCVLR",
                "Q8S5T1;Q7G6Z6",
                "Q8S5T1",
                1,  # position 1 is M, not C
                "mismatched",
                "Domain",
            ),
        ],
    )
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.sd04_dropped_coordinate_mismatch == 1


def test_non_cysteine_sd01_row_is_dropped_and_counted(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(
        sd01,
        [
            (
                "MTSGTVEPR",
                "S",
                "3",
                "1",
                "yes",
                "2",
                "1",
                "9",
                "Q6YZ09",
                "not cysteine",
            ),
        ],
    )
    _make_sd04(sd04, [])
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.sd01_dropped_non_cysteine == 1


def test_position_out_of_range_is_dropped(tmp_path: Path) -> None:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(
        sd01,
        [
            (
                "VTCGTVEPR",
                "C",
                "999",
                "1",
                "yes",
                "2",
                "1",
                "9",
                "Q6YZ09",
                "bad position",
            ),
        ],
    )
    _make_sd04(sd04, [])
    table = parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)
    assert table.total_verified_sites == 0
    assert table.sd01_dropped_coordinate_mismatch == 1


# ---------------------------------------------------------------------------
# Cross-species eval row construction
# ---------------------------------------------------------------------------


def _parsed_table(tmp_path: Path) -> PXD072089SiteTable:
    sd01 = tmp_path / "sd01.xlsx"
    sd04 = tmp_path / "sd04.xlsx"
    _make_sd01(
        sd01,
        [
            (
                "VTCGTVEPR",
                "C",
                "15",
                "1",
                "yes",
                "2",
                "12",
                "20",
                "Q6YZ09",
                "NADH protein",
            ),
        ],
    )
    _make_sd04(
        sd04,
        [
            (
                "AAVCEMPFATVASDDLGGVGGTCVLR",
                "Q8S5T1;Q7G6Z6",
                "Q8S5T1",
                5,
                "GR protein",
                "FAD/NAD(P)-binding",
            ),
        ],
    )
    return parse_pxd072089_sites(sd01, sd04, MINI_PROTEOME)


def test_eval_rows_mark_sites_positive_and_other_cysteines_unlabeled(
    tmp_path: Path,
) -> None:
    table = _parsed_table(tmp_path)
    rows = build_cross_species_eval_rows(
        table, MINI_PROTEOME, study_accession="PXD072089", source_sha256="deadbeef"
    )
    by_key = {(r["protein_accession"], r["cys_position_in_protein"]): r for r in rows}

    # (Q6YZ09, 15) is positive
    pos = by_key[("Q6YZ09", "15")]
    assert pos["label"] == "positive"
    assert pos["study_accession"] == "PXD072089"
    assert pos["source_sha256"] == "deadbeef"

    # (Q8S5T1, 5) is positive
    pos2 = by_key[("Q8S5T1", "5")]
    assert pos2["label"] == "positive"

    # (Q8S5T1, 24) is unlabeled (Cys but not a persulfidation site)
    unl = by_key[("Q8S5T1", "24")]
    assert unl["label"] == "unlabeled"
    assert unl["study_accession"] == ""

    # (A3BE48, 15) is unlabeled
    unl2 = by_key[("A3BE48", "15")]
    assert unl2["label"] == "unlabeled"


def test_eval_rows_cover_every_proteome_cysteine_exactly_once(
    tmp_path: Path,
) -> None:
    table = _parsed_table(tmp_path)
    rows = build_cross_species_eval_rows(table, MINI_PROTEOME)
    n_cys = sum(seq.count("C") for seq in MINI_PROTEOME.values())
    assert len(rows) == n_cys
    keys = [(r["protein_accession"], r["cys_position_in_protein"]) for r in rows]
    assert len(set(keys)) == len(keys)
    # Deterministic order: sorted by (accession, numeric position)
    assert keys == sorted(keys, key=lambda k: (k[0], int(k[1])))


def test_eval_rows_never_emit_non_cysteine_positions(tmp_path: Path) -> None:
    table = _parsed_table(tmp_path)
    rows = build_cross_species_eval_rows(table, MINI_PROTEOME)
    for row in rows:
        seq = MINI_PROTEOME[row["protein_accession"]]
        assert seq[int(row["cys_position_in_protein"]) - 1] == "C"


def test_eval_rows_positives_match_table_sites(tmp_path: Path) -> None:
    table = _parsed_table(tmp_path)
    rows = build_cross_species_eval_rows(table, MINI_PROTEOME)
    n_pos = sum(1 for r in rows if r["label"] == "positive")
    assert n_pos == len(table.sites)


# ---------------------------------------------------------------------------
# Real-data tests (skip unless registered downloads present)
# ---------------------------------------------------------------------------

REAL_SD01 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd01.xlsx")
REAL_SD04 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd04.xlsx")
REAL_SS = Path("data/raw/supplements/PXD072089/SS-all-peptides.tsv")
REAL_PROTEOME = Path("data/raw/references/rice_proteome_v1/uniprot_rice_v1.fasta")

_REAL_INPUTS_PRESENT = (
    REAL_SD01.is_file() and REAL_SD04.is_file() and REAL_PROTEOME.is_file()
)
_REAL_INPUTS_WITH_SS_PRESENT = _REAL_INPUTS_PRESENT and REAL_SS.is_file()


@pytest.mark.skipif(
    not _REAL_INPUTS_PRESENT,
    reason="registered PXD072089 downloads or rice proteome not present",
)
def test_real_pxd072089_sites_verify_against_rice_proteome() -> None:
    from plantpersulf.features.sequence import _load_proteome

    proteome = _load_proteome(REAL_PROTEOME)
    table = parse_pxd072089_sites(REAL_SD01, REAL_SD04, proteome)

    # SD01+SD04 alone: 897 coordinate-verified sites on 646 proteins
    assert table.total_verified_sites == 897
    assert table.total_verified_proteins == 646

    # Every parsed site must be Cys at the claimed position
    for site in table.sites:
        seq = proteome[site.protein_accession]
        assert seq[site.cys_position - 1] == "C"

    # Source breakdown sanity
    assert table.sd01_single_cys_parsed > 780
    assert table.sd04_parsed > 350
    assert table.sd01_multi_cys_skipped > 100

    # No coordinate mismatches in real data
    assert table.sd01_dropped_coordinate_mismatch == 0
    assert table.sd04_dropped_coordinate_mismatch <= 10

    # Missing accessions are expected (UniProt 2025-2026 cleanup)
    assert table.sd01_dropped_missing_accession > 400
    assert table.sd04_dropped_missing_accession > 200


@pytest.mark.skipif(
    not _REAL_INPUTS_WITH_SS_PRESENT,
    reason="registered PXD072089 downloads/SS-all-peptides/rice proteome not present",
)
def test_real_pxd072089_sites_with_ss_reach_929() -> None:
    from plantpersulf.features.sequence import _load_proteome

    proteome = _load_proteome(REAL_PROTEOME)
    table = parse_pxd072089_sites(
        REAL_SD01, REAL_SD04, proteome, ss_all_peptides_path=REAL_SS
    )

    # SD01+SD04+SS: 929 coordinate-verified sites (32 more than SD01+SD04 alone)
    assert table.total_verified_sites == 929

    for site in table.sites:
        seq = proteome[site.protein_accession]
        assert seq[site.cys_position - 1] == "C"

    # SS is a superset of SD01's single-Cys coverage — every combined-source
    # tag containing "ss" alone (no sd01/sd04) represents a genuinely new site.
    ss_only = [s for s in table.sites if s.source == "ss"]
    assert len(ss_only) == 32


@pytest.mark.skipif(
    not _REAL_INPUTS_PRESENT,
    reason="registered PXD072089 downloads or rice proteome not present",
)
def test_real_eval_rows_background_matches_rice_cysteine_count() -> None:
    from plantpersulf.features.sequence import _load_proteome

    proteome = _load_proteome(REAL_PROTEOME)
    table = parse_pxd072089_sites(REAL_SD01, REAL_SD04, proteome)
    rows = build_cross_species_eval_rows(table, proteome)

    n_cys = sum(seq.count("C") for seq in proteome.values())
    assert len(rows) == n_cys
    n_pos = sum(1 for r in rows if r["label"] == "positive")
    assert n_pos == len(table.sites)
    assert n_pos == 897
