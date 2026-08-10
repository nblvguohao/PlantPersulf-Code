"""kiae271 (SlLCD1-OE vs WT tomato leaf) differential persulfidation parser.

Zhang et al. 2024 (Plant Physiology, doi:10.1093/plphys/kiae271) reports a
proteome-wide 4D label-free differential persulfidation screen ("119
peptides") comparing SlLCD1-overexpression to WT tomato leaves. The main text
only tabulates one site (SlWRKY6 Cys396); the full table is not restated
anywhere in the paper's prose, but it IS the paper's own published
Supplementary Dataset S1 (mistitled "Persulfidation site analysis data of
SlWRKY6 by LC-MS/MS" — the sheet actually holds all 119 rows, not just
SlWRKY6's). Registered source:
``data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx``
(sha256 a437f1a9..., see ``data/registry/supplementary_sources.tsv``).

Every reported site must be a cysteine whose coordinate matches the tomato
reference proteome sequence and must clear a minimum localization-probability
bar; anything else is dropped and counted, never silently repaired or
guessed. ``regulation`` classifies each verified site by which condition(s)
it was quantified in (PU semantics — "unlabeled" is never a hard negative,
so a site quantified only in WT is not claimed as "H2S-suppressed", only as
"not detected as elevated under LCD1-OE"):

* ``lcd_gain``  — nonzero intensity in LCD1-OE only (H2S-associated gain)
* ``wt_only``   — nonzero intensity in WT only
* ``both``      — nonzero intensity in both conditions
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

SHEET_NAME = "Supplementary Dataset S1"
DEFAULT_LOCALIZATION_MIN = 0.75

REGULATION_LCD_GAIN = "lcd_gain"
REGULATION_WT_ONLY = "wt_only"
REGULATION_BOTH = "both"


@dataclass(frozen=True)
class Kiae271Site:
    protein_accession: str
    cys_position: int
    localization_prob: float
    regulation: str  # lcd_gain | wt_only | both
    intensity_lcd: float
    intensity_wt: float


@dataclass(frozen=True)
class Kiae271SiteTable:
    sites: tuple[Kiae271Site, ...]
    rows_total: int
    dropped_not_cysteine: int
    dropped_no_position_alignment: int
    dropped_missing_accession: int
    dropped_coordinate_mismatch: int
    dropped_low_localization: int
    dropped_no_quantified_intensity: int
    dropped_duplicate: int
    total_verified_sites: int
    total_verified_proteins: int


def _header_accession(field: str) -> str:
    """Same convention as ``plantpersulf.features.sequence._header_accession``:
    UniProt ``tr|ACC|NAME`` / ``sp|ACC|NAME`` takes the second pipe field."""
    parts = field.strip().split("|")
    if len(parts) >= 2 and parts[1]:
        return parts[1]
    return field.strip()


def _read_xlsx_sheet(path: Path, sheet_name: str) -> list[dict[str, Any]]:
    """Read one named sheet into a list of dicts — row 1 is a title, row 2 is
    the header, data starts at row 3 (same convention as the PXD072089
    parser's ``_read_xlsx``, extended with an explicit sheet selector since
    this workbook has five sheets, not one)."""
    import openpyxl  # type: ignore[import-untyped]

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    try:
        ws = wb[sheet_name]
        rows_iter = ws.iter_rows(values_only=True)
        try:
            next(rows_iter)  # title row
        except StopIteration:
            return []
        header = [str(c).strip() if c is not None else "" for c in next(rows_iter)]
        records: list[dict[str, Any]] = []
        for row in rows_iter:
            if all(v is None for v in row):
                continue
            record: dict[str, Any] = {}
            for i, value in enumerate(row):
                if i < len(header) and header[i]:
                    record[header[i]] = value
            records.append(record)
        return records
    finally:
        wb.close()


def _to_float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def parse_kiae271_sites(
    xlsx_path: Path,
    proteome: dict[str, str],
    localization_min: float = DEFAULT_LOCALIZATION_MIN,
    sheet_name: str = SHEET_NAME,
) -> Kiae271SiteTable:
    """Parse Supplementary Dataset S1 into a coordinate-verified site table.

    ``Positions within proteins`` and ``Proteins`` are parallel ``;``-joined
    lists (MaxQuant protein-group convention); the position matching
    ``Leading proteins`` is found by locating that accession's index in
    ``Proteins``, not by assuming index 0 — shared/razor peptides list
    multiple proteins and the leading one is not always first.
    """
    rows = _read_xlsx_sheet(xlsx_path, sheet_name)
    rows_total = len(rows)

    dropped_not_cysteine = 0
    dropped_no_position_alignment = 0
    dropped_missing_accession = 0
    dropped_coordinate_mismatch = 0
    dropped_low_localization = 0
    dropped_no_quantified_intensity = 0
    dropped_duplicate = 0

    seen: dict[tuple[str, int], Kiae271Site] = {}

    for row in rows:
        amino = str(row.get("Amino acid", "") or "").strip()
        if amino != "C":
            dropped_not_cysteine += 1
            continue

        proteins_field = str(row.get("Proteins", "") or "")
        positions_field = str(row.get("Positions within proteins", "") or "")
        leading = str(row.get("Leading proteins", "") or "").strip()
        proteins_list = [p.strip() for p in proteins_field.split(";") if p.strip()]
        positions_list = [p.strip() for p in positions_field.split(";") if p.strip()]
        if (
            not leading
            or leading not in proteins_list
            or len(proteins_list) != len(positions_list)
        ):
            dropped_no_position_alignment += 1
            continue
        idx = proteins_list.index(leading)
        pos_str = positions_list[idx]
        try:
            position = int(pos_str)
        except ValueError:
            dropped_no_position_alignment += 1
            continue

        accession = _header_accession(leading)
        sequence = proteome.get(accession)
        if sequence is None:
            dropped_missing_accession += 1
            continue
        if position < 1 or position > len(sequence) or sequence[position - 1] != "C":
            dropped_coordinate_mismatch += 1
            continue

        loc_prob = _to_float(row.get("Localization prob"))
        if loc_prob < localization_min:
            dropped_low_localization += 1
            continue

        intensity_lcd = _to_float(row.get("Intensity LCD"))
        intensity_wt = _to_float(row.get("Intensity WT"))
        if intensity_lcd <= 0 and intensity_wt <= 0:
            dropped_no_quantified_intensity += 1
            continue
        if intensity_lcd > 0 and intensity_wt > 0:
            regulation = REGULATION_BOTH
        elif intensity_lcd > 0:
            regulation = REGULATION_LCD_GAIN
        else:
            regulation = REGULATION_WT_ONLY

        key = (accession, position)
        if key in seen:
            dropped_duplicate += 1
            continue
        seen[key] = Kiae271Site(
            protein_accession=accession,
            cys_position=position,
            localization_prob=loc_prob,
            regulation=regulation,
            intensity_lcd=intensity_lcd,
            intensity_wt=intensity_wt,
        )

    sites = tuple(seen[k] for k in sorted(seen))
    proteins = {s.protein_accession for s in sites}

    return Kiae271SiteTable(
        sites=sites,
        rows_total=rows_total,
        dropped_not_cysteine=dropped_not_cysteine,
        dropped_no_position_alignment=dropped_no_position_alignment,
        dropped_missing_accession=dropped_missing_accession,
        dropped_coordinate_mismatch=dropped_coordinate_mismatch,
        dropped_low_localization=dropped_low_localization,
        dropped_no_quantified_intensity=dropped_no_quantified_intensity,
        dropped_duplicate=dropped_duplicate,
        total_verified_sites=len(sites),
        total_verified_proteins=len(proteins),
    )
