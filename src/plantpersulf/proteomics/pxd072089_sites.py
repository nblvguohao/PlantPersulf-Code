"""PXD072089 (Oryza sativa persulfidome) site-table parser.

Consumes three independent evidence tables for Xie et al. 2026 (PNAS,
doi:10.1073/pnas.2608150123) — the rice persulfidome mapped via NM-biotin
chemistry with DTT selective elution — and derives a coordinate-verified
site union:

* **SD01 single-Cys peptides** (PNAS Supplementary Dataset 1): for peptides
  containing exactly one cysteine, the ``Position`` column is the
  persulfidated position.

* **SD04 site-level table** (PNAS Supplementary Dataset 4): every row
  carries an explicit ``-SSH Cys position`` column that resolves which
  cysteine is persulfidated — including multi-Cys peptides.

* **SS-all-peptides.tsv** (submitter-deposited MaxQuant ``peptides.txt``
  export on PRIDE, PXD072089): the full pre-filter peptide table. For
  peptides with exactly one cysteine (``C Count == 1``), the site is the
  single Cys in the ``Sequence`` window given by ``Start position``/
  ``End position``. This is a superset of SD01's single-Cys coverage
  (proteins with deleted UniProt accessions are excluded from all three
  tables identically — the accession deletion happened upstream of every
  published/deposited output, so no table recovers those positions).

All three are merged (union on (accession, position)) and every site must
pass coordinate verification against the reference proteome — accession
exists, position in range, ``sequence[position - 1] == 'C'`` — before it
counts; every dropped row is counted, never silently repaired.

Chemistry context: NM-biotin alkylates both Cys-SH and Cys-SSH; DTT
selectively cleaves the -SSH disulfide, liberating only persulfidated
peptides for LC-MS/MS. Independent of the Seville tag-switch benchmark
(lab, chemistry, and species axes).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any

SOURCE_SD01 = "sd01"
SOURCE_SD04 = "sd04"
SOURCE_SS = "ss"
_SOURCE_ORDER = (SOURCE_SD01, SOURCE_SD04, SOURCE_SS)


def _combine_source(tags: set[str]) -> str:
    return "+".join(t for t in _SOURCE_ORDER if t in tags)


@dataclass(frozen=True)
class PXD072089Site:
    protein_accession: str
    cys_position: int
    source: str  # "+"-joined subset of {"sd01", "sd04", "ss"}, e.g. "sd01+ss"


@dataclass(frozen=True)
class PXD072089SiteTable:
    sites: tuple[PXD072089Site, ...]
    # SD01 drops
    sd01_rows_total: int
    sd01_dropped_non_cysteine: int
    sd01_dropped_missing_accession: int
    sd01_dropped_coordinate_mismatch: int
    sd01_single_cys_parsed: int
    sd01_multi_cys_skipped: int
    # SD04 drops
    sd04_rows_total: int
    sd04_dropped_missing_accession: int
    sd04_dropped_coordinate_mismatch: int
    sd04_parsed: int
    # SS-all-peptides drops (0 / unset when ss_path not provided)
    ss_rows_total: int
    ss_dropped_missing_accession: int
    ss_dropped_coordinate_mismatch: int
    ss_multi_cys_skipped: int
    ss_single_cys_parsed: int
    # Union
    total_verified_sites: int
    total_verified_proteins: int


SD01_REQUIRED_COLS = (
    "Sequence",
    "Amino acid",
    "Position",
    "Leading razor protein",
)

SD04_REQUIRED_COLS = (
    "Peptide sequence",
    "Protein_ID",
    "leading razor protein",
    "-SSH Cys position",
)

SS_REQUIRED_COLS = (
    "Sequence",
    "C Count",
    "Leading razor protein",
    "Start position",
)


def _read_xlsx(path: Path) -> list[dict[str, Any]]:
    """Read an xlsx into a list of dicts — first row is title, second row is header."""
    import openpyxl  # type: ignore[import-untyped]

    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows_iter = ws.iter_rows(values_only=True)
    # skip title row
    try:
        next(rows_iter)
    except StopIteration:
        wb.close()
        return []
    header = [str(c).strip() if c is not None else "" for c in next(rows_iter)]
    records: list[dict[str, Any]] = []
    for row in rows_iter:
        record: dict[str, Any] = {}
        for i, value in enumerate(row):
            if i < len(header):
                record[header[i]] = value
        records.append(record)
    wb.close()
    return records


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def parse_pxd072089_sites(
    sd01_path: Path,
    sd04_path: Path,
    proteome: dict[str, str],
    ss_all_peptides_path: Path | None = None,
    allowed_accessions: set[str] | None = None,
) -> PXD072089SiteTable:
    """Parse SD01 + SD04 (+ optional SS-all-peptides) into a coordinate-verified
    site union.

    Every site must pass coordinate verification against *proteome*: the
    accession must exist, the position must be in range, and
    ``sequence[position - 1]`` must be ``'C'``.
    """

    # --- SD01: single-Cys peptides only ---
    sd01 = _read_xlsx(sd01_path)
    sd01_rows_total = len(sd01)
    sd01_dropped_non_cysteine = 0
    sd01_dropped_missing_accession = 0
    sd01_dropped_coordinate_mismatch = 0
    sd01_single_cys_parsed = 0
    sd01_multi_cys_skipped = 0

    site_sources: dict[tuple[str, int], set[str]] = {}

    for row in sd01:
        amino = str(row.get("Amino acid", "") or "").strip()
        if amino != "C":
            sd01_dropped_non_cysteine += 1
            continue

        pos_str = str(row.get("Position", "") or "").strip()
        positions = [p.strip() for p in pos_str.split(",") if p.strip()]
        if len(positions) != 1:
            sd01_multi_cys_skipped += 1
            continue

        accession = str(row.get("Leading razor protein", "") or "").strip()
        if not accession:
            sd01_dropped_missing_accession += 1
            continue
        if allowed_accessions is not None and accession not in allowed_accessions:
            continue

        try:
            position = int(positions[0])
        except ValueError:
            sd01_dropped_coordinate_mismatch += 1
            continue

        sequence = proteome.get(accession)
        if sequence is None:
            sd01_dropped_missing_accession += 1
            continue
        if position < 1 or position > len(sequence) or sequence[position - 1] != "C":
            sd01_dropped_coordinate_mismatch += 1
            continue

        site_sources.setdefault((accession, position), set()).add(SOURCE_SD01)
        sd01_single_cys_parsed += 1

    # --- SD04: definitive site-level table ---
    sd04 = _read_xlsx(sd04_path)
    sd04_rows_total = len(sd04)
    sd04_dropped_missing_accession = 0
    sd04_dropped_coordinate_mismatch = 0
    sd04_parsed = 0

    for row in sd04:
        accession = str(row.get("leading razor protein", "") or "").strip()
        if not accession:
            sd04_dropped_missing_accession += 1
            continue
        if allowed_accessions is not None and accession not in allowed_accessions:
            continue

        pos_val = row.get("-SSH Cys position")
        try:
            position = int(str(pos_val))
        except (ValueError, TypeError):
            sd04_dropped_coordinate_mismatch += 1
            continue

        sequence = proteome.get(accession)
        if sequence is None:
            sd04_dropped_missing_accession += 1
            continue
        if position < 1 or position > len(sequence) or sequence[position - 1] != "C":
            sd04_dropped_coordinate_mismatch += 1
            continue

        site_sources.setdefault((accession, position), set()).add(SOURCE_SD04)
        sd04_parsed += 1

    # --- SS-all-peptides: full pre-filter MaxQuant peptide table (optional) ---
    ss_rows_total = 0
    ss_dropped_missing_accession = 0
    ss_dropped_coordinate_mismatch = 0
    ss_multi_cys_skipped = 0
    ss_single_cys_parsed = 0

    if ss_all_peptides_path is not None:
        ss = _read_tsv(ss_all_peptides_path)
        ss_rows_total = len(ss)
        for row in ss:
            try:
                n_cys = int(row.get("C Count", "0") or 0)
            except ValueError:
                n_cys = 0
            if n_cys == 0:
                continue
            if n_cys != 1:
                ss_multi_cys_skipped += 1
                continue

            accession = (row.get("Leading razor protein") or "").strip()
            if not accession:
                ss_dropped_missing_accession += 1
                continue
            if allowed_accessions is not None and accession not in allowed_accessions:
                continue

            sequence = proteome.get(accession)
            if sequence is None:
                ss_dropped_missing_accession += 1
                continue

            pep = row.get("Sequence") or ""
            try:
                start = int(row.get("Start position") or "")
            except ValueError:
                ss_dropped_coordinate_mismatch += 1
                continue

            cys_rel = [i for i, c in enumerate(pep) if c == "C"]
            if len(cys_rel) != 1:
                ss_dropped_coordinate_mismatch += 1
                continue
            position = start + cys_rel[0]

            if (
                position < 1
                or position > len(sequence)
                or sequence[position - 1] != "C"
            ):
                ss_dropped_coordinate_mismatch += 1
                continue

            site_sources.setdefault((accession, position), set()).add(SOURCE_SS)
            ss_single_cys_parsed += 1

    # --- Build sorted site list ---
    sites = tuple(
        PXD072089Site(
            protein_accession=acc,
            cys_position=pos,
            source=_combine_source(site_sources[(acc, pos)]),
        )
        for acc, pos in sorted(site_sources)
    )
    proteins = {s.protein_accession for s in sites}

    return PXD072089SiteTable(
        sites=sites,
        sd01_rows_total=sd01_rows_total,
        sd01_dropped_non_cysteine=sd01_dropped_non_cysteine,
        sd01_dropped_missing_accession=sd01_dropped_missing_accession,
        sd01_dropped_coordinate_mismatch=sd01_dropped_coordinate_mismatch,
        sd01_single_cys_parsed=sd01_single_cys_parsed,
        sd01_multi_cys_skipped=sd01_multi_cys_skipped,
        sd04_rows_total=sd04_rows_total,
        sd04_dropped_missing_accession=sd04_dropped_missing_accession,
        sd04_dropped_coordinate_mismatch=sd04_dropped_coordinate_mismatch,
        sd04_parsed=sd04_parsed,
        ss_rows_total=ss_rows_total,
        ss_dropped_missing_accession=ss_dropped_missing_accession,
        ss_dropped_coordinate_mismatch=ss_dropped_coordinate_mismatch,
        ss_multi_cys_skipped=ss_multi_cys_skipped,
        ss_single_cys_parsed=ss_single_cys_parsed,
        total_verified_sites=len(sites),
        total_verified_proteins=len(proteins),
    )


def build_cross_species_eval_rows(
    table: PXD072089SiteTable,
    proteome: dict[str, str],
    study_accession: str = "PXD072089",
    source_sha256: str = "",
) -> tuple[dict[str, str], ...]:
    """Benchmark-schema evaluation rows for the rice cross-species transfer track.

    Every parsed site becomes a ``positive`` row (with study/evidence/source
    provenance); every OTHER cysteine in the proteome becomes an ``unlabeled``
    row. PU semantics: unlabeled means "not detected", never a hard negative,
    and downstream consumers must not train on it as one. Rows are sorted by
    (accession, position) so the output is deterministic; every positive key
    is cysteine-verified by the parser, so the row count equals the proteome
    cysteine count exactly.
    """
    positive_keys = {(s.protein_accession, s.cys_position) for s in table.sites}
    rows: list[dict[str, str]] = []
    n_positive = 0
    for accession in sorted(proteome):
        sequence = proteome[accession]
        for position, residue in enumerate(sequence, start=1):
            if residue != "C":
                continue
            if (accession, position) in positive_keys:
                n_positive += 1
                rows.append(
                    {
                        "protein_accession": accession,
                        "cys_position_in_protein": str(position),
                        "label": "positive",
                        "study_accession": study_accession,
                        "evidence_level": "site_ms",
                        "source_sha256": source_sha256,
                    }
                )
            else:
                rows.append(
                    {
                        "protein_accession": accession,
                        "cys_position_in_protein": str(position),
                        "label": "unlabeled",
                        "study_accession": "",
                        "evidence_level": "",
                        "source_sha256": "",
                    }
                )
    if n_positive != len(table.sites):
        raise RuntimeError(
            "cross-species eval rows lost positives: "
            f"{n_positive} emitted vs {len(table.sites)} parsed"
        )
    return tuple(rows)
