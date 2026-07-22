"""PXD063170 (Magnaporthe oryzae S-sulfhydration) site-table parser.

Consumes the TSV export of Supplementary Data 1 from Hu et al. 2025
(Nat Commun 16, doi:10.1038/s41467-025-61582-8) — the CSE_OE vs WT
site-level S-sulfhydration table — and applies the same integrity rules as
the Arabidopsis benchmark: only confidently localized cysteine sites whose
coordinates verify against the reference proteome are kept; every dropped
row is counted, never silently repaired.

The TSV export is produced from the registered xlsx at fetch time (the
xlsx itself is the registered source of truth; see
``configs/supplementary_sources_v1.yaml``).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

REQUIRED_COLUMNS = (
    "Protein accession",
    "Position",
    "Amino acid",
    "Localization probability",
)


@dataclass(frozen=True)
class PXD063170Site:
    protein_accession: str
    cys_position: int
    localization_probability: float


@dataclass(frozen=True)
class PXD063170SiteTable:
    sites: tuple[PXD063170Site, ...]
    dropped_low_localization: int
    dropped_coordinate_mismatch: int
    dropped_missing_accession: int


def load_ensembl_fungi_proteome(path: Path) -> dict[str, str]:
    """Load an EnsemblFungi ``pep.all.fa`` file: header's first whitespace-
    separated token is the accession (e.g. ``>MGG_07573T0 pep chromosome:...``).
    """
    sequences: dict[str, str] = {}
    accession = ""
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if accession:
                sequences[accession] = "".join(chunks)
            accession = line[1:].strip().split()[0]
            chunks = []
        elif line.strip():
            chunks.append(line.strip())
    if accession:
        sequences[accession] = "".join(chunks)
    if not sequences:
        raise RuntimeError(f"proteome fasta is empty: {path}")
    return sequences


def parse_pxd063170_sites(
    source_tsv: Path,
    proteome: dict[str, str],
    min_localization: float = 0.75,
) -> PXD063170SiteTable:
    """Parse, filter, and coordinate-verify the PXD063170 site table.

    Rows failing the localization threshold, referencing proteins absent
    from the proteome, or whose position does not land on a cysteine in the
    sequence are dropped and counted. A row whose ``Amino acid`` field is
    not ``C`` is malformed (the table is defined as cysteine-only) and
    fails closed with a RuntimeError.
    """
    with source_tsv.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        if not all(col in fieldnames for col in REQUIRED_COLUMNS):
            raise RuntimeError(
                f"PXD063170 site table has invalid columns: {source_tsv}"
            )
        best: dict[tuple[str, int], float] = {}
        dropped_loc = 0
        dropped_mismatch = 0
        dropped_missing = 0
        for row in reader:
            amino = (row["Amino acid"] or "").strip()
            if amino != "C":
                raise RuntimeError(f"non-cysteine row in PXD063170 site table: {row}")
            accession = (row["Protein accession"] or "").strip()
            try:
                position = int((row["Position"] or "").strip())
                localization = float((row["Localization probability"] or "").strip())
            except ValueError as exc:
                raise RuntimeError(
                    f"malformed numeric field in PXD063170 site table: {row}"
                ) from exc
            if localization < min_localization:
                dropped_loc += 1
                continue
            sequence = proteome.get(accession)
            if sequence is None:
                dropped_missing += 1
                continue
            if (
                position < 1
                or position > len(sequence)
                or sequence[position - 1] != "C"
            ):
                dropped_mismatch += 1
                continue
            key = (accession, position)
            if localization > best.get(key, -1.0):
                best[key] = localization

    sites = tuple(
        PXD063170Site(
            protein_accession=accession,
            cys_position=position,
            localization_probability=best[(accession, position)],
        )
        for accession, position in sorted(best)
    )
    return PXD063170SiteTable(
        sites=sites,
        dropped_low_localization=dropped_loc,
        dropped_coordinate_mismatch=dropped_mismatch,
        dropped_missing_accession=dropped_missing,
    )


def build_cross_species_eval_rows(
    table: PXD063170SiteTable,
    proteome: dict[str, str],
    study_accession: str = "PXD063170",
    source_sha256: str = "",
) -> tuple[dict[str, str], ...]:
    """Benchmark-schema evaluation rows for the cross-species transfer track.

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
