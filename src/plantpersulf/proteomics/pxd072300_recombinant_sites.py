"""PXD072300 (recombinant rice protein persulfidation) control-site parser.

Consumes the 10 PEAKS DB protein-summary CSVs deposited on PRIDE for
Xie et al. 2026 (PNAS, doi:10.1073/pnas.2608150123): purified recombinant
rice proteins (MDHAR3, MDHAR4, MDHAR5, ALDP, FBA1, FBA3, PFP, TKT, RPI,
TAL) treated in vitro and assayed by the same NM-biotin + DTT persulfidation
workflow as the proteome-wide PXD072089 experiment, then identified by
PEAKS DB.

**Provenance class**: these are SAME-PAPER, SAME-LAB validations — an
orthogonal in-vitro assay (single purified protein) rather than an
independent-lab replication. They do NOT count toward Gate 2 condition 1
or condition 3 as an independent unit; they serve only as a same-paper
internal-consistency probe for the rice cross-species transfer track
(``scripts/validate_cross_species_rice.py``).

**File format quirks (reverse-engineered from the raw CSV, no header
documentation shipped)**:

* Each CSV is a PEAKS protein-summary export: a "protein" row (non-empty
  numeric ID in column 0) is followed by that protein's supporting
  "peptide" rows, until the next protein row.
* The *target* recombinant protein is not always the top (highest-score)
  protein — human keratin contaminants often outrank it. The correct
  block is located by matching the protein AC against a known label
  (e.g. ``MDAR4``, ``transketolase``) rather than by rank.
* The peptide row's ``Positions`` column (format ``"<n>,<before>,<after>"``)
  gives the position of the residue *preceding* the peptide, not the
  peptide's own start position — the peptide's first residue is at
  ``n + 1``. This was verified by aligning claimed peptide coordinates
  against the resolved UniProt sequence; using the field verbatim
  silently shifts every position by one residue.
* The ``Modification`` column marks the specific persulfidation-derived
  Carbamidomethyl[C] site as ``"<idx>,Carbamidomethyl[C]S;"`` (index
  within the peptide, 1-based); unmodified duplicate rows for the same
  peptide are also present (a background/lower-confidence identification)
  and are ignored.

Every extracted site is coordinate-verified against the rice reference
proteome using the resolved UniProt accession (identified by full-coverage
peptide-substring matching, not from the CSV's internal AC).
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

# Recombinant-protein label -> (internal AC used in the CSV, resolved
# UniProt accession). The UniProt accession was resolved by requiring
# every unique tryptic peptide from the protein's own block to occur as
# an exact substring of the candidate sequence (full coverage); ties
# between 100%-identical Japonica/Indica isoform accessions were broken
# by preferring a reviewed (sp|) entry or the more informative header.
RECOMBINANT_PROTEINS: dict[str, dict[str, str]] = {
    "MDHAR5": {"internal_ac": "MDAR5", "uniprot_accession": "Q84PW3"},
    "ALDP": {"internal_ac": "ALDP", "uniprot_accession": "Q40677"},
    "PFP": {"internal_ac": "PFP-ALPHA", "uniprot_accession": "Q6ZFT9"},
    "FBA1": {"internal_ac": "FBA1", "uniprot_accession": "P17784"},
    "MDHAR4": {"internal_ac": "MDAR4", "uniprot_accession": "Q6ZJ08"},
    "TKT": {"internal_ac": "transketolase", "uniprot_accession": ""},
    "RPI": {"internal_ac": "ribose-5-phosphate", "uniprot_accession": "Q7XVP0"},
    "FBA3": {"internal_ac": "FBA3", "uniprot_accession": "Q5N725"},
    "MDHAR3": {"internal_ac": "MDAR3", "uniprot_accession": "Q652L6"},
    "TAL": {"internal_ac": "Transaldolase", "uniprot_accession": "Q5JK10"},
}
# TKT/transketolase: no candidate covers 100% of its unique tryptic
# peptides (best coverage 60/84, ~71%) — the recombinant construct does
# not exactly match any single rice UniProt entry, so TKT is excluded
# from control registration (status=unmappable), not silently mapped to
# a partial-coverage guess.

_MOD_ENTRY_RE = re.compile(r"^(\d+),")
_POS_FIELD_RE = re.compile(r"^(\d+),")


@dataclass(frozen=True)
class RecombinantControlSite:
    gene: str
    uniprot_accession: str
    cys_position: int
    status: str  # "mapped" | "no_sites_detected" | "unmappable"


def _find_protein_block(
    rows: list[list[str]], internal_ac: str
) -> tuple[int, int] | None:
    """Return (start, end) row indices of the peptide rows following the
    protein row whose AC == internal_ac, or None if not found."""
    protein_row_indices = [
        i
        for i, r in enumerate(rows)
        if i >= 2 and r and r[0].strip().isdigit()
    ]
    for idx, row_idx in enumerate(protein_row_indices):
        if rows[row_idx][1] == internal_ac:
            start = row_idx + 1
            end = (
                protein_row_indices[idx + 1]
                if idx + 1 < len(protein_row_indices)
                else len(rows)
            )
            return start, end
    return None


def parse_recombinant_control_sites(
    csv_dir: Path,
    proteome: dict[str, str],
    proteins: dict[str, dict[str, str]] | None = None,
) -> tuple[RecombinantControlSite, ...]:
    """Parse PXD072300 PEAKS CSVs and return coordinate-verified
    persulfidation control sites.

    ``proteins`` defaults to the full :data:`RECOMBINANT_PROTEINS` registry;
    pass a smaller mapping to test a single gene in isolation. Genes whose
    recombinant construct's UniProt accession is unresolved (TKT) or whose
    block is not found are reported with ``status="unmappable"`` and no
    positions; every position on a resolved accession must independently
    verify as a cysteine at that coordinate.
    """
    sites: list[RecombinantControlSite] = []
    proteins = RECOMBINANT_PROTEINS if proteins is None else proteins

    for gene, info in proteins.items():
        uniprot_acc = info["uniprot_accession"]
        if not uniprot_acc:
            sites.append(
                RecombinantControlSite(
                    gene=gene,
                    uniprot_accession="",
                    cys_position=0,
                    status="unmappable",
                )
            )
            continue

        path = csv_dir / f"Persulfidation_{gene}.csv"
        if not path.is_file():
            sites.append(
                RecombinantControlSite(
                    gene=gene,
                    uniprot_accession=uniprot_acc,
                    cys_position=0,
                    status="unmappable",
                )
            )
            continue

        with path.open(encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.reader(handle))

        block = _find_protein_block(rows, info["internal_ac"])
        if block is None:
            sites.append(
                RecombinantControlSite(
                    gene=gene,
                    uniprot_accession=uniprot_acc,
                    cys_position=0,
                    status="unmappable",
                )
            )
            continue

        start, end = block
        sequence = proteome.get(uniprot_acc)
        if sequence is None:
            sites.append(
                RecombinantControlSite(
                    gene=gene,
                    uniprot_accession=uniprot_acc,
                    cys_position=0,
                    status="unmappable",
                )
            )
            continue

        verified_positions: set[int] = set()
        for row in rows[start:end]:
            if len(row) <= 11:
                continue
            pep = row[3]
            pos_field = row[11]
            mod_field = row[8] if len(row) > 8 else ""
            if not mod_field.strip():
                continue

            pos_match = _POS_FIELD_RE.match(pos_field)
            if not pos_match:
                continue
            # The field gives the position of the residue BEFORE the
            # peptide; the peptide itself starts one residue later.
            pep_start = int(pos_match.group(1)) + 1
            if sequence[pep_start - 1 : pep_start - 1 + len(pep)] != pep:
                continue  # claimed coordinate does not reproduce the peptide

            for mod_entry in mod_field.split(";"):
                mod_entry = mod_entry.strip()
                if not mod_entry or "Carbamidomethyl[C]" not in mod_entry:
                    continue
                mod_match = _MOD_ENTRY_RE.match(mod_entry)
                if not mod_match:
                    continue
                idx_in_pep = int(mod_match.group(1))
                if (
                    idx_in_pep < 1
                    or idx_in_pep > len(pep)
                    or pep[idx_in_pep - 1] != "C"
                ):
                    continue
                abs_pos = pep_start + idx_in_pep - 1
                if (
                    1 <= abs_pos <= len(sequence)
                    and sequence[abs_pos - 1] == "C"
                ):
                    verified_positions.add(abs_pos)

        if not verified_positions:
            # Accession and block both resolved; this experiment simply
            # detected zero persulfidated cysteines for this protein
            # (e.g. ALDP has no Cys-containing peptides at all in its
            # block) — distinct from an unresolved accession/block.
            sites.append(
                RecombinantControlSite(
                    gene=gene,
                    uniprot_accession=uniprot_acc,
                    cys_position=0,
                    status="no_sites_detected",
                )
            )
            continue

        for pos in sorted(verified_positions):
            sites.append(
                RecombinantControlSite(
                    gene=gene,
                    uniprot_accession=uniprot_acc,
                    cys_position=pos,
                    status="mapped",
                )
            )

    return tuple(sites)
