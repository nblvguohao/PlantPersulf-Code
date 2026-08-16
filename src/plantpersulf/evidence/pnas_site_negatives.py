"""Co-peptide negative extraction from a PNAS persulfidome supplement.

PXD072089 (Huang et al., *PNAS* 2025, ``10.1073/pnas.2608150123``) deposited
two complementary tables in its supplementary information:

- **Dataset S4** ("Persulfidated sites with corresponding domain"): one row
  per paper-confirmed ``-SSH`` site — peptide sequence, protein IDs, leading
  razor protein, and the modified Cys position (protein coordinates). A
  peptide with several modified Cys appears once per site, so the rows of a
  peptide enumerate its paper-confirmed modification set.
- **Dataset S1** ("Cys-SSH sites identified in rice..."): the MaxQuant
  peptide inventory — for every detected persulfidated peptide, the full Cys
  positions and total Cys count. ``Position`` lists **all** Cys of the
  peptide, not just the modified ones.

A co-peptide negative arises where a detected peptide carries more Cys than
the paper confirms modified: the confirmed positions are ``positive`` and
the remaining Cys of the same detected peptide are ``negative``, exactly as
the registered kiad070 tomato entries (paper-reported unambiguous
localisation on the same detected peptide, raw spectra not locally
available). ``classify_peptide_cys`` (copeptide_negatives module) is reused
for the three-state label.

This module is pure logic (no spreadsheet I/O); the scanner script reads the
xlsx supplements and builds :class:`Site` records. Coordinates: S4/S1 protein
positions live in the paper's own search space, which is the registered
UniProt rice space for 4 of the 5 partial peptides but off by one residue for
the 5th. Mapping is therefore done through the peptide's own Cys *rank*
(the in-peptide position of a modified Cys is invariant under coordinate
shifts), and the peptide is re-localised by unique string match in the
registered rice proteome for canonical coordinates.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import TypedDict


class Site(TypedDict):
    """One Dataset S4 row (paper-confirmed modified Cys)."""

    peptide: str
    protein_ids: list[str]
    leading_razor: str
    modified_positions: tuple[int, ...]
    domain: str


class CoPeptideGroup(TypedDict):
    """A peptide with partial modification -> co-peptide negative candidates."""

    peptide: str
    leading_razor: str
    protein_ids: list[str]
    modified_positions: tuple[int, ...]
    domains: set[str]
    cys_in_peptide: tuple[int, ...]
    modified_in_peptide: tuple[int, ...]
    negative_in_peptide: tuple[int, ...]


def parse_position_field(value: object) -> tuple[int, ...]:
    """Parse ``"413,417"`` / ``"413"`` / ``"144,146,153,157,160"`` -> ints."""
    if value is None:
        return ()
    text = str(value).strip()
    if not text:
        return ()
    parts = text.replace(",", ";").split(";")
    return tuple(int(part) for part in parts if part.strip())


def cys_in_peptide_positions(peptide: str) -> tuple[int, ...]:
    """1-based in-peptide positions of every Cys of ``peptide``."""
    return tuple(i + 1 for i, residue in enumerate(peptide) if residue == "C")


def map_modified_positions(
    modified_paper_positions: Sequence[int],
    inventory_positions: Sequence[int],
    cys_in_pep: Sequence[int],
    peptide: str,
) -> tuple[int, ...]:
    """Map paper-space modified positions to in-peptide positions.

    ``inventory_positions`` is the paper's full Cys inventory (all Cys of the
    peptide, protein coordinates, ascending). The in-peptide position of a
    modified Cys is found by its *rank* within the inventory — invariant
    under the coordinate shift between the paper's search space and the
    registered proteome. Fail-closed when a modified position is not in the
    inventory (the paper's coordinates are then inconsistent with S1).
    """
    ordered = sorted(inventory_positions)
    cys_list = sorted(cys_in_pep)
    if len(ordered) != len(cys_list):
        raise ValueError(
            f"inventory {ordered} and peptide Cys {cys_list} disagree "
            f"for {peptide!r}"
        )
    position_by_paper: dict[int, int] = {}
    for rank, paper_position in enumerate(ordered):
        position_by_paper[paper_position] = cys_list[rank]
    missing = [
        paper for paper in modified_paper_positions if paper not in position_by_paper
    ]
    if missing:
        raise ValueError(
            f"modified position(s) {missing} not in peptide inventory "
            f"{ordered} for {peptide!r}"
        )
    return tuple(position_by_paper[paper] for paper in modified_paper_positions)


def co_peptide_groups(
    sites: Sequence[Site],
    inventory: Mapping[tuple[str, str], tuple[int, ...]],
) -> list[CoPeptideGroup]:
    """Group sites into co-peptide candidates with partial modification.

    Only peptides with >= 2 Cys and fewer paper-confirmed modified Cys than
    total Cys produce a negative; fully-modified peptides are skipped (no
    co-peptide contrast). Each returned group carries the mapped in-peptide
    modified positions.
    """
    merged: dict[tuple[str, str], CoPeptideGroup] = {}
    position_sets: dict[tuple[str, str], set[int]] = {}
    domain_sets: dict[tuple[str, str], set[str]] = {}
    order: list[tuple[str, str]] = []
    for site in sites:
        key = (site["peptide"], site["leading_razor"])
        if key not in merged:
            merged[key] = {
                "peptide": site["peptide"],
                "leading_razor": site["leading_razor"],
                "protein_ids": list(site["protein_ids"]),
                "modified_positions": (),
                "domains": set(),
                "cys_in_peptide": (),
                "modified_in_peptide": (),
                "negative_in_peptide": (),
            }
            position_sets[key] = set()
            domain_sets[key] = set()
            order.append(key)
        position_sets[key].update(site["modified_positions"])
        if site["domain"]:
            domain_sets[key].add(site["domain"])

    results: list[CoPeptideGroup] = []
    for key in order:
        group = merged[key]
        group["modified_positions"] = tuple(sorted(position_sets[key]))
        group["domains"] = domain_sets[key]
        peptide = group["peptide"]
        cys = cys_in_peptide_positions(peptide)
        if len(cys) < 2:
            continue
        if len(group["modified_positions"]) >= len(cys):
            continue  # fully modified: no co-peptide negative
        inventory_positions = inventory.get(key)
        if inventory_positions is None:
            raise RuntimeError(
                f"no Dataset S1 inventory row for {key!r}"
            )
        group["cys_in_peptide"] = cys
        group["modified_in_peptide"] = map_modified_positions(
            group["modified_positions"], inventory_positions, cys, peptide
        )
        group["negative_in_peptide"] = tuple(
            position
            for position in cys
            if position not in group["modified_in_peptide"]
        )
        results.append(group)
    return results
