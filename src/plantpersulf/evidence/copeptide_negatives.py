"""Co-peptide negative evidence — the only class of explicit negatives the
PU contract permits without re-observing the sample.

``AGENTS.md``: *"Do not treat an unobserved persulfidation site as an
experimental negative. The primary task is positive-unlabeled unless
explicit negative evidence is registered."* A co-peptide negative is the
strongest available explicit negative: within ONE detected peptide — one
spectrum, one digestion, one enrichment, one protein — some Cys carry the
persulfidation modification while others do not, so detectability is matched
by construction.

Two conditions are required before an unmodified Cys is labelled ``negative``
(rather than ``undetermined``), following the methodology design:

1. at least one Cys of the peptide has a confidently localized modification
   (localization probability above threshold / unambiguous published MS/MS
   localisation) — otherwise the whole peptide is ``undetermined``;
2. the negative position is covered by site-determining ions in the spectrum
   (or, for paper-registered entries, the published localisation of the
   co-modified sites on the same peptide is unambiguous) — otherwise the
   unmodified Cys is ``undetermined``, because "not localized" is not "not
   modified".

The registry ``data/registry/copeptide_negatives_v1.tsv`` holds one row per
Cys of each registered peptide with its three-state label and full
provenance.
"""

from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path

REGISTRY_COLUMNS = (
    "site_id",
    "species",
    "protein_accession",
    "cys_position",
    "state",
    "peptide_sequence",
    "peptide_start",
    "peptide_end",
    "in_peptide_position",
    "positive_positions",
    "negative_positions",
    "undetermined_positions",
    "localization_confirmed",
    "site_determining_ion_coverage",
    "source_doi",
    "source_detail",
    "provenance",
)

STATES = ("positive", "negative", "undetermined")


def classify_peptide_cys(
    *,
    peptide: str,
    modified_in_peptide_positions: tuple[int, ...],
    localization_confirmed: bool,
    site_determining_ions_confirmed: bool,
    weak_modified_in_peptide_positions: tuple[int, ...] = (),
) -> dict[int, str]:
    """Classify every Cys of one peptide into positive/negative/undetermined.

    ``modified_in_peptide_positions`` are 1-based positions within
    ``peptide``. Unmodified Cys become ``negative`` only when BOTH
    localization is confirmed on the peptide and site-determining ions cover
    the position; otherwise ``undetermined``.

    ``weak_modified_in_peptide_positions`` are candidate positions with only
    weak localisation evidence (e.g. probability < threshold in a MaxQuant
    sites table): they are never ``negative``, because "not confidently
    localised" is not "not modified".
    """
    if not peptide:
        raise ValueError("peptide must not be empty")
    cys_positions = {
        i + 1 for i, residue in enumerate(peptide) if residue == "C"
    }
    modified_set = set(modified_in_peptide_positions)
    weak_set = set(weak_modified_in_peptide_positions)
    bad = modified_set - cys_positions
    if bad:
        raise ValueError(f"modified position(s) not Cys in peptide: {sorted(bad)}")
    bad_weak = weak_set - cys_positions
    if bad_weak:
        raise ValueError(
            f"weak-modified position(s) not Cys in peptide: {sorted(bad_weak)}"
        )
    overlap = modified_set & weak_set
    if overlap:
        raise ValueError(
            f"position(s) both strong and weak modified: {sorted(overlap)}"
        )

    states: dict[int, str] = {}
    for position in sorted(cys_positions):
        if position in weak_set:
            states[position] = "undetermined"
        elif position in modified_set:
            states[position] = "positive" if localization_confirmed else "undetermined"
        else:
            states[position] = (
                "negative"
                if localization_confirmed and site_determining_ions_confirmed
                else "undetermined"
            )
    return states


def _parse_positions(field: str) -> tuple[int, ...]:
    return tuple(int(value) for value in field.split(";") if value)


def _format_positions(positions: Iterable[int]) -> str:
    return ";".join(str(position) for position in sorted(positions))


def write_copeptide_registry(path: Path, rows: list[dict[str, str]]) -> None:
    """Write the registry TSV (fail-closed on schema drift)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(REGISTRY_COLUMNS),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def load_copeptide_registry(path: Path) -> list[dict[str, str]]:
    """Load and validate the registry (unique ids, schema, allowed states)."""
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != REGISTRY_COLUMNS:
            raise RuntimeError(f"copeptide registry has invalid columns: {path}")
        rows = [dict(record) for record in reader]

    seen: set[str] = set()
    for record in rows:
        if record["site_id"] in seen:
            raise RuntimeError(
                f"duplicate site_id in copeptide registry: {record['site_id']}"
            )
        seen.add(record["site_id"])
        if record["state"] not in STATES:
            raise RuntimeError(
                f"invalid state {record['state']!r} in copeptide registry: "
                f"{record['site_id']}"
            )
        # peptide coordinates must be consistent with the reported positions
        start = int(record["peptide_start"])
        end = int(record["peptide_end"])
        position = int(record["cys_position"])
        in_peptide = int(record["in_peptide_position"])
        if not start <= position <= end:
            raise RuntimeError(
                f"cys_position outside peptide range: {record['site_id']}"
            )
        if position - start + 1 != in_peptide:
            raise RuntimeError(
                f"in_peptide_position inconsistent with peptide_start: "
                f"{record['site_id']}"
            )
    return rows
