"""Co-peptide negative extraction from MaxQuant sites tables.

MaxQuant ``...Sites.txt`` output has one row per candidate localisation of a
modified residue. The probability-annotated column (e.g.
``Sulfide(C) Probabilities`` = ``VPSPTC(0.5)WC(0.5)SK``) carries the full
peptide sequence and every candidate position with its localisation
probability, and the ``Number of Sulfide(C)`` column carries the total
number of modified Cys in the detected peptide form. Together these make
the peptide's modification state fully determined (or provably
indeterminate), so the unmodified Cys of the same peptide are registered
co-peptide negatives (AGENTS.md explicit-negative class, detectability
matched by construction).

Classification follows the conservative rules of the methodology design:

- candidates with probability >= ``STRONG_PROBABILITY`` are the modified
  Cys (``positive``);
- candidates with lower probability are weak alternative localisations and
  are ``undetermined`` (never negative — "not localized" is not "not
  modified");
- Cys never appearing as candidates are ``negative`` only when the total
  modified count is determinate (either from the ``Number of ...`` column
  or, when every candidate is strong, from the candidate count) — a
  determinate count closes the possibility that a non-candidate Cys was
  modified but unobserved.

Peptide coordinates are re-localised by unique string match in the
registered reference proteome: MaxQuant coordinates live in its own search
database ID/coordinate space, which does not always agree with the
registered proteomes.
"""

from __future__ import annotations

import math
import re
from collections.abc import Sequence

from plantpersulf.evidence.copeptide_negatives import classify_peptide_cys

STRONG_PROBABILITY = 0.75

_PROBABILITY_ANNOTATION = re.compile(r"([A-Z])(?:\((\d+(?:\.\d+)?)\))?")


def parse_modified_sequence(modified_sequence: str) -> tuple[str, dict[int, float]]:
    """Parse a probability-annotated peptide like ``VPSPTC(0.5)WC(0.5)SK``.

    Returns ``(sequence, {in_peptide_position: probability})`` with 1-based
    positions; Cys without annotation carry no entry.
    """
    if not modified_sequence:
        raise ValueError("modified sequence must not be empty")
    sequence: list[str] = []
    candidates: dict[int, float] = {}
    covered_until = 0
    for match in _PROBABILITY_ANNOTATION.finditer(modified_sequence):
        if match.start() != covered_until:
            raise ValueError(
                f"probability annotation without residue in {modified_sequence!r}"
            )
        residue, probability = match.groups()
        sequence.append(residue)
        if probability is not None:
            candidates[len(sequence)] = float(probability)
        covered_until = match.end()
    if covered_until != len(modified_sequence):
        raise ValueError(f"unparsed residue in {modified_sequence!r}")
    return "".join(sequence), candidates


class PeptideEvidence:
    """One detected peptide, merged across every sites-table row sharing it."""

    def __init__(
        self,
        *,
        sequence: str,
        candidates: Sequence[tuple[int, float]],
        n_modified: int | None,
        source_rows: int,
    ) -> None:
        self.sequence = sequence
        self.candidates = tuple(candidates)
        self.n_modified = n_modified
        self.source_rows = source_rows

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, PeptideEvidence):
            return NotImplemented
        return (
            self.sequence,
            self.candidates,
            self.n_modified,
            self.source_rows,
        ) == (
            other.sequence,
            other.candidates,
            other.n_modified,
            other.source_rows,
        )

    def __repr__(self) -> str:
        return (
            f"PeptideEvidence(sequence={self.sequence!r}, "
            f"candidates={self.candidates!r}, n_modified={self.n_modified!r})"
        )


def _is_nan(value: object) -> bool:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return True
    return math.isnan(number)


def evidence_from_site_rows(
    rows: Sequence[dict[str, object]],
    *,
    probability_column: str,
    number_column: str,
) -> PeptideEvidence:
    """Merge sites-table rows of one peptide into a PeptideEvidence.

    Fail-closed on: inconsistent peptide sequences across rows,
    contradictory ``number_column`` values, and a determinate count that
    contradicts a fully-strong candidate set.
    """
    if not rows:
        raise ValueError("no rows to merge")
    sequence: str | None = None
    merged: dict[int, float] = {}
    numbers: list[int] = []
    for row in rows:
        row_sequence, row_candidates = parse_modified_sequence(
            str(row[probability_column])
        )
        if sequence is None:
            sequence = row_sequence
        elif row_sequence != sequence:
            raise RuntimeError(
                f"inconsistent peptide sequences across rows: "
                f"{row_sequence!r} vs {sequence!r}"
            )
        for position, probability in row_candidates.items():
            merged[position] = max(probability, merged.get(position, 0.0))
        raw_number = row.get(number_column)
        if raw_number is not None and not _is_nan(raw_number):
            numbers.append(int(float(raw_number)))  # type: ignore[arg-type]
    assert sequence is not None

    if numbers and any(number != numbers[0] for number in numbers):
        raise RuntimeError(f"contradictory {number_column} values: {numbers}")

    candidates = tuple(sorted(merged.items()))
    for position, _ in candidates:
        if sequence[position - 1] != "C":
            raise RuntimeError(
                f"candidate position {position} is not Cys in {sequence!r}"
            )
    all_strong = all(probability >= STRONG_PROBABILITY for _, probability in candidates)
    if numbers:
        n_modified = numbers[0]
        if all_strong and n_modified != len(candidates):
            raise RuntimeError(
                f"{number_column}={n_modified} contradicts {len(candidates)} "
                f"fully-strong candidates in {sequence!r}"
            )
    elif all_strong:
        n_modified = len(candidates)
    else:
        n_modified = None
    return PeptideEvidence(
        sequence=sequence,
        candidates=candidates,
        n_modified=n_modified,
        source_rows=len(rows),
    )


def classify_evidence(evidence: PeptideEvidence) -> dict[int, str]:
    """Three-state labels for every Cys of the peptide (see module docstring)."""
    strong = {
        position
        for position, probability in evidence.candidates
        if probability >= STRONG_PROBABILITY
    }
    weak = {
        position
        for position, probability in evidence.candidates
        if probability < STRONG_PROBABILITY
    }
    return classify_peptide_cys(
        peptide=evidence.sequence,
        modified_in_peptide_positions=tuple(sorted(strong)),
        localization_confirmed=bool(strong),
        site_determining_ions_confirmed=evidence.n_modified is not None,
        weak_modified_in_peptide_positions=tuple(sorted(weak)),
    )


def locate_peptide(peptide: str, protein_sequence: str) -> int | None:
    """1-based start of a UNIQUE occurrence of ``peptide``, else ``None``."""
    occurrences = protein_sequence.count(peptide)
    if occurrences != 1:
        return None
    return protein_sequence.index(peptide) + 1
