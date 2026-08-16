"""Shortcut-free biological context features for cysteine sites.

This module intentionally exposes only local sequence-derived biology. It
does not use protein length, structure confidence, study identity, or any
other dataset-specific shortcut.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass

from plantpersulf.features.sequence import (
    SequenceFeatureRow,
    _flanking_window,
    _hydrophobicity,
    _local_positive_charge_density,
)

BIOLOGY_FEATURE_NAMES: tuple[str, ...] = (
    "hydrophobicity",
    "protein_cys_density",
    "local_positive_charge_density",
    "local_negative_charge_density",
    "local_cys_density",
    "local_sequence_entropy",
)


@dataclass(frozen=True)
class SiteBiologyVector:
    """Named six-dimensional biological context for one cysteine site."""

    names: tuple[str, ...]
    values: tuple[float, ...]


def _density(window: str, residues: frozenset[str]) -> float:
    if not window:
        return 0.0
    return sum(residue in residues for residue in window) / len(window)


def _entropy(window: str) -> float:
    residues = [residue for residue in window if residue != "X"]
    if not residues:
        return 0.0
    counts = Counter(residues)
    length = len(residues)
    return -sum(
        (count / length) * math.log2(count / length) for count in counts.values()
    )


def build_site_biology_vector(row: SequenceFeatureRow) -> SiteBiologyVector:
    """Build the six core biology features from an existing sequence row."""
    values = (
        row.hydrophobicity,
        row.cys_density,
        row.local_positive_charge_density,
        _density(row.flanking_window, frozenset({"D", "E"})),
        _density(row.flanking_window, frozenset({"C"})),
        _entropy(row.flanking_window),
    )
    return SiteBiologyVector(names=BIOLOGY_FEATURE_NAMES, values=values)


def build_site_biology_vector_from_sequence(
    accession: str,
    position: int,
    sequence: str,
    radius: int = 10,
) -> SiteBiologyVector:
    """Build site biology directly from a cysteine coordinate and sequence."""
    if radius < 0:
        raise ValueError("radius must be non-negative")
    if position < 1 or position > len(sequence):
        raise ValueError("position must be within the sequence")
    if sequence[position - 1] != "C":
        raise ValueError("position must identify a cysteine")

    window = _flanking_window(sequence, position, radius)
    values = (
        _hydrophobicity(window),
        sequence.count("C") / len(sequence),
        _local_positive_charge_density(window),
        _density(window, frozenset({"D", "E"})),
        _density(window, frozenset({"C"})),
        _entropy(window),
    )
    return SiteBiologyVector(names=BIOLOGY_FEATURE_NAMES, values=values)
