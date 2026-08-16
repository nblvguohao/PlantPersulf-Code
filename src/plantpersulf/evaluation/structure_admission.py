"""Evidence gate for admitting an optional structure residual arm."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class StructureAdmissionDecision:
    admitted: bool
    delta: int
    reason: str


def admit_structure_residual(
    paired_structure_deltas: tuple[float, ...],
    paired_coverage_only_deltas: tuple[float, ...],
) -> StructureAdmissionDecision:
    """Admit only stable gains that cannot be explained by coverage alone."""
    if len(paired_structure_deltas) < 5:
        return StructureAdmissionDecision(False, 0, "insufficient_repeated_folds")
    if len(paired_structure_deltas) != len(paired_coverage_only_deltas):
        raise ValueError("paired delta lengths must match")
    structure_stable = min(paired_structure_deltas) > 0.0
    coverage_has_gain = max(paired_coverage_only_deltas) > 0.0
    if structure_stable and not coverage_has_gain:
        return StructureAdmissionDecision(True, 1, "matched_structure_gain")
    return StructureAdmissionDecision(False, 0, "coverage_or_instability_failure")
