"""Missing-neutral scalar structure residual scoring."""

from __future__ import annotations


def score_structure_residual(
    weights: tuple[float, ...],
    intercept: float,
    structure_features: list[list[float]],
    quality_mask: list[bool],
) -> tuple[float, ...]:
    """Score available structure rows while forcing unavailable rows to zero."""
    if len(structure_features) != len(quality_mask):
        raise ValueError("structure features and quality mask length mismatch")
    if len(weights) > 3:
        raise ValueError("structure residual is limited to three features")
    if any(len(row) != len(weights) for row in structure_features):
        raise ValueError("structure feature width mismatch")
    return tuple(
        intercept
        + sum(
            weight * value
            for weight, value in zip(weights, row, strict=True)
        )
        if available
        else 0.0
        for row, available in zip(structure_features, quality_mask, strict=True)
    )
