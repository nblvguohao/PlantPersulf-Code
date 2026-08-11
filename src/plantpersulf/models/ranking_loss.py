"""Pairwise and non-negative PU objectives for site ranking."""

from __future__ import annotations

from collections import defaultdict

import torch

ALLOWED_LABELS = frozenset({"positive", "unlabeled"})


def _validate_labels(labels: list[str]) -> None:
    forbidden = set(labels) - ALLOWED_LABELS
    if forbidden:
        raise ValueError(f"forbidden label(s): {sorted(forbidden)}")


def within_protein_pairs(
    labels: list[str], protein_ids: list[str]
) -> tuple[tuple[int, int], ...]:
    """Return positive-to-unlabeled index pairs within each protein."""
    _validate_labels(labels)
    if len(labels) != len(protein_ids):
        raise ValueError("labels and protein_ids length mismatch")
    unlabeled_by_protein: dict[str, list[int]] = defaultdict(list)
    for index, (label, protein) in enumerate(zip(labels, protein_ids, strict=True)):
        if label == "unlabeled":
            unlabeled_by_protein[protein].append(index)
    return tuple(
        (positive_index, unlabeled_index)
        for positive_index, (label, protein) in enumerate(
            zip(labels, protein_ids, strict=True)
        )
        if label == "positive"
        for unlabeled_index in unlabeled_by_protein[protein]
    )


def nnpu_logistic_risk(
    logits: torch.Tensor, labels: list[str], class_prior: float
) -> torch.Tensor:
    """Return the non-negative PU logistic risk for a batch."""
    _validate_labels(labels)
    if logits.ndim != 1 or logits.numel() != len(labels):
        raise ValueError("logits must be one-dimensional and match labels")
    if not 0.0 < class_prior < 1.0:
        raise ValueError("class_prior must be in (0, 1)")
    positive = torch.tensor(
        [label == "positive" for label in labels],
        dtype=torch.bool,
        device=logits.device,
    )
    unlabeled = ~positive
    if not positive.any() or not unlabeled.any():
        raise ValueError("nnPU requires positive and unlabeled rows")
    pos_positive = torch.nn.functional.softplus(-logits[positive]).mean()
    pos_other = torch.nn.functional.softplus(logits[positive]).mean()
    unl_other = torch.nn.functional.softplus(logits[unlabeled]).mean()
    residual_risk = unl_other - class_prior * pos_other
    return class_prior * pos_positive + torch.clamp(residual_risk, min=0.0)


def combined_pu_ranking_loss(
    logits: torch.Tensor,
    labels: list[str],
    protein_ids: list[str],
    class_prior: float,
    pairwise_weight: float,
    pairs: tuple[tuple[int, int], ...] | None = None,
) -> torch.Tensor:
    """Combine nnPU risk with protein-local positive-to-unlabeled ranking."""
    if pairwise_weight < 0.0:
        raise ValueError("pairwise_weight must be non-negative")
    risk = nnpu_logistic_risk(logits, labels, class_prior)
    resolved_pairs = (
        within_protein_pairs(labels, protein_ids) if pairs is None else pairs
    )
    if not resolved_pairs:
        return risk
    pair_loss = torch.stack(
        [
            torch.nn.functional.softplus(-(logits[positive] - logits[unlabeled]))
            for positive, unlabeled in resolved_pairs
        ]
    ).mean()
    return risk + pairwise_weight * pair_loss
