"""Pairwise and non-negative PU objectives for site ranking."""

from __future__ import annotations

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
    return tuple(
        (positive_index, unlabeled_index)
        for positive_index, positive_label in enumerate(labels)
        for unlabeled_index, unlabeled_label in enumerate(labels)
        if positive_label == "positive"
        and unlabeled_label == "unlabeled"
        and protein_ids[positive_index] == protein_ids[unlabeled_index]
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
) -> torch.Tensor:
    """Combine nnPU risk with protein-local positive-to-unlabeled ranking."""
    if pairwise_weight < 0.0:
        raise ValueError("pairwise_weight must be non-negative")
    risk = nnpu_logistic_risk(logits, labels, class_prior)
    pairs = within_protein_pairs(labels, protein_ids)
    if not pairs:
        return risk
    pair_loss = torch.stack(
        [
            torch.nn.functional.softplus(-(logits[positive] - logits[unlabeled]))
            for positive, unlabeled in pairs
        ]
    ).mean()
    return risk + pairwise_weight * pair_loss
