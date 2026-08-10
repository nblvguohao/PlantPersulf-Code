"""Numerical contract tests for the pairwise nnPU objective."""

import pytest
import torch

from plantpersulf.models.ranking_loss import (
    combined_pu_ranking_loss,
    within_protein_pairs,
)


def test_pairs_are_within_protein_and_positive_to_unlabeled() -> None:
    pairs = within_protein_pairs(
        ["positive", "unlabeled", "unlabeled", "positive"],
        ["group-a", "group-a", "group-b", "group-b"],
    )
    assert pairs == ((0, 1), (3, 2))


def test_forbidden_label_is_rejected() -> None:
    with pytest.raises(ValueError, match="forbidden"):
        combined_pu_ranking_loss(
            torch.tensor([0.0, 1.0]),
            ["positive", "negative"],
            ["group-a", "group-a"],
            class_prior=0.5,
            pairwise_weight=1.0,
        )
