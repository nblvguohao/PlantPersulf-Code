"""Structure contrast on the explicit co-peptide negative pool.

The co-peptide negatives (``copeptide_negatives_v1.tsv``) are MS-confirmed
modified (positive) and unmodified (negative) Cys within the SAME peptide —
same digestion, same enrichment, same local sequence context, so detection
is matched by construction. The frozen sequence model cannot separate them
(9 peptides, 8/9 unseparated; information-theoretically expected for +/-10
windows). This module tests whether structure features can: within each
peptide, is the modified Cys cleanly above or below its unmodified siblings
on each structure feature, and does that hold beyond permutation chance?

Direction is intentionally left free: the known-control analysis (n=12, LOO)
validated burial (``contact_number`` high) for *functional* sites across
proteins, but the BRG3 RING co-peptide case shows the modified Cys there is
the *less* buried of the cluster (C206 exposed, negative C209 buried). The
within-peptide contrast is the direct probe of which direction holds when
detection is matched.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

FEATURE_NAMES = (
    "plddt",
    "rsa_relative",
    "cys_count_8a",
    "nearest_sg_distance",
    "positive_residue_count_6a",
    "coulomb_potential_sg",
    "contact_number_10a",
)


def peptide_contrast(
    features: Mapping[int, Mapping[str, float]],
    positive_positions: Sequence[int],
    negative_positions: Sequence[int],
    feature_names: Sequence[str] = FEATURE_NAMES,
) -> dict[str, str]:
    """Per-feature separation direction within one co-peptide group.

    A feature separates the group cleanly when ALL positive Cys are above ALL
    negative Cys (``pos_higher``) or below all of them (``pos_lower``);
    otherwise ``mixed``. Positions must all exist in ``features``.
    """
    for position in list(positive_positions) + list(negative_positions):
        if position not in features:
            raise ValueError(f"position {position} missing from features")
    result: dict[str, str] = {}
    for feature in feature_names:
        pos_values = [float(features[p][feature]) for p in positive_positions]
        neg_values = [float(features[p][feature]) for p in negative_positions]
        if min(pos_values) > max(neg_values):
            result[feature] = "pos_higher"
        elif max(pos_values) < min(neg_values):
            result[feature] = "pos_lower"
        else:
            result[feature] = "mixed"
    return result


def aggregate_contrast(
    groups: Sequence[dict[str, Any]],
    feature_names: Sequence[str] = FEATURE_NAMES,
) -> dict[str, dict[str, int]]:
    """Count separation directions across all co-peptide groups.

    Each group: ``{"features": ..., "positive_positions": ...,
    "negative_positions": ...}``.
    """
    counts: dict[str, dict[str, int]] = {}
    for feature in feature_names:
        counts[feature] = {"pos_higher": 0, "pos_lower": 0, "mixed": 0}
        for group in groups:
            direction = peptide_contrast(
                group["features"],
                group["positive_positions"],
                group["negative_positions"],
                (feature,),
            )[feature]
            counts[feature][direction] += 1
    return counts


def contrast_permutation_null(
    groups: Sequence[dict[str, Any]],
    feature: str,
    n_perm: int,
    rng: np.random.RandomState,
) -> dict[str, list[int]]:
    """Permutation null for the aggregate separation counts.

    Within each peptide, the POS/NEG labels are shuffled among that peptide's
    own Cys (preserving each group's k_pos / k_neg), and the number of
    peptides that separate ``pos_higher`` / ``pos_lower`` is recounted.
    Returns two sorted null distributions, one per direction.
    """
    null_higher: list[int] = []
    null_lower: list[int] = []
    for _ in range(n_perm):
        n_higher = 0
        n_lower = 0
        for group in groups:
            positions = list(
                set(group["positive_positions"]) | set(group["negative_positions"])
            )
            k_pos = len(group["positive_positions"])
            rng.shuffle(positions)
            pos_set = set(positions[:k_pos])
            contrast = peptide_contrast(
                group["features"],
                [p for p in positions if p in pos_set],
                [p for p in positions if p not in pos_set],
                (feature,),
            )[feature]
            n_higher += contrast == "pos_higher"
            n_lower += contrast == "pos_lower"
        null_higher.append(n_higher)
        null_lower.append(n_lower)
    return {"pos_higher": sorted(null_higher), "pos_lower": sorted(null_lower)}
