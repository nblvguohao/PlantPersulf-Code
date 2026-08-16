"""Unit tests for the regime-routed MoE diagnostic helpers.

Covers ``evaluation.regime_moe`` — the pure-function layer behind the
regime-routing MoE design (``2026-08-15-regime-routing-moe-design.md``,
diagnostic A). Expert directions are written down as structural priors
(never learned from labeled data), so there is no LOO sign-learning bias:
the exchangeability permutation null is the honest reference.

Claim class ``diagnostic_only``; nothing here trains a model or touches a
frozen artifact.
"""

from __future__ import annotations

import numpy as np
import pytest

from plantpersulf.evaluation.regime_moe import (
    EXPERT_DIRECTIONS,
    moe_route_scores,
    routed_burden,
    routed_permutation_null,
)


def _rows(*values: tuple[float, float, float, float]) -> dict[int, dict[str, float]]:
    """Build feature rows (plddt, contact, rsa, nearest_sg) keyed by 1..n."""
    return {
        i + 1: {
            "plddt": p,
            "contact_number_10a": c,
            "rsa_relative": r,
            "nearest_sg_distance": s,
        }
        for i, (p, c, r, s) in enumerate(values)
    }


def test_expert_directions_prior_shape() -> None:
    """Three regime experts, each with the structural-prior direction."""
    assert set(EXPERT_DIRECTIONS) == {"folded", "linker", "disordered"}
    # folded: burial/contact is the significant feature (§9.4).
    assert EXPERT_DIRECTIONS["folded"]["contact_number_10a"] == 1.0
    # disordered: exposure/isolation prior (RSA high, nearest-Sγ far).
    assert EXPERT_DIRECTIONS["disordered"]["rsa_relative"] == 1.0
    assert EXPERT_DIRECTIONS["disordered"]["nearest_sg_distance"] == 1.0


def test_moe_route_scores_route_by_plddt() -> None:
    """A high-pLDDT (folded) site is scored with the folded prior; a low-pLDDT
    (disordered) site with the disordered prior."""
    rows = _rows(
        (92.0, 8.0, 0.3, 4.0),  # 1 folded, buried (contact high)
        (80.0, 6.0, 0.4, 5.0),  # 2 folded, less buried
        (30.0, 2.0, 0.9, 9.0),  # 3 disordered, exposed (RSA high, Sγ far)
        (25.0, 3.0, 0.2, 2.0),  # 4 disordered, buried-ish (RSA low, Sγ near)
    )
    scores = moe_route_scores(rows, EXPERT_DIRECTIONS)
    # Folded expert: higher contact → higher score.
    assert scores[1] > scores[2]
    # Disordered expert: higher RSA + farther Sγ → higher score.
    assert scores[3] > scores[4]


def test_moe_route_scores_z_scope_is_regime_local() -> None:
    """z-scores are computed within each expert's own regime bucket, so a
    site's score is unaffected by sites in a different regime."""
    rows = _rows(
        (92.0, 8.0, 0.3, 4.0),
        (80.0, 6.0, 0.4, 5.0),
        (30.0, 2.0, 0.9, 9.0),
    )
    base = moe_route_scores(rows, EXPERT_DIRECTIONS)
    # Add another folded site with extreme contact; the disordered site's
    # score must not change (its z-pool is regime-local).
    rows2 = _rows(
        (92.0, 8.0, 0.3, 4.0),
        (80.0, 6.0, 0.4, 5.0),
        (99.0, 50.0, 0.1, 6.0),
        (30.0, 2.0, 0.9, 9.0),
    )
    after = moe_route_scores(rows2, EXPERT_DIRECTIONS)
    assert after[4] == pytest.approx(base[3])


def test_routed_burden_counts_ranks_of_true_sites() -> None:
    records = [
        (
            _rows((92.0, 9.0, 0.1, 4.0), (80.0, 5.0, 0.4, 5.0)),
            1,  # true site = folded, most buried → rank 1
        ),
        (
            _rows((30.0, 1.0, 0.9, 9.0), (25.0, 2.0, 0.1, 2.0)),
            1,  # true site = disordered, most exposed → rank 1
        ),
    ]
    out = routed_burden(records, EXPERT_DIRECTIONS, "regime")
    assert out["n"] == 2
    assert out["top1"] == 2
    assert out["total_first_hit_burden"] == 2


def test_routed_burden_rejects_bad_scope() -> None:
    with pytest.raises(ValueError):
        routed_burden(
            [(_rows((90.0, 1.0, 0.1, 3.0)), 1)],
            EXPERT_DIRECTIONS,
            "nonsense",
        )


def test_routed_permutation_null_returns_sorted_counts() -> None:
    records = [
        (
            _rows((92.0, 9.0, 0.1, 4.0), (80.0, 5.0, 0.4, 5.0)),
            1,
        ),
        (
            _rows((30.0, 1.0, 0.9, 9.0), (25.0, 2.0, 0.1, 2.0)),
            1,
        ),
    ]
    rng = np.random.RandomState(20260815)
    null = routed_permutation_null(records, EXPERT_DIRECTIONS, "regime", 50, rng)
    assert len(null) == 50
    assert null == sorted(null)
    # Two 2-Cys proteins: each random true site ranks 1 or 2, total burden
    # therefore in [2, 4].
    assert min(null) >= 2 and max(null) <= 4
