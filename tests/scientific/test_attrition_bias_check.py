"""Tests for the PXD072089 attrition compositional-bias check.

Self-review-triggered (2026-07-23): does the ~55% UniProt-accession-deletion
attrition in the rice persulfidome dataset systematically remove or retain
redox/oxidoreductase-related proteins, in a way that could confound the
family-level (``cross_species_conservation.py``) or structural-level
(``cross_species_structural_context.py``) convergence findings?
"""

from __future__ import annotations

import math

from plantpersulf.evaluation.attrition_bias_check import (
    REDOX_KEYWORDS,
    attrition_bias_test,
    fisher_exact_two_sided,
    is_redox_related,
)


def test_is_redox_related_matches_known_keyword() -> None:
    assert is_redox_related(
        "Q6YZ09_ORYSJ NADH:ubiquinone reductase (non-electrogenic) "
        "OS=Oryza sativa subsp. japonica OX=39947 GN=Os08g0141400 PE=1 SV=1"
    )
    assert is_redox_related(
        "Q8S5T1_ORYSJ Glutathione reductase OS=Oryza sativa subsp. "
        "japonica OX=39947 GN=LOC_Os03g06740 PE=1 SV=1"
    )


def test_is_redox_related_rejects_unrelated_description() -> None:
    assert not is_redox_related(
        "A0A067Y2H7_ORYSJ Lectin OS=Oryza sativa subsp. japonica OX=39947 PE=2 SV=1"
    )
    assert not is_redox_related(
        "A0A0K2SR05_ORYSJ 50S ribosomal protein L16, chloroplastic (Fragment) "
        "OS=Oryza sativa subsp. japonica OX=39947 GN=rpl16 PE=3 SV=1"
    )


def test_is_redox_related_only_checks_functional_text_before_os() -> None:
    # "OX=" and organism text must not accidentally match a keyword.
    description = "Some protein OS=Oryza sativa OX=39947 GN=REDUCTASE_LIKE_LOCUS"
    assert not is_redox_related(description)


def test_redox_keywords_are_uppercase_and_nonempty() -> None:
    assert len(REDOX_KEYWORDS) > 0
    for kw in REDOX_KEYWORDS:
        assert kw == kw.upper()
        assert kw.strip() == kw


def test_fisher_exact_two_sided_symmetric_table_is_significant() -> None:
    # Strong, obviously non-null association: [[90, 10], [10, 90]]
    p = fisher_exact_two_sided(90, 10, 10, 90)
    assert p < 0.001


def test_fisher_exact_two_sided_balanced_table_is_not_significant() -> None:
    # Proportions nearly identical between groups -> large p-value.
    p = fisher_exact_two_sided(50, 50, 51, 49)
    assert p > 0.5


def test_fisher_exact_two_sided_matches_hand_computed_small_table() -> None:
    # Table [[1, 9], [11, 3]]: margins are (10, 14) x (12, 12), N=24.
    # Two-sided p = sum of pmf(k) for k in {0, 1, 9, 10} (the four tables
    # whose hypergeometric probability is <= pmf(observed=1)), hand-summed
    # via C(12,k)*C(12,10-k)/C(24,10): 66 + 2640 + 2640 + 66 = 5412,
    # divided by C(24,10) = 1961256 -> 0.00275920...
    p = fisher_exact_two_sided(1, 9, 11, 3)
    assert math.isclose(p, 5412 / 1961256, rel_tol=1e-6)


def test_fisher_exact_two_sided_is_symmetric_under_table_transpose() -> None:
    p1 = fisher_exact_two_sided(3, 7, 8, 2)
    p2 = fisher_exact_two_sided(8, 2, 3, 7)
    assert math.isclose(p1, p2, rel_tol=1e-9)


def test_fisher_exact_two_sided_p_value_never_exceeds_one() -> None:
    p = fisher_exact_two_sided(5, 5, 5, 5)
    assert 0.0 <= p <= 1.0


def test_attrition_bias_test_detects_no_difference() -> None:
    kept = ["Glutathione reductase OS=X"] * 20 + ["Lectin OS=X"] * 80
    dropped = ["Glutathione reductase OS=X"] * 20 + ["Lectin OS=X"] * 80
    result = attrition_bias_test(kept, dropped)
    assert result.kept_redox_fraction == result.dropped_redox_fraction
    assert result.p_value > 0.9


def test_attrition_bias_test_detects_strong_differential_loss() -> None:
    # Kept set is redox-enriched relative to dropped set.
    kept = ["NADH dehydrogenase OS=X"] * 80 + ["Lectin OS=X"] * 20
    dropped = ["NADH dehydrogenase OS=X"] * 10 + ["Lectin OS=X"] * 90
    result = attrition_bias_test(kept, dropped)
    assert result.kept_redox_fraction > result.dropped_redox_fraction
    assert result.p_value < 0.001
    assert result.odds_ratio is not None
    assert result.odds_ratio > 1.0


def test_attrition_bias_test_reports_raw_counts_not_just_fractions() -> None:
    kept = ["Glutathione reductase OS=X", "Lectin OS=X"]
    dropped = ["Lectin OS=X"]
    result = attrition_bias_test(kept, dropped)
    assert result.kept_redox == 1
    assert result.kept_non_redox == 1
    assert result.dropped_redox == 0
    assert result.dropped_non_redox == 1


def test_attrition_bias_test_handles_empty_dropped_group() -> None:
    kept = ["Glutathione reductase OS=X", "Lectin OS=X"]
    dropped: list[str] = []
    result = attrition_bias_test(kept, dropped)
    assert result.dropped_redox_fraction == 0.0
    assert result.p_value == 1.0
