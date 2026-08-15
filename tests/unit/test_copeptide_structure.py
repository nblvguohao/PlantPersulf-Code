"""Unit tests for copeptide_structure.py.

Tests the per-peptide structural contrast between MS-confirmed modified
(POS) and unmodified (NEG) Cys in the same peptide, plus the permutation
null that keeps each peptide's composition fixed.
"""

from __future__ import annotations

import numpy as np
import pytest

from plantpersulf.evaluation.copeptide_structure import (
    FEATURE_NAMES,
    aggregate_contrast,
    contrast_permutation_null,
    peptide_contrast,
)

_CONTACT = "contact_number_10a"


def _features(*contact_values):
    """Build {position: {contact_number_10a: value}} from positional args."""
    return {
        pos: {_CONTACT: float(value)}
        for pos, value in enumerate(contact_values, start=206)
    }


def _group(features, positive_positions, negative_positions):
    return {
        "features": features,
        "positive_positions": positive_positions,
        "negative_positions": negative_positions,
    }


class TestPeptideContrast:
    def test_pos_higher_when_all_positive_above_all_negative(self):
        features = _features(30.0, 20.0, 28.0)  # 206, 207, 208
        result = peptide_contrast(features, [206, 208], [207], (_CONTACT,))
        assert result == {_CONTACT: "pos_higher"}

    def test_pos_lower_when_all_positive_below_all_negative(self):
        features = _features(10.0, 30.0, 12.0)
        result = peptide_contrast(features, [206, 208], [207], (_CONTACT,))
        assert result == {_CONTACT: "pos_lower"}

    def test_mixed_when_sets_overlap(self):
        features = _features(5.0, 10.0, 15.0)
        result = peptide_contrast(features, [206, 208], [207], (_CONTACT,))
        assert result == {_CONTACT: "mixed"}

    def test_single_positive_single_negative(self):
        features = _features(22.0, 26.0)
        result = peptide_contrast(features, [206], [207], (_CONTACT,))
        assert result == {_CONTACT: "pos_lower"}

    def test_missing_position_raises(self):
        features = _features(22.0)
        with pytest.raises(ValueError, match="207 missing from features"):
            peptide_contrast(features, [206, 207], [], (_CONTACT,))

    def test_default_feature_names_full_set(self):
        features = {p: {f: 1.0 for f in FEATURE_NAMES} for p in (206, 207)}
        result = peptide_contrast(features, [206], [207])
        assert set(result) == set(FEATURE_NAMES)
        # all features identical -> no clean separation either direction
        assert set(result.values()) == {"mixed"}


class TestAggregateContrast:
    def test_counts_directions_across_groups(self):
        groups = [
            _group(_features(30.0, 20.0, 28.0), [206, 208], [207]),  # pos_higher
            _group(_features(10.0, 30.0, 12.0), [206, 208], [207]),  # pos_lower
            _group(_features(5.0, 10.0, 15.0), [206, 208], [207]),  # mixed
        ]
        result = aggregate_contrast(groups, (_CONTACT,))
        assert result[_CONTACT] == {"pos_higher": 1, "pos_lower": 1, "mixed": 1}

    def test_all_groups_have_all_direction_keys(self):
        features = {p: {f: 1.0 for f in FEATURE_NAMES} for p in (206, 207)}
        groups = [_group(features, [206], [207])]
        result = aggregate_contrast(groups)
        for feature in FEATURE_NAMES:
            assert set(result[feature]) == {"pos_higher", "pos_lower", "mixed"}


class TestPermutationNull:
    def test_deterministic_for_fixed_seed(self):
        groups = [
            _group(_features(30.0, 20.0, 28.0), [206, 208], [207]),
            _group(_features(10.0, 30.0, 12.0), [206, 208], [207]),
        ]
        rng = np.random.RandomState(7)
        first = contrast_permutation_null(groups, _CONTACT, 20, rng)
        rng = np.random.RandomState(7)
        second = contrast_permutation_null(groups, _CONTACT, 20, rng)
        assert first == second

    def test_null_counts_bounded_by_n_groups(self):
        groups = [
            _group(_features(30.0, 20.0), [206], [207]) for _ in range(5)
        ]
        null = contrast_permutation_null(
            groups, _CONTACT, 50, np.random.RandomState(0)
        )
        for direction in ("pos_higher", "pos_lower"):
            assert max(null[direction]) <= len(groups)
            assert min(null[direction]) >= 0
            assert null[direction] == sorted(null[direction])

    def test_null_changes_when_seed_changes(self):
        groups = [_group(_features(30.0, 20.0), [206], [207])]
        a = contrast_permutation_null(groups, _CONTACT, 100, np.random.RandomState(1))
        b = contrast_permutation_null(groups, _CONTACT, 100, np.random.RandomState(2))
        assert a != b

    def test_observed_direction_counts_dominated_in_constant_group(self):
        # A group with equal feature values can never separate cleanly.
        groups = [_group(_features(1.0, 1.0), [206], [207])]
        null = contrast_permutation_null(
            groups, _CONTACT, 30, np.random.RandomState(3)
        )
        assert null["pos_higher"] == [0] * 30
        assert null["pos_lower"] == [0] * 30
