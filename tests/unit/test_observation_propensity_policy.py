"""Fail-closed policy tests for observation-propensity correction."""

from pathlib import Path

import pytest

from plantpersulf.models.observation_propensity import resolve_observation_propensity


def test_missing_table_disables_propensity_without_substitution() -> None:
    status = resolve_observation_propensity(None, None)

    assert status.enabled is False
    assert status.reason == "complete_discovery_table_not_supplied"


def test_unregistered_table_is_rejected(tmp_path: Path) -> None:
    marker = tmp_path / "policy-marker.tsv"
    marker.write_text("policy-marker\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="registry"):
        resolve_observation_propensity(marker, None)
