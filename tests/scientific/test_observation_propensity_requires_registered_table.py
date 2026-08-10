"""Scientific guard that observation correction cannot silently substitute data."""

from plantpersulf.models.observation_propensity import resolve_observation_propensity


def test_missing_observation_table_has_no_source_checksum() -> None:
    status = resolve_observation_propensity(None, None)

    assert status.enabled is False
    assert status.source_sha256 is None
