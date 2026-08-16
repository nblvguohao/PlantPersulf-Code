"""Fail-closed observation-propensity policy for future SAR-PU work."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file


@dataclass(frozen=True)
class ObservationPropensityStatus:
    """Whether a registered discovery-observation table permits correction."""

    enabled: bool
    reason: str
    source_sha256: str | None


def resolve_observation_propensity(
    table_path: Path | None, registry_path: Path | None
) -> ObservationPropensityStatus:
    """Return a fail-closed status without substituting an observation table."""
    if table_path is None:
        return ObservationPropensityStatus(
            enabled=False,
            reason="complete_discovery_table_not_supplied",
            source_sha256=None,
        )
    if registry_path is None:
        raise RuntimeError("registry is required for an observation table")
    assert_registered_input(table_path, registry_path)
    return ObservationPropensityStatus(
        enabled=False,
        reason="observation_propensity_fitting_not_implemented",
        source_sha256=hash_file(table_path, "sha256"),
    )
