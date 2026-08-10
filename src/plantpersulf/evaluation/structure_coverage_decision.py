"""Phase S0, Task 6 — fail-closed six-condition structure-coverage decision.

All six conditions must be met for ``STRUCTURE_SIGNAL_STABLE``; any single
failure yields ``STRUCTURE_SIGNAL_UNSTABLE``. Gate 2 remains ``GATE2_STOP``
regardless of the outcome — this is a development-level diagnostic, not a
gate override.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Literal

REQUIRED_CONDITIONS = (
    "full_exceeds_sequence_both_studies",
    "cluster_ci_excludes_zero_both_studies",
    "all_seed_directions_positive_both_studies",
    "full_exceeds_coverage_only_both_studies",
    "top_cluster_removed_gain_positive_both_studies",
    "all_audits_pass",
)


@dataclass(frozen=True)
class StructureSignalDecision:
    status: Literal[
        "STRUCTURE_SIGNAL_STABLE",
        "STRUCTURE_SIGNAL_UNSTABLE",
    ]
    conditions: dict[str, bool]
    failed_conditions: tuple[str, ...]
    gate2_status: Literal["GATE2_STOP"] = "GATE2_STOP"


def decide_structure_signal(
    conditions: Mapping[str, bool],
) -> StructureSignalDecision:
    """Return ``STRUCTURE_SIGNAL_STABLE`` only when every required condition
    is ``True``. Rejects unknown or missing conditions."""
    seen = set(conditions)
    required = set(REQUIRED_CONDITIONS)
    extra = seen - required
    if extra:
        raise ValueError(f"unknown condition(s): {sorted(extra)}")
    missing = required - seen
    if missing:
        raise ValueError(f"missing condition(s): {sorted(missing)}")

    failed = tuple(c for c in REQUIRED_CONDITIONS if not conditions[c])
    status: Literal["STRUCTURE_SIGNAL_STABLE", "STRUCTURE_SIGNAL_UNSTABLE"] = (
        "STRUCTURE_SIGNAL_STABLE" if not failed else "STRUCTURE_SIGNAL_UNSTABLE"
    )
    return StructureSignalDecision(
        status=status,
        conditions={c: conditions[c] for c in REQUIRED_CONDITIONS},
        failed_conditions=failed,
    )
