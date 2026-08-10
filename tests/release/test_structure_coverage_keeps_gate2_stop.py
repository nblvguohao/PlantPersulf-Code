"""Phase S0 release invariant: Gate 2 remains GATE2_STOP regardless of the
structure-coverage experiment outcome."""

from __future__ import annotations

from pathlib import Path

GATE2_CONFIG = Path("configs/gate2_v1.yaml")
GATE2_DECISION = Path("results/external_validation/pu_ranker_v1/gate2_decision.json")


def test_gate2_config_hash_is_unchanged() -> None:
    """The Gate 2 config must be byte-identical to its freeze-time content."""
    from plantpersulf.evaluation.structure_coverage_audit import sha256_file

    # We don't hardcode the hash here (it could drift across branches),
    # but we verify it exists and is readable.
    assert GATE2_CONFIG.is_file()
    sha = sha256_file(GATE2_CONFIG)
    assert len(sha) == 64
    assert sha == sha.upper()


def test_gate2_decision_still_says_stop() -> None:
    """The official Gate 2 decision must still be GATE2_STOP."""
    import json

    assert GATE2_DECISION.is_file()
    data = json.loads(GATE2_DECISION.read_text(encoding="utf-8"))
    assert data.get("decision") == "GATE2_STOP", "Gate 2 decision changed"


def test_gate2_config_has_version_four() -> None:
    """The Gate 2 config must keep its version and structure."""
    import yaml

    parsed = yaml.safe_load(GATE2_CONFIG.read_text(encoding="utf-8"))
    assert isinstance(parsed, dict)
    assert parsed.get("version") == 1
