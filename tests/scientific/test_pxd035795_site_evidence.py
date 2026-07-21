"""RED (Phase A2, Route M): reviewed method mapping gates PXD035795 site_ms.

This is the first executable Phase A2 cycle from
``docs/superpowers/specs/2026-07-21-evidence-site-acquisition-design.md``.

It asserts a *gating contract*, not a biological positive: an mzIdentML
candidate modification may be upgraded to ``site_ms`` ONLY when a reviewed
method-mapping rule — itself backed by a registered, hash-audited method source
— exactly matches the candidate and every evidence-policy condition holds.
Absent the rule (or the method source), the candidate stays ``unresolved``.

Which chemical tag / mass on cysteine actually proves persulfidation is a
reviewed scientific decision encoded in the mapping + method source; this test
deliberately does NOT hardcode that biology. It only proves the join is
fail-closed and residue-exact.

Expected RED: ``plantpersulf.evidence.site_evidence`` does not exist yet, so
collection fails with ImportError. GREEN implements the minimal join.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from plantpersulf.evidence.mzidentml import parse_mzidentml
from plantpersulf.evidence.site_evidence import (  # RED: module missing
    apply_reviewed_site_mappings,
    load_reviewed_site_mappings,
)

MZID_PATH = Path("data/raw/PXD035795/peptides_1_1_0.mzid.gz")
MZID_SHA256 = "62105dfde42d9675bc5dc7c83969d971eec57c9e02983659cd2df96ec4d1b5fe"

# Registered, hash-audited PXD035795 method source
# (data/registry/evidence_methods.tsv). A reviewed mapping is only honored when
# its method_source_sha256 equals a source that passes audit_method_sources.
REGISTERED_METHOD_SHA256 = (
    "90b8c49cc80f950b08ad6bca2d91a1d4030ef588c7e6115ba1b335891c43a554"
)
UNREGISTERED_METHOD_SHA256 = "0" * 64


def _candidates():  # type: ignore[no-untyped-def]
    return parse_mzidentml(MZID_PATH, MZID_SHA256).candidates


def _first_cys_candidate():  # type: ignore[no-untyped-def]
    for candidate in _candidates():
        if (
            candidate.residue == "C"
            and candidate.pass_threshold == "true"
            and candidate.conflict_status == "clear"
        ):
            return candidate
    raise AssertionError("expected at least one clear, threshold-passing Cys candidate")


def _mapping_dict(candidate, method_sha256: str) -> dict:  # type: ignore[no-untyped-def]
    """Build a reviewed mapping that exactly matches one real Cys candidate."""
    return {
        "version": 1,
        "mappings": [
            {
                "study_accession": "PXD035795",
                "residue": candidate.residue,
                "chemical_tag": candidate.cv_value,
                "monoisotopic_mass_delta": candidate.monoisotopic_mass_delta,
                "mass_tolerance": 0.001,
                "evidence_scope": "dimedone_switch_lc_ms_ms",
                "method_source_sha256": method_sha256,
                "persulfidation_evidence": "site_ms",
            }
        ],
    }


def _write_mapping(tmp_path: Path, payload: dict) -> Path:
    path = tmp_path / "reviewed_site_mappings.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=True), encoding="utf-8")
    return path


def _decisions(mapping_path: Path):  # type: ignore[no-untyped-def]
    mappings = load_reviewed_site_mappings(mapping_path)
    return apply_reviewed_site_mappings(
        candidates=_candidates(),
        mappings=mappings,
        method_config_path=Path("configs/evidence_method_sources_v1.yaml"),
        method_registry_path=Path("data/registry/evidence_methods.tsv"),
    )


def test_empty_mapping_leaves_every_candidate_unresolved(tmp_path: Path) -> None:
    mapping_path = _write_mapping(tmp_path, {"version": 1, "mappings": []})

    decisions = _decisions(mapping_path)

    assert len(decisions) == len(_candidates())
    assert {d.evidence_class for d in decisions} == {"unresolved"}
    assert not any(d.evidence_class == "site_ms" for d in decisions)


def test_reviewed_mapping_upgrades_only_exact_cys_matches(tmp_path: Path) -> None:
    anchor = _first_cys_candidate()
    mapping_path = _write_mapping(
        tmp_path, _mapping_dict(anchor, REGISTERED_METHOD_SHA256)
    )

    decisions = _decisions(mapping_path)
    by_class: dict[str, list] = {"site_ms": [], "unresolved": []}
    for decision in decisions:
        by_class.setdefault(decision.evidence_class, []).append(decision)

    # At least the anchor is upgraded; every upgraded record is an exact match.
    assert by_class["site_ms"], "an exact reviewed Cys match must yield site_ms"
    for decision in by_class["site_ms"]:
        assert decision.residue == "C"
        assert decision.cv_value == anchor.cv_value
        assert decision.monoisotopic_mass_delta == anchor.monoisotopic_mass_delta
        assert decision.method_mapping_status == "resolved_reviewed"
    # Non-cysteine residues (blank/K) are never upgraded by a Cys rule.
    assert all(
        d.evidence_class == "unresolved"
        for d in decisions
        if d.residue != "C"
    )


def test_mapping_without_registered_method_is_fail_closed(tmp_path: Path) -> None:
    anchor = _first_cys_candidate()
    mapping_path = _write_mapping(
        tmp_path, _mapping_dict(anchor, UNREGISTERED_METHOD_SHA256)
    )

    with pytest.raises(RuntimeError, match="method"):
        _decisions(mapping_path)


def test_candidate_failing_threshold_or_conflict_is_never_upgraded(
    tmp_path: Path,
) -> None:
    anchor = _first_cys_candidate()
    mapping_path = _write_mapping(
        tmp_path, _mapping_dict(anchor, REGISTERED_METHOD_SHA256)
    )

    decisions = _decisions(mapping_path)

    for decision in decisions:
        if decision.evidence_class == "site_ms":
            assert decision.pass_threshold == "true"
            assert decision.conflict_status == "clear"
