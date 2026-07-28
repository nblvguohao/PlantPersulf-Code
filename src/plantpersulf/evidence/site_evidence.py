"""Phase A2 (Route M): reviewed method mapping → PXD035795 site_ms upgrade.

This is a fail-closed join. An mzIdentML candidate modification is upgraded to
``site_ms`` only when a reviewed mapping — itself backed by a registered,
hash-audited method source — exactly matches the candidate and every evidence
policy condition holds. Absent the mapping, an unregistered method source, or
any unmet condition, the candidate stays ``unresolved``.

The module encodes the *gating contract*, not the biology: which chemical tag
or mass on cysteine proves persulfidation is a reviewed scientific decision
supplied through the mapping plus its method source, never inferred here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.evidence.methods import audit_method_sources
from plantpersulf.evidence.mzidentml import MzidCandidate

SITE_EVIDENCE_CLASS = "site_ms"


@dataclass(frozen=True)
class ReviewedSiteMapping:
    study_accession: str
    residue: str
    chemical_tag: str
    monoisotopic_mass_delta: str
    mass_tolerance: float
    evidence_scope: str
    method_source_sha256: str
    persulfidation_evidence: str


@dataclass(frozen=True)
class SiteEvidenceDecision:
    candidate: MzidCandidate
    evidence_class: str
    method_mapping_status: str
    reason: str

    @property
    def residue(self) -> str:
        return self.candidate.residue

    @property
    def cv_value(self) -> str:
        return self.candidate.cv_value

    @property
    def monoisotopic_mass_delta(self) -> str:
        return self.candidate.monoisotopic_mass_delta

    @property
    def pass_threshold(self) -> str:
        return self.candidate.pass_threshold

    @property
    def conflict_status(self) -> str:
        return self.candidate.conflict_status


def _required_string(mapping: dict[str, Any], field: str) -> str:
    value = mapping.get(field)
    if not isinstance(value, str) or not value.strip():
        raise RuntimeError(f"reviewed site mapping has invalid {field}")
    return value


def load_reviewed_site_mappings(path: Path) -> tuple[ReviewedSiteMapping, ...]:
    """Load the exact reviewed mapping list; empty is valid."""
    try:
        loaded: object = yaml.safe_load(path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(f"invalid reviewed site mappings: {path}") from exc
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise RuntimeError("reviewed site mappings require version 1")
    raw_mappings = loaded.get("mappings")
    if not isinstance(raw_mappings, list):
        raise RuntimeError("reviewed site mappings require a mappings list")
    mappings: list[ReviewedSiteMapping] = []
    for raw in raw_mappings:
        if not isinstance(raw, dict):
            raise RuntimeError("reviewed site mapping must be a mapping")
        entry = cast(dict[str, Any], raw)
        tolerance = entry.get("mass_tolerance")
        if not isinstance(tolerance, (int, float)) or isinstance(
            tolerance, bool
        ) or tolerance < 0:
            raise RuntimeError("reviewed site mapping has invalid mass_tolerance")
        persulfidation = _required_string(entry, "persulfidation_evidence")
        if persulfidation != SITE_EVIDENCE_CLASS:
            raise RuntimeError(
                "reviewed site mapping may only assign site_ms evidence"
            )
        mappings.append(
            ReviewedSiteMapping(
                study_accession=_required_string(entry, "study_accession"),
                residue=_required_string(entry, "residue"),
                chemical_tag=_required_string(entry, "chemical_tag"),
                monoisotopic_mass_delta=_required_string(
                    entry, "monoisotopic_mass_delta"
                ),
                mass_tolerance=float(tolerance),
                evidence_scope=_required_string(entry, "evidence_scope"),
                method_source_sha256=_required_string(
                    entry, "method_source_sha256"
                ),
                persulfidation_evidence=persulfidation,
            )
        )
    return tuple(mappings)


def _mapping_matches(
    candidate: MzidCandidate,
    mapping: ReviewedSiteMapping,
) -> bool:
    if candidate.residue != mapping.residue:
        return False
    if candidate.cv_value != mapping.chemical_tag:
        return False
    try:
        observed = float(candidate.monoisotopic_mass_delta)
        expected = float(mapping.monoisotopic_mass_delta)
    except ValueError:
        return False
    return abs(observed - expected) <= mapping.mass_tolerance


def _policy_satisfied(candidate: MzidCandidate) -> bool:
    return (
        candidate.pass_threshold == "true"
        and candidate.conflict_status == "clear"
        and bool(candidate.source_locator)
        and bool(candidate.protein_accession)
        and bool(candidate.peptide_sequence)
        and bool(candidate.location)
    )


def apply_reviewed_site_mappings(
    candidates: tuple[MzidCandidate, ...],
    mappings: tuple[ReviewedSiteMapping, ...],
    method_config_path: Path = Path("configs/evidence_method_sources_v1.yaml"),
    method_registry_path: Path = Path("data/registry/evidence_methods.tsv"),
) -> tuple[SiteEvidenceDecision, ...]:
    """Upgrade candidates to site_ms only under a registered reviewed mapping."""
    registered = {
        (source.study_accession, source.sha256)
        for source in audit_method_sources(
            method_config_path, method_registry_path
        )
    }
    for mapping in mappings:
        if (
            mapping.study_accession,
            mapping.method_source_sha256,
        ) not in registered:
            raise RuntimeError(
                "reviewed site mapping references an unregistered method "
                f"source: {mapping.study_accession}"
            )

    decisions: list[SiteEvidenceDecision] = []
    for candidate in candidates:
        matched = next(
            (
                mapping
                for mapping in mappings
                if _mapping_matches(candidate, mapping)
            ),
            None,
        )
        if matched is None:
            decisions.append(
                SiteEvidenceDecision(
                    candidate=candidate,
                    evidence_class="unresolved",
                    method_mapping_status="absent",
                    reason="no_reviewed_mapping_matches_candidate",
                )
            )
            continue
        if not _policy_satisfied(candidate):
            decisions.append(
                SiteEvidenceDecision(
                    candidate=candidate,
                    evidence_class="unresolved",
                    method_mapping_status="policy_unmet",
                    reason="evidence_policy_conditions_not_satisfied",
                )
            )
            continue
        decisions.append(
            SiteEvidenceDecision(
                candidate=candidate,
                evidence_class=matched.persulfidation_evidence,
                method_mapping_status="resolved_reviewed",
                reason="reviewed_method_mapping_and_policy_satisfied",
            )
        )
    return tuple(decisions)
