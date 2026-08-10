"""Label-free novel-cysteine universe for the frozen tomato snapshot."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from plantpersulf.features.sequence import _load_proteome
from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites
from plantpersulf.proteomics.tomato_local_dataset import kiae271_excluded_keys
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file


@dataclass(frozen=True)
class TomatoCandidate:
    protein_accession: str
    cys_position: int
    site_key: str


@dataclass(frozen=True)
class TomatoCandidateRegistry:
    candidates: tuple[TomatoCandidate, ...]
    kiae271_sha256: str
    proteome_sha256: str
    excluded_known_positive_count: int
    excluded_ambiguous_count: int


def _site_key(accession: str, position: int) -> str:
    return f"{accession}:C{position}"


def build_tomato_candidate_registry(
    kiae271_xlsx: Path,
    proteome_path: Path,
    supplementary_registry: Path,
    model_input_registry: Path,
) -> TomatoCandidateRegistry:
    """Enumerate only unobserved, unambiguous Cys sites in the tomato snapshot."""
    assert_registered_input(kiae271_xlsx, supplementary_registry)
    assert_registered_input(proteome_path, model_input_registry)
    proteome = _load_proteome(proteome_path)
    positives = {
        (site.protein_accession, site.cys_position)
        for site in parse_kiae271_sites(kiae271_xlsx, proteome).sites
    }
    ambiguous = kiae271_excluded_keys(kiae271_xlsx, proteome)
    excluded = positives | ambiguous
    candidates = tuple(
        TomatoCandidate(
            protein_accession=accession,
            cys_position=position,
            site_key=_site_key(accession, position),
        )
        for accession in sorted(proteome)
        for position, residue in enumerate(proteome[accession], start=1)
        if residue == "C" and (accession, position) not in excluded
    )
    return TomatoCandidateRegistry(
        candidates=candidates,
        kiae271_sha256=hash_file(kiae271_xlsx, "sha256"),
        proteome_sha256=hash_file(proteome_path, "sha256"),
        excluded_known_positive_count=len(positives),
        excluded_ambiguous_count=len(ambiguous),
    )
