"""Shared frozen-bundle scoring for diagnostic site lists (Gate-0-safe).

Scores arbitrary ``(species, protein_accession, cys_position)`` lists with the
frozen candidate-release structure_ranker bundle using the registered
reference proteomes, mirroring the exact feature path of the release's
primary artifact (``top_k_candidates.tsv``): the structure branch is masked
for every scored site, because release v2's structure registry contained no
tomato accessions and every tomato candidate was therefore scored with the
structure branch masked.

This module performs NO fitting and mutates no frozen artifact; it is the
shared diagnostic-only scoring path used by the co-peptide-negative,
within-protein-ranking, and regime-hypothesis diagnostics. Coordinates are
fail-closed: the feature extractor raises if a position is not a Cys in the
registered reference proteome.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

from plantpersulf.models.structure_ranker import (
    BranchFeatures,
    StructureRankerBundle,
    score_structure_ranker_bundle,
)
from plantpersulf.proteomics.multispecies_v2_dataset import MultispeciesV2SiteRow
from plantpersulf.proteomics.multispecies_v2_sources import sequence_feature_map


@dataclass(frozen=True)
class ScoredSite:
    """One Cys site scored by the frozen bundle (structure branch masked)."""

    species: str
    protein_accession: str
    cys_position: int
    score: float
    uncertainty: float


def cys_positions(sequence: str) -> tuple[int, ...]:
    """1-based positions of every Cys in a protein sequence."""
    return tuple(i + 1 for i, residue in enumerate(sequence) if residue == "C")


def make_site_row(
    species: str, accession: str, position: int
) -> MultispeciesV2SiteRow:
    """A development-scope unlabeled row for one coordinate — the same row
    shape the frozen feature extractor consumes."""
    return MultispeciesV2SiteRow(
        species=species,
        protein_accession=accession,
        cys_position=position,
        label="unlabeled",
        study_accessions=(),
        global_protein_id=f"{species}|{accession}",
        cluster_id="",
        split="development",
        development_fold=None,
    )


def score_sites(
    sites: Sequence[tuple[str, str, int]],
    *,
    bundle: StructureRankerBundle,
    proteomes: dict[str, dict[str, str]],
) -> dict[tuple[str, str, int], ScoredSite]:
    """Score ``(species, accession, position)`` sites with the frozen bundle.

    Duplicate coordinates collapse (first occurrence wins); output order is
    the input order after deduplication. Raises ``RuntimeError`` when a
    coordinate is not a Cys in the registered reference proteome.
    """
    deduped: dict[tuple[str, str, int], None] = {}
    for site in sites:
        deduped.setdefault(site, None)
    ordered = list(deduped)

    rows = [make_site_row(*site) for site in ordered]
    raw_features = sequence_feature_map(rows, proteomes)
    ordered_keys = sorted(raw_features)
    sequence_columns: list[list[float]] = []
    for key in ordered_keys:
        values = raw_features[key]
        checked: list[float] = []
        for value in values:
            if value is None:
                raise RuntimeError(f"feature extraction returned None for {key}")
            checked.append(value)
        sequence_columns.append(checked)
    features = BranchFeatures(
        sequence=sequence_columns,
        esm=[[0.0] for _ in ordered_keys],
        structure=[[0.0, 0.0] for _ in ordered_keys],
        structure_mask=[False for _ in ordered_keys],
        study_ids=None,
    )
    output = score_structure_ranker_bundle(bundle, features, device_name="cpu")

    site_by_feature_key: dict[tuple[str, int], tuple[str, str, int]] = {
        (make_site_row(*site).global_protein_id, site[2]): site for site in ordered
    }
    result: dict[tuple[str, str, int], ScoredSite] = {}
    for key, score, uncertainty in zip(
        ordered_keys, output.scores, output.uncertainty, strict=True
    ):
        matched = site_by_feature_key.get(key)
        if matched is None:
            raise RuntimeError(f"feature key without matching site: {key}")
        result[matched] = ScoredSite(
            species=matched[0],
            protein_accession=matched[1],
            cys_position=matched[2],
            score=float(score),
            uncertainty=float(uncertainty),
        )
    return result
