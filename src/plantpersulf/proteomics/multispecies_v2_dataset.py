"""Split-bound positive-unlabelled rows for multispecies v2.

This module is intentionally independent from ``multispecies_dataset``:
the latter preserves the historical v1 PANTHER-based construction, while v2
uses the immutable global MMseqs2 cluster split as its sole partition source.
"""

from __future__ import annotations

import random
from collections.abc import Iterable
from dataclasses import dataclass

from plantpersulf.benchmark.multispecies_splits import (
    DEVELOPMENT_SPLIT,
    FrozenMultispeciesSplit,
    GlobalClusterRow,
    MultispeciesSiteRow,
    global_protein_id,
)


@dataclass(frozen=True)
class MultispeciesV2SiteRow:
    """One Cys coordinate with inherited global cluster and frozen split."""

    species: str
    protein_accession: str
    cys_position: int
    label: str
    study_accessions: tuple[str, ...]
    global_protein_id: str
    cluster_id: str
    split: str
    development_fold: int | None


@dataclass(frozen=True)
class V2DevelopmentFold:
    """The only rows a development-fold fitter may receive."""

    fit_rows: tuple[MultispeciesV2SiteRow, ...]
    validation_rows: tuple[MultispeciesV2SiteRow, ...]


def build_split_bound_v2_rows(
    *,
    positives: Iterable[MultispeciesSiteRow],
    all_cysteines: Iterable[tuple[str, str, int]],
    clusters: Iterable[GlobalClusterRow],
    frozen_split: FrozenMultispeciesSplit,
) -> tuple[MultispeciesV2SiteRow, ...]:
    """Return development-only Cys rows assigned to the immutable split.

    ``all_cysteines`` must come from SHA-verified reference proteomes.  The
    caller therefore supplies coordinates rather than sequences, keeping this
    policy layer unable to invent or repair biological input. Frozen-test rows
    are deliberately omitted: test labels are only available through the
    separately unlocked scoring boundary.
    """
    cluster_by_protein: dict[str, GlobalClusterRow] = {}
    for cluster_row in clusters:
        if cluster_row.global_protein_id in cluster_by_protein:
            raise RuntimeError(
                f"duplicate global cluster: {cluster_row.global_protein_id}"
            )
        cluster_by_protein[cluster_row.global_protein_id] = cluster_row
    split_by_protein = {}
    for split_row in frozen_split.rows:
        if split_row.global_protein_id in split_by_protein:
            raise RuntimeError(
                f"duplicate frozen split row: {split_row.global_protein_id}"
            )
        split_by_protein[split_row.global_protein_id] = split_row
    positive_studies: dict[tuple[str, str, int], set[str]] = {}
    for positive_row in positives:
        if positive_row.label != "positive":
            raise ValueError("v2 positive input must use label='positive'")
        assigned_split = split_by_protein.get(positive_row.global_protein_id)
        if assigned_split is None:
            raise RuntimeError(
                "positive protein missing frozen split: "
                f"{positive_row.global_protein_id}"
            )
        if assigned_split.split != DEVELOPMENT_SPLIT:
            raise RuntimeError(
                f"frozen test positive supplied to development builder: "
                f"{positive_row.global_protein_id}"
            )
        key = (
            positive_row.species,
            positive_row.protein_accession,
            positive_row.cys_position,
        )
        positive_studies.setdefault(key, set()).add(positive_row.study_accession)

    cysteine_keys: set[tuple[str, str, int]] = set()
    for species, accession, position in all_cysteines:
        key = (species, accession, position)
        if position < 1:
            raise ValueError("cysteine positions must be one-based positive integers")
        if key in cysteine_keys:
            raise RuntimeError(f"duplicate reference cysteine coordinate: {key}")
        cysteine_keys.add(key)
    missing_positive = set(positive_studies).difference(cysteine_keys)
    if missing_positive:
        missing_key = sorted(missing_positive)[0]
        raise RuntimeError(
            f"positive coordinate missing from reference cysteines: {missing_key}"
        )

    result: list[MultispeciesV2SiteRow] = []
    for species, accession, position in sorted(cysteine_keys):
        protein_id = global_protein_id(species, accession)
        cluster = cluster_by_protein.get(protein_id)
        if cluster is None:
            raise RuntimeError(f"missing global cluster for protein: {protein_id}")
        split = split_by_protein.get(protein_id)
        if split is None:
            raise RuntimeError(f"missing frozen split for protein: {protein_id}")
        if split.cluster_id != cluster.cluster_id:
            raise RuntimeError(
                f"global cluster and frozen split disagree: {protein_id}"
            )
        if split.split != DEVELOPMENT_SPLIT:
            continue
        key = (species, accession, position)
        studies = tuple(sorted(positive_studies.get(key, set())))
        result.append(
            MultispeciesV2SiteRow(
                species=species,
                protein_accession=accession,
                cys_position=position,
                label="positive" if studies else "unlabeled",
                study_accessions=studies,
                global_protein_id=protein_id,
                cluster_id=cluster.cluster_id,
                split=split.split,
                development_fold=split.development_fold,
            )
        )
    return tuple(result)


def sample_v2_unlabeled_panels(
    rows: Iterable[MultispeciesV2SiteRow], *, per_positive: int, seed: int
) -> tuple[MultispeciesV2SiteRow, ...]:
    """Keep every positive and sample unlabelled Cys within each frozen split."""
    if per_positive < 1:
        raise ValueError("per_positive must be >= 1")
    by_split: dict[str, list[MultispeciesV2SiteRow]] = {}
    seen: set[tuple[str, str, int]] = set()
    for row in rows:
        key = (row.species, row.protein_accession, row.cys_position)
        if key in seen:
            raise RuntimeError(f"duplicate split-bound site coordinate: {key}")
        seen.add(key)
        by_split.setdefault(row.split, []).append(row)

    selected: list[MultispeciesV2SiteRow] = []
    for index, (_split, partition) in enumerate(sorted(by_split.items())):
        positives = [row for row in partition if row.label == "positive"]
        unlabeled = [row for row in partition if row.label == "unlabeled"]
        sample_size = min(len(unlabeled), per_positive * len(positives))
        sampled = random.Random(seed + index).sample(unlabeled, sample_size)
        selected.extend(positives)
        selected.extend(sampled)
    return tuple(
        sorted(
            selected,
            key=lambda row: (
                row.split,
                row.species,
                row.protein_accession,
                row.cys_position,
            ),
        )
    )


def prepare_v2_development_fold(
    rows: Iterable[MultispeciesV2SiteRow],
    *,
    validation_fold: int,
    n_development_folds: int,
) -> V2DevelopmentFold:
    """Return one fit/validation pair while excluding every frozen-test row."""
    if not 0 <= validation_fold < n_development_folds:
        raise ValueError("validation_fold outside configured development folds")
    fit_rows: list[MultispeciesV2SiteRow] = []
    validation_rows: list[MultispeciesV2SiteRow] = []
    protein_assignment: dict[str, tuple[str, int | None]] = {}
    for row in rows:
        assignment = (row.split, row.development_fold)
        previous = protein_assignment.setdefault(row.global_protein_id, assignment)
        if previous != assignment:
            raise RuntimeError(
                f"protein crosses v2 partitions: {row.global_protein_id}"
            )
        if row.split != DEVELOPMENT_SPLIT:
            continue
        if row.development_fold is None:
            raise RuntimeError(
                f"development protein lacks fold: {row.global_protein_id}"
            )
        if row.development_fold == validation_fold:
            validation_rows.append(row)
        else:
            fit_rows.append(row)
    return V2DevelopmentFold(tuple(fit_rows), tuple(validation_rows))
