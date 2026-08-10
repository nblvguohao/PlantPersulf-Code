"""Multi-species joint persulfidation-site dataset (multi-species training).

Merges the four registered site-level evidence sources into one PU dataset
so a ranker can learn the *shared* chemistry of persulfidation (thiolate
reactivity, local environment) across species instead of being limited to
the ~390 same-lab Arabidopsis positives:

* **arabidopsis** — ``benchmark_v1/sites.tsv`` positives (PXD006140 +
  PXD024061, Seville tag-switch).
* **tomato** — kiae271 differential persulfidome (99 sites, Zhang lab,
  Plant Physiology 2024; ``plantpersulf.proteomics.kiae271_sites``).
* **rice** — PXD072089 (929 sites, Xie lab, PNAS 2026; NM-biotin;
  ``plantpersulf.proteomics.pxd072089_sites``).
* **magnaporthe** — PXD063170 (1,482 sites, Chen lab, Nat Commun 2025;
  ``plantpersulf.proteomics.pxd063170_sites``).

Cross-species homology grouping uses PANTHER **family-level** IDs
(``PTHR#####``, subfamily suffixes dropped — subfamilies are too
fine-grained to be comparable across distantly related species), so a
homology family's members never span CV folds. Conservation proxies v1
are computed from the PANTHER annotation layer alone (species breadth,
member count), never from the label set — label-derived "conservation"
would leak training information.

Literature basis (multi-species PTM training): the soybean multi-organism
phosphosite study (CD-HIT redundancy removal; cluster conflict resolved as
positive since other species' "negatives" are likely undiscovered
positives — identical to PU semantics), GPS 6.0's general-model-plus-
transfer architecture over 185 species, and the consistent finding that
evolutionary conservation is among the strongest PTM-site features
(PSSM-based features beating BLOSUM and amino-acid composition in
S-sulfenylation prediction).
"""

from __future__ import annotations

import csv
import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

SPECIES_ARABIDOPSIS = "arabidopsis"
SPECIES_TOMATO = "tomato"
SPECIES_RICE = "rice"
SPECIES_MAGNAPORTHE = "magnaporthe"
ALL_SPECIES = (
    SPECIES_ARABIDOPSIS,
    SPECIES_TOMATO,
    SPECIES_RICE,
    SPECIES_MAGNAPORTHE,
)


@dataclass(frozen=True)
class MultispeciesSite:
    protein_accession: str
    cys_position: int
    species: str
    study_accession: str
    panther_family: str  # family-level PANTHER id, e.g. "PTHR10782"


# ---------------------------------------------------------------------------
# PANTHER family-level parsing
# ---------------------------------------------------------------------------


def load_panther_family_ids(path: Path) -> dict[str, str]:
    """Parse a UniProt ``fields=accession,xref_panther`` TSV stream into
    ``accession -> family-level PTHR id``.

    Subfamily suffixes (``:SF##``) are dropped; a protein annotated with
    several families gets the lexicographically first family id
    (deterministic). Proteins with no family annotation are absent from
    the mapping.
    """
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            acc = (row.get("Entry") or "").strip()
            raw = row.get("PANTHER") or ""
            families = sorted(
                {
                    t.split(":")[0]
                    for t in raw.split(";")
                    if t.strip() and t.strip().startswith("PTHR")
                }
            )
            if acc and families:
                mapping[acc] = families[0]
    return mapping


def family_stats_from_panther(
    panther_by_species: dict[str, dict[str, str]],
) -> dict[str, tuple[int, int]]:
    """``family -> (species_breadth, member_count)`` from the PANTHER
    annotation layer only (all annotated members across all species) —
    never from the label set."""
    stats: dict[str, tuple[set[str], int]] = {}
    for species, mapping in panther_by_species.items():
        for family in mapping.values():
            species_set, count = stats.get(family, (set(), 0))
            species_set.add(species)
            stats[family] = (species_set, count + 1)
    return {
        family: (len(species_set), count)
        for family, (species_set, count) in stats.items()
    }


# ---------------------------------------------------------------------------
# merge
# ---------------------------------------------------------------------------


def build_multispecies_sites(
    arabidopsis: list[MultispeciesSite] = (),
    tomato: list[MultispeciesSite] = (),
    rice: list[MultispeciesSite] = (),
    magnaporthe: list[MultispeciesSite] = (),
) -> tuple[MultispeciesSite, ...]:
    """Merge per-species site lists into one deduplicated, sorted tuple.

    Duplicate (accession, position) keys collapse to the first occurrence
    (identical sites reported by more than one study are not double
    counted). Sorting is by (species, accession, position) for
    determinism.
    """
    seen: dict[tuple[str, int], MultispeciesSite] = {}
    for site in (*arabidopsis, *tomato, *rice, *magnaporthe):
        key = (site.protein_accession, site.cys_position)
        if key not in seen:
            seen[key] = site
    order = {s: i for i, s in enumerate(ALL_SPECIES)}
    return tuple(
        sorted(seen.values(), key=lambda s: (order[s.species], s.protein_accession, s.cys_position))
    )


# ---------------------------------------------------------------------------
# family-grouped CV
# ---------------------------------------------------------------------------


NO_FAMILY = "__nofamily__"


def family_grouped_folds(
    sites: list[MultispeciesSite] | tuple[MultispeciesSite, ...],
    n_folds: int,
    seed: int,
) -> list[tuple[list[int], list[int]]]:
    """Stratified-by-family CV fold assignment: a PANTHER family's members
    (across all species) always land in exactly one fold, so no homology
    family spans the train/test boundary. Proteins WITHOUT a PANTHER
    family annotation are grouped per-protein (each becomes its own
    singleton group) so they spread across folds instead of pooling into
    one giant '__nofamily__' fold. Groups are distributed largest-first
    into the currently smallest fold (deterministic given ``seed``).
    Returns ``(train_idx, test_idx)`` per fold."""
    if n_folds < 2:
        raise ValueError("n_folds must be >= 2")
    family_to_idx: dict[str, list[int]] = defaultdict(list)
    for i, site in enumerate(sites):
        group = (
            site.panther_family
            if site.panther_family != NO_FAMILY
            else f"__protein__{site.protein_accession}"
        )
        family_to_idx[group].append(i)

    families = sorted(family_to_idx)
    rng = random.Random(seed)
    rng.shuffle(families)
    families.sort(key=lambda f: -len(family_to_idx[f]))

    fold_families: list[list[str]] = [[] for _ in range(n_folds)]
    for family in families:
        smallest = min(range(n_folds), key=lambda k: len(fold_families[k]))
        fold_families[smallest].append(family)

    folds: list[tuple[list[int], list[int]]] = []
    for test_families in fold_families:
        test_idx = sorted(
            i for f in test_families for i in family_to_idx[f]
        )
        test_set = set(test_idx)
        train_idx = [i for i in range(len(sites)) if i not in test_set]
        folds.append((train_idx, test_idx))
    return folds
