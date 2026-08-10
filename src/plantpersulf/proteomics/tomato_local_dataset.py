"""Within-tomato PU dataset assembly for the tomato-local model (v1).

Everything the project has trained on so far is Arabidopsis (390 positives,
two same-lab Seville tag-switch studies). The kiae271 panel
(``plantpersulf.proteomics.kiae271_sites``: 99 coordinate-verified Solanum
lycopersicum Cys persulfidation sites over 88 proteins, Zhang et al. 2024
Plant Physiology, doi:10.1093/plphys/kiae271, Supplementary Dataset S1) has
previously only been used as an *external* cross-species control panel. This
module uses it the other way round — as target-species supervision for a
model that is both trained and evaluated inside tomato.

Three things this module exists to guarantee:

1. **The 20 unverified Dataset S1 rows are neither positive nor negative.**
   A row that failed coordinate verification, position alignment, or the
   localization-probability bar is *ambiguous evidence*, not background. It
   would be a silent data repair to recycle it into the unlabeled pool, so
   ``kiae271_excluded_keys`` enumerates every (accession, position) the
   supplementary table touches and the builder removes all of them that are
   not verified positives.

2. **Two explicitly named arenas, because they answer different questions.**

   * ``ARENA_PROTEOME`` — background = all other Cys anywhere in the tomato
     reference proteome (subsampled). This is the realistic
     candidate-prioritisation task ("where in the tomato proteome are the
     persulfidation sites?"), but it is confounded by mass-spectrometry
     protein-level detection bias: the 88 panel proteins are, among other
     things, simply proteins abundant enough for a 4D label-free screen to
     see. A model can score well here by ranking *proteins*.
   * ``ARENA_PANEL`` — background = only the *other* cysteines of those same
     88 detected proteins. Protein-level abundance is held constant by
     construction, so the only thing left to learn is which Cys within an
     already-detected protein carries the modification. This is the harder
     and mechanistically more meaningful question, and it is the arena where
     a per-residue reactivity feature such as
     ``local_positive_charge_density`` has to earn its keep.

3. **Grouped folds whose unit is a homology cluster, not a site.** Cysteines
   of one protein are highly correlated (shared ``cys_density``, shared
   ``protein_length``, overlapping flanking windows for adjacent Cys), and
   paralogs share flanking-window composition. Splitting at the site level
   would leak. ``assign_grouped_folds`` keeps every member of an MMseqs2
   cluster on one side of every fold boundary.

The subsample convention (keep every positive, keep ``ratio`` x n_positives
unlabeled rows under a fixed seed, then sort into a stable order) is the same
one used by ``scripts/run_experiment.py::_subsample_unlabeled``.
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from itertools import product
from pathlib import Path

from plantpersulf.proteomics.kiae271_sites import (
    DEFAULT_LOCALIZATION_MIN,
    SHEET_NAME,
    _header_accession,
    _read_xlsx_sheet,
    parse_kiae271_sites,
)

ARENA_PROTEOME = "proteome"
"""Background = all other Cys in the tomato reference proteome (subsampled)."""

ARENA_PANEL = "panel"
"""Background = only the other Cys of the 88 kiae271-detected proteins."""

ARENAS = (ARENA_PROTEOME, ARENA_PANEL)

LABEL_POSITIVE = "positive"
LABEL_UNLABELED = "unlabeled"


@dataclass(frozen=True)
class TomatoPuRow:
    protein_accession: str
    cys_position: int
    label: str  # positive | unlabeled
    cluster_id: str = ""
    fold: int = -1


def kiae271_excluded_keys(
    xlsx_path: Path,
    proteome: dict[str, str],
    localization_min: float = DEFAULT_LOCALIZATION_MIN,
    sheet_name: str = SHEET_NAME,
) -> frozenset[tuple[str, int]]:
    """Every (accession, Cys position) Dataset S1 touches but does not verify.

    A supplementary row lists ``Proteins`` and ``Positions within proteins``
    as parallel ``;``-joined lists. When the two lists align, the parallel
    pairing is used. When they do not align (8 of the 119 rows), the pairing
    is genuinely unknown, so the *cartesian product* of the listed proteins
    and positions is excluded — deliberately over-inclusive, because the cost
    of dropping a handful of Cys from a 250k-cysteine background is nil while
    the cost of scoring a real-but-unverifiable site as background is a
    corrupted negative class.

    Only positions that really are cysteines in the reference proteome are
    returned (anything else could never have entered the Cys background).
    Verified positives are subtracted: they are positives, not exclusions.
    """
    rows = _read_xlsx_sheet(xlsx_path, sheet_name)
    touched: set[tuple[str, int]] = set()
    for row in rows:
        proteins_field = str(row.get("Proteins", "") or "")
        positions_field = str(row.get("Positions within proteins", "") or "")
        proteins_list = [p.strip() for p in proteins_field.split(";") if p.strip()]
        positions_list = [p.strip() for p in positions_field.split(";") if p.strip()]
        if not proteins_list or not positions_list:
            continue
        if len(proteins_list) == len(positions_list):
            pairs = list(zip(proteins_list, positions_list, strict=True))
        else:
            pairs = list(product(proteins_list, positions_list))
        for protein_field, pos_str in pairs:
            accession = _header_accession(protein_field)
            try:
                position = int(pos_str)
            except ValueError:
                continue
            sequence = proteome.get(accession)
            if sequence is None:
                continue
            if 1 <= position <= len(sequence) and sequence[position - 1] == "C":
                touched.add((accession, position))

    table = parse_kiae271_sites(
        xlsx_path, proteome, localization_min=localization_min, sheet_name=sheet_name
    )
    verified = {(s.protein_accession, s.cys_position) for s in table.sites}
    return frozenset(touched - verified)


def _all_cysteine_keys(
    proteome: dict[str, str],
    accessions: list[str],
) -> list[tuple[str, int]]:
    keys: list[tuple[str, int]] = []
    for accession in accessions:
        sequence = proteome[accession]
        for index, residue in enumerate(sequence):
            if residue == "C":
                keys.append((accession, index + 1))
    return keys


def build_tomato_pu_rows(
    xlsx_path: Path,
    proteome: dict[str, str],
    arena: str,
    ratio: int,
    seed: int,
    localization_min: float = DEFAULT_LOCALIZATION_MIN,
) -> list[TomatoPuRow]:
    """Assemble the positive + unlabeled row set for one evaluation arena.

    Keeps every verified positive; draws ``ratio`` x n_positives unlabeled
    rows (or all of them, whichever is fewer) with ``random.Random(seed)``;
    returns them in a stable (accession, position) order so downstream feature
    extraction is reproducible.
    """
    if arena not in ARENAS:
        raise ValueError(f"unknown arena: {arena!r} (expected one of {ARENAS})")

    table = parse_kiae271_sites(
        xlsx_path, proteome, localization_min=localization_min
    )
    positive_keys = {(s.protein_accession, s.cys_position) for s in table.sites}
    excluded = kiae271_excluded_keys(
        xlsx_path, proteome, localization_min=localization_min
    )
    panel_accessions = sorted({s.protein_accession for s in table.sites})

    candidates = (
        panel_accessions if arena == ARENA_PANEL else sorted(proteome)
    )
    background = [
        key
        for key in _all_cysteine_keys(proteome, candidates)
        if key not in positive_keys and key not in excluded
    ]

    keep_n = min(len(background), max(1, ratio * max(1, len(positive_keys))))
    rng = random.Random(seed)
    sampled = rng.sample(background, keep_n)

    rows = [
        TomatoPuRow(accession, position, LABEL_POSITIVE)
        for accession, position in sorted(positive_keys)
    ] + [
        TomatoPuRow(accession, position, LABEL_UNLABELED)
        for accession, position in sampled
    ]
    rows.sort(key=lambda r: (r.protein_accession, r.cys_position))
    return rows


def cluster_id_for(accession: str, cluster_map: dict[str, str]) -> str:
    """Homology cluster of a protein, or its own singleton cluster.

    Same ``__singleton__`` convention as
    ``scripts/validate_cross_species_rice.py`` — an unclustered protein is
    still a valid grouping unit, it is simply a group of one.
    """
    cluster = cluster_map.get(accession)
    return cluster if cluster else f"__singleton__{accession}"


def assign_grouped_folds(
    rows: list[TomatoPuRow],
    cluster_map: dict[str, str],
    n_folds: int,
    seed: int,
) -> list[TomatoPuRow]:
    """Deal homology clusters into ``n_folds`` grouped cross-validation folds.

    Positive-bearing clusters are shuffled and dealt round-robin *first*, so
    every fold receives positives (average precision is undefined on a fold
    with none); background-only clusters are then dealt round-robin on top.
    Every member of a cluster lands in exactly one fold, which is what makes
    the split leakage-safe at the homology level, not merely the protein
    level.
    """
    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")

    labelled = [
        TomatoPuRow(
            protein_accession=row.protein_accession,
            cys_position=row.cys_position,
            label=row.label,
            cluster_id=cluster_id_for(row.protein_accession, cluster_map),
            fold=-1,
        )
        for row in rows
    ]

    positive_clusters = sorted(
        {r.cluster_id for r in labelled if r.label == LABEL_POSITIVE}
    )
    other_clusters = sorted(
        {r.cluster_id for r in labelled} - set(positive_clusters)
    )
    if len(positive_clusters) < n_folds:
        raise ValueError(
            f"only {len(positive_clusters)} positive-bearing clusters for "
            f"{n_folds} folds — some fold would have no positive"
        )

    rng = random.Random(seed)
    rng.shuffle(positive_clusters)
    rng.shuffle(other_clusters)

    fold_of: dict[str, int] = {}
    for i, cluster in enumerate(positive_clusters):
        fold_of[cluster] = i % n_folds
    for i, cluster in enumerate(other_clusters):
        fold_of[cluster] = i % n_folds

    return [
        TomatoPuRow(
            protein_accession=row.protein_accession,
            cys_position=row.cys_position,
            label=row.label,
            cluster_id=row.cluster_id,
            fold=fold_of[row.cluster_id],
        )
        for row in labelled
    ]


def read_cluster_map(path: Path) -> dict[str, str]:
    """Read a ``protein_accession``/``cluster_id`` TSV (the schema produced by
    MMseqs2 easy-cluster and used by ``configs/splits/split_config_v2.yaml``)."""
    import csv

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != ("protein_accession", "cluster_id"):
            raise RuntimeError(f"cluster file has invalid columns: {path}")
        return {row["protein_accession"]: row["cluster_id"] for row in reader}
