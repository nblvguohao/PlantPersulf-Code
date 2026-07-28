"""Cross-species structural context of persulfidation targeting.

Companion to ``cross_species_conservation.py`` (family-level convergence):
this module asks whether persulfidated cysteines occupy a *systematically
different structural context* than other cysteines in the same proteins,
and whether that difference — if any — points the same direction across
three independent species. This is the mechanistic layer for the
"convergence is structural/functional, not sequence-level" hypothesis that
explains the weak cross-species *sequence* transfer already measured
(§5.3 of the evidence audit, AP ~1.3x base rate).

**PU semantics preserved**: the comparison set for a persulfidated
protein's other cysteines is labelled ``"unlabeled"`` (never
``"negative"``) — a non-persulfidated call is "not detected", not
"confirmed absent", exactly as everywhere else in this project
(``evaluation/metrics.py``).

**Two independent structural measurements, both reported (2026-07-23
self-review finding)**: the original ``contact_number_proxy``
(``features/structure.py``, explicitly documented there as a
solvent-accessibility PROXY, not rigorous SASA — a Cα contact count) is
kept as-is; a genuine Shrake-Rupley SASA calculation
(``features/sasa.py``, whole-residue and side-chain-sulfur/SG-specific,
full heavy-atom geometry) is added as an independent cross-check. Both
pathways are exposed here and both must be reported — neither silently
replaces the other, so a reader can see whether the proxy-based finding
replicates under a rigorous calculation.

**Test**: reuses ``evaluation/permutation.py``'s ``permutation_test``
engine with a mean-difference metric — the same from-scratch,
scipy-free permutation-test infrastructure already used for the
seed-ensemble AP significance tests elsewhere in this project.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from plantpersulf.evaluation.permutation import PermutationResult, permutation_test
from plantpersulf.features.sasa import cys_residue_sasa, parse_all_atoms
from plantpersulf.features.structure import (
    CysStructureFeature,
    extract_cys_structure_features,
)

Scored = list[tuple[float, str]]


def mean_difference_metric(scored: Scored) -> float | None:
    """mean(score | label=='positive') - mean(score | label=='unlabeled').

    Returns ``None`` if either group is empty (undefined difference), so
    ``permutation_test`` treats it as a non-contributing permutation
    (matching its existing ``None``-handling contract).
    """
    positive = [s for s, y in scored if y == "positive"]
    unlabeled = [s for s, y in scored if y == "unlabeled"]
    if not positive or not unlabeled:
        return None
    return sum(positive) / len(positive) - sum(unlabeled) / len(unlabeled)


def _absolute_mean_difference_metric(scored: Scored) -> float | None:
    """Two-sided wrapper: ``permutation_test`` only tests "observed >=
    permuted" (appropriate for one-sided metrics like AP); the structural
    preference direction is not assumed here, so the p-value is computed
    on ``|mean_diff|`` and the signed value is reported separately."""
    diff = mean_difference_metric(scored)
    return None if diff is None else abs(diff)


@dataclass(frozen=True)
class StructuralContextRow:
    protein_accession: str
    cys_position: int
    label: str  # "positive" | "unlabeled"
    feature: CysStructureFeature


def build_structural_context_rows(
    persulfidated_keys: set[tuple[str, int]],
    proteome: dict[str, str],
    structure_dir: Path,
    accession_to_structure_file: dict[str, Path],
) -> list[StructuralContextRow]:
    """For every protein that (a) contains >=1 persulfidated cysteine and
    (b) has a downloaded AlphaFold structure, extract structure features
    for EVERY cysteine in that protein.

    A cysteine at a key in ``persulfidated_keys`` is labelled
    ``"positive"``; every other cysteine in the same protein is
    ``"unlabeled"``. Proteins with no entry in
    ``accession_to_structure_file`` (no downloaded structure — isoform,
    404, or simply not yet fetched) are skipped entirely; their cysteines
    contribute no rows (never imputed).

    Row order is fully deterministic (sorted by accession, then position)
    regardless of ``persulfidated_keys``' set-iteration order — Python's
    string hash is randomized per process, so iterating an unsorted set
    would make the downstream permutation test's exact p-value
    non-reproducible across runs even at a fixed seed (2026-07-23
    self-review finding: two runs of the same analysis differed in the
    third decimal place of a p-value because of this).
    """
    proteins_with_positive: dict[str, set[int]] = {}
    for accession, position in persulfidated_keys:
        proteins_with_positive.setdefault(accession, set()).add(position)

    rows: list[StructuralContextRow] = []
    for accession in sorted(proteins_with_positive):
        positive_positions = proteins_with_positive[accession]
        structure_path = accession_to_structure_file.get(accession)
        sequence = proteome.get(accession)
        if structure_path is None or sequence is None or not structure_path.is_file():
            continue

        cys_positions = [
            i + 1 for i, residue in enumerate(sequence) if residue == "C"
        ]
        if not cys_positions:
            continue

        pdb_text = structure_path.read_text(encoding="utf-8")
        features = extract_cys_structure_features(pdb_text, accession, cys_positions)

        for feature in features:
            if not feature.has_structure:
                continue
            label = (
                "positive"
                if feature.cys_position in positive_positions
                else "unlabeled"
            )
            rows.append(
                StructuralContextRow(
                    protein_accession=accession,
                    cys_position=feature.cys_position,
                    label=label,
                    feature=feature,
                )
            )
    return rows


@dataclass(frozen=True)
class StructuralContextResult:
    species: str
    n_proteins: int
    n_positive_cys: int
    n_unlabeled_cys: int
    mean_contact_positive: float
    mean_contact_unlabeled: float
    mean_diff: float
    permutation: PermutationResult


def structural_context_test(
    species: str,
    rows: list[StructuralContextRow],
    n_perm: int = 1000,
    seed: int = 0,
) -> StructuralContextResult:
    """Permutation test: is the mean accessibility-proxy contact number of
    persulfidated cysteines different from other cysteines in the same
    proteins, more than expected by chance?

    A *lower* ``contact_number_proxy`` indicates fewer nearby residues
    (more surface-exposed); a *higher* value indicates burial. The sign of
    ``mean_diff`` therefore indicates the direction of any structural
    preference — this function reports it, it does not assume a direction.
    """
    scored: Scored = [
        (row.feature.contact_number_proxy, row.label)
        for row in rows
        if row.feature.contact_number_proxy is not None
    ]
    positive_vals = [s for s, y in scored if y == "positive"]
    unlabeled_vals = [s for s, y in scored if y == "unlabeled"]
    mean_pos = sum(positive_vals) / len(positive_vals) if positive_vals else 0.0
    mean_unl = sum(unlabeled_vals) / len(unlabeled_vals) if unlabeled_vals else 0.0

    perm = permutation_test(
        scored, _absolute_mean_difference_metric, n_perm=n_perm, seed=seed
    )

    return StructuralContextResult(
        species=species,
        n_proteins=len({row.protein_accession for row in rows}),
        n_positive_cys=len(positive_vals),
        n_unlabeled_cys=len(unlabeled_vals),
        mean_contact_positive=mean_pos,
        mean_contact_unlabeled=mean_unl,
        mean_diff=mean_pos - mean_unl,
        permutation=perm,
    )


# ---------------------------------------------------------------------------
# Real Shrake-Rupley SASA pathway (features/sasa.py) — an independent
# cross-check on the contact_number_proxy finding above, added 2026-07-23.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class StructuralContextRowSasa:
    protein_accession: str
    cys_position: int
    label: str  # "positive" | "unlabeled"
    residue_sasa: float
    sg_sasa: float | None  # None if the SG atom is missing from the record


def build_structural_context_rows_sasa(
    persulfidated_keys: set[tuple[str, int]],
    proteome: dict[str, str],
    accession_to_structure_file: dict[str, Path],
    n_sphere_points: int = 92,
) -> list[StructuralContextRowSasa]:
    """SASA analogue of ``build_structural_context_rows``: for every protein
    with >=1 persulfidated cysteine and a downloaded structure, compute real
    Shrake-Rupley SASA (whole-residue and SG-specific) for EVERY cysteine in
    that protein, using the full heavy-atom geometry (not just Cα).

    Chain mapping: AlphaFold DB single-chain models use chain ``"A"`` and
    the UniProt sequence position directly as the PDB residue number — the
    same assumption ``features/structure.py`` already makes for the
    contact-proxy pathway, so results from the two pathways are directly
    comparable position-for-position.

    Row order is fully deterministic (sorted by accession) for the same
    reason as ``build_structural_context_rows`` — see that function's
    docstring.
    """
    proteins_with_positive: dict[str, set[int]] = {}
    for accession, position in persulfidated_keys:
        proteins_with_positive.setdefault(accession, set()).add(position)

    rows: list[StructuralContextRowSasa] = []
    for accession in sorted(proteins_with_positive):
        positive_positions = proteins_with_positive[accession]
        structure_path = accession_to_structure_file.get(accession)
        sequence = proteome.get(accession)
        if structure_path is None or sequence is None or not structure_path.is_file():
            continue

        cys_positions = [
            i + 1 for i, residue in enumerate(sequence) if residue == "C"
        ]
        if not cys_positions:
            continue

        pdb_text = structure_path.read_text(encoding="utf-8")
        atoms = parse_all_atoms(pdb_text)

        for position in cys_positions:
            result = cys_residue_sasa(
                atoms, chain_id="A", res_seq=position, n_sphere_points=n_sphere_points
            )
            if result is None:
                continue
            label = "positive" if position in positive_positions else "unlabeled"
            rows.append(
                StructuralContextRowSasa(
                    protein_accession=accession,
                    cys_position=position,
                    label=label,
                    residue_sasa=result.residue_sasa,
                    sg_sasa=result.sg_sasa,
                )
            )
    return rows


@dataclass(frozen=True)
class StructuralContextResultSasa:
    species: str
    metric_name: str  # "residue_sasa" | "sg_sasa"
    n_proteins: int
    n_positive_cys: int
    n_unlabeled_cys: int
    mean_sasa_positive: float
    mean_sasa_unlabeled: float
    mean_diff: float
    permutation: PermutationResult


def sasa_structural_context_test(
    species: str,
    rows: list[StructuralContextRowSasa],
    metric: str = "residue_sasa",
    n_perm: int = 1000,
    seed: int = 0,
) -> StructuralContextResultSasa:
    """Permutation test on real SASA (Å²): is persulfidated-cysteine
    solvent accessibility different from other cysteines in the same
    proteins, more than expected by chance?

    Unlike ``contact_number_proxy``, SASA has an unambiguous physical
    direction: LOWER SASA means MORE buried, HIGHER SASA means MORE
    solvent-exposed. ``metric`` selects whole-residue SASA (sum over all
    atoms in the cysteine) or SG-only SASA (the reactive thiol/persulfide
    sulfur specifically, mechanistically the more direct measurement for a
    thiol-modifying PTM).
    """
    if metric == "residue_sasa":
        scored: Scored = [(row.residue_sasa, row.label) for row in rows]
    elif metric == "sg_sasa":
        scored = [
            (row.sg_sasa, row.label) for row in rows if row.sg_sasa is not None
        ]
    else:
        raise ValueError(f"unknown metric: {metric!r}")

    positive_vals = [s for s, y in scored if y == "positive"]
    unlabeled_vals = [s for s, y in scored if y == "unlabeled"]
    mean_pos = sum(positive_vals) / len(positive_vals) if positive_vals else 0.0
    mean_unl = sum(unlabeled_vals) / len(unlabeled_vals) if unlabeled_vals else 0.0

    perm = permutation_test(
        scored, _absolute_mean_difference_metric, n_perm=n_perm, seed=seed
    )

    return StructuralContextResultSasa(
        species=species,
        metric_name=metric,
        n_proteins=len({row.protein_accession for row in rows}),
        n_positive_cys=len(positive_vals),
        n_unlabeled_cys=len(unlabeled_vals),
        mean_sasa_positive=mean_pos,
        mean_sasa_unlabeled=mean_unl,
        mean_diff=mean_pos - mean_unl,
        permutation=perm,
    )
