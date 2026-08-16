"""Cross-species conservation of persulfidation targeting via shared PANTHER
ortholog subfamilies.

Framing shift (2026-07-23): rather than asking whether a model trained on
one species' sites *predicts* another species' sites (the Gate 2 / transfer
question, structurally blocked by study non-independence — see
``docs/phase_z_evidence_audit.md``), this module asks a different,
data-available question: **is the set of PANTHER ortholog families
containing a persulfidated protein in one species enriched for also
containing a persulfidated protein in an independent species, beyond
chance?** This uses exactly the same independent-lab/chemistry/species
datasets (Arabidopsis benchmark_v1, PXD072089 rice, kiae271 tomato,
PXD063170 Magnaporthe) that Gate 2 rejected as *predictive* evidence, but as
*associational* evidence instead — a well-established comparative-genomics
methodology (ortholog-group co-annotation enrichment), not a downgrade in
rigor.

**Orthology proxy**: PANTHER subfamily IDs (``PTHR#####:SF#``), fetched in
bulk from UniProt's ``xref_panther`` field, are used as a fine-grained
ortholog-group proxy — two proteins from different species sharing a
subfamily ID are treated as belonging to the same ortholog group. This is
coarser than a curated 1:1 orthology call (e.g. reciprocal-best-hit) but
requires no new tooling and has near-complete proteome coverage for the
plant species (100% for Arabidopsis/rice at the family level; the
subfamily level — used here, since it is the orthology-informative tier —
covers 77.2% of Arabidopsis and 61.6% of rice; proteins without a
subfamily call are excluded from the comparison, not assumed
non-orthologous).

**Magnaporthe accession bridge**: PXD063170's site table uses EnsemblFungi
gene IDs (``MGG_#####T0``), not UniProt accessions. UniProt's Magnaporthe
oryzae 70-15 entries (taxon 242507) carry the bare gene ID (``MGG_#####``,
no transcript suffix) in the ``Gene Names (ORF)`` field, giving a
99.2%-coverage bridge from MG8 accessions to UniProt PANTHER calls
(12,494/12,593 MG8 proteins resolved).

**Test**: for a species pair, restrict to PANTHER subfamilies present in
both proteomes (``shared_families``). Build the 2x2 contingency table of
"family contains a persulfidated protein in species A" x "... in species
B" and run a one-sided Fisher's exact test (hypergeometric tail) for
enrichment of co-persulfidation. No scipy dependency — the hypergeometric
tail is computed directly via ``math.comb``, consistent with this
project's existing from-scratch statistics (``permutation.py``,
``effect_size.py``).
"""

from __future__ import annotations

import csv
import random
from collections import Counter
from dataclasses import dataclass
from math import comb
from pathlib import Path


@dataclass(frozen=True)
class SpeciesPantherMap:
    """accession -> frozenset of PANTHER subfamily IDs (``PTHR#####:SF#``
    only; family-level-only calls without a subfamily are dropped)."""

    family_by_accession: dict[str, frozenset[str]]

    def all_families(self) -> frozenset[str]:
        families: set[str] = set()
        for fams in self.family_by_accession.values():
            families |= fams
        return frozenset(families)

    def persulfidated_families(
        self, persulfidated_accessions: set[str]
    ) -> frozenset[str]:
        families: set[str] = set()
        for acc in persulfidated_accessions:
            families |= self.family_by_accession.get(acc, frozenset())
        return frozenset(families)


def _parse_subfamilies(raw: str) -> frozenset[str]:
    return frozenset(t for t in raw.split(";") if t and ":SF" in t)


def load_panther_annotations(path: Path) -> SpeciesPantherMap:
    """Parse a UniProt ``fields=accession,xref_panther`` TSV stream.

    Keyed by UniProt accession (column ``Entry``).
    """
    mapping: dict[str, frozenset[str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            acc = row["Entry"]
            mapping[acc] = _parse_subfamilies(row.get("PANTHER", ""))
    return SpeciesPantherMap(family_by_accession=mapping)


def load_panther_annotations_by_orf_gene(
    path: Path,
    accession_to_gene: dict[str, str],
) -> SpeciesPantherMap:
    """Parse a UniProt ``fields=accession,gene_orf,xref_panther`` TSV stream,
    then re-key from UniProt accession to the caller's accession space via
    a ``caller_accession -> ORF gene name`` bridge.

    Used for Magnaporthe: PXD063170 site-table accessions (``MGG_#####T0``)
    are EnsemblFungi gene IDs, not UniProt accessions; UniProt's own
    Magnaporthe entries carry the bare gene ID (``MGG_#####``) in the
    ``Gene Names (ORF)`` field, which bridges the two accession spaces.
    """
    gene_to_families: dict[str, frozenset[str]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = (row.get("Gene Names (ORF)") or "").strip()
            if not gene:
                continue
            gene_to_families[gene] = _parse_subfamilies(row.get("PANTHER", ""))

    mapping: dict[str, frozenset[str]] = {}
    for caller_acc, gene in accession_to_gene.items():
        families = gene_to_families.get(gene)
        if families is not None:
            mapping[caller_acc] = families
    return SpeciesPantherMap(family_by_accession=mapping)


def mg8_accession_to_gene(accession: str) -> str:
    """``MGG_07573T0`` -> ``MGG_07573`` (strip the Ensembl Fungi transcript
    suffix to recover the bare gene ID UniProt indexes Magnaporthe by)."""
    if "T" in accession:
        gene, _, suffix = accession.rpartition("T")
        if suffix.isdigit():
            return gene
    return accession


@dataclass(frozen=True)
class ConservationTestResult:
    species_a: str
    species_b: str
    shared_families: int
    persulfidated_families_a: int
    persulfidated_families_b: int
    co_persulfidated_families: int
    expected_co_persulfidated_under_independence: float
    odds_ratio: float | None
    p_value: float


def _hypergeometric_upper_tail(
    successes_drawn: int,
    total_successes: int,
    total_draws: int,
    population: int,
) -> float:
    """P(X >= successes_drawn) for X ~ Hypergeometric(population,
    total_successes, total_draws) — the one-sided Fisher exact tail for
    enrichment. Computed exactly via ``math.comb``."""
    upper = min(total_draws, total_successes)
    denom = comb(population, total_draws)
    total = 0
    for x in range(successes_drawn, upper + 1):
        total += comb(total_successes, x) * comb(
            population - total_successes, total_draws - x
        )
    return total / denom


def cross_species_conservation_test(
    species_a: str,
    map_a: SpeciesPantherMap,
    persulfidated_a: set[str],
    species_b: str,
    map_b: SpeciesPantherMap,
    persulfidated_b: set[str],
) -> ConservationTestResult:
    """One-sided Fisher's exact test for enrichment of co-persulfidation
    among PANTHER subfamilies shared by both proteomes.

    Contingency table over the shared-family universe:
    ``a`` = families persulfidated in both species, ``b`` = persulfidated
    in A only, ``c`` = persulfidated in B only, ``d`` = persulfidated in
    neither. Tests whether ``a`` exceeds the count expected if
    persulfidation status in A and B were independent given the marginals.
    """
    shared = map_a.all_families() & map_b.all_families()
    fam_a = map_a.persulfidated_families(persulfidated_a) & shared
    fam_b = map_b.persulfidated_families(persulfidated_b) & shared

    both = fam_a & fam_b
    a = len(both)
    b = len(fam_a) - a
    c = len(fam_b) - a
    d = len(shared) - a - b - c

    n = len(shared)
    expected = (len(fam_a) * len(fam_b) / n) if n > 0 else 0.0
    odds_ratio = (a * d) / (b * c) if b > 0 and c > 0 else None

    p_value = (
        _hypergeometric_upper_tail(
            successes_drawn=a,
            total_successes=len(fam_b),
            total_draws=len(fam_a),
            population=n,
        )
        if n > 0
        else 1.0
    )

    return ConservationTestResult(
        species_a=species_a,
        species_b=species_b,
        shared_families=n,
        persulfidated_families_a=len(fam_a),
        persulfidated_families_b=len(fam_b),
        co_persulfidated_families=a,
        expected_co_persulfidated_under_independence=expected,
        odds_ratio=odds_ratio,
        p_value=p_value,
    )


@dataclass(frozen=True)
class CorrectedPValue:
    p_value: float
    bonferroni_p_value: float
    benjamini_hochberg_q_value: float
    significant_bonferroni: bool
    significant_bh: bool


def bonferroni_and_bh_correction(
    p_values: list[float],
    alpha: float = 0.05,
) -> list[CorrectedPValue]:
    """Multiple-testing correction for a small family of p-values from
    independent tests (here: the three pairwise species-conservation
    tests).

    Returns results in the SAME order as the input ``p_values`` list.

    - Bonferroni: ``p * m``, capped at 1.0 — controls family-wise error
      rate, appropriate given only ``m=3`` tests here.
    - Benjamini-Hochberg: standard step-up FDR procedure, included as a
      less conservative cross-check; with ``m=3`` the two procedures can
      diverge only modestly.

    No scipy/statsmodels dependency, consistent with this project's
    existing from-scratch statistics (``permutation.py``,
    ``effect_size.py``, and the hypergeometric tail above).
    """
    m = len(p_values)
    if m == 0:
        return []

    bonferroni = [min(1.0, p * m) for p in p_values]

    # Benjamini-Hochberg: sort ascending, compute q = p * m / rank, then
    # enforce monotonicity by taking a running minimum from the largest
    # p-value down to the smallest.
    order = sorted(range(m), key=lambda i: p_values[i])
    bh_sorted = [0.0] * m
    running_min = 1.0
    for rank_from_largest, idx in enumerate(reversed(order)):
        rank = m - rank_from_largest  # 1-indexed rank ascending
        q = p_values[idx] * m / rank
        running_min = min(running_min, q)
        bh_sorted[idx] = running_min
    bh_sorted = [min(1.0, q) for q in bh_sorted]

    return [
        CorrectedPValue(
            p_value=p_values[i],
            bonferroni_p_value=bonferroni[i],
            benjamini_hochberg_q_value=bh_sorted[i],
            significant_bonferroni=bonferroni[i] < alpha,
            significant_bh=bh_sorted[i] < alpha,
        )
        for i in range(m)
    ]


# ---------------------------------------------------------------------------
# N-way conservation spectrum
#
# The pairwise Fisher test above answers "do species A and B co-target the
# same ortholog subfamilies?". With four datasets the biologically pointed
# question is the n-way one: how many subfamilies carry a persulfidated
# protein in *every* species, and is that more than independence predicts?
# Pairwise enrichment does not imply an n-way excess (three pairwise
# associations are compatible with an empty triple intersection), so the top
# cell needs its own null rather than an inference from the pairwise table.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ConservationSpectrum:
    """Distribution of ortholog subfamilies over "persulfidated in k of the
    S species", restricted to the subfamilies present in every proteome.

    ``observed_by_species_count`` lists only the k values that actually
    occur; ``expected_by_species_count`` covers every k in ``0..S`` so the
    observed/expected comparison is defined at the top cell even when it is
    empty. Expectations are the Poisson-binomial means implied by each
    species' own marginal targeting rate over the shared universe, i.e. the
    independence null with all marginals held fixed.
    """

    species: tuple[str, ...]
    universe_size: int
    persulfidated_families_by_species: dict[str, int]
    observed_by_species_count: dict[int, int]
    expected_by_species_count: dict[int, float]
    conserved_in_all: tuple[str, ...]


def _shared_universe(species_maps: dict[str, SpeciesPantherMap]) -> frozenset[str]:
    universe: frozenset[str] | None = None
    for species_map in species_maps.values():
        families = species_map.all_families()
        universe = families if universe is None else (universe & families)
    return universe if universe is not None else frozenset()


def _poisson_binomial_pmf(probabilities: list[float]) -> list[float]:
    """Exact PMF of the number of successes among independent Bernoulli
    trials with heterogeneous probabilities, by direct convolution."""
    pmf = [1.0]
    for p in probabilities:
        nxt = [0.0] * (len(pmf) + 1)
        for k, mass in enumerate(pmf):
            nxt[k] += mass * (1.0 - p)
            nxt[k + 1] += mass * p
        pmf = nxt
    return pmf


def conservation_spectrum(
    species_maps: dict[str, SpeciesPantherMap],
    persulfidated: dict[str, set[str]],
) -> ConservationSpectrum:
    """Count shared ortholog subfamilies by how many species persulfidate them.

    ``species_maps`` insertion order is preserved in ``species`` so callers
    control the reporting order; the returned family lists are sorted, so
    the whole result is deterministic.
    """
    species = tuple(species_maps)
    universe = _shared_universe(species_maps)
    universe_size = len(universe)

    per_species: dict[str, frozenset[str]] = {
        name: species_maps[name].persulfidated_families(persulfidated[name]) & universe
        for name in species
    }

    counts: Counter[str] = Counter()
    for families in per_species.values():
        counts.update(families)

    histogram: Counter[int] = Counter()
    for family in universe:
        histogram[counts.get(family, 0)] += 1

    if universe_size > 0:
        rates = [len(per_species[name]) / universe_size for name in species]
        pmf = _poisson_binomial_pmf(rates)
        expected = {k: universe_size * mass for k, mass in enumerate(pmf)}
    else:
        expected = {k: 0.0 for k in range(len(species) + 1)}

    conserved_in_all = tuple(
        sorted(family for family in universe if counts.get(family, 0) == len(species))
    )

    return ConservationSpectrum(
        species=species,
        universe_size=universe_size,
        persulfidated_families_by_species={
            name: len(per_species[name]) for name in species
        },
        observed_by_species_count=dict(sorted(histogram.items())),
        expected_by_species_count=expected,
        conserved_in_all=conserved_in_all,
    )


@dataclass(frozen=True)
class ConservationSpectrumPermutationResult:
    species_count: int
    observed: int
    expected_under_independence: float
    mean_null: float
    p_value: float
    n_perm: int
    seed: int


def conservation_spectrum_permutation_test(
    species_maps: dict[str, SpeciesPantherMap],
    persulfidated: dict[str, set[str]],
    *,
    species_count: int | None = None,
    n_perm: int = 1000,
    seed: int = 0,
) -> ConservationSpectrumPermutationResult:
    """Permutation null for "families persulfidated in at least ``species_count``
    species", defaulting to all of them.

    Each replicate re-draws every species' persulfidated-family set uniformly
    without replacement from the shared universe, holding that species' family
    count fixed. This preserves all marginals and destroys only the
    cross-species alignment, which is exactly the hypothesis under test. The
    p-value is add-one smoothed, so it is bounded in ``[1/(n_perm+1), 1]`` and
    never reported as exactly zero.
    """
    species = tuple(species_maps)
    threshold = len(species) if species_count is None else species_count
    universe = _shared_universe(species_maps)
    universe_size = len(universe)

    per_species_families = [
        species_maps[name].persulfidated_families(persulfidated[name]) & universe
        for name in species
    ]
    per_species_counts = [len(families) for families in per_species_families]

    observed_counts: Counter[str] = Counter()
    for families in per_species_families:
        observed_counts.update(families)
    observed = sum(1 for c in observed_counts.values() if c >= threshold)

    spectrum = conservation_spectrum(species_maps, persulfidated)
    expected = sum(
        mass for k, mass in spectrum.expected_by_species_count.items() if k >= threshold
    )

    rng = random.Random(seed)
    indices = list(range(universe_size))
    null_stats: list[int] = []
    for _ in range(n_perm):
        drawn: Counter[int] = Counter()
        for k in per_species_counts:
            drawn.update(rng.sample(indices, k) if k else ())
        null_stats.append(sum(1 for c in drawn.values() if c >= threshold))

    at_least_observed = sum(1 for stat in null_stats if stat >= observed)
    p_value = (1 + at_least_observed) / (n_perm + 1)
    mean_null = (sum(null_stats) / n_perm) if n_perm else 0.0

    return ConservationSpectrumPermutationResult(
        species_count=threshold,
        observed=observed,
        expected_under_independence=expected,
        mean_null=mean_null,
        p_value=p_value,
        n_perm=n_perm,
        seed=seed,
    )


@dataclass(frozen=True)
class DetectabilityFloor:
    """Smallest co-persulfidation count that would have been significant.

    A non-significant pairwise result means nothing on its own when one of
    the two datasets is shallow: if independence already predicts fewer than
    one shared family, no attainable observation is significant and "n.s."
    is a power statement. Reporting the floor next to the observation makes
    that explicit instead of leaving the reader to infer it.
    """

    alpha: float
    observed: int
    expected_under_independence: float
    minimum_significant_count: int | None
    minimum_significant_fold_enrichment: float | None


def conservation_detectability_floor(
    result: ConservationTestResult,
    alpha: float = 0.05,
) -> DetectabilityFloor:
    """Invert the one-sided hypergeometric tail: the smallest overlap count
    reaching ``p < alpha`` under ``result``'s marginals, or ``None`` when even
    the maximum attainable overlap cannot."""
    n = result.shared_families
    fam_a = result.persulfidated_families_a
    fam_b = result.persulfidated_families_b
    expected = result.expected_co_persulfidated_under_independence

    minimum: int | None = None
    if n > 0:
        # The upper tail is non-increasing in the observed count, so the first
        # count that clears alpha is the floor.
        for count in range(0, min(fam_a, fam_b) + 1):
            tail = _hypergeometric_upper_tail(
                successes_drawn=count,
                total_successes=fam_b,
                total_draws=fam_a,
                population=n,
            )
            if tail < alpha:
                minimum = count
                break

    fold = minimum / expected if minimum is not None and expected > 0 else None
    return DetectabilityFloor(
        alpha=alpha,
        observed=result.co_persulfidated_families,
        expected_under_independence=expected,
        minimum_significant_count=minimum,
        minimum_significant_fold_enrichment=fold,
    )
