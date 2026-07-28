"""Cross-species conservation of persulfidation targeting via shared PANTHER
ortholog subfamilies.

Framing shift (2026-07-23): rather than asking whether a model trained on
one species' sites *predicts* another species' sites (the Gate 2 / transfer
question, structurally blocked by study non-independence — see
``docs/phase_z_evidence_audit.md``), this module asks a different,
data-available question: **is the set of PANTHER ortholog families
containing a persulfidated protein in one species enriched for also
containing a persulfidated protein in an independent species, beyond
chance?** This uses exactly the same three independent-lab/chemistry/
species datasets (Arabidopsis benchmark_v1, PXD072089 rice, PXD063170
Magnaporthe) that Gate 2 rejected as *predictive* evidence, but as
*associational* evidence instead — a well-established comparative-genomics
methodology (ortholog-group co-annotation enrichment), not a downgrade in
rigor.

**Orthology proxy**: PANTHER subfamily IDs (``PTHR#####:SF#``), fetched in
bulk from UniProt's ``xref_panther`` field, are used as a fine-grained
ortholog-group proxy — two proteins from different species sharing a
subfamily ID are treated as belonging to the same ortholog group. This is
coarser than a curated 1:1 orthology call (e.g. reciprocal-best-hit) but
requires no new tooling and has near-complete proteome coverage for all
three species (100% for Arabidopsis/rice at the family level; the
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
