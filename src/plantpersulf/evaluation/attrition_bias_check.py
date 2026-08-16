"""Attrition-bias check for PXD072089 (2026-07-23, self-review response).

The rice persulfidome dataset (PXD072089) lost ~55% of its originally
reported UniProt accessions to a 2025-2026 TrEMBL cleanup (documented in
``docs/phase_z_evidence_audit.md`` §5.3) — accessions were deleted, not
merged, so those peptides/sites cannot be coordinate-verified against the
current reference proteome and are dropped.

This module tests a specific, concrete worry raised in self-review: is
that attrition **compositionally biased** with respect to redox/
oxidoreductase functional annotation — the functional class driving the
§5.5.1 family-level enrichment and §5.5.2 structural-context findings? If
proteins with deleted accessions are disproportionately *non*-redox
compared to surviving proteins, the verified set could be artificially
redox-enriched as an artifact of differential data loss, not a real
biological signal. If the two groups have statistically indistinguishable
redox-keyword composition, the attrition is not a plausible confound for
those findings.

**Method**: classify each protein's UniProt-style functional description
(the free-text field carried by SD01/SD04 even for later-deleted
accessions — the description was captured at the time of the original MS
search, before the accession was deleted) as redox-related via a fixed,
documented keyword list. Build a 2x2 contingency table (kept vs. dropped)
x (redox-related vs. not) and run a two-sided Fisher's exact test — exact
hypergeometric enumeration via ``math.comb``, no scipy, consistent with
this project's other from-scratch statistics
(``cross_species_conservation.py``, ``permutation.py``).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb

# Fixed, documented keyword list for "redox/oxidoreductase-class enzyme".
# Not exhaustive — a reasonable, auditable proxy for the functional class
# implicated by the §5.5.1 triple-conserved-family identities (NADH
# dehydrogenase, sorbitol dehydrogenase, aminotransferase, etc.), applied
# uniformly to both the kept and dropped groups so any bias in matching is
# shared, not differential.
REDOX_KEYWORDS: tuple[str, ...] = (
    "DEHYDROGENASE",
    "OXIDASE",
    "OXIDOREDUCTASE",
    "REDUCTASE",
    "PEROXIDASE",
    "PEROXIREDOXIN",
    "CATALASE",
    "THIOREDOXIN",
    "GLUTATHIONE",
    "FERREDOXIN",
    "CYTOCHROME",
    "SUPEROXIDE DISMUTASE",
    "NAD(P)",
    "FLAVIN",
    "AMINOTRANSFERASE",
)


def is_redox_related(description: str) -> bool:
    """Keyword match on the functional-annotation portion of a UniProt-style
    description (the text before ``OS=``), case-insensitive substring
    match against ``REDOX_KEYWORDS``."""
    functional_text = description.split("OS=")[0].upper()
    return any(kw in functional_text for kw in REDOX_KEYWORDS)


def _hypergeom_pmf(k: int, successes: int, draws: int, population: int) -> float:
    return (
        comb(successes, k)
        * comb(population - successes, draws - k)
        / comb(population, draws)
    )


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher's exact test p-value for the 2x2 table
    ``[[a, b], [c, d]]``: sums the hypergeometric probability of every
    table sharing the same margins whose probability is no greater than
    the observed table's — the standard two-sided definition. Exact via
    ``math.comb``, no scipy.
    """
    row1 = a + b
    col1 = a + c
    population = row1 + c + d
    draws = row1
    successes = col1

    observed_p = _hypergeom_pmf(a, successes, draws, population)
    lo = max(0, draws - (population - successes))
    hi = min(draws, successes)

    total = 0.0
    # Small relative tolerance guards against floating-point noise placing
    # the observed table itself just outside its own inclusion threshold.
    threshold = observed_p * (1.0 + 1e-9)
    for k in range(lo, hi + 1):
        p = _hypergeom_pmf(k, successes, draws, population)
        if p <= threshold:
            total += p
    return min(1.0, total)


@dataclass(frozen=True)
class AttritionBiasResult:
    kept_redox: int
    kept_non_redox: int
    dropped_redox: int
    dropped_non_redox: int
    kept_redox_fraction: float
    dropped_redox_fraction: float
    odds_ratio: float | None
    p_value: float


def attrition_bias_test(
    kept_descriptions: list[str],
    dropped_descriptions: list[str],
) -> AttritionBiasResult:
    """Is redox-keyword composition different between kept (accession
    survives in the current reference proteome) and dropped (accession
    deleted) proteins?
    """
    kept_redox = sum(1 for d in kept_descriptions if is_redox_related(d))
    kept_non_redox = len(kept_descriptions) - kept_redox
    dropped_redox = sum(1 for d in dropped_descriptions if is_redox_related(d))
    dropped_non_redox = len(dropped_descriptions) - dropped_redox

    odds_ratio = (
        (kept_redox * dropped_non_redox) / (kept_non_redox * dropped_redox)
        if kept_non_redox > 0 and dropped_redox > 0
        else None
    )
    p_value = fisher_exact_two_sided(
        kept_redox, kept_non_redox, dropped_redox, dropped_non_redox
    )

    return AttritionBiasResult(
        kept_redox=kept_redox,
        kept_non_redox=kept_non_redox,
        dropped_redox=dropped_redox,
        dropped_non_redox=dropped_non_redox,
        kept_redox_fraction=(
            kept_redox / len(kept_descriptions) if kept_descriptions else 0.0
        ),
        dropped_redox_fraction=(
            dropped_redox / len(dropped_descriptions) if dropped_descriptions else 0.0
        ),
        odds_ratio=odds_ratio,
        p_value=p_value,
    )
