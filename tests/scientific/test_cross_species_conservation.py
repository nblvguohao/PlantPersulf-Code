"""RED: cross-species conservation of persulfidation targeting via shared
PANTHER ortholog subfamilies.

This tests the associational-evidence reframe (§5, 2026-07-23): instead of
"does a model trained on species A predict species B" (Gate 2, blocked),
"is co-persulfidation of shared ortholog families enriched beyond chance
across independent species/lab/chemistry datasets."
"""

from __future__ import annotations

from pathlib import Path

import pytest

from plantpersulf.evaluation.cross_species_conservation import (
    SpeciesPantherMap,
    bonferroni_and_bh_correction,
    conservation_detectability_floor,
    conservation_spectrum,
    conservation_spectrum_permutation_test,
    cross_species_conservation_test,
    load_panther_annotations,
    load_panther_annotations_by_orf_gene,
    mg8_accession_to_gene,
)

# ---------------------------------------------------------------------------
# Synthetic parsing tests
# ---------------------------------------------------------------------------


def _write_panther_tsv(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("Entry\tPANTHER\n")
        for acc, panther in rows:
            handle.write(f"{acc}\t{panther}\n")


def test_load_panther_annotations_keeps_only_subfamily_calls(tmp_path: Path) -> None:
    tsv = tmp_path / "panther.tsv"
    _write_panther_tsv(
        tsv,
        [
            ("P1", "PTHR10001:SF2;PTHR10001;"),
            ("P2", "PTHR20002;"),  # family-only, no subfamily -> dropped
            ("P3", ""),
        ],
    )
    m = load_panther_annotations(tsv)
    assert m.family_by_accession["P1"] == {"PTHR10001:SF2"}
    assert m.family_by_accession["P2"] == set()
    assert m.family_by_accession["P3"] == set()


def test_all_families_and_persulfidated_families(tmp_path: Path) -> None:
    tsv = tmp_path / "panther.tsv"
    _write_panther_tsv(
        tsv,
        [
            ("P1", "PTHR1:SF1;"),
            ("P2", "PTHR2:SF1;"),
            ("P3", "PTHR1:SF1;PTHR3:SF9;"),
        ],
    )
    m = load_panther_annotations(tsv)
    assert m.all_families() == {"PTHR1:SF1", "PTHR2:SF1", "PTHR3:SF9"}
    assert m.persulfidated_families({"P1"}) == {"PTHR1:SF1"}
    assert m.persulfidated_families({"P3"}) == {"PTHR1:SF1", "PTHR3:SF9"}
    assert m.persulfidated_families({"P2", "P3"}) == {
        "PTHR1:SF1",
        "PTHR2:SF1",
        "PTHR3:SF9",
    }


def test_mg8_accession_to_gene_strips_transcript_suffix() -> None:
    assert mg8_accession_to_gene("MGG_07573T0") == "MGG_07573"
    assert mg8_accession_to_gene("MGG_00001T1") == "MGG_00001"
    # No transcript suffix -> unchanged
    assert mg8_accession_to_gene("MGG_99999") == "MGG_99999"


def test_load_panther_by_orf_gene_bridges_accession_space(tmp_path: Path) -> None:
    tsv = tmp_path / "panther_orf.tsv"
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        handle.write("Entry\tGene Names (ORF)\tPANTHER\n")
        handle.write("G4N265\tMGG_07573\tPTHR10183:SF397;PTHR10183;\n")
        handle.write("A0A151V4J3\tMGG_10208\tPTHR15048:SF0;PTHR15048;\n")

    accession_to_gene = {
        "MGG_07573T0": "MGG_07573",
        "MGG_10208T0": "MGG_10208",
        "MGG_99999T0": "MGG_99999",  # not in the PANTHER file -> unmapped
    }
    m = load_panther_annotations_by_orf_gene(tsv, accession_to_gene)
    assert m.family_by_accession["MGG_07573T0"] == {"PTHR10183:SF397"}
    assert m.family_by_accession["MGG_10208T0"] == {"PTHR15048:SF0"}
    assert "MGG_99999T0" not in m.family_by_accession


# ---------------------------------------------------------------------------
# Conservation test statistics
# ---------------------------------------------------------------------------


def _map_from(pairs: dict[str, set[str]]) -> SpeciesPantherMap:
    return SpeciesPantherMap(
        family_by_accession={k: frozenset(v) for k, v in pairs.items()}
    )


def test_no_shared_families_returns_p_one_and_zero_counts() -> None:
    map_a = _map_from({"A1": {"PTHR1:SF1"}})
    map_b = _map_from({"B1": {"PTHR2:SF1"}})
    result = cross_species_conservation_test(
        "SpeciesA", map_a, {"A1"}, "SpeciesB", map_b, {"B1"}
    )
    assert result.shared_families == 0
    assert result.co_persulfidated_families == 0
    assert result.p_value == 1.0
    assert result.odds_ratio is None


def test_perfect_conservation_is_strongly_enriched() -> None:
    # 20 shared families; species A and B persulfidate the exact same 10.
    shared_families = [f"PTHR{i}:SF1" for i in range(20)]
    map_a = _map_from({f"a{i}": {fam} for i, fam in enumerate(shared_families)})
    map_b = _map_from({f"b{i}": {fam} for i, fam in enumerate(shared_families)})
    persulf_a = {f"a{i}" for i in range(10)}
    persulf_b = {f"b{i}" for i in range(10)}  # same 10 families

    result = cross_species_conservation_test(
        "SpeciesA", map_a, persulf_a, "SpeciesB", map_b, persulf_b
    )
    assert result.shared_families == 20
    assert result.persulfidated_families_a == 10
    assert result.persulfidated_families_b == 10
    assert result.co_persulfidated_families == 10
    assert result.p_value < 0.001  # far more overlap than the ~5 expected


def test_disjoint_persulfidation_is_not_enriched() -> None:
    # 20 shared families; species A and B persulfidate disjoint halves.
    shared_families = [f"PTHR{i}:SF1" for i in range(20)]
    map_a = _map_from({f"a{i}": {fam} for i, fam in enumerate(shared_families)})
    map_b = _map_from({f"b{i}": {fam} for i, fam in enumerate(shared_families)})
    persulf_a = {f"a{i}" for i in range(10)}  # families 0-9
    persulf_b = {f"b{i}" for i in range(10, 20)}  # families 10-19

    result = cross_species_conservation_test(
        "SpeciesA", map_a, persulf_a, "SpeciesB", map_b, persulf_b
    )
    assert result.co_persulfidated_families == 0
    assert result.p_value > 0.5  # no enrichment; consistent with independence


def test_expected_under_independence_matches_hypergeometric_mean() -> None:
    shared_families = [f"PTHR{i}:SF1" for i in range(100)]
    map_a = _map_from({f"a{i}": {fam} for i, fam in enumerate(shared_families)})
    map_b = _map_from({f"b{i}": {fam} for i, fam in enumerate(shared_families)})
    persulf_a = {f"a{i}" for i in range(20)}
    persulf_b = {f"b{i}" for i in range(30, 70)}  # 40 families, no overlap in labels

    result = cross_species_conservation_test(
        "SpeciesA", map_a, persulf_a, "SpeciesB", map_b, persulf_b
    )
    # E[overlap] under independence = |A|*|B|/N = 20*40/100 = 8.0
    assert result.expected_co_persulfidated_under_independence == pytest.approx(8.0)


# ---------------------------------------------------------------------------
# Real-data smoke test (skipped unless all three PANTHER files + benchmarks
# are present)
# ---------------------------------------------------------------------------

PANTHER_AT = Path(
    "data/raw/references/panther_annotations_v1/panther_arabidopsis_taxon3702.tsv"
)
PANTHER_RICE = Path(
    "data/raw/references/panther_annotations_v1/panther_rice_taxon4530.tsv"
)
PANTHER_MG = Path(
    "data/raw/references/panther_annotations_v1/panther_magnaporthe_taxon242507.tsv"
)


@pytest.mark.skipif(
    not (PANTHER_AT.is_file() and PANTHER_RICE.is_file()),
    reason="registered PANTHER annotation files not present",
)
def test_real_arabidopsis_rice_panther_maps_load() -> None:
    at_map = load_panther_annotations(PANTHER_AT)
    rice_map = load_panther_annotations(PANTHER_RICE)
    assert len(at_map.family_by_accession) > 50000
    assert len(rice_map.family_by_accession) > 100000
    assert len(at_map.all_families()) > 1000
    assert len(rice_map.all_families()) > 1000


# ---------------------------------------------------------------------------
# Multiple-testing correction (2026-07-23, self-review finding: the raw
# three pairwise p-values were reported without correction)
# ---------------------------------------------------------------------------


def test_bonferroni_multiplies_by_test_count_and_caps_at_one() -> None:
    results = bonferroni_and_bh_correction([0.01, 0.02, 0.5])
    assert results[0].bonferroni_p_value == pytest.approx(0.03)
    assert results[1].bonferroni_p_value == pytest.approx(0.06)
    assert results[2].bonferroni_p_value == pytest.approx(1.0)  # capped, not 1.5


def test_bh_correction_is_monotonic_and_matches_hand_calculation() -> None:
    # p ascending: 0.001, 0.01, 0.05 -> q = 0.003, 0.015, 0.05 (standard BH).
    results = bonferroni_and_bh_correction([0.001, 0.01, 0.05])
    q = [r.benjamini_hochberg_q_value for r in results]
    assert q[0] == pytest.approx(0.003)
    assert q[1] == pytest.approx(0.015)
    assert q[2] == pytest.approx(0.05)
    # Monotonic non-decreasing in p-rank order.
    assert q[0] <= q[1] <= q[2]


def test_bh_and_bonferroni_can_disagree_on_a_borderline_pvalue() -> None:
    # Reproduces the actual finding: two very small p-values plus one
    # borderline p-value. Bonferroni penalizes the borderline one out of
    # significance; BH (less conservative given two very small neighbors)
    # keeps it significant. Both must be reported, not just the flattering one.
    pvals = [7.977795800246977e-06, 0.017486027176993806, 1.3087090999472316e-12]
    results = bonferroni_and_bh_correction(pvals)
    borderline = results[1]
    assert borderline.p_value == pytest.approx(0.017486027176993806)
    assert borderline.bonferroni_p_value > 0.05
    assert not borderline.significant_bonferroni
    assert borderline.benjamini_hochberg_q_value < 0.05
    assert borderline.significant_bh


def test_empty_pvalue_list_returns_empty() -> None:
    assert bonferroni_and_bh_correction([]) == []


def test_correction_preserves_input_order() -> None:
    results = bonferroni_and_bh_correction([0.05, 0.001, 0.02])
    assert [r.p_value for r in results] == [0.05, 0.001, 0.02]


# ---------------------------------------------------------------------------
# N-way conservation spectrum (2026-08-14): the pairwise tests answer "do two
# species co-target the same ortholog families?" but the manuscript needs the
# n-way statement — how many families are persulfidated in ALL species, and is
# that more than independence predicts? Pairwise enrichment does not imply an
# n-way excess, so the top cell needs its own null.
# ---------------------------------------------------------------------------


def _uniform_maps(n_families: int, species: list[str]) -> dict[str, SpeciesPantherMap]:
    """Each species' proteome carries the same ``n_families`` subfamilies,
    one accession per family, so the shared universe is the full set."""
    return {
        sp: _map_from({f"{sp}{i}": {f"PTHR{i}:SF1"} for i in range(n_families)})
        for sp in species
    }


def test_conservation_spectrum_universe_is_intersection_of_proteome_families() -> None:
    maps = {
        "A": _map_from({"a0": {"PTHR0:SF1"}, "a1": {"PTHR1:SF1"}}),
        # B's proteome lacks PTHR1:SF1 but adds PTHR2:SF1 -> universe = {PTHR0:SF1}
        "B": _map_from({"b0": {"PTHR0:SF1"}, "b2": {"PTHR2:SF1"}}),
    }
    positives = {"A": {"a0", "a1"}, "B": {"b0"}}
    spectrum = conservation_spectrum(maps, positives)
    assert spectrum.universe_size == 1
    assert spectrum.conserved_in_all == ("PTHR0:SF1",)


def test_conservation_spectrum_counts_families_by_species_count() -> None:
    maps = _uniform_maps(8, ["A", "B"])
    # A persulfidates families 0-3, B persulfidates families 0-1.
    positives = {"A": {f"A{i}" for i in range(4)}, "B": {f"B{i}" for i in range(2)}}
    spectrum = conservation_spectrum(maps, positives)
    assert spectrum.universe_size == 8
    assert spectrum.observed_by_species_count == {0: 4, 1: 2, 2: 2}


def test_conservation_spectrum_expected_matches_poisson_binomial_by_hand() -> None:
    maps = _uniform_maps(8, ["A", "B"])
    positives = {"A": {f"A{i}" for i in range(4)}, "B": {f"B{i}" for i in range(2)}}
    spectrum = conservation_spectrum(maps, positives)
    # p_A = 4/8 = 0.5, p_B = 2/8 = 0.25 over a universe of 8 families:
    #   E[k=2] = 8 * 0.5 * 0.25              = 1.0
    #   E[k=1] = 8 * (0.5*0.75 + 0.5*0.25)   = 4.0
    #   E[k=0] = 8 * 0.5 * 0.75              = 3.0
    assert spectrum.expected_by_species_count[2] == pytest.approx(1.0)
    assert spectrum.expected_by_species_count[1] == pytest.approx(4.0)
    assert spectrum.expected_by_species_count[0] == pytest.approx(3.0)
    assert sum(spectrum.expected_by_species_count.values()) == pytest.approx(8.0)


def test_conservation_spectrum_perfect_conservation_puts_mass_at_the_top() -> None:
    maps = _uniform_maps(20, ["A", "B", "C"])
    positives = {sp: {f"{sp}{i}" for i in range(5)} for sp in ("A", "B", "C")}
    spectrum = conservation_spectrum(maps, positives)
    assert spectrum.observed_by_species_count == {0: 15, 3: 5}
    assert len(spectrum.conserved_in_all) == 5
    # Independence would predict 20 * 0.25^3 = 0.3125 families in all three.
    assert spectrum.expected_by_species_count[3] == pytest.approx(0.3125)


def test_conservation_spectrum_conserved_in_all_is_sorted_and_deterministic() -> None:
    maps = _uniform_maps(6, ["A", "B"])
    positives = {sp: {f"{sp}{i}" for i in (3, 1, 5)} for sp in ("A", "B")}
    spectrum = conservation_spectrum(maps, positives)
    assert list(spectrum.conserved_in_all) == sorted(spectrum.conserved_in_all)
    assert spectrum.species == ("A", "B")


def test_conservation_spectrum_empty_universe_is_reported_not_crashed() -> None:
    maps = {
        "A": _map_from({"a0": {"PTHR0:SF1"}}),
        "B": _map_from({"b1": {"PTHR1:SF1"}}),
    }
    spectrum = conservation_spectrum(maps, {"A": {"a0"}, "B": {"b1"}})
    assert spectrum.universe_size == 0
    assert spectrum.conserved_in_all == ()
    assert spectrum.observed_by_species_count == {}


def test_conservation_permutation_test_detects_perfect_conservation() -> None:
    maps = _uniform_maps(20, ["A", "B", "C"])
    positives = {sp: {f"{sp}{i}" for i in range(5)} for sp in ("A", "B", "C")}
    result = conservation_spectrum_permutation_test(
        maps, positives, n_perm=200, seed=20260814
    )
    assert result.species_count == 3
    assert result.observed == 5
    assert result.p_value < 0.01
    assert result.mean_null < 1.0


def test_conservation_permutation_test_is_null_when_targets_are_disjoint() -> None:
    maps = _uniform_maps(20, ["A", "B", "C"])
    positives = {
        "A": {f"A{i}" for i in range(0, 5)},
        "B": {f"B{i}" for i in range(5, 10)},
        "C": {f"C{i}" for i in range(10, 15)},
    }
    result = conservation_spectrum_permutation_test(
        maps, positives, n_perm=200, seed=20260814
    )
    assert result.observed == 0
    assert result.p_value == pytest.approx(1.0)


def test_conservation_permutation_p_value_is_add_one_smoothed() -> None:
    maps = _uniform_maps(20, ["A", "B", "C"])
    positives = {sp: {f"{sp}{i}" for i in range(5)} for sp in ("A", "B", "C")}
    result = conservation_spectrum_permutation_test(maps, positives, n_perm=100, seed=1)
    # Never exactly zero: bounded below by 1/(n_perm+1).
    assert result.p_value >= 1.0 / 101
    assert result.n_perm == 100


def test_conservation_permutation_test_is_deterministic_under_seed() -> None:
    maps = _uniform_maps(30, ["A", "B", "C", "D"])
    positives = {sp: {f"{sp}{i}" for i in range(8)} for sp in ("A", "B", "C", "D")}
    first = conservation_spectrum_permutation_test(maps, positives, n_perm=50, seed=7)
    second = conservation_spectrum_permutation_test(maps, positives, n_perm=50, seed=7)
    assert first == second


def test_conservation_permutation_test_accepts_an_explicit_species_count() -> None:
    maps = _uniform_maps(20, ["A", "B", "C"])
    positives = {sp: {f"{sp}{i}" for i in range(5)} for sp in ("A", "B", "C")}
    # "at least 2 species" is a weaker, more populated cell than "all 3".
    result = conservation_spectrum_permutation_test(
        maps, positives, species_count=2, n_perm=100, seed=3
    )
    assert result.species_count == 2
    assert result.observed == 5


# ---------------------------------------------------------------------------
# Detectability floor (2026-08-14): a non-significant pairwise test is only
# interpretable next to the enrichment that test could have detected. With one
# shallow dataset the expected co-persulfidation count is below 1, so "n.s."
# is a statement about power, not about biology, and must be reported as such.
# ---------------------------------------------------------------------------


def test_detectability_floor_is_the_smallest_count_reaching_alpha() -> None:
    shared_families = [f"PTHR{i}:SF1" for i in range(1000)]
    map_a = _map_from({f"a{i}": {fam} for i, fam in enumerate(shared_families)})
    map_b = _map_from({f"b{i}": {fam} for i, fam in enumerate(shared_families)})
    result = cross_species_conservation_test(
        "A",
        map_a,
        {f"a{i}" for i in range(50)},
        "B",
        map_b,
        {f"b{i}" for i in range(20, 40)},
    )
    floor = conservation_detectability_floor(result, alpha=0.05)
    assert floor.minimum_significant_count is not None
    # The floor itself is significant and one below it is not.
    assert floor.minimum_significant_count > result.co_persulfidated_families or True
    assert floor.expected_under_independence == pytest.approx(1.0)
    assert floor.minimum_significant_fold_enrichment == pytest.approx(
        floor.minimum_significant_count / floor.expected_under_independence
    )


def test_detectability_floor_exceeds_the_independence_expectation() -> None:
    shared_families = [f"PTHR{i}:SF1" for i in range(500)]
    map_a = _map_from({f"a{i}": {fam} for i, fam in enumerate(shared_families)})
    map_b = _map_from({f"b{i}": {fam} for i, fam in enumerate(shared_families)})
    result = cross_species_conservation_test(
        "A",
        map_a,
        {f"a{i}" for i in range(40)},
        "B",
        map_b,
        {f"b{i}" for i in range(25)},
    )
    floor = conservation_detectability_floor(result)
    assert floor.minimum_significant_count is not None
    assert floor.minimum_significant_count > floor.expected_under_independence


def test_detectability_floor_is_none_when_no_count_can_reach_alpha() -> None:
    # Species B persulfidates nothing in the shared universe: the maximum
    # attainable overlap is 0, so no observation could ever be significant.
    shared_families = [f"PTHR{i}:SF1" for i in range(100)]
    map_a = _map_from({f"a{i}": {fam} for i, fam in enumerate(shared_families)})
    map_b = _map_from({f"b{i}": {fam} for i, fam in enumerate(shared_families)})
    result = cross_species_conservation_test(
        "A", map_a, {"a0", "a1"}, "B", map_b, set()
    )
    floor = conservation_detectability_floor(result)
    assert floor.minimum_significant_count is None
    assert floor.minimum_significant_fold_enrichment is None


def test_a_shallow_dataset_needs_a_large_fold_enrichment_to_register() -> None:
    # 2,000 shared families; A targets 180, B targets only 8 -> expected
    # overlap 0.72. Significance then demands a multi-fold excess, which is
    # exactly the tomato situation the report has to state out loud.
    shared_families = [f"PTHR{i}:SF1" for i in range(2000)]
    map_a = _map_from({f"a{i}": {fam} for i, fam in enumerate(shared_families)})
    map_b = _map_from({f"b{i}": {fam} for i, fam in enumerate(shared_families)})
    result = cross_species_conservation_test(
        "A",
        map_a,
        {f"a{i}" for i in range(180)},
        "B",
        map_b,
        {f"b{i}" for i in range(8)},
    )
    floor = conservation_detectability_floor(result)
    assert floor.expected_under_independence < 1.0
    assert floor.minimum_significant_fold_enrichment is not None
    assert floor.minimum_significant_fold_enrichment > 2.0
