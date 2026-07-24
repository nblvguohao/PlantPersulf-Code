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
