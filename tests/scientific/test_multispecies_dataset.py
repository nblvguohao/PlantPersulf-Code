"""Multi-species joint-training dataset construction (RED).

Merges the four registered site-level evidence sources (Arabidopsis
benchmark positives, tomato kiae271, rice PXD072089, Magnaporthe
PXD063170) into one PU dataset with PANTHER family-level cross-species
homology grouping for leakage-safe CV.
"""

from __future__ import annotations

from pathlib import Path

from plantpersulf.proteomics.multispecies_dataset import (
    MultispeciesSite,
    build_multispecies_sites,
    family_grouped_folds,
    family_stats_from_panther,
    load_panther_family_ids,
)

# --- PANTHER parsing --------------------------------------------------------


def test_panther_family_ids_take_family_level_first_sorted(tmp_path: Path) -> None:
    p = tmp_path / "panther.tsv"
    p.write_text(
        "Entry\tPANTHER\n"
        "P1\tPTHR10000:SF3;PTHR10000;\n"
        "P2\tPTHR20000:SF1;PTHR30000;\n"
        "P3\t\n",
        encoding="utf-8",
    )
    fams = load_panther_family_ids(p)
    assert fams["P1"] == "PTHR10000"
    # two families -> deterministic sorted-first pick
    assert fams["P2"] == "PTHR20000"
    assert "P3" not in fams


# --- merge ------------------------------------------------------------------


def _site(acc: str, pos: int, species: str, family: str) -> MultispeciesSite:
    return MultispeciesSite(
        protein_accession=acc,
        cys_position=pos,
        species=species,
        study_accession=f"STUDY_{species}",
        panther_family=family,
    )


def test_build_multispecies_sites_merges_and_assigns_families() -> None:
    arab = [
        _site("AT1", 5, "arabidopsis", "PTHR10000"),
        _site("AT2", 7, "arabidopsis", "PTHR20000"),
    ]
    tomato = [_site("SL1", 9, "tomato", "PTHR10000")]
    rice = [_site("OS1", 3, "rice", "PTHR30000")]
    magna = [_site("MG1", 11, "magnaporthe", "PTHR30000")]
    sites = build_multispecies_sites(
        arabidopsis=arab, tomato=tomato, rice=rice, magnaporthe=magna
    )
    assert len(sites) == 5
    species = {s.species for s in sites}
    assert species == {"arabidopsis", "tomato", "rice", "magnaporthe"}
    by_key = {(s.protein_accession, s.cys_position): s for s in sites}
    assert by_key[("SL1", 9)].panther_family == "PTHR10000"
    # same family across species keeps the same family id
    assert by_key[("OS1", 3)].panther_family == by_key[("MG1", 11)].panther_family


def test_duplicate_accession_position_collapses() -> None:
    a = _site("AT1", 5, "arabidopsis", "PTHR10000")
    b = _site("AT1", 5, "arabidopsis", "PTHR10000")  # duplicate key
    sites = build_multispecies_sites(arabidopsis=[a, b])
    assert len(sites) == 1


def test_same_accession_position_in_different_species_is_not_collapsed() -> None:
    """Deduplicating without species would erase one source species' evidence."""
    arab = _site("P_SHARED", 5, "arabidopsis", "PTHR10000")
    tomato = _site("P_SHARED", 5, "tomato", "PTHR10000")

    sites = build_multispecies_sites(arabidopsis=[arab], tomato=[tomato])

    observed = {
        (site.species, site.protein_accession, site.cys_position) for site in sites
    }
    assert observed == {
        ("arabidopsis", "P_SHARED", 5),
        ("tomato", "P_SHARED", 5),
    }


def test_duplicate_site_preserves_all_study_accessions() -> None:
    """Keeping only the first study would erase the site provenance chain."""
    first = MultispeciesSite("AT1", 5, "arabidopsis", "PXD_A", "F1")
    second = MultispeciesSite("AT1", 5, "arabidopsis", "PXD_B", "F1")

    merged = build_multispecies_sites(arabidopsis=[first, second])

    assert merged[0].study_accessions == ("PXD_A", "PXD_B")


# --- family-grouped CV -------------------------------------------------------


def test_family_grouped_folds_never_split_a_family() -> None:
    sites = [_site(f"P{i}", 5, "arabidopsis", f"F{i % 3}") for i in range(9)]
    folds = family_grouped_folds(sites, n_folds=3, seed=42)
    assert len(folds) == 3
    family_of = {s.protein_accession: s.panther_family for s in sites}
    for train_idx, test_idx in folds:
        test_fams = {family_of[sites[i].protein_accession] for i in test_idx}
        train_fams = {family_of[sites[i].protein_accession] for i in train_idx}
        assert test_fams.isdisjoint(train_fams)


def test_family_grouped_folds_cover_each_site_exactly_once() -> None:
    sites = [_site(f"P{i}", 5, "arabidopsis", f"F{i % 5}") for i in range(10)]
    folds = family_grouped_folds(sites, n_folds=5, seed=0)
    seen: set[int] = set()
    for _, test_idx in folds:
        overlap = seen & set(test_idx)
        assert not overlap
        seen |= set(test_idx)
    assert seen == set(range(10))


def test_family_grouped_folds_are_deterministic_for_fixed_seed() -> None:
    sites = [_site(f"P{i}", 5, "arabidopsis", f"F{i % 4}") for i in range(8)]
    a = family_grouped_folds(sites, n_folds=4, seed=7)
    b = family_grouped_folds(sites, n_folds=4, seed=7)
    assert a == b


# --- conservation features (v1) ---------------------------------------------


def test_family_stats_species_breadth_and_member_count() -> None:
    """Family-level conservation proxies come from the PANTHER annotation
    layer (all proteome members, all species), never from the label set —
    computing them from positive sites only would leak label information."""
    panther_by_species = {
        "arabidopsis": {"AT1": "F1", "AT2": "F2"},
        "tomato": {"SL1": "F1"},
        "rice": {"OS1": "F1", "OS2": "F3"},
        "magnaporthe": {"MG1": "F2"},
    }
    stats = family_stats_from_panther(panther_by_species)
    # F1 present in 3 species, 3 members; F2 in 2 species, 2 members
    assert stats["F1"] == (3, 3)
    assert stats["F2"] == (2, 2)
    assert stats["F3"] == (1, 1)
