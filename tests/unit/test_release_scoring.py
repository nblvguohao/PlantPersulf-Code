"""Unit tests for the shared frozen-bundle scoring helpers."""

from __future__ import annotations

import pytest

from plantpersulf.evaluation.release_scoring import (
    ScoredSite,
    cys_positions,
    make_site_row,
    score_sites,
)


def test_cys_positions_returns_1based_positions() -> None:
    assert cys_positions("MSCSAC") == (3, 6)


def test_cys_positions_empty_sequence() -> None:
    assert cys_positions("MKLV") == ()


def test_cys_positions_leading_and_trailing_cys() -> None:
    assert cys_positions("CAGCC") == (1, 4, 5)


def test_make_site_row_builds_development_unlabeled_row() -> None:
    row = make_site_row("tomato", "A0A3Q7EW23", 206)
    assert row.species == "tomato"
    assert row.protein_accession == "A0A3Q7EW23"
    assert row.cys_position == 206
    assert row.label == "unlabeled"
    assert row.global_protein_id == "tomato|A0A3Q7EW23"
    assert row.split == "development"
    assert row.study_accessions == ()


def test_score_sites_assembles_results_in_input_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The structure branch is masked and every site gets score + uncertainty."""

    class FakeOutput:
        scores = [0.42, 0.17]
        uncertainty = [0.01, 0.02]

    captured: list[dict[str, object]] = []

    def fake_scorer(bundle: object, features: object, *, device_name: str | None):
        captured.append(
            {
                "sequence": features.sequence,
                "structure_mask": features.structure_mask,
                "device": device_name,
            }
        )
        return FakeOutput()

    import plantpersulf.evaluation.release_scoring as mod

    monkeypatch.setattr(mod, "score_structure_ranker_bundle", fake_scorer)

    sites = [
        ("tomato", "A0A3Q7EW23", 206),
        ("tomato", "A0A3Q7EW23", 209),
        ("tomato", "A0A3Q7EW23", 206),  # duplicate collapses
    ]
    proteomes = {"tomato": {"A0A3Q7EW23": "C" * 250}}
    bundle = object()
    result = score_sites(sites, bundle=bundle, proteomes=proteomes)

    assert len(captured) == 1
    assert captured[0]["device"] == "cpu"
    assert captured[0]["structure_mask"] == [False, False]
    # sequence features: 3 columns (hydrophobicity, cys_density, charge density)
    assert len(captured[0]["sequence"][0]) == 3

    assert set(result) == {("tomato", "A0A3Q7EW23", 206), ("tomato", "A0A3Q7EW23", 209)}
    assert result[("tomato", "A0A3Q7EW23", 206)].score == 0.42
    assert result[("tomato", "A0A3Q7EW23", 209)].uncertainty == 0.02
    assert isinstance(result[("tomato", "A0A3Q7EW23", 206)], ScoredSite)


def test_score_sites_fails_closed_on_non_cys_coordinate(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The registered-proteome Cys contract must reject a non-Cys position."""
    import plantpersulf.evaluation.release_scoring as mod

    monkeypatch.setattr(
        mod,
        "score_structure_ranker_bundle",
        lambda bundle, features, *, device_name=None: None,
    )
    sites = [("tomato", "A0A3Q7EW23", 2)]  # S, not C (position 2)
    proteomes = {"tomato": {"A0A3Q7EW23": "MS" + "C" * 248}}
    with pytest.raises(RuntimeError, match="unverified v2 Cys coordinate"):
        score_sites(sites, bundle=object(), proteomes=proteomes)
