"""Unit tests for the accuracy-campaign runner's pure helpers."""

from __future__ import annotations

import json

import pytest

from scripts.accuracy_campaign.run_literature_ranker_scores import (
    load_expected_panel_sha256s,
    per_species_average_precision,
    write_ranker_scores_tsv,
)


def test_load_expected_panel_sha256s_returns_seed_map(tmp_path) -> None:
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "runs": [
                    {"seed": 0, "model": "structure_ranker", "panel_sha256": "aaa"},
                    {"seed": 0, "model": "sul_bertgru", "panel_sha256": "aaa"},
                    {"seed": 1, "model": "structure_ranker", "panel_sha256": "bbb"},
                ]
            }
        ),
        encoding="utf-8",
    )
    assert load_expected_panel_sha256s(summary) == {0: "aaa", 1: "bbb"}


def test_load_expected_panel_sha256s_rejects_cross_model_disagreement(
    tmp_path,
) -> None:
    summary = tmp_path / "summary.json"
    summary.write_text(
        json.dumps(
            {
                "runs": [
                    {"seed": 0, "model": "structure_ranker", "panel_sha256": "aaa"},
                    {"seed": 0, "model": "sul_bertgru", "panel_sha256": "bbb"},
                ]
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError):
        load_expected_panel_sha256s(summary)


def test_write_ranker_scores_tsv_writes_aligned_rows(tmp_path) -> None:
    path = tmp_path / "structure_ranker_seed0_test.tsv"
    write_ranker_scores_tsv(
        path,
        partition="test",
        rows=(("arabidopsis|A", 10, "positive"), ("tomato|B", 5, "unlabeled")),
        scores=[0.8, 0.2],
        uncertainty=[0.05, 0.01],
    )
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == (
        "partition\tglobal_protein_id\tcys_position\tlabel\tscore\tuncertainty"
    )
    assert lines[1] == "test\tarabidopsis|A\t10\tpositive\t0.8\t0.05"
    assert lines[2] == "test\ttomato|B\t5\tunlabeled\t0.2\t0.01"


def test_write_ranker_scores_tsv_rejects_misalignment(tmp_path) -> None:
    with pytest.raises(ValueError):
        write_ranker_scores_tsv(
            tmp_path / "x.tsv",
            partition="test",
            rows=(("a", 1, "positive"),),
            scores=[0.5, 0.6],
            uncertainty=[0.1],
        )


def test_per_species_average_precision_separates_species() -> None:
    rows = [
        ("arabidopsis|A", 1, "positive"),
        ("arabidopsis|A", 2, "unlabeled"),
        ("tomato|B", 3, "positive"),
        ("tomato|B", 4, "unlabeled"),
    ]
    ap = per_species_average_precision(rows, [0.9, 0.1, 0.5, 0.4])
    assert ap == {"arabidopsis": pytest.approx(1.0), "tomato": pytest.approx(1.0)}


def test_per_species_average_precision_skips_species_without_positives() -> None:
    rows = [
        ("arabidopsis|A", 1, "positive"),
        ("arabidopsis|A", 2, "unlabeled"),
        ("tomato|B", 3, "unlabeled"),
    ]
    ap = per_species_average_precision(rows, [0.9, 0.1, 0.5])
    assert set(ap) == {"arabidopsis"}
