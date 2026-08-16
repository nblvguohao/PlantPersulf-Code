"""Unit tests for deterministic rank uncertainty summaries."""

from plantpersulf.evaluation.rank_ensemble import summarize_rank_runs


def test_rank_summary_uses_run_distribution() -> None:
    summary = summarize_rank_runs(
        ("site-a", "site-b"),
        ((0.9, 0.1), (0.8, 0.2), (0.7, 0.3)),
    )

    assert summary[0].site_key == "site-a"
    assert summary[0].median_rank == 1.0
