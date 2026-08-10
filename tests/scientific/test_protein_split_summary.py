"""Analysis pipeline for the protein-split display track.

The frozen protein-split experiment (`pu_ranker_protein_split_v1`) is the
literature-comparable display track (Sul-BertGRU-comparable regime: random
20% protein holdout per seed, 10 seeds, no homology control). Its numbers are
structurally excluded from Gate 2 (see
tests/release/test_gate2_ignores_within_dataset_split_metrics.py); this test
locks in the summary-generation contract: fail-closed input validation,
per-ablation 10-seed statistics, paired-delta CI via the same
``paired_delta_ci`` used by the Gate-2 structure-gain measurement, and the
Route A/B handling of the release-vs-baseline delta.
"""

from __future__ import annotations

import json
from pathlib import Path

from scripts.analyze_protein_split import (
    ABLATIONS,
    N_SEEDS_EXPECTED,
    ROUTE_B_WORDING,
    build_summary,
    load_metrics,
    parse_model_name,
    per_ablation_summary,
)

CONFIG_TEXT = """\
version: 1

experiment:
  name: pu_ranker_protein_split_v1
  benchmark_labels: data/processed/benchmark_v1/sites.tsv
  reference_proteome: data/raw/references/arabidopsis_ref_proteome_v1.fasta

splits:
  mode: random_protein
  test_ratio: 0.2

features:
  - name: sequence
    module: plantpersulf.features.sequence

models:
  - structure_ranker

ablations:
  - name: sequence_only
    use_esm: false
    use_structure: false
    use_study_context: false
  - name: seq_structure
    use_esm: false
    use_study_context: false

evaluation:
  seeds: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

subsample:
  unlabeled_per_positive: 20
  seed: 12345

output:
  directory: results/experiments/pu_ranker_protein_split_v1

limitation: |
  RANDOM PROTEIN-LEVEL SPLIT (within-dataset, no homology control) — NOT
  cross-study validation. Reported only to be comparable with published
  cysteine-PTM predictors; must never support a cross-study claim.
"""


def _synthetic_metrics(
    tmp_path: Path,
    *,
    n_ablations: int = 7,
    drop_seed: int | None = None,
    tag: str = "structure_ranker",
) -> Path:
    path = tmp_path / "metrics.tsv"
    lines = ["model\tseed\tval_ap\ttest_ap\ttest_recall_10\ttest_recall_50\ttest_mrr"]
    for ab in ABLATIONS[:n_ablations]:
        for seed in range(N_SEEDS_EXPECTED):
            if drop_seed is not None and seed == drop_seed:
                continue
            # deterministic synthetic AP: seq_structure highest, others lower
            base = 0.10 if ab == "seq_structure" else 0.05
            ap = base + seed * 0.001
            lines.append(
                f"protein_split_seed{seed}|{tag}:{ab}\t{seed}\t0.5\t"
                f"{ap:.6f}\t0.1\t0.2\t0.3"
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _config_file(tmp_path: Path) -> Path:
    path = tmp_path / "config.yaml"
    path.write_text(CONFIG_TEXT, encoding="utf-8")
    return path


def test_parse_model_name_shapes() -> None:
    seed, model, ab = parse_model_name(
        "protein_split_seed7|structure_ranker:seq_structure"
    )
    assert (seed, model, ab) == (7, "structure_ranker", "seq_structure")


def test_parse_model_name_rejects_wrong_tag() -> None:
    import pytest

    with pytest.raises(ValueError, match="unexpected model tag"):
        parse_model_name("cluster_v1|structure_ranker:full")


def test_load_metrics_accepts_full_grid(tmp_path: Path) -> None:
    rows = load_metrics(_synthetic_metrics(tmp_path))
    assert len(rows) == len(ABLATIONS) * N_SEEDS_EXPECTED


def test_load_metrics_fails_closed_on_missing_seed(tmp_path: Path) -> None:
    import pytest

    rows_path = _synthetic_metrics(tmp_path, drop_seed=4)
    with pytest.raises(ValueError, match="seeds"):
        load_metrics(rows_path)


def test_per_ablation_summary_shapes(tmp_path: Path) -> None:
    rows = load_metrics(_synthetic_metrics(tmp_path))
    summary = per_ablation_summary(rows)
    assert set(summary) == set(ABLATIONS)
    assert "test_ap_mean" in summary["seq_structure"]
    # seq_structure synthetic APs are strictly above sequence_only's
    assert (
        summary["seq_structure"]["test_ap_mean"]
        > summary["sequence_only"]["test_ap_mean"]
    )


def test_build_summary_structure_gain_and_gap_record(tmp_path: Path) -> None:
    out = tmp_path / "out"
    summary = build_summary(
        _synthetic_metrics(tmp_path),
        _config_file(tmp_path),
        baseline_path=None,
        out_dir=out,
    )
    assert summary["admissible_for_gate2"] is False
    assert summary["display_role"] == "literature_comparable"
    assert "limitation" in summary and "RANDOM PROTEIN-LEVEL" in summary["limitation"]
    assert summary["base_rate"] == 0.047619047619047616
    sg = summary["structure_gain"]
    assert sg["n_deltas"] == N_SEEDS_EXPECTED
    assert sg["point"] > 0.03  # synthetic: 0.10 vs 0.05 base
    gap = summary["release_vs_baseline_delta"]
    assert gap["present"] is False
    assert gap["locked_wording"] == ROUTE_B_WORDING
    # files written
    parsed = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert parsed["experiment"] == "pu_ranker_protein_split_v1"
    assert (out / "summary.md").exists()


def test_build_summary_with_baseline_delta(tmp_path: Path) -> None:
    out = tmp_path / "out2"
    baseline = _synthetic_metrics(tmp_path, n_ablations=1, tag="pu_logistic").rename(
        tmp_path / "baseline.tsv"
    )
    # rewrite baseline rows with a constant AP
    lines = ["model\tseed\tval_ap\ttest_ap\ttest_recall_10\ttest_recall_50\ttest_mrr"]
    for seed in range(N_SEEDS_EXPECTED):
        lines.append(
            f"protein_split_seed{seed}|pu_logistic\t{seed}\t0.5\t0.040000\t0.1\t0.2\t0.3"
        )
    baseline.write_text("\n".join(lines) + "\n", encoding="utf-8")
    summary = build_summary(
        _synthetic_metrics(tmp_path),
        _config_file(tmp_path),
        baseline_path=baseline,
        out_dir=out,
    )
    d = summary["release_vs_baseline_delta"]
    assert d["present"] is True
    assert d["n_deltas"] == N_SEEDS_EXPECTED
    assert d["point"] > 0.05  # 0.10+ vs 0.04
    assert d["lower"] > 0.0
