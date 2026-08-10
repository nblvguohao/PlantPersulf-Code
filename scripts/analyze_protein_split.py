"""Standalone analysis for the random protein-split display track.

The frozen protein-split experiment (`pu_ranker_protein_split_v1`) runs under
a Sul-BertGRU-comparable regime (random 20% protein holdout per seed, 10
seeds = 10 repetitions, no homology control). `score_release.py` cannot
express protein splits (it is structurally leave-study-out-only), so this
script is the canonical summary generator for that track.

What it computes:
- per-ablation summary over the 10 seeds (mean ± std of test_ap,
  test_recall_10, test_recall_50, test_mrr) with the subsampled base rate
  reported alongside;
- structure gain: 10 paired per-seed deltas (seq_structure − sequence_only)
  → ``paired_delta_ci`` (same function/defaults as the Gate-2 structure-gain
  measurement, so the two tracks are statistically comparable by
  construction);
- release-vs-baseline effect delta: when ``--baseline-metrics`` is given,
  10 paired per-seed deltas (seq_structure − pu_logistic) → ``paired_delta_ci``;
  when absent, a machine-readable gap record plus the locked limitation
  wording (Route B fallback) is emitted instead.

Fail-closed by design: any row/column/tag invariant violation aborts the
script — no partial summary is ever produced.

The output summary.json carries ``admissible_for_gate2: false`` and the
frozen config's limitation text verbatim. This track is display-only; Gate 2
admissibility is locked by
tests/release/test_gate2_ignores_within_dataset_split_metrics.py.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from plantpersulf.evaluation.effect_size import paired_delta_ci  # noqa: E402

METRICS_COLUMNS = (
    "model",
    "seed",
    "val_ap",
    "test_ap",
    "test_recall_10",
    "test_recall_50",
    "test_mrr",
)
ABLATIONS = (
    "sequence_only",
    "seq_esm",
    "seq_structure",
    "full",
    "no_plddt",
    "no_accessibility",
    "no_study_context",
)
N_SEEDS_EXPECTED = 10
MODEL_TAG = "protein_split_seed"
RELEASE_ABLATION = "seq_structure"
ABLATED_ABLATION = "sequence_only"
BASELINE_MODEL = "pu_logistic"

ROUTE_B_WORDING = (
    "The release-vs-baseline effect delta is not quantifiable on this track "
    "because the baseline arm (pu_logistic) was not part of the frozen "
    "protein-split run; the track reports absolute AP vs base rate only. "
    "This must not be framed as an effect claim."
)


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass(frozen=True)
class ParsedRow:
    seed: int
    model_part: str
    ablation: str
    test_ap: float
    test_recall_10: float
    test_recall_50: float
    test_mrr: float


def parse_model_name(model: str) -> tuple[int, str, str]:
    """Parse ``protein_split_seed<k>|structure_ranker:<ablation>``.

    Baseline rows carry no ablation suffix
    (``protein_split_seed<k>|pu_logistic``) — those parse with ablation="".
    """
    head, rest = model.split("|", 1)
    if not head.startswith(MODEL_TAG):
        raise ValueError(f"unexpected model tag: {model!r}")
    seed = int(head[len(MODEL_TAG) :])
    if ":" in rest:
        model_part, ablation = rest.split(":", 1)
    else:
        model_part, ablation = rest, ""
    return seed, model_part, ablation


def load_metrics(
    path: Path,
    *,
    expect_full_ablation_grid: bool = True,
) -> list[ParsedRow]:
    """Load and validate a protein-split metrics.tsv (fail-closed).

    ``expect_full_ablation_grid=True`` enforces the 7 ablations × 10 seeds
    invariant of the release experiment; baseline runs (single model, no
    ablation suffix) must pass ``False``.
    """
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != METRICS_COLUMNS:
            raise ValueError(
                f"unexpected metrics columns: {tuple(reader.fieldnames or ())}"
            )
        rows = [dict(r) for r in reader]
    parsed: list[ParsedRow] = []
    for r in rows:
        seed, model_part, ablation = parse_model_name(r["model"])
        parsed.append(
            ParsedRow(
                seed=seed,
                model_part=model_part,
                ablation=ablation,
                test_ap=float(r["test_ap"]),
                test_recall_10=float(r["test_recall_10"]),
                test_recall_50=float(r["test_recall_50"]),
                test_mrr=float(r["test_mrr"]),
            )
        )
    if expect_full_ablation_grid:
        # invariant: every ablation x every seed present, exactly once
        seen: dict[str, set[int]] = {}
        for p in parsed:
            seen.setdefault(p.ablation, set()).add(p.seed)
        for ab in ABLATIONS:
            seeds = seen.get(ab, set())
            if seeds != set(range(N_SEEDS_EXPECTED)):
                raise ValueError(
                    f"ablation {ab!r} has seeds {sorted(seeds)}, expected 0..9"
                )
    return parsed


def _seed_ap(rows: list[ParsedRow], ablation: str) -> dict[int, float]:
    return {p.seed: p.test_ap for p in rows if p.ablation == ablation}


def per_ablation_summary(rows: list[ParsedRow]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for ab in ABLATIONS:
        sel = [p for p in rows if p.ablation == ab]
        out[ab] = {
            "test_ap_mean": statistics.mean(p.test_ap for p in sel),
            "test_ap_std": statistics.stdev(p.test_ap for p in sel)
            if len(sel) > 1
            else 0.0,
            "test_ap_min": min(p.test_ap for p in sel),
            "test_ap_max": max(p.test_ap for p in sel),
            "recall10_mean": statistics.mean(p.test_recall_10 for p in sel),
            "recall50_mean": statistics.mean(p.test_recall_50 for p in sel),
            "mrr_mean": statistics.mean(p.test_mrr for p in sel),
        }
    return out


def paired_delta_ci_result(
    deltas: list[float],
) -> dict[str, Any]:
    result = paired_delta_ci(deltas, n_boot=2000, alpha=0.05, seed=0)
    return {
        "present": True,
        "point": result.point,
        "lower": result.lower,
        "upper": result.upper,
        "n_boot": result.n_boot,
        "n_deltas": len(deltas),
    }


def base_rate(subsample_ratio: int) -> float:
    """Subsampled base rate: positives / (positives + ratio * positives)."""
    return 1.0 / (1.0 + subsample_ratio)


def build_summary(
    metrics_path: Path,
    config_path: Path,
    baseline_path: Path | None,
    out_dir: Path,
    *,
    subsample_ratio: int = 20,
) -> dict[str, Any]:
    rows = load_metrics(metrics_path)
    cfg = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    limitation = str(cfg.get("limitation", ""))

    summary: dict[str, Any] = {
        "experiment": cfg["experiment"]["name"],
        "admissible_for_gate2": False,
        "display_role": "literature_comparable",
        "config_sha256": _sha256(config_path),
        "limitation": limitation,
        "metrics_sha256": _sha256(metrics_path),
        "n_seeds": N_SEEDS_EXPECTED,
        "base_rate": base_rate(subsample_ratio),
        "per_ablation": per_ablation_summary(rows),
    }

    release = _seed_ap(rows, RELEASE_ABLATION)
    ablated = _seed_ap(rows, ABLATED_ABLATION)
    structure_deltas = [release[s] - ablated[s] for s in sorted(release)]
    summary["structure_gain"] = paired_delta_ci_result(structure_deltas)

    if baseline_path is not None and baseline_path.is_file():
        base_rows = load_metrics(baseline_path, expect_full_ablation_grid=False)
        base = {p.seed: p.test_ap for p in base_rows if p.model_part == BASELINE_MODEL}
        missing = sorted(set(range(N_SEEDS_EXPECTED)) - set(base))
        if missing:
            raise ValueError(f"baseline metrics missing seeds {missing}")
        effect_deltas = [release[s] - base[s] for s in sorted(release)]
        summary["release_vs_baseline_delta"] = paired_delta_ci_result(effect_deltas)
    else:
        summary["release_vs_baseline_delta"] = {
            "present": False,
            "reason": "baseline arm not run in frozen config",
            "locked_wording": ROUTE_B_WORDING,
        }

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (out_dir / "summary.md").write_text(render_markdown(summary), encoding="utf-8")
    return summary


def render_markdown(summary: dict[str, Any]) -> str:
    lines = [
        f"# Protein-split display track summary — {summary['experiment']}",
        "",
        f"- seeds: {summary['n_seeds']} random protein splits (Sul-BertGRU-"
        "comparable regime, no homology control)",
        f"- subsampled base rate: {summary['base_rate']:.4f}",
        f"- admissible for Gate 2: **{summary['admissible_for_gate2']}**",
        "",
        "| ablation | test_ap mean±std | min-max | recall@10 | recall@50 | MRR |",
        "|---|---|---|---|---|---|",
    ]
    for ab in ABLATIONS:
        s = summary["per_ablation"][ab]
        lines.append(
            f"| {ab} | {s['test_ap_mean']:.4f}±{s['test_ap_std']:.4f} "
            f"| {s['test_ap_min']:.4f}-{s['test_ap_max']:.4f} "
            f"| {s['recall10_mean']:.4f} | {s['recall50_mean']:.4f} "
            f"| {s['mrr_mean']:.4f} |"
        )
    lines += ["", "## structure_gain (seq_structure − sequence_only)"]
    sg = summary["structure_gain"]
    lines.append(
        f"- point {sg['point']:.4f} [95% CI {sg['lower']:.4f}, "
        f"{sg['upper']:.4f}] (n_boot={sg['n_boot']}, n_deltas={sg['n_deltas']})"
    )
    d = summary["release_vs_baseline_delta"]
    if d.get("present", True):
        lines += ["", "## release_vs_baseline_delta (seq_structure − pu_logistic)"]
        lines.append(
            f"- point {d['point']:.4f} [95% CI {d['lower']:.4f}, "
            f"{d['upper']:.4f}] (n_boot={d['n_boot']})"
        )
    else:
        lines += ["", "## release_vs_baseline_delta — GAP"]
        lines.append(f"- reason: {d['reason']}")
        lines.append(f"- locked wording: {d['locked_wording']}")
    lines += [
        "",
        "## limitation (verbatim)",
        "",
        "> " + summary["limitation"].replace("\n", "\n> "),
    ]
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--metrics", required=True, type=Path)
    p.add_argument(
        "--config",
        required=True,
        type=Path,
        help="frozen experiment config (limitation + identity)",
    )
    p.add_argument("--baseline-metrics", type=Path, default=None)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--subsample-ratio", type=int, default=20)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        summary = build_summary(
            args.metrics,
            args.config,
            args.baseline_metrics,
            args.out_dir,
            subsample_ratio=args.subsample_ratio,
        )
    except ValueError as exc:
        print(f"FAIL-CLOSED: {exc}", file=sys.stderr)
        return 1
    print(
        f"wrote summary.json/md to {args.out_dir} (experiment {summary['experiment']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
