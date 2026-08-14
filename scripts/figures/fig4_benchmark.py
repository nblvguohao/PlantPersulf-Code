"""Figure 4 — Within-dataset literature-comparable benchmark.

Contract:
- Core conclusion: under the split protocol customary in the cysteine-PTM
  tool literature, the frozen gated-fusion PU ranker leads five baselines
  on every seed; the per-species pattern tracks training-positive and
  structure-feature availability.
- Evidence chain: (a) per-seed macro AP across six models; (b) per-species
  AP of the ranker against each species' base rate.
- Archetype: quantitative grid, hero = (a).
- Data: literature_random_protein/summary.json (60 runs = 6 models x 10
  seeds). Claim discipline: all text says within-dataset.
- Export: double column, editable text, 600 dpi.
"""

import json
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt

from _style import DOUBLE_COL, PALETTE, REPO, apply_style, panel_label, save

SUMMARY = (
    REPO
    / "results/experiments/multispecies_v2_global_clusters_v11"
    / "literature_random_protein/summary.json"
)

MODEL_LABELS = {
    "structure_ranker": "structure_ranker (frozen)",
    "esm_linear_head": "ESM-2 linear head",
    "xgboost": "XGBoost",
    "pu_logistic": "PU logistic",
    "random_forest": "random forest",
    "sul_bertgru": "Sul-BertGRU (retrained)",
}
SPECIES = ["arabidopsis", "rice", "tomato"]
SPECIES_LABELS = {"arabidopsis": "Arabidopsis", "rice": "rice", "tomato": "tomato"}


def main() -> None:
    apply_style()
    runs = json.loads(SUMMARY.read_text(encoding="utf-8"))["runs"]

    macro = defaultdict(list)
    per_species = defaultdict(lambda: defaultdict(list))
    base_rates = {}
    for run in runs:
        model = run["model"]
        macro[model].append(run["three_crop_macro_average_precision"])
        for sp in SPECIES:
            cell = run["three_crop_primary_test"][sp]
            per_species[model][sp].append(cell["average_precision"])
            base_rates.setdefault(sp, cell["base_rate"])

    order = sorted(macro, key=lambda m: -np.mean(macro[m]))

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.36),
        gridspec_kw={"width_ratios": [1.45, 1.0], "wspace": 0.30},
    )

    # (a) per-seed macro AP
    rng = np.random.default_rng(20260813)
    for i, model in enumerate(order):
        vals = np.asarray(macro[model])
        color = PALETTE["blue_main"] if model == "structure_ranker" else PALETTE["neutral_mid"]
        x = i + rng.uniform(-0.16, 0.16, size=len(vals))
        ax1.scatter(x, vals, s=9, color=color, alpha=0.75, linewidths=0, zorder=3)
        mean, sd = vals.mean(), vals.std()
        ax1.errorbar(
            i, mean, yerr=sd, fmt="_", ms=0, color=PALETTE["neutral_black"],
            elinewidth=1.0, capsize=3, zorder=4,
        )
        ax1.plot([i - 0.25, i + 0.25], [mean, mean], color=PALETTE["neutral_black"],
                 lw=1.6, zorder=4)
    ax1.set_xticks(range(len(order)))
    ax1.set_xticklabels([MODEL_LABELS[m] for m in order], rotation=18, ha="right")
    ax1.set_ylabel("Macro average precision (3 crop species)")
    ax1.set_xlabel("")
    ax1.set_title("within-dataset, random protein holdout, 10 seeds", fontsize=6.5,
                  color=PALETTE["neutral_dark"], loc="left")
    panel_label(ax1, "a")

    # (b) per-species AP of the frozen ranker vs base rate
    vals = per_species["structure_ranker"]
    for i, sp in enumerate(SPECIES):
        ap = np.asarray(vals[sp])
        x = i + rng.uniform(-0.14, 0.14, size=len(ap))
        ax2.scatter(x, ap, s=9, color=PALETTE["teal"], alpha=0.75, linewidths=0, zorder=3)
        ax2.plot([i - 0.25, i + 0.25], [ap.mean(), ap.mean()],
                 color=PALETTE["neutral_black"], lw=1.6, zorder=4)
        ax2.hlines(
            base_rates[sp], i - 0.32, i + 0.32,
            color=PALETTE["red_strong"], linestyles="--", linewidth=1.0, zorder=2,
        )
    ax2.set_xticks(range(len(SPECIES)))
    ax2.set_xticklabels([SPECIES_LABELS[s] for s in SPECIES])
    ax2.set_ylabel("Average precision")
    ax2.set_title("frozen ranker, per species; dashed = positive base rate",
                  fontsize=6.5, color=PALETTE["neutral_dark"], loc="left")
    panel_label(ax2, "b")

    for p in save(fig, "fig4_benchmark"):
        print(p)


if __name__ == "__main__":
    main()
