"""Figure 7 — Preregistered power landscape.

Contract:
- Core conclusion: power was tabulated honestly over the full design grid
  before any blind data existed; at a fixed assay budget the recommended
  K = 200, 2-controls design balances power against feasibility.
- Evidence chain: three heatmaps (control ratios 1 / 2 / 5) of simulated
  power over K x assumed OR at the 2% background confirmation rate; the
  recommended design cell is marked.
- Archetype: quantitative grid (heatmap row + shared colorbar).
- Data: power_table.json cells (144; 10,000 simulations per cell, seed
  20260813).
- Export: double column, editable text, 600 dpi.
"""

import json

import numpy as np
import matplotlib.pyplot as plt

from _style import DOUBLE_COL, PALETTE, REPO, apply_style, save

POWER = REPO / "results/candidates/multispecies_v2_candidate_release_v1/power_table.json"
KS = [50, 100, 200, 300]
ORS = [1.5, 2, 3, 5]
RATIOS = [1, 2, 5]
P0 = 0.02


def main() -> None:
    apply_style()
    cells = json.loads(POWER.read_text(encoding="utf-8"))["cells"]
    assert len(cells) == 144

    fig, axes = plt.subplots(
        1, 3, figsize=(DOUBLE_COL, DOUBLE_COL * 0.30),
        gridspec_kw={"wspace": 0.12},
    )
    im = None
    for ax, ratio in zip(axes, RATIOS):
        grid = np.full((len(KS), len(ORS)), np.nan)
        for c in cells:
            if c["controls_per_candidate"] == ratio and abs(c["p0_control_rate"] - P0) < 1e-9:
                grid[KS.index(c["k"]), ORS.index(c["assumed_odds_ratio"])] = c["power"]
        assert not np.isnan(grid).any()
        im = ax.imshow(grid, cmap="magma", vmin=0, vmax=1, aspect="auto")
        for i, k in enumerate(KS):
            for j, o in enumerate(ORS):
                val = grid[i, j]
                color = "white" if val < 0.55 else PALETTE["neutral_black"]
                ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                        fontsize=6.3, color=color)
        if ratio == 2:
            ax.add_patch(
                plt.Rectangle((ORS.index(3) - 0.5, KS.index(200) - 0.5), 1, 1,
                              fill=False, edgecolor="#22D7E6", linewidth=1.6)
            )
            ax.annotate("recommended", xy=(ORS.index(3), KS.index(200) - 0.5),
                        xytext=(0.1, -0.55), fontsize=6.2, color="#22D7E6",
                        arrowprops=dict(arrowstyle="-", color="#22D7E6", lw=0.8))
        ax.set_xticks(range(len(ORS)))
        ax.set_xticklabels([str(o) for o in ORS])
        ax.set_yticks(range(len(KS)))
        if ratio == 1:
            ax.set_yticklabels(KS)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("assumed odds ratio")
        ax.set_title(f"{ratio} control{'s' if ratio > 1 else ''} per candidate",
                     fontsize=6.5, color=PALETTE["neutral_dark"])
        ax.set_frame_on(False)
    axes[0].set_ylabel("candidate count K")

    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label(
        f"power (one-sided Fisher, α = 0.05)\nbackground rate {P0:.0%}, "
        "10,000 sims per cell",
        fontsize=6.2,
    )
    cbar.ax.tick_params(labelsize=6)

    for p in save(fig, "fig7_power_landscape"):
        print(p)


if __name__ == "__main__":
    main()
