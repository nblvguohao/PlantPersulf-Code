"""Figure 6 — Pre-blind calibration against published sites.

Contract:
- Core conclusion: published sites the frozen model never saw rank
  moderately above the candidate-table median but do not concentrate at
  the list top — expectations for the blind cohort are modest by design.
- Evidence chain: (a) percentile ranks of the 20 never-trained kiae271
  tomato sites (79/99 panel overlap annotated); (b) percentile ranks of
  the 7 independent-laboratory Arabidopsis controls.
- Archetype: quantitative grid, two strip panels sharing the 0-100
  percentile axis.
- Data: results/known_controls/*_release_bundle_recovery_v1.rows.tsv.
- Export: double column, editable text, 600 dpi.
"""

import csv

import numpy as np
import matplotlib.pyplot as plt

from _style import DOUBLE_COL, PALETTE, REPO, apply_style, panel_label, save

CONTROLS = REPO / "results/known_controls"


def main() -> None:
    apply_style()
    with open(CONTROLS / "kiae271_release_bundle_recovery_v1.rows.tsv",
              encoding="utf-8") as fh:
        kiae = list(csv.DictReader(fh, delimiter="\t"))
    with open(CONTROLS / "arabidopsis_release_bundle_recovery_v1.rows.tsv",
              encoding="utf-8") as fh:
        arab = list(csv.DictReader(fh, delimiter="\t"))

    assert len(kiae) == 99 and len(arab) == 7
    clean = [r for r in kiae if r["in_candidate_table"] == "True"]
    assert len(clean) == 20
    clean_pct = np.array([float(r["percentile_rank"]) for r in clean])
    assert abs(clean_pct.mean() - 60.5) < 0.5, clean_pct.mean()
    arab_pct = np.array([float(r["percentile_rank"]) for r in arab])
    arab_genes = [r["gene"] for r in arab]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.33),
        gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.26},
    )

    # (a) never-trained tomato sites
    rng = np.random.default_rng(20260813)
    order = np.argsort(clean_pct)
    ys = np.arange(len(clean_pct))
    colors = np.where(clean_pct >= 50, PALETTE["blue_main"], PALETTE["neutral_mid"])
    ax1.hlines(ys, 0, clean_pct[order], color=colors[order], linewidth=0.9, alpha=0.8)
    ax1.scatter(clean_pct[order], ys, s=14, color=colors[order], linewidths=0, zorder=3)
    ax1.axvline(50, color=PALETTE["red_strong"], linestyle="--", linewidth=1.0)
    ax1.text(50.8, len(clean_pct) - 0.6, "table median", fontsize=6.2,
             color=PALETTE["red_strong"], va="top")
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("Percentile within frozen tomato candidate table")
    ax1.set_ylabel("Never-trained kiae271 site (ranked)")
    ax1.set_yticks([])
    n_above = int((clean_pct > 50).sum())
    ax1.set_title("20 kiae271 sites the model never saw", fontsize=6.5,
                  color=PALETTE["neutral_dark"], loc="left")
    ax1.text(
        97, 1.2,
        f"{n_above}/20 above median (mean {clean_pct.mean():.1f})\n"
        "79/99 published sites sit in the training panel\n"
        "0/99 in the candidate-table top 2,000",
        ha="right", va="bottom", fontsize=6.0, color=PALETTE["neutral_dark"],
        linespacing=1.45,
    )
    panel_label(ax1, "a")

    # (b) Arabidopsis independent-laboratory controls
    order = np.argsort(arab_pct)
    ys = np.arange(len(arab_pct))
    ax2.hlines(ys, 0, arab_pct[order], color=PALETTE["green_3"], linewidth=0.9, alpha=0.85)
    ax2.scatter(arab_pct[order], ys, s=16, color=PALETTE["green_3"], linewidths=0, zorder=3)
    ax2.axvline(50, color=PALETTE["red_strong"], linestyle="--", linewidth=1.0)
    ax2.set_yticks(ys)
    ax2.set_yticklabels([arab_genes[i] for i in order], fontsize=6)
    ax2.set_xlim(0, 100)
    ax2.set_xlabel("Percentile within Arabidopsis candidate table")
    ax2.set_title("7 independent-laboratory Arabidopsis controls",
                  fontsize=6.5, color=PALETTE["neutral_dark"], loc="left")
    ax2.text(
        97, 0.2,
        f"{int((arab_pct > 50).sum())}/7 above median "
        f"(mean {arab_pct.mean():.1f})",
        ha="right", va="bottom", fontsize=6.0, color=PALETTE["neutral_dark"],
    )
    panel_label(ax2, "b")

    for p in save(fig, "fig6_preblind_calibration"):
        print(p)


if __name__ == "__main__":
    main()
