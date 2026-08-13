"""Figure 3 — Frozen tomato candidate release.

Contract:
- Core conclusion: the frozen release scores every non-panel tomato
  cysteine site with an extremely flat top, and matched controls sit
  within the frozen 0.25 SD caliper.
- Evidence chain: (a) score distribution of all 179,736 candidates with
  the top-200 band marked, inset showing the flat top; (b) matched-control
  distance distribution against the caliper.
- Archetype: quantitative grid.
- Data: top_k_candidates.tsv (179,736 rows), matched_controls.tsv (966).
- Export: double column, editable text, 600 dpi.
"""

import csv

import numpy as np
import matplotlib.pyplot as plt

from _style import DOUBLE_COL, PALETTE, REPO, apply_style, panel_label, save

RELEASE = REPO / "results/candidates/multispecies_v2_candidate_release_v1"


def main() -> None:
    apply_style()
    with open(RELEASE / "top_k_candidates.tsv", encoding="utf-8") as fh:
        scores = np.array(
            [float(r["score"]) for r in csv.DictReader(fh, delimiter="\t")]
        )
    with open(RELEASE / "matched_controls.tsv", encoding="utf-8") as fh:
        dists = np.array(
            [float(r["sequence_distance_sd"]) for r in csv.DictReader(fh, delimiter="\t")]
        )
    assert len(scores) == 179736, len(scores)
    assert len(dists) == 966, len(dists)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.34),
        gridspec_kw={"wspace": 0.28},
    )

    # (a) full score distribution + flat-top inset
    ax1.hist(scores, bins=200, color=PALETTE["blue_secondary"], edgecolor="none")
    top200 = np.sort(scores)[-200:]
    ax1.axvspan(top200.min(), scores.max(), color=PALETTE["red_strong"], alpha=0.25)
    ax1.set_yscale("log")
    ax1.set_xlabel("Frozen model score")
    ax1.set_ylabel("Sites (log scale)")
    ax1.set_title(f"tomato candidate table, n = {len(scores):,}", fontsize=6.5,
                  color=PALETTE["neutral_dark"], loc="left")
    ax1.annotate("top 200", xy=((top200.min() + scores.max()) / 2, 1.6),
                 ha="center", fontsize=6.5, color=PALETTE["red_strong"])

    inset = ax1.inset_axes([0.06, 0.42, 0.52, 0.52])
    inset.plot(np.arange(1, 201), top200[::-1], color=PALETTE["blue_main"], lw=1.2)
    inset.set_xlabel("rank", fontsize=6)
    inset.set_ylabel("score", fontsize=6)
    inset.set_title("top 200, flat head", fontsize=6, color=PALETTE["neutral_dark"])
    inset.tick_params(labelsize=5.5)
    panel_label(ax1, "a")

    # (b) matched-control distances
    ax2.hist(dists, bins=40, color=PALETTE["teal"], edgecolor="none")
    ax2.axvline(0.25, color=PALETTE["red_strong"], linestyle="--", linewidth=1.0)
    ax2.text(0.255, ax2.get_ylim()[1] * 0.9, "caliper 0.25 SD", fontsize=6.5,
             color=PALETTE["red_strong"])
    ax2.set_xlabel("Control-candidate distance (z-scored SD)")
    ax2.set_ylabel("Matched pairs")
    ax2.set_title(f"matched controls, n = {len(dists)}", fontsize=6.5,
                  color=PALETTE["neutral_dark"], loc="left")
    panel_label(ax2, "b")

    for p in save(fig, "fig3_frozen_release"):
        print(p)


if __name__ == "__main__":
    main()
