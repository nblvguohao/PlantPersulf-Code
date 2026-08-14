"""Figure 1 — Cross-kingdom conservation of persulfidation targeting.

Contract:
- Core conclusion: co-persulfidation of orthologous PANTHER subfamilies is
  enriched beyond chance in two of three plant/fungal pairwise comparisons
  after correction, and eight subfamilies (mostly redox/central-carbon
  metabolism enzymes) are targeted in Arabidopsis, rice and Magnaporthe
  alike; tomato's pairwise cells are underpowered, not negative, and the
  figure shows the detectability floor alongside each one so that is not
  left to the reader to infer.
- Evidence chain: (a) fold enrichment per species pair with Bonferroni/BH
  significance and detectability floor; (b) observed vs independence-
  expected count of subfamilies persulfidated in >=k species; (c) identity
  of the eight subfamilies persulfidated in >=3 species.
- Archetype: dumbbell/dot plot + grouped bar + text list.
- Data: results/cross_species_conservation/conservation_v2.json.
- Export: double column, editable text, 600 dpi.
"""

import json

import matplotlib.pyplot as plt

from _style import DOUBLE_COL, PALETTE, REPO, apply_style, panel_label, save

DATA = REPO / "results/cross_species_conservation/conservation_v2.json"

# Table order (matches manuscript Table 1): no-tomato pairs first, then the
# three underpowered tomato pairs.
PAIR_ORDER = [
    ("Arabidopsis", "Rice"),
    ("Rice", "Magnaporthe"),
    ("Arabidopsis", "Magnaporthe"),
    ("Arabidopsis", "Tomato"),
    ("Rice", "Tomato"),
    ("Tomato", "Magnaporthe"),
]
PAIR_LABELS = {
    ("Arabidopsis", "Rice"): "Arabidopsis x Rice",
    ("Rice", "Magnaporthe"): "Rice x Magnaporthe",
    ("Arabidopsis", "Magnaporthe"): "Arabidopsis x Magnaporthe",
    ("Arabidopsis", "Tomato"): "Arabidopsis x Tomato",
    ("Rice", "Tomato"): "Rice x Tomato",
    ("Tomato", "Magnaporthe"): "Tomato x Magnaporthe",
}


def _pair_key(row: dict) -> tuple:
    return (row["species_a"], row["species_b"])


def main() -> None:
    apply_style()
    doc = json.loads(DATA.read_text(encoding="utf-8"))
    pairwise = {_pair_key(r): r for r in doc["pairwise_tests"]}
    assert set(pairwise) | {(b, a) for a, b in pairwise} >= set(PAIR_ORDER)
    spectrum = doc["conservation_spectrum"]
    assert spectrum["universe_size"] == 2012

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.62))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1.0], hspace=0.62, wspace=0.30)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    # --- (a) pairwise fold enrichment, dumbbell + detectability floor ------
    ys = list(range(len(PAIR_ORDER)))[::-1]
    ax_a.axvline(1.0, color=PALETTE["neutral_mid"], linewidth=0.8, linestyle="-")
    for y, pair in zip(ys, PAIR_ORDER, strict=True):
        row = pairwise.get(pair) or pairwise[(pair[1], pair[0])]
        fold = row["fold_enrichment"] or 0.0
        floor_fold = row["minimum_significant_fold_enrichment"]
        bonf_sig = row["significant_bonferroni_0.05"]
        bh_sig = row["benjamini_hochberg_q_value"] < 0.05
        ax_a.plot([1.0, fold], [y, y], color=PALETTE["neutral_light"],
                  linewidth=1.2, zorder=1)
        if floor_fold is not None:
            ax_a.plot(
                [floor_fold], [y], marker="|", markersize=9, markeredgewidth=1.3,
                color=PALETTE["red_strong"], zorder=2,
            )
        if bonf_sig:
            color, filled = PALETTE["blue_main"], True
        elif bh_sig:
            color, filled = PALETTE["teal"], True
        else:
            color, filled = PALETTE["neutral_mid"], False
        ax_a.scatter(
            [fold], [y], s=42,
            facecolor=color if filled else "white",
            edgecolor=color, linewidths=1.2, zorder=3,
        )
        p = row["p_value"]
        p_str = f"p={p:.1e}" if p < 0.001 else f"p={p:.2f}"
        ax_a.text(
            max(fold, floor_fold or 0) + 0.18, y, p_str,
            fontsize=6.0, va="center", color=PALETTE["neutral_dark"],
        )
    ax_a.set_yticks(ys)
    ax_a.set_yticklabels([PAIR_LABELS[p] for p in PAIR_ORDER])
    ax_a.set_xlabel(
        "fold enrichment of co-persulfidated subfamilies (independence = 1x)")
    ax_a.set_xlim(-0.15, 5.4)
    ax_a.set_title(
        "pairwise conservation (dot: observed; red tick: detectability floor)",
        fontsize=6.5, color=PALETTE["neutral_dark"], loc="left",
    )
    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="",
                   markerfacecolor=PALETTE["blue_main"],
                   markeredgecolor=PALETTE["blue_main"], label="Bonferroni sig. (m=6)"),
        plt.Line2D([0], [0], marker="o", linestyle="", markerfacecolor=PALETTE["teal"],
                   markeredgecolor=PALETTE["teal"], label="BH sig. only"),
        plt.Line2D([0], [0], marker="o", linestyle="", markerfacecolor="white",
                   markeredgecolor=PALETTE["neutral_mid"], label="not significant"),
    ]
    ax_a.legend(handles=handles, loc="upper right", fontsize=5.8, handletextpad=0.4,
                borderaxespad=0.1)
    panel_label(ax_a, "a")

    # --- (b) n-way spectrum: observed vs expected for k>=2, k>=3 -----------
    perms = {p["species_count"]: p for p in spectrum["permutation_tests"]}
    ks = [2, 3]
    x = list(range(len(ks)))
    width = 0.32
    observed = [perms[k]["observed"] for k in ks]
    expected = [perms[k]["expected_under_independence"] for k in ks]
    ax_b.bar([i - width / 2 for i in x], observed, width=width,
              color=PALETTE["blue_main"], label="observed")
    ax_b.bar([i + width / 2 for i in x], expected, width=width,
              color=PALETTE["neutral_light"], label="expected (independence)")
    ax_b.set_yscale("log")
    ax_b.set_ylim(0.3, 500)
    for i, k in enumerate(ks):
        p = perms[k]["p_value"]
        ax_b.text(i, max(observed[i], expected[i]) * 1.6, f"p={p:.1e}",
                   ha="center", fontsize=6.0, color=PALETTE["neutral_dark"])
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([f"$\\geq${k} species" for k in ks])
    ax_b.set_ylabel("shared subfamilies (log scale)")
    ax_b.legend(fontsize=5.6, loc="upper right", handlelength=1.2, handletextpad=0.4)
    ax_b.set_title(
        "n-way spectrum (universe: 2,012 subfamilies\nshared by all 4 proteomes)",
        fontsize=6.3, color=PALETTE["neutral_dark"], loc="left",
    )
    k4 = perms[4]
    ax_b.text(
        0.5, -0.22,
        f"k$\\geq$4: 0 vs {k4['expected_under_independence']:.3f} expected "
        f"(p={k4['p_value']:.2f}) -- no power given tomato",
        transform=ax_b.transAxes, fontsize=5.3, color=PALETTE["neutral_dark"],
        va="top", ha="center", clip_on=False,
    )
    panel_label(ax_b, "b")

    # --- (c) identity of the 8 subfamilies persulfidated in >=3 species ----
    names = spectrum["conserved_in_all_species"]  # k=4 cell (empty by construction)
    assert names["count"] == 0
    # The manuscript's k>=3 set is exactly the three-species (Arabidopsis x
    # Rice x Magnaporthe) intersection retained from the v1 report.
    triple = doc["triple_conserved_families"]
    assert triple["count"] == 8
    enzyme_names = list(triple["annotated_names"].values())
    ax_c.axis("off")
    ax_c.set_title(
        "8 subfamilies persulfidated in Arabidopsis, rice\nand Magnaporthe alike",
        fontsize=6.5, color=PALETTE["neutral_dark"], loc="left",
    )
    for i, name in enumerate(enzyme_names):
        marker_color = (PALETTE["neutral_mid"]
                        if ("RCC1" in name or "isomerase" in name)
                        else PALETTE["blue_main"])
        ax_c.scatter([0.0], [1 - i / (len(enzyme_names) - 1 + 0.001) * 0.92], s=18,
                     color=marker_color, transform=ax_c.transAxes, clip_on=False)
        ax_c.text(
            0.06, 1 - i / (len(enzyme_names) - 1 + 0.001) * 0.92, name,
            transform=ax_c.transAxes, fontsize=6.0, va="center",
            color=PALETTE["neutral_black"],
        )
    ax_c.text(
        0.0, -0.06, "grey: the two non-redox exceptions (PPIase D-related, RCC1)",
        transform=ax_c.transAxes, fontsize=5.5, color=PALETTE["neutral_mid"],
    )
    panel_label(ax_c, "c")

    for p in save(fig, "fig1_conservation"):
        print(p)


if __name__ == "__main__":
    main()
