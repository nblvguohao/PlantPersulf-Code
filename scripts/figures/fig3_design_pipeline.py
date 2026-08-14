"""Figure 3 — Study design and gate pipeline (schematic-led composite).

Contract:
- Core conclusion: every stage of the study is frozen or preregistered
  before any blind data exist; the outcome is reported regardless of
  direction.
- Evidence chain: six sequential stages, each annotated with its concrete
  artifact; a bottom band records the binding no-backfitting rule and the
  preregistered decision matrix.
- Archetype: schematic-led composite (single-row pipeline + rule band).
- Export: double column (174 mm), SVG/PDF/TIFF 600 dpi, editable text.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from _style import DOUBLE_COL, PALETTE, apply_style, save

STAGES = [
    ("1. Dataset", "4 species, 5 studies\n389,609 sites\n2,334 positives\nfull provenance"),
    ("2. Benchmark", "6 models, 10 seeds\nwithin-dataset splits\nliterature-comparable\nclaim class only"),
    ("3. Freeze", "weights + features +\ncandidates + controls\nSHA256-registered\nco-signed 2026-08-13"),
    ("4. Calibration", "mock-blind diagnostic\n20 unseen tomato sites\n7 Arabidopsis controls\nexpectations set pre-blind"),
    ("5. Preregistration", "SAP + analysis plan\nOSF-registered DOI\nprimary endpoint fixed\nsuccess criteria fixed"),
    ("6. Blind cohort", "593 sites assayed blind\none-sided Fisher exact\nreported regardless\nof outcome"),
]

DECISIONS = [
    ("enrichment met + mechanism", "upgraded submission target"),
    ("enrichment met, no mechanism", "resource and validation report"),
    ("enrichment not met", "negative result, first-class finding"),
]


def main() -> None:
    apply_style()
    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.42))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 42)
    ax.axis("off")

    n = len(STAGES)
    box_w, box_h = 13.6, 17
    y_top = 33
    gap = (100 - n * box_w) / (n + 1)
    centers = []
    for i, (title, body) in enumerate(STAGES):
        x = gap + i * (box_w + gap)
        y = y_top - box_h
        fill = PALETTE["blue_main"] if i in (2, 4) else "white"
        edge = PALETTE["blue_main"]
        text_c = "white" if i in (2, 4) else PALETTE["neutral_black"]
        ax.add_patch(
            FancyBboxPatch(
                (x, y), box_w, box_h,
                boxstyle="round,pad=0.6",
                facecolor=fill, edgecolor=edge, linewidth=1.0,
            )
        )
        ax.text(x + box_w / 2, y + box_h - 1.6, title, ha="center", va="top",
                fontsize=7.5, fontweight="bold", color=text_c)
        ax.text(x + box_w / 2, y + box_h - 4.6, body, ha="center", va="top",
                fontsize=6.2, color=text_c, linespacing=1.35)
        centers.append((x, x + box_w, y + box_h / 2))
        if i:
            ax.add_patch(
                FancyArrowPatch(
                    (centers[i - 1][1] + 0.4, centers[i - 1][2]),
                    (x - 0.4, y + box_h / 2),
                    arrowstyle="-|>", mutation_scale=8,
                    color=PALETTE["neutral_dark"], linewidth=1.0,
                )
            )

    # Binding rule band
    ax.add_patch(
        FancyBboxPatch(
            (gap, 7.6), 100 - 2 * gap, 5.4,
            boxstyle="round,pad=0.4",
            facecolor="#F2F2F2", edgecolor=PALETTE["neutral_mid"], linewidth=0.8,
        )
    )
    ax.text(
        50, 10.3,
        "Binding rule — no model, feature, candidate-table, threshold, or "
        "control-set change after unblinding;\nevery artifact hash-verified; "
        "negative results and exclusions preserved and archived.",
        ha="center", va="center", fontsize=6.3, color=PALETTE["neutral_black"],
        linespacing=1.5,
    )

    # Decision matrix
    ax.text(gap, 5.6, "Preregistered decision matrix:", fontsize=6.6,
            fontweight="bold", color=PALETTE["neutral_black"], va="center")
    d_w = (100 - 2 * gap - 2 * 2.0) / 3
    for i, (cond, outcome) in enumerate(DECISIONS):
        x = gap + i * (d_w + 2.0)
        ax.add_patch(
            FancyBboxPatch(
                (x, 0.2), d_w, 4.6,
                boxstyle="round,pad=0.3",
                facecolor="white", edgecolor=PALETTE["teal"], linewidth=0.8,
            )
        )
        ax.text(x + d_w / 2, 4.3, cond, ha="center", va="top", fontsize=5.9,
                color=PALETTE["neutral_dark"], linespacing=1.3)
        ax.text(x + d_w / 2, 1.1, outcome, ha="center", va="center",
                fontsize=5.9, style="italic", color=PALETTE["teal"],
                linespacing=1.3)

    for p in save(fig, "fig3_design_pipeline"):
        print(p)


if __name__ == "__main__":
    main()
