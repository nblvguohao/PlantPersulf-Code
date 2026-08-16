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

Design source: user-provided AI mockup v1 (figure_prompts/
fig3_design_ai_mockup_v1.png, 2026-08-14), replicated faithfully: white
canvas; one row of six stage boxes with thin navy (#045599) outline
(stages 3 and 5 solid navy with white text); navy arrows between stages;
a full-width light-grey band carrying the binding rule in dark grey; three
navy-outlined decision cells at the bottom. Intentional deviations from
the mockup: (1) uniform inter-box gaps with an arrow added between stages
1 and 2 (missing in the mockup); (2) all text re-set with the exact
manuscript strings at Plant Physiology-compliant 6-8 pt.
"""

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from _style import DOUBLE_COL, apply_style, save

# Mockup palette (fig3_design_ai_mockup_v1.png).
NAVY = "#045599"
BAND_GREY = "#EDEDED"
BAND_TEXT = "#2C2E34"

# Mockup canvas geometry, in mockup pixels (data units below are px).
W, H = 2172, 724

STAGES = [
    ("1. Dataset",
     "4 species, 5 studies\n389,609 sites\n2,334 positives\nfull provenance"),
    ("2. Benchmark",
     "6 models, 10 seeds\nwithin-dataset splits\nliterature-comparable\n"
     "claim class only"),
    ("3. Freeze",
     "weights + features +\ncandidates + controls\nSHA256-registered\n"
     "co-signed 2026-08-13"),
    ("4. Calibration",
     "mock-blind diagnostic\n20 unseen tomato sites\n7 Arabidopsis controls\n"
     "expectations set pre-blind"),
    ("5. Preregistration",
     "SAP + analysis plan\nOSF-registered DOI\nprimary endpoint fixed\n"
     "success criteria fixed"),
    ("6. Blind cohort",
     "593 sites assayed blind\none-sided Fisher exact\nreported regardless\n"
     "of outcome"),
]

BINDING_LINES = [
    "Binding rule — no model, feature, candidate-table, threshold, or "
    "control-set change after unblinding;",
    "every artifact hash-verified; negative results and exclusions "
    "preserved and archived.",
]

DECISIONS = [
    ("enrichment met + mechanism", "upgraded submission target"),
    ("enrichment met, no mechanism", "resource and validation report"),
    ("enrichment not met", "negative result, first-class finding"),
]

# Mockup box row: y 20..378, height 358. Box widths sized proportionally
# to each box's measured widest text line (verified programmatically; the
# mockup's own widths left no room for box 4's longest line).
BOX_TOP, BOX_H = 20, 358
BOX_WIDTHS = [276, 297, 305, 359, 333, 325]
LEFT = 21
GAP = (W - LEFT - 21 - sum(BOX_WIDTHS)) / (len(BOX_WIDTHS) - 1)

BAND_X, BAND_Y, BAND_H = 24, 426, 92
CELL_XS = [21, 755, 1472]
CELL_W = 675
CELL_Y, CELL_H = 559, 140


def _y(mock_y: float) -> float:
    """Convert mockup pixel y (top-down) to matplotlib axis y."""
    return H - mock_y


def main() -> None:
    apply_style()
    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * H / W))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, H)
    ax.axis("off")

    # --- stage boxes ---------------------------------------------------
    x = float(LEFT)
    rights: list[float] = []
    for i, ((title, body), w) in enumerate(zip(STAGES, BOX_WIDTHS, strict=True)):
        solid = i in (2, 4)
        ax.add_patch(
            FancyBboxPatch(
                (x, _y(BOX_TOP + BOX_H)), w, BOX_H,
                boxstyle="round,pad=0,rounding_size=10",
                facecolor=NAVY if solid else "white",
                edgecolor=NAVY,
                linewidth=0.8,
            )
        )
        text_c = "white" if solid else NAVY
        ax.text(
            x + w / 2, _y(BOX_TOP + 0.10 * BOX_H), title,
            ha="center", va="top", fontsize=8, fontweight="bold", color=text_c,
        )
        body_lines = body.split("\n")
        for j, line in enumerate(body_lines):
            ax.text(
                x + w / 2, _y(BOX_TOP + (0.32 + 0.14 * j) * BOX_H), line,
                ha="center", va="top", fontsize=6.2, color=text_c,
                linespacing=1.3,
            )
        rights.append(x + w)
        x += w + GAP

    # --- arrows --------------------------------------------------------
    y_mid = _y(BOX_TOP + BOX_H / 2)
    for i in range(len(rights) - 1):
        left_next = rights[i + 1] - BOX_WIDTHS[i + 1]
        ax.add_patch(
            FancyArrowPatch(
                (rights[i] + 0.12 * GAP, y_mid),
                (left_next - 0.12 * GAP, y_mid),
                arrowstyle="-|>", mutation_scale=10,
                color=NAVY, linewidth=0.8,
            )
        )

    # --- binding-rule band ----------------------------------------------
    ax.add_patch(
        FancyBboxPatch(
            (BAND_X, _y(BAND_Y + BAND_H)), W - 2 * BAND_X, BAND_H,
            boxstyle="round,pad=0,rounding_size=10",
            facecolor=BAND_GREY, edgecolor="none",
        )
    )
    for j, line in enumerate(BINDING_LINES):
        ax.text(
            W / 2, _y(BAND_Y + BAND_H / 2) + (0.5 - j) * 18, line,
            ha="center", va="center", fontsize=6.3, color=BAND_TEXT,
        )

    # --- decision matrix -------------------------------------------------
    for cx, (cond, outcome) in zip(CELL_XS, DECISIONS, strict=True):
        ax.add_patch(
            FancyBboxPatch(
                (cx, _y(CELL_Y + CELL_H)), CELL_W, CELL_H,
                boxstyle="round,pad=0,rounding_size=10",
                facecolor="white", edgecolor=NAVY, linewidth=0.8,
            )
        )
        ax.text(
            cx + CELL_W / 2, _y(CELL_Y + 0.22 * CELL_H), cond,
            ha="center", va="top", fontsize=6.4, color=NAVY,
        )
        ax.text(
            cx + CELL_W / 2, _y(CELL_Y + 0.55 * CELL_H), outcome,
            ha="center", va="top", fontsize=6.3, style="italic",
            color=NAVY,
        )

    for p in save(fig, "fig3_design_pipeline"):
        print(p)


if __name__ == "__main__":
    main()
