"""Shared publication style for the Plant Physiology manuscript figures.

Journal/export contract (Plant Physiology): single column 85 mm, double
column 174 mm; Arial/Helvetica 6-8 pt; line art >= 600 dpi; editable text
in vector exports. All figure scripts import from here so the contract is
applied once.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt

MM = 1 / 25.4
SINGLE_COL = 85 * MM
DOUBLE_COL = 174 * MM

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "green_3": "#8BCF8B",
    "red_strong": "#B64342",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "neutral_black": "#272727",
}

REPO = Path(__file__).resolve().parents[2]
FIGDIR = REPO / "manuscripts" / "plant_physiology" / "figures"


def apply_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
            "xtick.labelsize": 6.5,
            "ytick.labelsize": 6.5,
            "axes.labelsize": 7,
        }
    )


def save(fig, name: str, dpi: int = 600) -> list[str]:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    base = FIGDIR / name
    saved = []
    for ext in ("svg", "pdf", "tiff"):
        path = f"{base}.{ext}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        saved.append(path)
    qa = f"{base}_qa.png"
    fig.savefig(qa, dpi=200, bbox_inches="tight")
    saved.append(qa)
    plt.close(fig)
    return saved


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        ha="left",
        va="bottom",
    )
