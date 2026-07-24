#!/usr/bin/env python
"""Publication figures for the cross-species conservation reframe (§5.5,
docs/phase_z_evidence_audit.md) — targeting Redox Biology / New Phytologist.

Core conclusion: persulfidation targeting converges across three
independent species (Arabidopsis, rice, Magnaporthe) at the ortholog-family
level and at the structural-burial level, even though naive sequence-based
cross-species prediction is weak (§5.2-5.4).

Three panels, each carrying a distinct piece of evidence:
    A. Family-level co-persulfidation enrichment (odds ratio, three
       species pairs), annotated with raw / Bonferroni / BH-corrected
       significance.
    B. Structural context, contact-number accessibility proxy pathway
       (positive vs unlabeled cysteines, three species).
    C. Structural context, real Shrake-Rupley SASA pathway (whole-residue
       and SG-specific), the independent cross-check on panel B.

Usage::

    python scripts/make_conservation_figures.py --output-dir figures
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# ---------------------------------------------------------------------------
# Mandatory publication style (editable SVG text, Nature-style axes)
# ---------------------------------------------------------------------------
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 9
plt.rcParams["axes.spines.right"] = False
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.linewidth"] = 0.9
plt.rcParams["legend.frameon"] = False

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "neutral_light": "#CFCECE",
    "neutral_mid": "#767676",
    "neutral_dark": "#4D4D4D",
    "red_strong": "#B64342",
    "green_ok": "#3F8F5A",
}

SPECIES_COLORS = {
    "Arabidopsis": PALETTE["blue_main"],
    "Rice": PALETTE["teal"],
    "Magnaporthe": PALETTE["violet"],
}

CONSERVATION_JSON = Path("results/cross_species_conservation/conservation_v1.json")
STRUCTURAL_JSON = Path(
    "results/cross_species_conservation/structural_context_v1.json"
)


def _sig_marker(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


def _finalize(
    fig: plt.Figure, out_path: Path, dpi: int = 600, use_tight_layout: bool = True
) -> list[str]:
    if use_tight_layout:
        fig.tight_layout(pad=1.5)
    os.makedirs(out_path.parent, exist_ok=True)
    saved = []
    for ext in ("svg", "pdf", "png"):
        p = out_path.with_suffix(f".{ext}")
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        saved.append(str(p))
    plt.close(fig)
    return saved


# ---------------------------------------------------------------------------
# Panel A: family-level enrichment (odds ratio + 3-way significance)
# ---------------------------------------------------------------------------


def panel_a_family_enrichment(ax: plt.Axes, conservation: dict) -> None:
    pairs = conservation["pairwise_tests"]
    labels = [f"{p['species_a']} ×\n{p['species_b']}" for p in pairs]
    odds_ratios = [p["odds_ratio"] for p in pairs]
    fold = [p["fold_enrichment"] for p in pairs]
    raw_p = [p["p_value"] for p in pairs]
    bonf_p = [p["bonferroni_p_value"] for p in pairs]
    bh_q = [p["benjamini_hochberg_q_value"] for p in pairs]
    bonf_sig = [p["significant_bonferroni_0.05"] for p in pairs]

    y = np.arange(len(pairs))[::-1]
    colors = [
        PALETTE["blue_main"] if s else PALETTE["neutral_mid"] for s in bonf_sig
    ]
    ax.barh(y, odds_ratios, color=colors, edgecolor="black", linewidth=0.8, height=0.5)
    ax.axvline(1.0, color=PALETTE["neutral_dark"], linestyle="--", linewidth=1.0)

    x_max = max(odds_ratios) * 1.85
    stats_zip = zip(y, odds_ratios, fold, raw_p, bonf_p, bh_q, strict=True)
    for yi, orv, f, rp, bp, bq in stats_zip:
        ax.text(
            orv + x_max * 0.03,
            yi,
            f"OR={orv:.2f} ({f:.2f}×)\n"
            f"raw p={rp:.1e}\nBonf={bp:.1e}\nBH q={bq:.1e}",
            va="center",
            ha="left",
            fontsize=5.8,
            color=PALETTE["neutral_dark"],
            linespacing=1.4,
        )

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_ylim(-0.6, len(pairs) - 0.4)
    ax.set_xlabel("Odds ratio (co-persulfidated ortholog families)")
    ax.set_xlim(0, x_max)
    ax.set_title(
        "A  Family-level co-persulfidation enrichment\n"
        "(PANTHER ortholog subfamilies, Fisher exact)",
        fontsize=8.5,
        loc="left",
    )

    legend_handles = [
        plt.Rectangle(
            (0, 0), 1, 1, color=PALETTE["blue_main"], label="Bonferroni sig."
        ),
        plt.Rectangle(
            (0, 0), 1, 1, color=PALETTE["neutral_mid"], label="Bonferroni n.s."
        ),
    ]
    ax.legend(
        handles=legend_handles, loc="upper center", bbox_to_anchor=(0.5, -0.18),
        ncol=2, fontsize=6.5,
    )


# ---------------------------------------------------------------------------
# Panel B: structural context, contact-number proxy pathway
# ---------------------------------------------------------------------------


def panel_b_contact_proxy(ax: plt.Axes, structural: dict) -> None:
    rows = structural["contact_proxy_pathway"]["per_species"]
    species = [r["species"] for r in rows]
    pos = [r["mean_contact_positive"] for r in rows]
    unl = [r["mean_contact_unlabeled"] for r in rows]
    pvals = [r["p_value_two_sided"] for r in rows]

    x = np.arange(len(species))
    w = 0.32
    ax.bar(
        x - w / 2, pos, width=w, label="Persulfidated (positive)",
        color=[SPECIES_COLORS[s] for s in species], edgecolor="black", linewidth=0.8,
    )
    ax.bar(
        x + w / 2, unl, width=w, label="Other Cys (unlabeled)",
        color=PALETTE["neutral_light"], edgecolor="black", linewidth=0.8,
    )
    y_max = max(max(pos), max(unl)) * 1.28
    for xi, p in zip(x, pvals, strict=True):
        y_top = max(pos[xi], unl[xi]) + y_max * 0.03
        ax.text(xi, y_top, _sig_marker(p), ha="center", fontsize=8, fontweight="bold")
        ax.text(xi, y_top + y_max * 0.07, f"p={p:.3f}", ha="center", fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(species, fontsize=8)
    ax.set_ylim(0, y_max)
    ax.set_ylabel("Contact-number accessibility proxy\n(higher = more buried)")
    ax.set_title(
        "B  Structural context: accessibility PROXY\n"
        "(Cα contact count, not rigorous SASA)",
        fontsize=8.5,
        loc="left",
    )
    ax.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=6.5,
    )


# ---------------------------------------------------------------------------
# Panel C: structural context, real Shrake-Rupley SASA (cross-check)
# ---------------------------------------------------------------------------


def panel_c_real_sasa(ax: plt.Axes, structural: dict) -> None:
    rows = structural["sasa_pathway"]["per_species_per_metric"]
    residue_rows = {r["species"]: r for r in rows if r["metric"] == "residue_sasa"}
    species = list(residue_rows.keys())

    pos = [residue_rows[s]["mean_sasa_positive_A2"] for s in species]
    unl = [residue_rows[s]["mean_sasa_unlabeled_A2"] for s in species]
    pvals = [residue_rows[s]["p_value_two_sided"] for s in species]

    x = np.arange(len(species))
    w = 0.32
    ax.bar(
        x - w / 2, pos, width=w, label="Persulfidated (positive)",
        color=[SPECIES_COLORS[s] for s in species], edgecolor="black", linewidth=0.8,
    )
    ax.bar(
        x + w / 2, unl, width=w, label="Other Cys (unlabeled)",
        color=PALETTE["neutral_light"], edgecolor="black", linewidth=0.8,
    )
    y_max = max(max(pos), max(unl)) * 1.28
    for xi, p in zip(x, pvals, strict=True):
        y_top = max(pos[xi], unl[xi]) + y_max * 0.03
        ax.text(xi, y_top, _sig_marker(p), ha="center", fontsize=8, fontweight="bold")
        ax.text(xi, y_top + y_max * 0.07, f"p={p:.3f}", ha="center", fontsize=6)

    ax.set_xticks(x)
    ax.set_xticklabels(species, fontsize=8)
    ax.set_ylim(0, y_max)
    ax.set_ylabel("Whole-residue SASA (Å²)\n(lower = more buried)")
    ax.set_title(
        "C  Structural context: real Shrake-Rupley SASA\n"
        "(independent cross-check on panel B)",
        fontsize=8.5,
        loc="left",
    )
    ax.legend(
        loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=6.5,
    )


def build_figure(output_dir: Path) -> list[str]:
    with CONSERVATION_JSON.open(encoding="utf-8") as f:
        conservation = json.load(f)
    with STRUCTURAL_JSON.open(encoding="utf-8") as f:
        structural = json.load(f)

    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.6))
    panel_a_family_enrichment(axes[0], conservation)
    panel_b_contact_proxy(axes[1], structural)
    panel_c_real_sasa(axes[2], structural)

    fig.suptitle(
        "Cross-species convergence of H$_2$S-mediated persulfidation "
        "targeting at the family and structural level",
        fontsize=10,
        y=1.02,
    )
    fig.subplots_adjust(bottom=0.28, wspace=0.55, top=0.82)
    return _finalize(
        fig,
        output_dir / "fig1_cross_species_conservation",
        dpi=600,
        use_tight_layout=False,
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="build cross-species conservation figures")
    p.add_argument("--output-dir", type=Path, default=Path("figures"))
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    saved = build_figure(args.output_dir)
    print("wrote:")
    for s in saved:
        print(" ", s)
