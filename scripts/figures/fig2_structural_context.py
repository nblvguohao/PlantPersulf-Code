"""Figure 2 — A model-confidence confound in the structural context of
persulfidation targeting.

Contract:
- Core conclusion: an apparent burial/exposure signal, inconsistent in
  direction across species, tracked AlphaFold per-residue confidence
  (pLDDT) rather than a shared structural mechanism; restricting every
  comparison to confident residues removes the signal in all four species.
  This is a self-correction of the project's own earlier interpretation,
  reported because it is true, not because it is flattering.
- Evidence chain: (a) naive mean difference (persulfidated vs other
  cysteines of the same protein), three metrics x four species; (b) the
  pLDDT difference driving it; (c) the same three metrics x four species
  after restricting to pLDDT>=70 residues only.
- Archetype: small-multiple grouped bars, before/after.
- Data: results/cross_species_conservation/structural_context_v2.json.
- Export: double column, editable text, 600 dpi.
"""

import json

import matplotlib.pyplot as plt

from _style import DOUBLE_COL, PALETTE, REPO, apply_style, panel_label, save

DATA = REPO / "results/cross_species_conservation/structural_context_v2.json"
SPECIES = ["Arabidopsis", "Rice", "Tomato", "Magnaporthe"]
SPECIES_SHORT = ["At", "Os", "Sl", "Mo"]
SPECIES_COLOR = {
    "Arabidopsis": PALETTE["blue_main"],
    "Rice": PALETTE["green_3"],
    "Tomato": PALETTE["red_strong"],
    "Magnaporthe": PALETTE["violet"],
}
METRICS = [
    ("contact_proxy", "contact-number proxy", "higher = more buried"),
    ("residue_sasa", "residue SASA (Å$^2$)", "lower = more buried"),
    ("sg_sasa", "S$\\gamma$ SASA (Å$^2$)", "lower = more buried"),
]


def _index(rows: list[dict], key: str) -> dict:
    return {r["species"]: r for r in rows}


def _index_metric(rows: list[dict], metric: str) -> dict:
    return {r["species"]: r for r in rows if r["metric"] == metric}


def _bar_panel(ax, values: dict, ps: dict, ylabel: str, note: str, sig_thresh=0.05,
                row_note: str | None = None) -> None:
    x = range(len(SPECIES))
    heights = [values[sp] for sp in SPECIES]
    colors = [SPECIES_COLOR[sp] for sp in SPECIES]
    ymin, ymax = min(0.0, *heights), max(0.0, *heights)
    pad = (ymax - ymin + 1e-9) * 0.22
    ax.bar(x, heights, color=colors, width=0.62)
    ax.axhline(0, color=PALETTE["neutral_dark"], linewidth=0.7)
    ax.set_ylim(ymin - pad, ymax + pad)
    for i, sp in enumerate(SPECIES):
        p = ps[sp]
        marker = "*" if p < sig_thresh else "ns"
        y = heights[i]
        va = "bottom" if y >= 0 else "top"
        yy = y + pad * 0.28 if y >= 0 else y - pad * 0.28
        ax.text(i, yy, marker, ha="center", va=va,
                fontsize=6.6 if marker == "*" else 5.2,
                color=PALETTE["neutral_dark"], clip_on=False)
    ax.set_xticks(list(x))
    ax.set_xticklabels(SPECIES_SHORT, fontsize=6.2)
    ax.set_ylabel(ylabel, fontsize=6.0)
    title = note if row_note is None else f"{row_note}\n{note}"
    ax.set_title(title, fontsize=5.4, color=PALETTE["neutral_mid"], loc="left", pad=3,
                 linespacing=1.3)
    ax.tick_params(axis="y", labelsize=5.8)


def main() -> None:
    apply_style()
    doc = json.loads(DATA.read_text(encoding="utf-8"))
    proxy = _index(doc["contact_proxy_pathway"]["per_species"], "species")
    sasa = doc["sasa_pathway"]["per_species_per_metric"]
    sasa_residue = _index_metric(sasa, "residue_sasa")
    sasa_sg = _index_metric(sasa, "sg_sasa")
    plddt = _index(
        doc["model_confidence_control"]["plddt_at_persulfidated_vs_other_cysteines"],
        "species",
    )
    conf_proxy = _index(
        doc["model_confidence_control"]["contact_proxy_confident_only"], "species")
    conf_sasa = doc["model_confidence_control"]["sasa_confident_only"]
    conf_residue = _index_metric(conf_sasa, "residue_sasa")
    conf_sg = _index_metric(conf_sasa, "sg_sasa")
    min_plddt = doc["model_confidence_control"]["min_plddt"]

    assert abs(proxy["Tomato"]["mean_diff"] - (-4.049)) < 0.01
    assert abs(plddt["Tomato"]["mean_diff"] - (-7.497)) < 0.01
    assert abs(conf_residue["Tomato"]["mean_diff_A2"] - 2.167) < 0.01

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.80))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 0.66, 1.0],
                          hspace=0.62, wspace=0.55,
                          top=0.94, bottom=0.10, left=0.07, right=0.98)

    # --- (a) naive accessibility signal, three metrics -----------------
    naive_sources = [
        (proxy, "mean_diff", "p_value_two_sided"),
        (sasa_residue, "mean_diff_A2", "p_value_two_sided"),
        (sasa_sg, "mean_diff_A2", "p_value_two_sided"),
    ]
    for col, ((_key, label, note), (rows, vkey, pkey)) in enumerate(
        zip(METRICS, naive_sources, strict=True)
    ):
        ax = fig.add_subplot(gs[0, col])
        values = {sp: rows[sp][vkey] for sp in SPECIES}
        ps = {sp: rows[sp][pkey] for sp in SPECIES}
        row_note = "naive, same protein" if col == 0 else None
        _bar_panel(ax, values, ps, label, note, row_note=row_note)
        if col == 0:
            panel_label(ax, "a")

    # --- (b) pLDDT difference driving the naive signal ------------------
    ax_b = fig.add_subplot(gs[1, :])
    values = {sp: plddt[sp]["mean_diff"] for sp in SPECIES}
    ps = {sp: plddt[sp]["p_value_two_sided"] for sp in SPECIES}
    _bar_panel(ax_b, values, ps, "$\\Delta$ pLDDT (persulfidated $-$ other)",
               "AlphaFold confidence difference: its sign tracks panel a's sign"
               " in every species")
    panel_label(ax_b, "b")

    # --- (c) confidence-restricted retest, three metrics -----------------
    conf_sources = [
        (conf_proxy, "mean_diff", "p_value_two_sided"),
        (conf_residue, "mean_diff_A2", "p_value_two_sided"),
        (conf_sg, "mean_diff_A2", "p_value_two_sided"),
    ]
    for col, ((_key, label, note), (rows, vkey, pkey)) in enumerate(
        zip(METRICS, conf_sources, strict=True)
    ):
        ax = fig.add_subplot(gs[2, col])
        values = {sp: rows[sp][vkey] for sp in SPECIES}
        ps = {sp: rows[sp][pkey] for sp in SPECIES}
        row_note = f"pLDDT $\\geq$ {min_plddt:g} only" if col == 0 else None
        _bar_panel(ax, values, ps, label, note, row_note=row_note)
        if col == 0:
            panel_label(ax, "c")

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=SPECIES_COLOR[sp]) for sp in SPECIES
    ]
    fig.legend(handles, ["Arabidopsis (At)", "Rice (Os)", "Tomato (Sl)",
                         "Magnaporthe (Mo)"],
               loc="lower center", ncol=4, fontsize=6.0, frameon=False,
               bbox_to_anchor=(0.5, 0.0))

    for p in save(fig, "fig2_structural_context"):
        print(p)


if __name__ == "__main__":
    main()
