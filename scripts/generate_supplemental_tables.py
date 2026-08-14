"""Generate Supplemental Table S1 (release training panel composition).

S1 documents, per species, the composition of the frozen release training
panel (multispecies-v2-candidate-release-v1): rows, positives, positive
rate, and source datasets. The panel is the union of the ten
literature-random-track seed panels (development split only); counts are
pinned to the frozen release manifest (SHA256 30e432fb...) and asserted
in tests/scientific/test_supplemental_tables.py.

Outputs (manuscripts/plant_physiology/supplements/):
- supplemental_table_s1.tsv — machine-readable copy
- supplemental_table_s1.tex — standalone LaTeX table (booktabs)
"""

import csv
import glob
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SEED_GLOB = (
    REPO
    / "results/experiments/multispecies_v2_global_clusters_v11"
    / "literature_random_protein/sul_bertgru_seed*/shared_panel.tsv"
)
OUTDIR = REPO / "manuscripts" / "plant_physiology" / "supplements"

SPECIES_ORDER = ["arabidopsis", "rice", "tomato", "magnaporthe"]

SPECIES_INFO = {
    "arabidopsis": {
        "display": "Arabidopsis (\\textit{A. thaliana})",
        "sources": "PXD006140, PXD024061",
        "published_sites": 390,
    },
    "rice": {
        "display": "Rice (\\textit{O. sativa})",
        "sources": "PXD072089",
        "published_sites": 929,
    },
    "tomato": {
        "display": "Tomato (\\textit{S. lycopersicum})",
        "sources": "kiae271",
        "published_sites": 99,
    },
    "magnaporthe": {
        "display": "\\textit{Magnaporthe oryzae}",
        "sources": "PXD063170",
        "published_sites": 1482,
    },
}


def composition_from_seed_panels(seed_glob: str) -> dict:
    """Union the seed panels by (protein, position) and tabulate per species."""
    rows: dict[tuple[str, str], str] = {}
    for path in sorted(glob.glob(seed_glob)):
        with open(path, encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                key = (r["global_protein_id"], r["cys_position"])
                rows.setdefault(key, r["adapter_label"])
    species: dict[str, dict[str, int]] = defaultdict(
        lambda: {"rows": 0, "positives": 0}
    )
    for (protein_id, _), label in rows.items():
        sp = protein_id.split("|", 1)[0]
        species[sp]["rows"] += 1
        if label == "positive":
            species[sp]["positives"] += 1
    return {
        "n_rows": len(rows),
        "n_positives": sum(s["positives"] for s in species.values()),
        "species": {sp: dict(species.get(sp, {"rows": 0, "positives": 0})) for sp in SPECIES_ORDER},
    }


def _latex_table(comp: dict) -> str:
    lines = [
        "\\documentclass[11pt]{article}",
        "\\usepackage[margin=2.2cm]{geometry}",
        "\\usepackage{booktabs}",
        "\\begin{document}",
        "\\begin{table}[h]",
        "\\caption{Supplemental Table S1. Composition of the frozen release "
        "training panel (release \\texttt{multispecies-v2-candidate-release-v1}, "
        "SHA256 \\texttt{30e432fb...}) by species. The panel is the union of "
        "the ten literature-random-track seed panels (development split only). "
        "Panel counts reflect coordinate-verified rows; published-site counts "
        "are given for reference.}",
        "\\label{tab:supplemental_s1}",
        "\\begin{tabular}{@{}lrrrrl@{}}",
        "\\toprule",
        "Species & Panel rows & Positives & Rate (\\%) & "
        "Published sites & Source datasets \\\\",
        "\\midrule",
    ]
    for sp in SPECIES_ORDER:
        s = comp["species"][sp]
        info = SPECIES_INFO[sp]
        rate = 100.0 * s["positives"] / s["rows"]
        lines.append(
            f"{info['display']} & {s['rows']:,} & {s['positives']:,} & "
            f"{rate:.2f} & {info['published_sites']:,} & "
            f"\\texttt{{{info['sources']}}} \\\\"
        )
    lines += [
        "\\midrule",
        f"Total & {comp['n_rows']:,} & {comp['n_positives']:,} & "
        f"{100.0 * comp['n_positives'] / comp['n_rows']:.2f} & --- & --- \\\\",
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}",
        "\\end{document}",
        "",
    ]
    return "\n".join(lines)


def _tsv(comp: dict) -> str:
    header = ["species", "panel_rows", "positives", "rate_percent",
              "published_sites", "source_datasets"]
    lines = ["\t".join(header)]
    for sp in SPECIES_ORDER:
        s = comp["species"][sp]
        info = SPECIES_INFO[sp]
        lines.append(
            "\t".join([
                sp, str(s["rows"]), str(s["positives"]),
                f"{100.0 * s['positives'] / s['rows']:.4f}",
                str(info["published_sites"]),
                info["sources"],
            ])
        )
    return "\n".join(lines) + "\n"


def main() -> None:
    comp = composition_from_seed_panels(str(SEED_GLOB))
    OUTDIR.mkdir(parents=True, exist_ok=True)
    tsv_path = OUTDIR / "supplemental_table_s1.tsv"
    tsv_path.write_text(_tsv(comp), encoding="utf-8")
    tex_path = OUTDIR / "supplemental_table_s1.tex"
    tex_path.write_text(_latex_table(comp), encoding="utf-8")
    print(f"wrote {tsv_path}")
    print(f"wrote {tex_path}")
    print(f"union rows={comp['n_rows']:,} positives={comp['n_positives']:,}")


if __name__ == "__main__":
    main()
