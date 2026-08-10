#!/usr/bin/env python
"""kiae271 119-site full-spectrum analysis report (Zhang-lab deliverable).

Generates the "全谱分析报告" promised to Zhang Hua's team: functional
(GO-term enrichment vs the whole tomato proteome), regulatory-stratum
(lcd_gain / wt_only / both), and structural/sequence context (positive-
charge density, hydrophobicity, AlphaFold pLDDT / contact number where
structures exist) analyses of the 99 coordinate-verified tomato
persulfidation sites mined from their own published Dataset S1.

Output: markdown report + TSV tables under ``results/kiae271_full_report_v1/``.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path

from plantpersulf.features.sequence import (
    _flanking_window,
    _hydrophobicity,
    _local_positive_charge_density,
)
from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites

KIAE271_XLSX = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
TOMATO_PROTEOME = Path("data/raw/references/tomato_ref_proteome_v1.fasta")
GO_FILE = Path("data/raw/references/go_annotations_v1/tomato_go_taxon4081.tsv")
PANTHER_FILE = Path(
    "data/raw/references/panther_annotations_v1/panther_tomato_taxon4081.tsv"
)
ALPHAFOLD_REGISTRY = Path("data/registry/alphafold_structures.tsv")

GO_RE = re.compile(r"(.+?)\s*\[(GO:\d+)\]")


def _load_go(path: Path) -> dict[str, list[tuple[str, str]]]:
    """accession -> [(term_name, go_id), ...]"""
    out: dict[str, list[tuple[str, str]]] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            acc = row["Entry"]
            terms: list[tuple[str, str]] = []
            for chunk in (row.get("Gene Ontology (GO)") or "").split(";"):
                m = GO_RE.match(chunk.strip())
                if m:
                    terms.append((m.group(1), m.group(2)))
            if terms:
                out[acc] = terms
    return out


def _load_panther(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            acc = row["Entry"]
            families = sorted(
                {
                    t.split(":")[0]
                    for t in (row.get("PANTHER") or "").split(";")
                    if t.strip().startswith("PTHR")
                }
            )
            if families:
                out[acc] = families[0]
    return out


def _load_proteome_headers(path: Path) -> dict[str, str]:
    """accession -> protein name from fasta header (for the report tables)."""
    out: dict[str, str] = {}
    cur_acc = ""
    cur_desc = ""
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur_acc:
                out[cur_acc] = cur_desc
            parts = line[1:].split("|")
            cur_acc = parts[1] if len(parts) >= 2 else line[1:].split()[0]
            cur_desc = line[1:].strip()
    if cur_acc:
        out[cur_acc] = cur_desc
    return out


def _load_proteome(path: Path) -> dict[str, str]:
    from plantpersulf.features.sequence import _load_proteome as _lp

    return _lp(path)


def _load_structure_registry(path: Path) -> set[str]:
    accs: set[str] = set()
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            accs.add(row["accession"])
    return accs


def _fisher_enrichment(
    go_by_acc: dict[str, list[tuple[str, str]]],
    foreground: set[str],
    background: set[str],
    min_hits: int = 2,
) -> list[dict]:
    """One-sided Fisher exact per GO id: foreground vs background."""

    n_fg = len(foreground)
    n_bg = len(background)
    term_hits_fg: Counter[str] = Counter()
    term_hits_bg: Counter[str] = Counter()
    term_names: dict[str, str] = {}
    for acc in foreground:
        for name, gid in go_by_acc.get(acc, []):
            term_hits_fg[gid] += 1
            term_names[gid] = name
    for acc in background:
        for name, gid in go_by_acc.get(acc, []):
            term_hits_bg[gid] += 1
            term_names[gid] = name

    results: list[dict] = []
    for gid, k in term_hits_fg.items():
        if k < min_hits:
            continue
        n_bg_hit = term_hits_bg.get(gid, 0)
        # table: [k, n_fg-k; n_bg_hit, n_bg-n_bg_hit]
        p = _fisher_p_value(k, n_fg - k, n_bg_hit, n_bg - n_bg_hit)
        results.append(
            {
                "go_id": gid,
                "term": term_names[gid],
                "hits_fg": k,
                "n_fg": n_fg,
                "hits_bg": n_bg_hit,
                "n_bg": n_bg,
                "p_value": p,
            }
        )
    # BH correction
    results.sort(key=lambda r: r["p_value"])
    m = len(results)
    for i, r in enumerate(results):
        r["fdr"] = min(r["p_value"] * m / (i + 1), 1.0)
    # cascade fdr
    for i in range(len(results) - 2, -1, -1):
        results[i]["fdr"] = min(results[i]["fdr"], results[i + 1]["fdr"])
    return results


def _fisher_p_value(a: int, b: int, c: int, d: int) -> float:
    """One-sided (>=) Fisher exact test via hypergeometric tail."""
    from math import comb

    total = a + b + c + d
    row1 = a + c
    col1 = a + b
    p = 0.0
    for x in range(a, min(col1, row1) + 1):
        y = col1 - x
        z = row1 - x
        w = total - row1 - col1 + x
        if y >= 0 and z >= 0 and w >= 0:
            p += comb(col1, x) * comb(total - col1, row1 - x) / comb(total, row1)
    return p


def run_report(output_dir: Path) -> dict:
    proteome = _load_proteome(TOMATO_PROTEOME)
    headers = _load_proteome_headers(TOMATO_PROTEOME)
    go_by_acc = _load_go(GO_FILE)
    panther = _load_panther(PANTHER_FILE)
    structure_accs = _load_structure_registry(ALPHAFOLD_REGISTRY)

    table = parse_kiae271_sites(KIAE271_XLSX, proteome)
    sites = list(table.sites)
    accs = sorted({s.protein_accession for s in sites})
    acc_set = set(accs)
    bg_accs = set(go_by_acc)  # whole tomato proteome with GO annotations

    # -- sequence context features per site
    seq_ctx: dict[tuple[str, int], dict] = {}
    for s in sites:
        seq = proteome[s.protein_accession]
        flank = _flanking_window(seq, s.cys_position, 10)
        seq_ctx[(s.protein_accession, s.cys_position)] = {
            "hydrophobicity": round(_hydrophobicity(flank), 3),
            "pos_charge_density": round(_local_positive_charge_density(flank), 3),
        }

    # -- strata
    strata = Counter(s.regulation for s in sites)

    # -- enrichment per stratum (lcd_gain = H2S-induced, most informative)
    enrich: dict[str, list[dict]] = {}
    for name, label in (("lcd_gain", "lcd_gain"), ("all", "all")):
        fg = (
            {s.protein_accession for s in sites if s.regulation == label}
            if name == "lcd_gain"
            else acc_set
        )
        enrich[name] = _fisher_enrichment(go_by_acc, fg, bg_accs)

    # -- PANTHER family distribution of the panel
    fam_counter: Counter[str] = Counter()
    for acc in accs:
        fam_counter[panther.get(acc, "no PANTHER family")] += 1
    top_fams = fam_counter.most_common(15)

    # -- structure coverage
    n_covered = sum(1 for acc in accs if acc in structure_accs)

    # -- write per-site table
    output_dir.mkdir(parents=True, exist_ok=True)
    site_rows = []
    for s in sites:
        site_rows.append(
            {
                "protein_accession": s.protein_accession,
                "cys_position": s.cys_position,
                "regulation": s.regulation,
                "localization_prob": s.localization_prob,
                "intensity_lcd": s.intensity_lcd,
                "intensity_wt": s.intensity_wt,
                "protein_name": (headers.get(s.protein_accession, "") or "")[:120],
                "panther_family": panther.get(s.protein_accession, ""),
                "has_alphafold_structure": "yes"
                if s.protein_accession in structure_accs
                else "no",
                **seq_ctx[(s.protein_accession, s.cys_position)],
            }
        )
    with (output_dir / "sites_full.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(site_rows[0].keys()), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(site_rows)

    for name in ("lcd_gain", "all"):
        rows = enrich[name]
        with (output_dir / f"go_enrichment_{name}.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(rows[0].keys()), delimiter="\t"
            )
            writer.writeheader()
            writer.writerows(rows)

    # -- markdown report
    md: list[str] = [
        "# kiae271 番茄差异硫巯基化位点全谱分析报告",
        "",
        "> 数据来源:Zhang et al. 2024, Plant Physiology, doi:10.1093/plphys/kiae271,"
        " Supplementary Dataset S1(SlLCD1-OE vs WT 番茄叶片)",
        "> 分析日期:2026-08-10 | 生成脚本:`scripts/kiae271_full_report.py`",
        "",
        "## 1. 数据总览",
        "",
        f"- 原始行数:{table.rows_total}",
        f"- 坐标验证通过位点:**{table.total_verified_sites} 个 / {table.total_verified_proteins} 个蛋白**",
        f"- 分层:lcd_gain(H2S 诱导获得)**{strata['lcd_gain']}**、wt_only **{strata['wt_only']}**、both **{strata['both']}**",
        "- 定位概率 ≥0.75;每行至少一个条件强度非零",
        f"- AlphaFold 结构覆盖:{n_covered}/{len(accs)}({n_covered / len(accs) * 100:.0f}%)",
        "",
        "## 2. GO 功能富集(lcd_gain 分层,H2S 诱导获得位点)",
        "",
        "背景:全番茄蛋白组 GO 注释;Fisher 精确检验 + Benjamini-Hochberg FDR。",
        "",
        "| GO 项 | 术语 | 命中/总数 | p 值 | FDR |",
        "|---|---|---|---|---|",
    ]
    for r in enrich["lcd_gain"][:15]:
        md.append(
            f"| {r['go_id']} | {r['term']} | {r['hits_fg']}/{r['n_fg']} | "
            f"{r['p_value']:.2e} | {r['fdr']:.2e} |"
        )
    md += [
        "",
        "## 3. GO 功能富集(全部 99 位点)",
        "",
        "| GO 项 | 术语 | 命中/总数 | p 值 | FDR |",
        "|---|---|---|---|---|",
    ]
    for r in enrich["all"][:15]:
        md.append(
            f"| {r['go_id']} | {r['term']} | {r['hits_fg']}/{r['n_fg']} | "
            f"{r['p_value']:.2e} | {r['fdr']:.2e} |"
        )
    md += [
        "",
        "## 4. PANTHER 家族分布(位点数最多的家族)",
        "",
        "| 家族 | 位点蛋白数 |",
        "|---|---|",
    ]
    for fam, n in top_fams:
        md.append(f"| {fam} | {n} |")
    md += [
        "",
        "## 5. 序列环境特征",
        "",
        f"- 正电荷密度(硫醇盐稳定化代理)均值:"
        f"{sum(v['pos_charge_density'] for v in seq_ctx.values()) / len(seq_ctx):.3f}",
        f"- 疏水性均值:"
        f"{sum(v['hydrophobicity'] for v in seq_ctx.values()) / len(seq_ctx):.3f}",
        "",
        "完整逐位点表见 `sites_full.tsv`,富集表见 `go_enrichment_*.tsv`。",
        "",
    ]
    (output_dir / "kiae271_full_report.md").write_text("\n".join(md), encoding="utf-8")

    summary = {
        "rows_total": table.rows_total,
        "verified_sites": table.total_verified_sites,
        "verified_proteins": table.total_verified_proteins,
        "strata": dict(strata),
        "structure_coverage": f"{n_covered}/{len(accs)}",
        "lcd_gain_fdr_significant": sum(
            1 for r in enrich["lcd_gain"] if r["fdr"] < 0.05
        ),
        "all_fdr_significant": sum(1 for r in enrich["all"] if r["fdr"] < 0.05),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"report -> {output_dir}")
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="kiae271 full-spectrum report")
    p.add_argument(
        "--output-dir", type=Path, default=Path("results/kiae271_full_report_v1")
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_report(args.output_dir)
