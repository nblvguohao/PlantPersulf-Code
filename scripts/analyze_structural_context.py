#!/usr/bin/env python
"""Cross-species structural context of persulfidation targeting (2026-07-23).

Mechanistic companion to ``analyze_cross_species_conservation.py``: within
each species, do persulfidated cysteines sit in a systematically different
structural context (AlphaFold contact-number accessibility proxy) than
other cysteines in the SAME proteins? And does the direction of any such
preference agree across three independent species?

Requires the bulk AlphaFold structure download
(``scripts/download_alphafold_structures_bulk.py``) to have populated
``data/registry/alphafold_structures.tsv``.

Usage::

    python scripts/analyze_structural_context.py \
        --output results/cross_species_conservation/structural_context_v1.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

BENCHMARK = Path("data/processed/benchmark_v1/sites.tsv")
ARABIDOPSIS_PROTEOME = Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta")
RICE_SD01 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd01.xlsx")
RICE_SD04 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd04.xlsx")
RICE_SS = Path("data/raw/supplements/PXD072089/SS-all-peptides.tsv")
RICE_PROTEOME = Path("data/raw/references/rice_proteome_v1/uniprot_rice_v1.fasta")
MG_SITE_TSV = Path("data/raw/supplements/PXD063170/PXD063170_sites_moesm3.tsv")
MG_PROTEOME = Path("data/raw/supplements/PXD063170/Magnaporthe_oryzae.MG8.pep.all.fa")
PANTHER_MG = Path(
    "data/raw/references/panther_annotations_v1/panther_magnaporthe_taxon242507.tsv"
)
ALPHAFOLD_REGISTRY = Path("data/registry/alphafold_structures.tsv")


def _arabidopsis_keys() -> set[tuple[str, int]]:
    keys: set[tuple[str, int]] = set()
    with BENCHMARK.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["label"] == "positive":
                pos = int(row["cys_position_in_protein"])
                keys.add((row["protein_accession"], pos))
    return keys


def _rice_keys() -> tuple[set[tuple[str, int]], dict[str, str]]:
    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.proteomics.pxd072089_sites import parse_pxd072089_sites

    proteome = _load_proteome(RICE_PROTEOME)
    table = parse_pxd072089_sites(
        RICE_SD01, RICE_SD04, proteome, ss_all_peptides_path=RICE_SS
    )
    keys = {(s.protein_accession, s.cys_position) for s in table.sites}
    return keys, proteome


def _magnaporthe_keys() -> tuple[set[tuple[str, int]], dict[str, str]]:
    from plantpersulf.evaluation.cross_species_conservation import (
        mg8_accession_to_gene,
    )
    from plantpersulf.proteomics.pxd063170_sites import (
        load_ensembl_fungi_proteome,
        parse_pxd063170_sites,
    )

    mg_proteome_by_gene_acc = load_ensembl_fungi_proteome(MG_PROTEOME)
    table = parse_pxd063170_sites(
        MG_SITE_TSV, mg_proteome_by_gene_acc, min_localization=0.75
    )

    gene_to_uniprot: dict[str, str] = {}
    with PANTHER_MG.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            gene = (row.get("Gene Names (ORF)") or "").strip()
            if gene:
                gene_to_uniprot[gene] = row["Entry"]

    keys: set[tuple[str, int]] = set()
    uniprot_proteome: dict[str, str] = {}
    for site in table.sites:
        gene = mg8_accession_to_gene(site.protein_accession)
        uniprot_acc = gene_to_uniprot.get(gene)
        if uniprot_acc is None:
            continue
        keys.add((uniprot_acc, site.cys_position))
        uniprot_proteome[uniprot_acc] = mg_proteome_by_gene_acc[site.protein_accession]
    return keys, uniprot_proteome


def _structure_file_map() -> dict[str, Path]:
    if not ALPHAFOLD_REGISTRY.is_file():
        raise RuntimeError(
            f"AlphaFold structure registry missing: {ALPHAFOLD_REGISTRY} — "
            "run scripts/download_alphafold_structures_bulk.py first"
        )
    mapping: dict[str, Path] = {}
    with ALPHAFOLD_REGISTRY.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            local_path = (ALPHAFOLD_REGISTRY.parent / row["local_path"]).resolve()
            mapping[row["accession"]] = local_path
    return mapping


def run_structural_context_analysis(
    output_path: Path,
    n_perm: int = 1000,
    seed: int = 0,
    n_sphere_points: int = 92,
) -> dict[str, Any]:
    from plantpersulf.evaluation.cross_species_structural_context import (
        build_structural_context_rows,
        build_structural_context_rows_sasa,
        sasa_structural_context_test,
        structural_context_test,
    )
    from plantpersulf.features.sequence import _load_proteome

    structure_files = _structure_file_map()
    print(f"AlphaFold structures registered: {len(structure_files)}")

    at_keys = _arabidopsis_keys()
    at_proteome = _load_proteome(ARABIDOPSIS_PROTEOME)
    rice_keys, rice_proteome = _rice_keys()
    mg_keys, mg_proteome = _magnaporthe_keys()

    species_specs = [
        ("Arabidopsis", at_keys, at_proteome),
        ("Rice", rice_keys, rice_proteome),
        ("Magnaporthe", mg_keys, mg_proteome),
    ]

    results: list[dict[str, Any]] = []
    sasa_results: list[dict[str, Any]] = []
    for species, keys, proteome in species_specs:
        # --- pathway 1: contact_number_proxy (existing, Cα-only) ---
        rows = build_structural_context_rows(
            persulfidated_keys=keys,
            proteome=proteome,
            structure_dir=ALPHAFOLD_REGISTRY.parent,
            accession_to_structure_file=structure_files,
        )
        if rows:
            result = structural_context_test(species, rows, n_perm=n_perm, seed=seed)
            results.append(
                {
                    "species": result.species,
                    "n_proteins": result.n_proteins,
                    "n_positive_cys": result.n_positive_cys,
                    "n_unlabeled_cys": result.n_unlabeled_cys,
                    "mean_contact_positive": result.mean_contact_positive,
                    "mean_contact_unlabeled": result.mean_contact_unlabeled,
                    "mean_diff": result.mean_diff,
                    "direction": (
                        "positive_more_exposed"
                        if result.mean_diff < 0
                        else "positive_more_buried"
                    ),
                    "p_value_two_sided": result.permutation.p_value,
                    "n_perm": result.permutation.n_perm,
                }
            )
            print(
                f"  [proxy] {species}: proteins={result.n_proteins} "
                f"positive={result.n_positive_cys} "
                f"unlabeled={result.n_unlabeled_cys} "
                f"mean_diff={result.mean_diff:+.3f} "
                f"p={result.permutation.p_value:.4f}"
            )
        else:
            print(f"[proxy] {species}: no rows (no structures resolved)")

        # --- pathway 2: real Shrake-Rupley SASA (independent cross-check,
        # 2026-07-23 self-review addition) ---
        sasa_rows = build_structural_context_rows_sasa(
            persulfidated_keys=keys,
            proteome=proteome,
            accession_to_structure_file=structure_files,
            n_sphere_points=n_sphere_points,
        )
        if not sasa_rows:
            print(f"[SASA] {species}: no rows (no structures resolved)")
            continue
        for metric in ("residue_sasa", "sg_sasa"):
            sasa_result = sasa_structural_context_test(
                species, sasa_rows, metric=metric, n_perm=n_perm, seed=seed
            )
            sasa_results.append(
                {
                    "species": sasa_result.species,
                    "metric": sasa_result.metric_name,
                    "n_proteins": sasa_result.n_proteins,
                    "n_positive_cys": sasa_result.n_positive_cys,
                    "n_unlabeled_cys": sasa_result.n_unlabeled_cys,
                    "mean_sasa_positive_A2": sasa_result.mean_sasa_positive,
                    "mean_sasa_unlabeled_A2": sasa_result.mean_sasa_unlabeled,
                    "mean_diff_A2": sasa_result.mean_diff,
                    "direction": (
                        "positive_more_buried"
                        if sasa_result.mean_diff < 0
                        else "positive_more_exposed"
                    ),
                    "p_value_two_sided": sasa_result.permutation.p_value,
                    "n_perm": sasa_result.permutation.n_perm,
                }
            )
            print(
                f"  [SASA:{metric}] {species}: "
                f"positive={sasa_result.mean_sasa_positive:.2f}A2 "
                f"unlabeled={sasa_result.mean_sasa_unlabeled:.2f}A2 "
                f"diff={sasa_result.mean_diff:+.3f} "
                f"p={sasa_result.permutation.p_value:.4f}"
            )

    directions = {r["species"]: r["direction"] for r in results}
    consistent = len(set(directions.values())) == 1 if directions else False
    print(f"\n[proxy] direction consistent across species: {consistent} ({directions})")

    sasa_directions_by_metric: dict[str, dict[str, str]] = {}
    for r in sasa_results:
        sasa_directions_by_metric.setdefault(r["metric"], {})[r["species"]] = r[
            "direction"
        ]
    sasa_consistency = {
        metric: len(set(d.values())) == 1
        for metric, d in sasa_directions_by_metric.items()
    }
    for metric, dirs in sasa_directions_by_metric.items():
        print(
            f"[SASA:{metric}] direction consistent: {sasa_consistency[metric]} ({dirs})"
        )

    summary: dict[str, Any] = {
        "framing": (
            "mechanistic layer: is the structural accessibility context of "
            "persulfidated cysteines (vs other cysteines in the SAME "
            "proteins) consistent in direction across three independent "
            "species? Two INDEPENDENT measurement pathways are reported: "
            "(1) features/structure.py's contact_number_proxy (documented "
            "accessibility PROXY, not rigorous SASA — Cα-only geometry); "
            "(2) features/sasa.py's genuine Shrake-Rupley SASA on full "
            "heavy-atom geometry (whole-residue and SG-specific), added "
            "2026-07-23 as an independent cross-check per self-review. "
            "Neither pathway silently replaces the other."
        ),
        "pu_semantics": (
            "non-persulfidated cysteines are labelled 'unlabeled' "
            "(not detected), never 'negative' (not confirmed absent)"
        ),
        "contact_proxy_pathway": {
            "per_species": results,
            "direction_consistent_across_species": consistent,
            "directions": directions,
        },
        "sasa_pathway": {
            "per_species_per_metric": sasa_results,
            "direction_consistent_across_species_by_metric": sasa_consistency,
            "directions_by_metric": sasa_directions_by_metric,
            "n_sphere_points": n_sphere_points,
        },
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(f"\nwrote summary -> {output_path}")
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="cross-species structural context of persulfidation targeting"
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("results/cross_species_conservation/structural_context_v1.json"),
    )
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-sphere-points", type=int, default=92)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_structural_context_analysis(
        output_path=args.output,
        n_perm=args.n_perm,
        seed=args.seed,
        n_sphere_points=args.n_sphere_points,
    )
