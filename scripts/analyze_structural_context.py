#!/usr/bin/env python
"""Cross-species structural context of persulfidation targeting (2026-07-23).

Mechanistic companion to ``analyze_cross_species_conservation.py``: within
each species, do persulfidated cysteines sit in a systematically different
structural context (AlphaFold contact-number accessibility proxy) than
other cysteines in the SAME proteins? And does the direction of any such
preference agree across four independent species?

Requires the bulk AlphaFold structure download
(``scripts/download_alphafold_structures_bulk.py``) to have populated
``data/registry/alphafold_structures.tsv``.

Usage::

    python scripts/analyze_structural_context.py \
        --output results/cross_species_conservation/structural_context_v2.json
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
KIAE271_XLSX = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
TOMATO_PROTEOME = Path("data/raw/references/tomato_ref_proteome_v1.fasta")
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


def _tomato_keys() -> tuple[set[tuple[str, int]], dict[str, str]]:
    """kiae271 (Zhang et al. 2024) coordinate-verified tomato sites.

    Structures resolve for the kiae271 proteins only, which is fine *here*
    and not fine in the ranking model. This test compares persulfidated
    cysteines against the other cysteines of the SAME protein, so structure
    availability is constant within every comparison unit and cannot be
    confounded with label status. A proteome-wide ranker has no such
    protection: there, giving structures only to proteins that carry a
    published site would make "has a structure" a direct proxy for "is a
    known target", which is why the frozen release registry carries no
    tomato accessions at all (see the release exclusion ledger, EX-003).
    """
    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites

    proteome = _load_proteome(TOMATO_PROTEOME)
    table = parse_kiae271_sites(KIAE271_XLSX, proteome)
    keys = {(s.protein_accession, s.cys_position) for s in table.sites}
    return keys, proteome


def _tomato_keys_by_regulation() -> tuple[
    dict[str, set[tuple[str, int]]], dict[str, str]
]:
    """kiae271 sites split by regulation class.

    kiae271 is a *differential* dataset (SlLCD1-OE vs WT), unlike the three
    steady-state persulfidomes: ``lcd_gain`` sites appear when H2S rises,
    ``wt_only`` sites disappear, and ``both`` are detected in either
    condition. Pooling gain and loss sites could manufacture a structural
    difference that belongs to neither class, so each class is also tested
    on its own.
    """
    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites

    proteome = _load_proteome(TOMATO_PROTEOME)
    table = parse_kiae271_sites(KIAE271_XLSX, proteome)
    strata: dict[str, set[tuple[str, int]]] = {}
    for site in table.sites:
        strata.setdefault(site.regulation, set()).add(
            (site.protein_accession, site.cys_position)
        )
    return strata, proteome


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
    min_plddt: float = 70.0,
) -> dict[str, Any]:
    from plantpersulf.evaluation.cross_species_structural_context import (
        build_structural_context_rows,
        build_structural_context_rows_sasa,
        filter_rows_by_plddt,
        filter_sasa_rows_by_plddt,
        plddt_by_key,
        plddt_context_test,
        sasa_structural_context_test,
        structural_context_test,
    )
    from plantpersulf.features.sequence import _load_proteome

    structure_files = _structure_file_map()
    print(f"AlphaFold structures registered: {len(structure_files)}")

    at_keys = _arabidopsis_keys()
    at_proteome = _load_proteome(ARABIDOPSIS_PROTEOME)
    rice_keys, rice_proteome = _rice_keys()
    tomato_keys, tomato_proteome = _tomato_keys()
    mg_keys, mg_proteome = _magnaporthe_keys()

    species_specs = [
        ("Arabidopsis", at_keys, at_proteome),
        ("Rice", rice_keys, rice_proteome),
        ("Tomato", tomato_keys, tomato_proteome),
        ("Magnaporthe", mg_keys, mg_proteome),
    ]

    # kiae271 is the only differential dataset in the set, so its pooled
    # result is also run per regulation class. These go through the same
    # code path as the species, then are reported separately so they never
    # enter the cross-species direction-consistency call.
    tomato_strata, _ = _tomato_keys_by_regulation()
    stratum_specs = [
        (f"Tomato[{regulation}]", tomato_strata[regulation], tomato_proteome)
        for regulation in sorted(tomato_strata)
    ]
    sensitivity_names = {name for name, _, _ in stratum_specs}

    results: list[dict[str, Any]] = []
    sasa_results: list[dict[str, Any]] = []
    plddt_results: list[dict[str, Any]] = []
    confident_results: list[dict[str, Any]] = []
    confident_sasa_results: list[dict[str, Any]] = []
    for species, keys, proteome in [*species_specs, *stratum_specs]:
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
            # --- model-confidence control ---
            plddt_result = plddt_context_test(
                species, rows, n_perm=n_perm, seed=seed
            )
            plddt_results.append(
                {
                    "species": plddt_result.species,
                    "n_positive_cys": plddt_result.n_positive_cys,
                    "n_unlabeled_cys": plddt_result.n_unlabeled_cys,
                    "mean_plddt_positive": plddt_result.mean_plddt_positive,
                    "mean_plddt_unlabeled": plddt_result.mean_plddt_unlabeled,
                    "mean_diff": plddt_result.mean_diff,
                    "p_value_two_sided": plddt_result.permutation.p_value,
                    "n_perm": plddt_result.permutation.n_perm,
                }
            )
            print(
                f"  [pLDDT] {species}: positive={plddt_result.mean_plddt_positive:.1f} "
                f"unlabeled={plddt_result.mean_plddt_unlabeled:.1f} "
                f"diff={plddt_result.mean_diff:+.2f} "
                f"p={plddt_result.permutation.p_value:.4f}"
            )

            confident_rows = filter_rows_by_plddt(rows, min_plddt=min_plddt)
            if confident_rows:
                confident = structural_context_test(
                    species, confident_rows, n_perm=n_perm, seed=seed
                )
                confident_results.append(
                    {
                        "species": confident.species,
                        "min_plddt": min_plddt,
                        "n_proteins": confident.n_proteins,
                        "n_positive_cys": confident.n_positive_cys,
                        "n_unlabeled_cys": confident.n_unlabeled_cys,
                        "mean_contact_positive": confident.mean_contact_positive,
                        "mean_contact_unlabeled": confident.mean_contact_unlabeled,
                        "mean_diff": confident.mean_diff,
                        "direction": (
                            "positive_more_exposed"
                            if confident.mean_diff < 0
                            else "positive_more_buried"
                        ),
                        "p_value_two_sided": confident.permutation.p_value,
                        "n_perm": confident.permutation.n_perm,
                    }
                )
                print(
                    f"  [proxy>=pLDDT{min_plddt:g}] {species}: "
                    f"positive={confident.n_positive_cys} "
                    f"unlabeled={confident.n_unlabeled_cys} "
                    f"mean_diff={confident.mean_diff:+.3f} "
                    f"p={confident.permutation.p_value:.4f}"
                )
            else:
                print(f"  [proxy>=pLDDT{min_plddt:g}] {species}: no confident rows")
        else:
            print(f"[proxy] {species}: no rows (no structures resolved)")
            rows = []

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

        confident_sasa_rows = filter_sasa_rows_by_plddt(
            sasa_rows, plddt_by_key(rows), min_plddt=min_plddt
        )
        for metric in ("residue_sasa", "sg_sasa"):
            if not confident_sasa_rows:
                print(f"  [SASA>=pLDDT{min_plddt:g}:{metric}] {species}: no rows")
                continue
            confident_sasa = sasa_structural_context_test(
                species, confident_sasa_rows, metric=metric, n_perm=n_perm, seed=seed
            )
            confident_sasa_results.append(
                {
                    "species": confident_sasa.species,
                    "metric": confident_sasa.metric_name,
                    "min_plddt": min_plddt,
                    "n_proteins": confident_sasa.n_proteins,
                    "n_positive_cys": confident_sasa.n_positive_cys,
                    "n_unlabeled_cys": confident_sasa.n_unlabeled_cys,
                    "mean_sasa_positive_A2": confident_sasa.mean_sasa_positive,
                    "mean_sasa_unlabeled_A2": confident_sasa.mean_sasa_unlabeled,
                    "mean_diff_A2": confident_sasa.mean_diff,
                    "direction": (
                        "positive_more_buried"
                        if confident_sasa.mean_diff < 0
                        else "positive_more_exposed"
                    ),
                    "p_value_two_sided": confident_sasa.permutation.p_value,
                    "n_perm": confident_sasa.permutation.n_perm,
                }
            )
            print(
                f"  [SASA>=pLDDT{min_plddt:g}:{metric}] {species}: "
                f"positive={confident_sasa.mean_sasa_positive:.2f}A2 "
                f"unlabeled={confident_sasa.mean_sasa_unlabeled:.2f}A2 "
                f"diff={confident_sasa.mean_diff:+.3f} "
                f"p={confident_sasa.permutation.p_value:.4f}"
            )

    sensitivity_plddt = [r for r in plddt_results if r["species"] in sensitivity_names]
    sensitivity_confident = [
        r for r in confident_results if r["species"] in sensitivity_names
    ]
    sensitivity_confident_sasa = [
        r for r in confident_sasa_results if r["species"] in sensitivity_names
    ]
    plddt_results = [r for r in plddt_results if r["species"] not in sensitivity_names]
    confident_results = [
        r for r in confident_results if r["species"] not in sensitivity_names
    ]
    confident_sasa_results = [
        r for r in confident_sasa_results if r["species"] not in sensitivity_names
    ]

    sensitivity_results = [r for r in results if r["species"] in sensitivity_names]
    sensitivity_sasa = [r for r in sasa_results if r["species"] in sensitivity_names]
    results = [r for r in results if r["species"] not in sensitivity_names]
    sasa_results = [r for r in sasa_results if r["species"] not in sensitivity_names]

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
            "proteins) consistent in direction across four independent "
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
        "model_confidence_control": {
            "question": (
                "Is the accessibility difference partly a model-confidence "
                "difference? Very-low-pLDDT regions are predicted as extended "
                "chain and read as exposed regardless of the real structure, "
                "and mass spectrometry independently favours flexible regions."
            ),
            "min_plddt": min_plddt,
            "plddt_at_persulfidated_vs_other_cysteines": plddt_results,
            "contact_proxy_confident_only": confident_results,
            "sasa_confident_only": confident_sasa_results,
            "confident_directions": {
                r["species"]: r["direction"] for r in confident_results
            },
            "confident_direction_consistent_across_species": (
                len({r["direction"] for r in confident_results}) == 1
                if confident_results
                else False
            ),
        },
        "tomato_regulation_sensitivity": {
            "question": (
                "kiae271 is a differential dataset (SlLCD1-OE vs WT). Does the "
                "pooled tomato direction survive when H2S-gained (lcd_gain), "
                "WT-only (wt_only) and condition-independent (both) sites are "
                "tested separately, or is it an artefact of pooling classes "
                "with opposite biology?"
            ),
            "contact_proxy_pathway": sensitivity_results,
            "sasa_pathway": sensitivity_sasa,
            "plddt_control": sensitivity_plddt,
            "contact_proxy_confident_only": sensitivity_confident,
            "sasa_confident_only": sensitivity_confident_sasa,
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
        default=Path("results/cross_species_conservation/structural_context_v2.json"),
    )
    p.add_argument("--n-perm", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-sphere-points", type=int, default=92)
    p.add_argument("--min-plddt", type=float, default=70.0)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_structural_context_analysis(
        output_path=args.output,
        n_perm=args.n_perm,
        seed=args.seed,
        n_sphere_points=args.n_sphere_points,
        min_plddt=args.min_plddt,
    )
