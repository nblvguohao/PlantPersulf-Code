#!/usr/bin/env python
"""Cross-species conservation of persulfidation targeting (2026-07-23 reframe).

Rather than "does a model trained on species A predict species B" (Gate 2,
structurally blocked by study non-independence), this asks: **is
co-persulfidation of shared PANTHER ortholog subfamilies enriched beyond
chance across three independent lab/chemistry/species datasets** —
Arabidopsis (benchmark_v1: PXD006140+PXD024061, Seville, tag-switch),
rice (PXD072089, Xie/Nanjing Forestry, NM-biotin+DTT), and Magnaporthe
(PXD063170, Chen/Huazhong Agri, IAA-PEO-biotin)?

This is associational evidence, not predictive evidence — a standard
comparative-genomics enrichment methodology, using exactly the datasets
Gate 2 rejected for prediction, evaluated with a different, appropriate
statistical question.

Usage::

    python scripts/analyze_cross_species_conservation.py \
        --output results/cross_species_conservation/conservation_v1.json
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

BENCHMARK = Path("data/processed/benchmark_v1/sites.tsv")
RICE_SD01 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd01.xlsx")
RICE_SD04 = Path("data/raw/supplements/PXD072089/pnas.2608150123.sd04.xlsx")
RICE_SS = Path("data/raw/supplements/PXD072089/SS-all-peptides.tsv")
RICE_PROTEOME = Path("data/raw/references/rice_proteome_v1/uniprot_rice_v1.fasta")
MG_SITE_TSV = Path("data/raw/supplements/PXD063170/PXD063170_sites_moesm3.tsv")
MG_PROTEOME = Path("data/raw/supplements/PXD063170/Magnaporthe_oryzae.MG8.pep.all.fa")
PANTHER_AT = Path(
    "data/raw/references/panther_annotations_v1/panther_arabidopsis_taxon3702.tsv"
)
PANTHER_RICE = Path(
    "data/raw/references/panther_annotations_v1/panther_rice_taxon4530.tsv"
)
PANTHER_MG = Path(
    "data/raw/references/panther_annotations_v1/panther_magnaporthe_taxon242507.tsv"
)

# Manually curated (2026-07-23, via UniProt PANTHER EntryName lookup on one
# representative Arabidopsis accession per family) — recorded here as a
# static annotation layer, not recomputed at run time, so the script stays
# network-free after the registered inputs are in place.
TRIPLE_CONSERVED_FAMILY_NAMES: dict[str, str] = {
    "PTHR11071:SF561": "Peptidyl-prolyl cis-trans isomerase D-related",
    "PTHR11556:SF1": "Fructose-bisphosphatase",
    "PTHR11773:SF1": "Glycine dehydrogenase (decarboxylating), mitochondrial",
    "PTHR21152:SF24": "Alanine--glyoxylate aminotransferase 1",
    "PTHR43060:SF17": "L-threonate dehydrogenase (aldolase superfamily)",
    "PTHR43105:SF13": "NADH-ubiquinone oxidoreductase 75 kDa subunit, mitochondrial",
    "PTHR43161:SF9": "Sorbitol dehydrogenase",
    "PTHR45982:SF1": "Regulator of chromosome condensation (RCC1)",
}


def _arabidopsis_positives() -> set[str]:
    positives: set[str] = set()
    with BENCHMARK.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["label"] == "positive":
                positives.add(row["protein_accession"])
    return positives


def _rice_positives() -> tuple[set[str], int]:
    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.proteomics.pxd072089_sites import parse_pxd072089_sites

    proteome = _load_proteome(RICE_PROTEOME)
    table = parse_pxd072089_sites(
        RICE_SD01, RICE_SD04, proteome, ss_all_peptides_path=RICE_SS
    )
    accessions = {s.protein_accession for s in table.sites}
    return accessions, table.total_verified_sites


def _magnaporthe_positives() -> tuple[set[str], int]:
    from plantpersulf.proteomics.pxd063170_sites import (
        load_ensembl_fungi_proteome,
        parse_pxd063170_sites,
    )

    proteome = load_ensembl_fungi_proteome(MG_PROTEOME)
    table = parse_pxd063170_sites(MG_SITE_TSV, proteome, min_localization=0.75)
    accessions = {s.protein_accession for s in table.sites}
    return accessions, len(table.sites)


def run_conservation_analysis(output_path: Path) -> dict[str, Any]:
    from plantpersulf.evaluation.cross_species_conservation import (
        SpeciesPantherMap,
        bonferroni_and_bh_correction,
        cross_species_conservation_test,
        load_panther_annotations,
        load_panther_annotations_by_orf_gene,
        mg8_accession_to_gene,
    )
    from plantpersulf.proteomics.pxd063170_sites import load_ensembl_fungi_proteome

    required = (
        BENCHMARK,
        RICE_SD01,
        RICE_SD04,
        RICE_SS,
        RICE_PROTEOME,
        MG_SITE_TSV,
        MG_PROTEOME,
        PANTHER_AT,
        PANTHER_RICE,
        PANTHER_MG,
    )
    for path in required:
        if not path.is_file():
            raise RuntimeError(f"required input missing: {path}")

    at_positives = _arabidopsis_positives()
    rice_positives, rice_n_sites = _rice_positives()
    mg_positives, mg_n_sites = _magnaporthe_positives()
    print(f"Arabidopsis: {len(at_positives)} persulfidated proteins (390 sites)")
    print(f"Rice: {len(rice_positives)} persulfidated proteins ({rice_n_sites} sites)")
    print(
        f"Magnaporthe: {len(mg_positives)} persulfidated proteins ({mg_n_sites} sites)"
    )

    at_map = load_panther_annotations(PANTHER_AT)
    rice_map = load_panther_annotations(PANTHER_RICE)

    mg_proteome = load_ensembl_fungi_proteome(MG_PROTEOME)
    mg_accession_to_gene = {acc: mg8_accession_to_gene(acc) for acc in mg_proteome}
    mg_map = load_panther_annotations_by_orf_gene(PANTHER_MG, mg_accession_to_gene)
    print(
        f"MG8 accessions bridged to PANTHER: "
        f"{len(mg_map.family_by_accession)}/{len(mg_proteome)}"
    )

    PairSpec = tuple[str, SpeciesPantherMap, set[str], str, SpeciesPantherMap, set[str]]
    pairs: list[PairSpec] = [
        ("Arabidopsis", at_map, at_positives, "Rice", rice_map, rice_positives),
        (
            "Arabidopsis",
            at_map,
            at_positives,
            "Magnaporthe",
            mg_map,
            mg_positives,
        ),
        ("Rice", rice_map, rice_positives, "Magnaporthe", mg_map, mg_positives),
    ]

    results: list[dict[str, Any]] = []
    for species_a, map_a, pos_a, species_b, map_b, pos_b in pairs:
        result = cross_species_conservation_test(
            species_a, map_a, pos_a, species_b, map_b, pos_b
        )
        results.append(
            {
                "species_a": result.species_a,
                "species_b": result.species_b,
                "shared_families": result.shared_families,
                "persulfidated_families_a": result.persulfidated_families_a,
                "persulfidated_families_b": result.persulfidated_families_b,
                "co_persulfidated_families": result.co_persulfidated_families,
                "expected_co_persulfidated_under_independence": (
                    result.expected_co_persulfidated_under_independence
                ),
                "fold_enrichment": (
                    result.co_persulfidated_families
                    / result.expected_co_persulfidated_under_independence
                    if result.expected_co_persulfidated_under_independence > 0
                    else None
                ),
                "odds_ratio": result.odds_ratio,
                "p_value": result.p_value,
            }
        )
        print(
            f"  {species_a} x {species_b}: shared_families={result.shared_families} "
            f"co_persulfidated={result.co_persulfidated_families} "
            f"expected={result.expected_co_persulfidated_under_independence:.2f} "
            f"OR={result.odds_ratio} p={result.p_value:.3e}"
        )

    # Multiple-testing correction across the 3 pairwise tests (2026-07-23
    # self-review finding: raw p-values were previously reported
    # uncorrected). Both Bonferroni (conservative, family-wise error rate)
    # and Benjamini-Hochberg (less conservative, FDR) are reported — they
    # can disagree on a borderline p-value and both must be shown, not just
    # whichever is more flattering.
    corrected = bonferroni_and_bh_correction([r["p_value"] for r in results])
    for r, c in zip(results, corrected, strict=True):
        r["bonferroni_p_value"] = c.bonferroni_p_value
        r["benjamini_hochberg_q_value"] = c.benjamini_hochberg_q_value
        r["significant_bonferroni_0.05"] = c.significant_bonferroni
        r["significant_bh_0.05"] = c.significant_bh
        print(
            f"  corrected: {r['species_a']} x {r['species_b']}: "
            f"bonferroni={c.bonferroni_p_value:.3e} (sig={c.significant_bonferroni}) "
            f"BH_q={c.benjamini_hochberg_q_value:.3e} (sig={c.significant_bh})"
        )

    at_fam = at_map.persulfidated_families(at_positives)
    rice_fam = rice_map.persulfidated_families(rice_positives)
    mg_fam = mg_map.persulfidated_families(mg_positives)
    triple_conserved = at_fam & rice_fam & mg_fam
    print(f"\nFamilies co-persulfidated in ALL THREE species: {len(triple_conserved)}")
    for fam in sorted(triple_conserved):
        print(f"  {fam}")

    summary: dict[str, Any] = {
        "framing": (
            "associational (co-persulfidation enrichment across shared "
            "PANTHER ortholog subfamilies), NOT predictive — "
            "does not bear on Gate 2 conditions"
        ),
        "species_positives": {
            "Arabidopsis": {
                "n_proteins": len(at_positives),
                "n_sites": 390,
                "source": "benchmark_v1 (PXD006140+PXD024061)",
            },
            "Rice": {
                "n_proteins": len(rice_positives),
                "n_sites": rice_n_sites,
                "source": "PXD072089 (SD01+SD04+SS-all-peptides)",
            },
            "Magnaporthe": {
                "n_proteins": len(mg_positives),
                "n_sites": mg_n_sites,
                "source": "PXD063170",
            },
        },
        "ortholog_bridge_coverage": {
            "magnaporthe_mg8_to_panther": (
                f"{len(mg_map.family_by_accession)}/{len(mg_proteome)}"
            ),
        },
        "pairwise_tests": results,
        "triple_conserved_families": {
            "count": len(triple_conserved),
            "panther_subfamily_ids": sorted(triple_conserved),
            "annotated_names": {
                fam: TRIPLE_CONSERVED_FAMILY_NAMES.get(fam, "(not annotated)")
                for fam in sorted(triple_conserved)
            },
            "note": (
                "Families with >=1 persulfidated protein in all three "
                "independent species; names are a static, manually curated "
                "annotation snapshot (2026-07-23), not recomputed at run "
                "time — see TRIPLE_CONSERVED_FAMILY_NAMES."
            ),
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
        description="cross-species persulfidation-targeting conservation analysis"
    )
    p.add_argument(
        "--output",
        type=Path,
        default=Path("results/cross_species_conservation/conservation_v1.json"),
    )
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_conservation_analysis(output_path=args.output)
