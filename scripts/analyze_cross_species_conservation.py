#!/usr/bin/env python
"""Cross-species conservation of persulfidation targeting (2026-07-23 reframe).

Rather than "does a model trained on species A predict species B" (Gate 2,
structurally blocked by study non-independence), this asks: **is
co-persulfidation of shared PANTHER ortholog subfamilies enriched beyond
chance across four independent lab/chemistry/species datasets** —
Arabidopsis (benchmark_v1: PXD006140+PXD024061, Seville, tag-switch),
rice (PXD072089, Xie/Nanjing Forestry, NM-biotin+DTT), tomato (kiae271,
Zhang/South China Agri, tag-switch), and Magnaporthe (PXD063170,
Chen/Huazhong Agri, IAA-PEO-biotin)?

Beyond the pairwise tests this also reports the n-way conservation
spectrum: how many shared ortholog subfamilies carry a persulfidated
protein in k of the four species, against an independence null that holds
every species' own targeting rate fixed, with a permutation test on the
all-species cell.

This is associational evidence, not predictive evidence — a standard
comparative-genomics enrichment methodology, using exactly the datasets
Gate 2 rejected for prediction, evaluated with a different, appropriate
statistical question.

Usage::

    python scripts/analyze_cross_species_conservation.py \
        --output results/cross_species_conservation/conservation_v2.json
"""

from __future__ import annotations

import argparse
import csv
import json
from itertools import combinations
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
KIAE271_XLSX = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
TOMATO_PROTEOME = Path("data/raw/references/tomato_ref_proteome_v1.fasta")
PANTHER_TOMATO = Path(
    "data/raw/references/panther_annotations_v1/panther_tomato_taxon4081.tsv"
)

# Manually curated (2026-07-23, via UniProt PANTHER EntryName lookup on one
# representative Arabidopsis accession per family) — recorded here as a
# static annotation layer, not recomputed at run time, so the script stays
# network-free after the registered inputs are in place.
FAMILY_NAMES: dict[str, str] = {
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


def _tomato_positives() -> tuple[set[str], int]:
    """kiae271 (Zhang et al. 2024) coordinate-verified tomato sites.

    This is the same parse the release pipeline uses, so the conservation
    layer and the ranking layer see one tomato site set, not two.
    """
    from plantpersulf.features.sequence import _load_proteome
    from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites

    proteome = _load_proteome(TOMATO_PROTEOME)
    table = parse_kiae271_sites(KIAE271_XLSX, proteome)
    accessions = {s.protein_accession for s in table.sites}
    return accessions, len(table.sites)


def _magnaporthe_positives() -> tuple[set[str], int]:
    from plantpersulf.proteomics.pxd063170_sites import (
        load_ensembl_fungi_proteome,
        parse_pxd063170_sites,
    )

    proteome = load_ensembl_fungi_proteome(MG_PROTEOME)
    table = parse_pxd063170_sites(MG_SITE_TSV, proteome, min_localization=0.75)
    accessions = {s.protein_accession for s in table.sites}
    return accessions, len(table.sites)


def run_conservation_analysis(
    output_path: Path,
    n_perm: int = 10000,
    seed: int = 20260814,
) -> dict[str, Any]:
    from plantpersulf.evaluation.cross_species_conservation import (
        bonferroni_and_bh_correction,
        conservation_detectability_floor,
        conservation_spectrum,
        conservation_spectrum_permutation_test,
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
        KIAE271_XLSX,
        TOMATO_PROTEOME,
        PANTHER_AT,
        PANTHER_RICE,
        PANTHER_MG,
        PANTHER_TOMATO,
    )
    for path in required:
        if not path.is_file():
            raise RuntimeError(f"required input missing: {path}")

    at_positives = _arabidopsis_positives()
    rice_positives, rice_n_sites = _rice_positives()
    tomato_positives, tomato_n_sites = _tomato_positives()
    mg_positives, mg_n_sites = _magnaporthe_positives()
    print(f"Arabidopsis: {len(at_positives)} persulfidated proteins (390 sites)")
    print(f"Rice: {len(rice_positives)} persulfidated proteins ({rice_n_sites} sites)")
    print(
        f"Tomato: {len(tomato_positives)} persulfidated proteins "
        f"({tomato_n_sites} sites)"
    )
    print(
        f"Magnaporthe: {len(mg_positives)} persulfidated proteins ({mg_n_sites} sites)"
    )

    at_map = load_panther_annotations(PANTHER_AT)
    rice_map = load_panther_annotations(PANTHER_RICE)
    tomato_map = load_panther_annotations(PANTHER_TOMATO)

    mg_proteome = load_ensembl_fungi_proteome(MG_PROTEOME)
    mg_accession_to_gene = {acc: mg8_accession_to_gene(acc) for acc in mg_proteome}
    mg_map = load_panther_annotations_by_orf_gene(PANTHER_MG, mg_accession_to_gene)
    print(
        f"MG8 accessions bridged to PANTHER: "
        f"{len(mg_map.family_by_accession)}/{len(mg_proteome)}"
    )

    # Ordered species registry — every pairwise test, the n-way spectrum and
    # the JSON report are derived from this one list, so adding a fifth
    # dataset later means editing exactly one place.
    species_maps = {
        "Arabidopsis": at_map,
        "Rice": rice_map,
        "Tomato": tomato_map,
        "Magnaporthe": mg_map,
    }
    species_positives = {
        "Arabidopsis": at_positives,
        "Rice": rice_positives,
        "Tomato": tomato_positives,
        "Magnaporthe": mg_positives,
    }
    names = list(species_maps)
    pairs = [
        (
            a,
            species_maps[a],
            species_positives[a],
            b,
            species_maps[b],
            species_positives[b],
        )
        for a, b in combinations(names, 2)
    ]

    results: list[dict[str, Any]] = []
    for species_a, map_a, pos_a, species_b, map_b, pos_b in pairs:
        result = cross_species_conservation_test(
            species_a, map_a, pos_a, species_b, map_b, pos_b
        )
        floor = conservation_detectability_floor(result)
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
                # A non-significant pair is only interpretable next to the
                # enrichment the test could have detected at all.
                "minimum_significant_co_persulfidated": (
                    floor.minimum_significant_count
                ),
                "minimum_significant_fold_enrichment": (
                    floor.minimum_significant_fold_enrichment
                ),
            }
        )
        print(
            f"  {species_a} x {species_b}: shared_families={result.shared_families} "
            f"co_persulfidated={result.co_persulfidated_families} "
            f"expected={result.expected_co_persulfidated_under_independence:.2f} "
            f"OR={result.odds_ratio} p={result.p_value:.3e} "
            f"detectable_at>={floor.minimum_significant_count}"
        )

    # Multiple-testing correction across the pairwise tests (2026-07-23
    # self-review finding: raw p-values were previously reported
    # uncorrected). Both Bonferroni (conservative, family-wise error rate)
    # and Benjamini-Hochberg (less conservative, FDR) are reported — they
    # can disagree on a borderline p-value and both must be shown, not just
    # whichever is more flattering. The family grew from 3 to 6 tests when
    # tomato joined, so both corrections tighten accordingly.
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

    # N-way spectrum: how many shared subfamilies are persulfidated in k of
    # the four species, against the independence null with every species'
    # own targeting rate held fixed. Pairwise enrichment does not imply an
    # n-way excess, so the top cell gets its own permutation null.
    spectrum = conservation_spectrum(species_maps, species_positives)
    # Every cell from 2 upward gets its own null. The all-species cell is the
    # headline question, but with one shallow dataset it can be empty by
    # construction, and reporting only that cell would hide the k = 3 excess
    # that carries the actual signal.
    permutations = [
        conservation_spectrum_permutation_test(
            species_maps,
            species_positives,
            species_count=k,
            n_perm=n_perm,
            seed=seed,
        )
        for k in range(2, len(names) + 1)
    ]
    print(
        f"\nShared subfamilies present in all {len(names)} proteomes: "
        f"{spectrum.universe_size}"
    )
    for k in sorted(spectrum.expected_by_species_count):
        observed_k = spectrum.observed_by_species_count.get(k, 0)
        print(
            f"  persulfidated in {k} species: observed {observed_k}, "
            f"expected {spectrum.expected_by_species_count[k]:.2f}"
        )
    for perm in permutations:
        print(
            f"  >= {perm.species_count} species: observed {perm.observed}, "
            f"expected {perm.expected_under_independence:.3f}, "
            f"permutation mean {perm.mean_null:.3f}, p={perm.p_value:.3e} "
            f"({perm.n_perm} permutations, seed {perm.seed})"
        )
    for fam in spectrum.conserved_in_all:
        print(f"    {fam}: {FAMILY_NAMES.get(fam, '(not annotated)')}")

    # The three-species intersection is retained so the v1 report and this
    # one can be compared directly; adding a species can only shrink it.
    at_fam = at_map.persulfidated_families(at_positives)
    rice_fam = rice_map.persulfidated_families(rice_positives)
    mg_fam = mg_map.persulfidated_families(mg_positives)
    triple_conserved = at_fam & rice_fam & mg_fam

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
            "Tomato": {
                "n_proteins": len(tomato_positives),
                "n_sites": tomato_n_sites,
                "source": "kiae271 Supplementary Dataset S1 (coordinate-verified)",
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
        "conservation_spectrum": {
            "species": list(spectrum.species),
            "universe_size": spectrum.universe_size,
            "persulfidated_families_by_species": (
                spectrum.persulfidated_families_by_species
            ),
            "observed_by_species_count": {
                str(k): v for k, v in spectrum.observed_by_species_count.items()
            },
            "expected_by_species_count": {
                str(k): v for k, v in spectrum.expected_by_species_count.items()
            },
            "conserved_in_all_species": {
                "count": len(spectrum.conserved_in_all),
                "panther_subfamily_ids": list(spectrum.conserved_in_all),
                "annotated_names": {
                    fam: FAMILY_NAMES.get(fam, "(not annotated)")
                    for fam in spectrum.conserved_in_all
                },
            },
            "permutation_tests": [
                {
                    "statistic": (
                        f"number of shared subfamilies persulfidated in at "
                        f"least {perm.species_count} species"
                    ),
                    "species_count": perm.species_count,
                    "observed": perm.observed,
                    "expected_under_independence": (perm.expected_under_independence),
                    "mean_null": perm.mean_null,
                    "p_value": perm.p_value,
                    "n_perm": perm.n_perm,
                    "seed": perm.seed,
                }
                for perm in permutations
            ],
            "permutation_note": (
                "Each replicate re-draws every species' persulfidated-family "
                "set uniformly without replacement from the shared universe, "
                "holding that species' family count fixed; only the "
                "cross-species alignment is destroyed. Add-one smoothed, so "
                "the p-value is bounded below by 1/(n_perm+1). Every cell from "
                "2 upward is reported, not only the one with the smallest "
                "p-value."
            ),
        },
        "triple_conserved_families": {
            "count": len(triple_conserved),
            "species": ["Arabidopsis", "Rice", "Magnaporthe"],
            "panther_subfamily_ids": sorted(triple_conserved),
            "annotated_names": {
                fam: FAMILY_NAMES.get(fam, "(not annotated)")
                for fam in sorted(triple_conserved)
            },
            "note": (
                "Retained from the three-species v1 report for direct "
                "comparison; names are a static, manually curated annotation "
                "snapshot (2026-07-23), not recomputed at run time — see "
                "FAMILY_NAMES."
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
        default=Path("results/cross_species_conservation/conservation_v2.json"),
    )
    p.add_argument("--n-perm", type=int, default=10000)
    p.add_argument("--seed", type=int, default=20260814)
    return p


if __name__ == "__main__":
    args = build_parser().parse_args()
    run_conservation_analysis(
        output_path=args.output, n_perm=args.n_perm, seed=args.seed
    )
