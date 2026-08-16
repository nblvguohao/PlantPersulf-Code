"""Real-input exclusions for the tomato candidate registry."""

from pathlib import Path

from plantpersulf.features.sequence import _load_proteome
from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites
from plantpersulf.proteomics.tomato_candidate_registry import (
    build_tomato_candidate_registry,
)
from plantpersulf.provenance.audit import assert_registered_input


def test_real_registry_excludes_every_known_positive_and_keeps_only_cys() -> None:
    xlsx = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
    fasta = Path("data/raw/references/tomato_ref_proteome_v1.fasta")
    supplementary_registry = Path("data/registry/supplementary_sources.tsv")
    model_input_registry = Path("data/registry/model_inputs.tsv")
    assert_registered_input(xlsx, supplementary_registry)
    assert_registered_input(fasta, model_input_registry)
    proteome = _load_proteome(fasta)
    positive_keys = {
        (site.protein_accession, site.cys_position)
        for site in parse_kiae271_sites(xlsx, proteome).sites
    }
    registry = build_tomato_candidate_registry(
        xlsx, fasta, supplementary_registry, model_input_registry
    )
    keys = {(row.protein_accession, row.cys_position) for row in registry.candidates}

    assert keys.isdisjoint(positive_keys)
    assert all(proteome[accession][position - 1] == "C" for accession, position in keys)
