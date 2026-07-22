"""RED: same-species independent-lab control registration (AtG6PD6 Cys159).

The controls framework was built for cross-species (tomato) controls that are
absent from the benchmark. AtG6PD6 Cys159 (PXD043969, Li lab NWAFU, New
Phytol 2023, doi:10.1111/nph.19188) is a same-species, independent-lab,
biochemically validated persulfidation site: it IS present in the benchmark
as an unlabeled row (PU pool) and must therefore be excluded from the
scorer's training sample before scoring, and it must never be a training
positive.

Expected RED: ``plantpersulf.evaluation.known_controls`` does not exist yet.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from plantpersulf.evaluation.known_controls import (  # RED: module missing
    REGISTERED_CONTROLS,
    filter_out_keys,
    scores_against_benchmark_proteome,
)

BENCHMARK = Path("data/processed/benchmark_v1/sites.tsv")
PROTEOME = Path("data/raw/references/arabidopsis_ref_proteome_v1.fasta")
AT_G6PD6 = ("Q9FJI5", 159)


def _control(lineage: str):
    matches = [c for c in REGISTERED_CONTROLS if c.mechanism_lineage_id == lineage]
    assert len(matches) == 1
    return matches[0]


def test_atg6pd6_control_is_registered_with_same_species_semantics() -> None:
    ctrl = _control("ATG6PD6_H2S_G6PD_SALT")
    assert (ctrl.uniprot_accession, ctrl.cys_position) == AT_G6PD6
    assert ctrl.gene == "AtG6PD6"
    assert ctrl.doi == "10.1111/nph.19188"
    assert ctrl.status == "mapped"
    assert ctrl.control_species == "Arabidopsis thaliana"
    assert ctrl.in_benchmark_as == "unlabeled"
    assert scores_against_benchmark_proteome(ctrl)


def test_tomato_controls_keep_cross_species_semantics() -> None:
    for lineage in ("SLWRKY6_H2S_PHOSPHORYLATION", "SLERFD2_H2S_ETHYLENE"):
        ctrl = _control(lineage)
        assert ctrl.control_species == "Solanum lycopersicum"
        assert ctrl.in_benchmark_as == "absent"
        assert not scores_against_benchmark_proteome(ctrl)


def test_new_control_is_a_distinct_lineage() -> None:
    lineages = [c.mechanism_lineage_id for c in REGISTERED_CONTROLS]
    assert len(set(lineages)) == len(lineages)


def test_filter_out_keys_removes_exactly_the_control_rows() -> None:
    rows = [
        {"protein_accession": "Q9FJI5", "cys_position_in_protein": "159"},
        {"protein_accession": "Q9FJI5", "cys_position_in_protein": "31"},
        {"protein_accession": "Q43727", "cys_position_in_protein": "159"},
    ]
    kept = filter_out_keys(rows, {AT_G6PD6})
    assert [(r["protein_accession"], r["cys_position_in_protein"]) for r in kept] == [
        ("Q9FJI5", "31"),
        ("Q43727", "159"),
    ]


@pytest.mark.skipif(not PROTEOME.is_file(), reason="reference proteome missing")
def test_atg6pd6_cys159_is_cysteine_in_reference_proteome() -> None:
    from plantpersulf.features.sequence import _load_proteome

    proteome = _load_proteome(PROTEOME)
    sequence = proteome["Q9FJI5"]
    assert sequence[158] == "C"


@pytest.mark.skipif(not BENCHMARK.is_file(), reason="benchmark missing")
def test_atg6pd6_cys159_is_unlabeled_never_positive_in_benchmark() -> None:
    labels = []
    with BENCHMARK.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["protein_accession"] == "Q9FJI5":
                labels.append((int(row["cys_position_in_protein"]), row["label"]))
    by_pos = dict(labels)
    assert by_pos[159] == "unlabeled"
    # And no G6PD6 cysteine is a training positive (self-recovery guard).
    assert set(by_pos.values()) == {"unlabeled"}
