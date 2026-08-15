"""Scientific tests: co-peptide negative registry entries are real Cys at
real coordinates in the registered reference proteomes (tomato v1,
arabidopsis v2), their peptides match uniquely, and their three-state
labels are consistent with the registered peptide patterns.

These tests use only registered real data (no invented biology).
"""

from __future__ import annotations

from pathlib import Path

from plantpersulf.evidence.copeptide_negatives import load_copeptide_registry
from plantpersulf.features.sequence import _load_proteome

_REPO_ROOT = Path(__file__).resolve().parents[2]
REGISTRY = _REPO_ROOT / "data" / "registry" / "copeptide_negatives_v1.tsv"
PROTEOMES = {
    "tomato": (
        _REPO_ROOT / "data" / "raw" / "references" / "tomato_ref_proteome_v1.fasta"
    ),
    "arabidopsis": (
        _REPO_ROOT / "data" / "raw" / "references" / "arabidopsis_ref_proteome_v2.fasta"
    ),
}


def _registered_rows() -> list[dict[str, str]]:
    return load_copeptide_registry(REGISTRY)


def test_registry_is_nonempty_and_has_both_negatives() -> None:
    rows = _registered_rows()
    assert len(rows) >= 4
    assert {r["state"] for r in rows} <= {"positive", "negative", "undetermined"}
    negatives = [r for r in rows if r["state"] == "negative"]
    assert {r["site_id"] for r in negatives} >= {"BRG3_C209_NEG", "RNF144B_C127_NEG"}


def test_every_entry_is_a_true_cys_at_its_peptide_coordinates() -> None:
    proteomes = {name: _load_proteome(path) for name, path in PROTEOMES.items()}
    for record in _registered_rows():
        proteome = proteomes[record["species"]]
        sequence = proteome[record["protein_accession"]]
        position = int(record["cys_position"])
        assert sequence[position - 1] == "C", record["site_id"]
        start = int(record["peptide_start"])
        end = int(record["peptide_end"])
        assert (
            sequence[start - 1 : end] == record["peptide_sequence"]
        ), record["site_id"]
        # the peptide matches exactly once in its reference proteome
        assert sequence.count(record["peptide_sequence"]) == 1, record["site_id"]
        # in-peptide position is the same Cys as the protein position
        in_peptide = int(record["in_peptide_position"])
        assert record["peptide_sequence"][in_peptide - 1] == "C", record["site_id"]


def test_every_registered_peptide_has_at_least_one_positive() -> None:
    """A registered peptide must carry a confidently localised modification —
    otherwise no unmodified Cys of it may be called negative."""
    rows = _registered_rows()
    peptides = {
        (r["species"], r["protein_accession"], r["peptide_sequence"]) for r in rows
    }
    for key in peptides:
        peptide_rows = [
            r
            for r in rows
            if (r["species"], r["protein_accession"], r["peptide_sequence"]) == key
        ]
        assert any(r["state"] == "positive" for r in peptide_rows), key
        assert peptide_rows[0]["localization_confirmed"] == "True"


def test_arabidopsis_entries_are_from_pxd024061_and_consistent() -> None:
    rows = _registered_rows()
    arabidopsis = [r for r in rows if r["species"] == "arabidopsis"]
    assert len(arabidopsis) >= 12
    for record in arabidopsis:
        assert record["source_doi"] == "10.3390/antiox10040508"
        assert record["site_id"].startswith("PXD024061_")
        assert record["site_determining_ion_coverage"] == "True"
    # every arabidopsis peptide contributes at least one negative
    peptides = {r["peptide_sequence"] for r in arabidopsis}
    for peptide in peptides:
        peptide_rows = [r for r in arabidopsis if r["peptide_sequence"] == peptide]
        assert any(r["state"] == "negative" for r in peptide_rows), peptide


def test_brg3_peptide_modification_pattern_matches_labels() -> None:
    rows = _registered_rows()
    brg3 = [
        r
        for r in rows
        if r["protein_accession"] == "A0A3Q7EW23"
        and r["peptide_sequence"] == "SSCMICLPCR"
    ]
    assert {r["cys_position"] for r in brg3} == {"206", "209", "212"}
    by_state = {
        state: {r["cys_position"] for r in brg3 if r["state"] == state}
        for state in {"positive", "negative", "undetermined"}
    }
    assert by_state["positive"] == {"206", "212"}
    assert by_state["negative"] == {"209"}
    # the peptide embeds all three Cys pairwise 3 residues apart
    peptide = brg3[0]["peptide_sequence"]
    cys_positions = [i + 1 for i, aa in enumerate(peptide) if aa == "C"]
    assert cys_positions == [3, 6, 9]


def test_rnf144b_peptide_modification_pattern_matches_labels() -> None:
    rows = _registered_rows()
    rnf = [
        r
        for r in rows
        if r["protein_accession"] == "A0A3Q7GXU6"
        and r["peptide_sequence"] == "FYCPYKDCSAMLVNDSDEIVR"
    ]
    assert {r["cys_position"] for r in rnf} == {"122", "127"}
    by_state = {
        state: {r["cys_position"] for r in rnf if r["state"] == state}
        for state in {"positive", "negative", "undetermined"}
    }
    assert by_state["positive"] == {"122"}
    assert by_state["negative"] == {"127"}


def test_negative_entries_never_overlap_positive_positions() -> None:
    for record in _registered_rows():
        if record["state"] != "negative":
            continue
        positives = {p for p in record["positive_positions"].split(";") if p}
        assert record["cys_position"] not in positives
        assert record["cys_position"] in {
            p for p in record["negative_positions"].split(";") if p
        }


def test_every_entry_has_registered_source_provenance() -> None:
    for record in _registered_rows():
        assert record["source_doi"]
        assert record["source_detail"]
        assert record["provenance"]
