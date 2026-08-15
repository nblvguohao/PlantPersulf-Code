"""Unit tests for MaxQuant sites-table → co-peptide negative extraction.

Covers the parser of the ``Sulfide(C) Probabilities``-style column, the
per-peptide merge across sites-table rows, the three-state classification
with weak-candidate handling, and unique-peptide localisation in a
reference proteome. RED phase of the PXD024061 expansion.
"""

from __future__ import annotations

import pytest

from plantpersulf.evidence.maxquant_sites_negatives import (
    PeptideEvidence,
    classify_evidence,
    evidence_from_site_rows,
    locate_peptide,
    parse_modified_sequence,
)


def test_parse_modified_sequence_unmodified_peptide() -> None:
    sequence, candidates = parse_modified_sequence("SSCMICLPCR")
    assert sequence == "SSCMICLPCR"
    assert candidates == {}


def test_parse_modified_sequence_single_candidate_with_probability() -> None:
    sequence, candidates = parse_modified_sequence("VEAAMVNARIC(1)KTVR")
    assert sequence == "VEAAMVNARIC(1)KTVR".replace("(1)", "")
    assert candidates == {11: 1.0}  # C is the 11th residue


def test_parse_modified_sequence_multiple_candidates() -> None:
    sequence, candidates = parse_modified_sequence("VPSPTC(0.5)WC(0.5)SK")
    assert sequence == "VPSPTCWCSK"
    assert candidates == {6: 0.5, 8: 0.5}


def test_parse_modified_sequence_mixed_probabilities() -> None:
    sequence, candidates = parse_modified_sequence(
        "LLC(0.825)NLC(0.148)KVKSGAC(0.027)IR"
    )
    assert sequence == "LLCNLCKVKSGACIR"
    assert candidates == {3: 0.825, 6: 0.148, 13: 0.027}  # C at 3, 6, 13


def test_parse_modified_sequence_rejects_garbage() -> None:
    with pytest.raises(ValueError):
        parse_modified_sequence("")
    with pytest.raises(ValueError):
        parse_modified_sequence("(0.5)VPSPTCWCSK")  # probability with no residue


def test_evidence_from_site_rows_merges_ambiguous_rows() -> None:
    rows = [
        {
            "Proteins": "A0A1I9LQH7",
            "Mod. peptide IDs": "104915",
            "Number of Sulfide(C)": 1.0,
            "Sulfide(C) Probabilities": "VPSPTC(0.5)WC(0.5)SK",
        },
        {
            "Proteins": "A0A1I9LQH7",
            "Mod. peptide IDs": "104915",
            "Number of Sulfide(C)": 1.0,
            "Sulfide(C) Probabilities": "VPSPTC(0.5)WC(0.5)SK",
        },
    ]
    evidence = evidence_from_site_rows(
        rows, probability_column="Sulfide(C) Probabilities",
        number_column="Number of Sulfide(C)",
    )
    assert evidence.sequence == "VPSPTCWCSK"
    assert evidence.candidates == ((6, 0.5), (8, 0.5))
    assert evidence.n_modified == 1
    assert evidence.source_rows == 2


def test_evidence_from_site_rows_infers_count_from_strong_candidates() -> None:
    rows = [
        {
            "Mod. peptide IDs": "47943",
            "Number of Sulfide(C)": float("nan"),
            "Sulfide(C) Probabilities": "IVDC(1)KPNVILTC(1)NAVKR",
        }
    ]
    evidence = evidence_from_site_rows(
        rows, probability_column="Sulfide(C) Probabilities",
        number_column="Number of Sulfide(C)",
    )
    assert evidence.n_modified == 2


def test_evidence_weak_candidate_without_count_is_indeterminate() -> None:
    rows = [
        {
            "Mod. peptide IDs": "49489",
            "Number of Sulfide(C)": float("nan"),
            "Sulfide(C) Probabilities": "KCTLC(0.5)KEC(0.5)VR",
        }
    ]
    evidence = evidence_from_site_rows(
        rows, probability_column="Sulfide(C) Probabilities",
        number_column="Number of Sulfide(C)",
    )
    assert evidence.n_modified is None


def test_evidence_from_site_rows_rejects_contradictory_count() -> None:
    rows = [
        {
            "Mod. peptide IDs": "1",
            "Number of Sulfide(C)": 1.0,
            "Sulfide(C) Probabilities": "VLAC(1)VVC(1)GR",
        }
    ]
    with pytest.raises(RuntimeError):
        evidence_from_site_rows(
            rows, probability_column="Sulfide(C) Probabilities",
            number_column="Number of Sulfide(C)",
        )


def test_evidence_from_site_rows_rejects_candidate_on_non_cys() -> None:
    rows = [
        {
            "Mod. peptide IDs": "1",
            "Sulfide(C) Probabilities": "VPSPT(0.5)CWCSK",
        }
    ]
    with pytest.raises(RuntimeError, match="not Cys"):
        evidence_from_site_rows(
            rows, probability_column="Sulfide(C) Probabilities",
            number_column="Number of Sulfide(C)",
        )


def test_evidence_from_site_rows_rejects_inconsistent_sequences() -> None:
    rows = [
        {
            "Mod. peptide IDs": "1",
            "Sulfide(C) Probabilities": "VLAC(1)VVC(1)GR",
        },
        {
            "Mod. peptide IDs": "1",
            "Sulfide(C) Probabilities": "VLAC(1)VVC(1)GK",
        },
    ]
    with pytest.raises(RuntimeError):
        evidence_from_site_rows(
            rows, probability_column="Sulfide(C) Probabilities",
            number_column="Number of Sulfide(C)",
        )


def test_classify_evidence_strong_positive_rest_negative() -> None:
    evidence = PeptideEvidence(
        sequence="SDEVKACIVTCGGLCPGINTVIR",
        candidates=((7, 1.0),),
        n_modified=1,
        source_rows=1,
    )
    states = classify_evidence(evidence)
    # Cys at in-peptide 7 (modified), 11 and 15 (unmodified)
    assert states[7] == "positive"
    assert states[11] == "negative"
    assert states[15] == "negative"


def test_classify_evidence_ambiguous_candidates_are_undetermined() -> None:
    evidence = PeptideEvidence(
        sequence="VPSPTCWCSK",
        candidates=((6, 0.5), (8, 0.5)),
        n_modified=1,
        source_rows=1,
    )
    states = classify_evidence(evidence)
    assert states == {6: "undetermined", 8: "undetermined"}


def test_classify_evidence_weak_candidate_never_negative() -> None:
    evidence = PeptideEvidence(
        sequence="LLCNLCKVKSGACIR",
        candidates=((3, 0.825), (6, 0.148), (13, 0.027)),
        n_modified=1,
        source_rows=1,
    )
    states = classify_evidence(evidence)
    assert states[3] == "positive"
    assert states[6] == "undetermined"
    assert states[13] == "undetermined"


def test_classify_evidence_indeterminate_count_blocks_negatives() -> None:
    evidence = PeptideEvidence(
        sequence="KCTLCKECVR",
        candidates=((5, 0.5), (8, 0.5)),
        n_modified=None,
        source_rows=1,
    )
    states = classify_evidence(evidence)
    assert set(states.values()) == {"undetermined"}


def test_locate_peptide_unique_match() -> None:
    assert locate_peptide("SSCMICLPCR", "MSSCMICLPCRR") == 2
    assert locate_peptide("SSCMICLPCR", "SSCMICLPCR") == 1


def test_locate_peptide_absent_or_ambiguous_returns_none() -> None:
    assert locate_peptide("SSCMICLPCR", "SSCMICLPCK") is None
    assert locate_peptide("CYS", "CYSXYCYS") is None
