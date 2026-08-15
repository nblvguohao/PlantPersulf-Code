"""Unit tests for co-peptide negative evidence classification and registry."""

from __future__ import annotations

from plantpersulf.evidence.copeptide_negatives import (
    classify_peptide_cys,
    load_copeptide_registry,
    write_copeptide_registry,
)

# --- three-state classification ---------------------------------------------


def test_classify_modified_and_unmodified_cys() -> None:
    """SSCMICLPCR with C206/C212 (in-peptide 3/9) modified, C209 (6) not."""
    states = classify_peptide_cys(
        peptide="SSCMICLPCR",
        modified_in_peptide_positions=(3, 9),
        localization_confirmed=True,
        site_determining_ions_confirmed=True,
    )
    assert states == {3: "positive", 6: "negative", 9: "positive"}


def test_classify_returns_undetermined_when_localization_not_confirmed() -> None:
    states = classify_peptide_cys(
        peptide="SSCMICLPCR",
        modified_in_peptide_positions=(3, 9),
        localization_confirmed=False,
        site_determining_ions_confirmed=False,
    )
    assert states == {3: "undetermined", 6: "undetermined", 9: "undetermined"}


def test_classify_unmodified_cys_undetermined_without_site_determining_ions() -> None:
    """The 'unobserved != unmodified' rule: without site-determining ion
    coverage the unmodified Cys must be undetermined, not negative."""
    states = classify_peptide_cys(
        peptide="SSCMICLPCR",
        modified_in_peptide_positions=(3, 9),
        localization_confirmed=True,
        site_determining_ions_confirmed=False,
    )
    assert states[6] == "undetermined"
    assert states[3] == "positive"


def test_classify_rejects_modified_position_that_is_not_cys() -> None:
    import pytest

    with pytest.raises(ValueError, match="modified position"):
        classify_peptide_cys(
            peptide="SSCMICLPCR",
            modified_in_peptide_positions=(2, 9),  # position 2 is S
            localization_confirmed=True,
            site_determining_ions_confirmed=True,
        )


def test_classify_weak_candidate_is_undetermined_even_with_confirmed_site() -> None:
    """A weak alternative localisation (probability < threshold) is never a
    negative: 'not confidently localised' is not 'not modified'."""
    states = classify_peptide_cys(
        peptide="LLCNLCKVKSGACIR",
        modified_in_peptide_positions=(3,),
        localization_confirmed=True,
        site_determining_ions_confirmed=True,
        weak_modified_in_peptide_positions=(6, 13),  # Cys at 3, 6, 13
    )
    assert states == {3: "positive", 6: "undetermined", 13: "undetermined"}


def test_classify_weak_candidate_does_not_block_other_negatives() -> None:
    states = classify_peptide_cys(
        peptide="CKTCLKCV",
        modified_in_peptide_positions=(7,),  # Cys at 1, 4, 7
        localization_confirmed=True,
        site_determining_ions_confirmed=True,
        weak_modified_in_peptide_positions=(4,),
    )
    # Cys1 is a non-candidate with determinate count -> negative; Cys4 weak;
    # Cys7 strong.
    assert states == {1: "negative", 4: "undetermined", 7: "positive"}


def test_classify_rejects_weak_position_that_is_not_cys() -> None:
    import pytest

    with pytest.raises(ValueError, match="weak-modified position"):
        classify_peptide_cys(
            peptide="SSCMICLPCR",
            modified_in_peptide_positions=(3, 9),
            localization_confirmed=True,
            site_determining_ions_confirmed=True,
            weak_modified_in_peptide_positions=(2,),  # position 2 is S
        )


def test_classify_rejects_position_in_both_strong_and_weak() -> None:
    import pytest

    with pytest.raises(ValueError, match="both strong and weak"):
        classify_peptide_cys(
            peptide="SSCMICLPCR",
            modified_in_peptide_positions=(3, 9),
            localization_confirmed=True,
            site_determining_ions_confirmed=True,
            weak_modified_in_peptide_positions=(3,),
        )


def test_classify_rejects_empty_peptide() -> None:
    import pytest

    with pytest.raises(ValueError, match="peptide"):
        classify_peptide_cys(
            peptide="",
            modified_in_peptide_positions=(),
            localization_confirmed=True,
            site_determining_ions_confirmed=True,
        )


# --- registry round-trip ----------------------------------------------------


def test_registry_round_trip_preserves_rows(tmp_path) -> None:
    rows = [
        {
            "site_id": "BRG3_C209_NEG",
            "species": "tomato",
            "protein_accession": "A0A3Q7EW23",
            "cys_position": 209,
            "state": "negative",
            "peptide_sequence": "SSCMICLPCR",
            "peptide_start": 204,
            "peptide_end": 213,
            "in_peptide_position": 6,
            "positive_positions": "206;212",
            "negative_positions": "209",
            "undetermined_positions": "",
            "localization_confirmed": "True",
            "site_determining_ion_coverage": "True",
            "source_doi": "10.1093/plphys/kiad070",
            "source_detail": "Fig. 9B",
            "provenance": "paper-reported MS/MS peptide",
        }
    ]
    path = tmp_path / "copeptide_negatives_v1.tsv"
    write_copeptide_registry(path, rows)
    loaded = load_copeptide_registry(path)
    assert len(loaded) == 1
    assert loaded[0]["site_id"] == "BRG3_C209_NEG"
    assert loaded[0]["state"] == "negative"
    assert loaded[0]["positive_positions"] == "206;212"


def test_registry_load_rejects_duplicate_site_ids(tmp_path) -> None:
    import pytest

    rows = [
        {
            "site_id": "BRG3_C209_NEG",
            "species": "tomato",
            "protein_accession": "A0A3Q7EW23",
            "cys_position": 209,
            "state": "negative",
            "peptide_sequence": "SSCMICLPCR",
            "peptide_start": 204,
            "peptide_end": 213,
            "in_peptide_position": 6,
            "positive_positions": "206;212",
            "negative_positions": "209",
            "undetermined_positions": "",
            "localization_confirmed": "True",
            "site_determining_ion_coverage": "True",
            "source_doi": "10.1093/plphys/kiad070",
            "source_detail": "Fig. 9B",
            "provenance": "paper-reported MS/MS peptide",
        }
    ]
    path = tmp_path / "copeptide_negatives_v1.tsv"
    write_copeptide_registry(path, rows + [dict(rows[0])])
    with pytest.raises(RuntimeError, match="duplicate site_id"):
        load_copeptide_registry(path)


def test_registry_load_rejects_missing_column(tmp_path) -> None:
    import pytest

    path = tmp_path / "copeptide_negatives_v1.tsv"
    path.write_text(
        "site_id\tspecies\tcys_position\nBRG3_C209_NEG\ttomato\t209\n",
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="columns"):
        load_copeptide_registry(path)


def test_registry_load_rejects_invalid_state(tmp_path) -> None:
    import pytest

    rows = [
        {
            "site_id": "BRG3_C209_NEG",
            "species": "tomato",
            "protein_accession": "A0A3Q7EW23",
            "cys_position": 209,
            "state": "bogus",
            "peptide_sequence": "SSCMICLPCR",
            "peptide_start": 204,
            "peptide_end": 213,
            "in_peptide_position": 6,
            "positive_positions": "206;212",
            "negative_positions": "209",
            "undetermined_positions": "",
            "localization_confirmed": "True",
            "site_determining_ion_coverage": "True",
            "source_doi": "10.1093/plphys/kiad070",
            "source_detail": "Fig. 9B",
            "provenance": "paper-reported MS/MS peptide",
        }
    ]
    path = tmp_path / "copeptide_negatives_v1.tsv"
    write_copeptide_registry(path, rows)
    with pytest.raises(RuntimeError, match="state"):
        load_copeptide_registry(path)
