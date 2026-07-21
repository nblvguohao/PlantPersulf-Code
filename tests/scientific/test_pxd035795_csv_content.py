from pathlib import Path

PEPTIDE_PATH = Path("data/raw/PXD035795/peptide.csv")
PROTEIN_PATH = Path("data/raw/PXD035795/proteins.csv")
PEPTIDE_SHA256 = (
    "63f2df50a2748d5614a8cfad16cf19c4d7dfac4825d67136382b1fc1fb62ab47"
)
PROTEIN_SHA256 = (
    "cd816699fd3427634c347bb0eb52ee6392c47218bc2bcee5cd55ef367e2115ba"
)


def test_real_peptide_csv_preserves_all_rows_and_has_no_inferred_site() -> None:
    from plantpersulf.evidence.csv_content import parse_peptide_csv

    audit = parse_peptide_csv(PEPTIDE_PATH, PEPTIDE_SHA256)

    assert len(audit.header) == 19
    assert len(audit.rows) == 9326
    assert audit.record_type == "peptide_csv"
    assert audit.rows[0].row_number == 2
    assert audit.rows[0].raw_fields[0] == "ESTLGFVDLLR"
    assert audit.rows[0].raw_fields[4] == "1.2568E9"
    assert audit.rows[0].raw_fields[17] == "O03042|RBL_ARATH"
    assert audit.rows[0].raw_fields[18] == ""
    assert all(
        record.site_localization_status == "not_available"
        for record in audit.rows
    )
    assert all(
        record.evidence_class == "identification_only"
        for record in audit.rows
    )


def test_real_protein_csv_is_never_promoted_to_site_level() -> None:
    from plantpersulf.evidence.csv_content import parse_protein_csv

    audit = parse_protein_csv(PROTEIN_PATH, PROTEIN_SHA256)

    assert len(audit.header) == 28
    assert len(audit.rows) == 1461
    assert audit.record_type == "protein_csv"
    assert audit.rows[0].raw_fields[2] == "O03042|RBL_ARATH"
    assert audit.rows[0].raw_fields[8] == "1.394E7"
    assert audit.rows[0].declared_modifications == (
        "NBF_C",
        "NBF_K",
        "Amidation",
        "Sulphone",
        "Oxidation (HW)",
        "23 more",
    )
    assert all(
        record.site_localization_status == "not_applicable"
        for record in audit.rows
    )
    assert not any(record.evidence_class == "site_ms" for record in audit.rows)


def test_csv_quantification_strings_are_not_converted_or_imputed() -> None:
    from plantpersulf.evidence.csv_content import (
        parse_peptide_csv,
        parse_protein_csv,
    )

    peptide = parse_peptide_csv(PEPTIDE_PATH, PEPTIDE_SHA256)
    protein = parse_protein_csv(PROTEIN_PATH, PROTEIN_SHA256)

    assert all(isinstance(value, str) for value in peptide.rows[0].raw_fields)
    assert all(isinstance(value, str) for value in protein.rows[0].raw_fields)
    assert any("" in record.raw_fields for record in peptide.rows)
    assert any("" in record.raw_fields for record in protein.rows)
