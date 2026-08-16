from pathlib import Path

MZID_PATH = Path("data/raw/PXD035795/peptides_1_1_0.mzid.gz")
MZID_SHA256 = "62105dfde42d9675bc5dc7c83969d971eec57c9e02983659cd2df96ec4d1b5fe"


def _audit():  # type: ignore[no-untyped-def]
    from plantpersulf.evidence.mzidentml import parse_mzidentml

    return parse_mzidentml(MZID_PATH, MZID_SHA256)


def test_real_mzid_preserves_registered_structure() -> None:
    audit = _audit()

    assert audit.namespace == "http://psidev.info/psi/pi/mzIdentML/1.1"
    assert audit.version == "1.1.0"
    assert audit.structural_counts == {
        "DBSequence": 5,
        "Peptide": 11409,
        "PeptideEvidence": 95,
        "SpectrumIdentificationResult": 661,
        "SpectrumIdentificationItem": 699,
        "Modification": 124,
    }


def test_real_mzid_emits_only_unresolved_identified_modifications() -> None:
    audit = _audit()

    assert len(audit.candidates) == 25
    assert {row.evidence_class for row in audit.candidates} == {"unresolved"}
    assert {row.method_mapping_status for row in audit.candidates} == {"absent_v1"}
    assert {row.cv_accession for row in audit.candidates} == {"MS:1001460"}
    assert {row.cv_name for row in audit.candidates} == {"unknown modification"}
    assert {row.cv_value for row in audit.candidates} == {"DCP", "NBF_C", "NBF_N"}
    assert {row.monoisotopic_mass_delta for row in audit.candidates} == {
        "163.0012",
        "168.0786",
    }
    assert all(row.pass_threshold == "true" for row in audit.candidates)
    assert not any(row.evidence_class == "site_ms" for row in audit.candidates)


def test_real_mzid_resolves_exact_references_and_preserves_locators() -> None:
    audit = _audit()

    assert all(row.peptide_id.startswith("PEPTIDE_") for row in audit.candidates)
    assert all(
        row.peptide_evidence_id.startswith("PEPTIDEEVIDENCE_")
        for row in audit.candidates
    )
    assert all(row.db_sequence_id.startswith("DBSEQUENCE_") for row in audit.candidates)
    assert all(row.protein_accession for row in audit.candidates)
    assert all(row.source_sha256 == MZID_SHA256 for row in audit.candidates)
    assert all(
        row.spectrum_result_id in row.source_locator
        and row.spectrum_identification_item_id in row.source_locator
        for row in audit.candidates
    )
    assert all(row.conflict_status in {"clear", "conflict"} for row in audit.candidates)
    assert all(conflict.source_sha256 == MZID_SHA256 for conflict in audit.conflicts)
