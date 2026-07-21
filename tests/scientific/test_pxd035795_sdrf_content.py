from pathlib import Path

SDRF_PATH = Path("data/raw/PXD035795/SDRF.txt")
FILES_PATH = Path("data/registry/files.tsv")
SDRF_SHA256 = "335934a0c9b307c99cbdd1c6be05d9da3813d4ee0d9645006f3cbf528a9c5fd1"


def _audit():  # type: ignore[no-untyped-def]
    from plantpersulf.evidence.sdrf import parse_sdrf

    return parse_sdrf(SDRF_PATH, SDRF_SHA256, FILES_PATH)


def test_real_sdrf_preserves_duplicate_headers_and_rows() -> None:
    audit = _audit()

    assert len(audit.header) == 26
    assert len(audit.rows) == 6
    assert audit.duplicate_header_positions[
        "comment[modification parameters]"
    ] == (19, 20, 21, 22, 23)
    assert audit.rows[0].source_name == "Sample 1"
    assert audit.rows[0].assay_name == "Run 1"
    assert audit.rows[0].replicate == "1"
    assert audit.rows[0].data_file_raw == "NPC_rep1.raw, 0h%25201.mgf"
    assert audit.rows[0].modifications == (
        "NT= NBF_N;  MT=variable; PP=Anywhere; AC=.; TA=N; MM=163.00",
        "NT= NBF_K;  MT=variable; PP=Anywhere; AC=.; TA=K; MM=163.00",
        "NT= NBF_C;  MT=variable; PP=Anywhere; AC=.; TA=C; MM=163.00",
        "NT= DCP;  MT=variable; PP=Anywhere; AC=.; TA=C; MM=168.08",
        "NT= dcp-ac;  MT=variable; PP=Anywhere; AC=.; TA=C; MM=196.08",
    )


def test_real_sdrf_maps_only_exact_registry_filenames() -> None:
    audit = _audit()

    assert len(audit.file_mappings) == 12
    exact = [row for row in audit.file_mappings if row.mapping_status == "exact"]
    conflicted = [
        row for row in audit.file_mappings if row.mapping_status == "conflict"
    ]
    assert [row.data_file for row in exact] == [
        "NPC_rep1.raw",
        "NPC_rep2.raw",
        "NPC_rep3.raw",
        "APC_rep1.raw",
        "APC_rep2.raw",
        "APC_rep3.raw",
    ]
    assert [row.data_file for row in conflicted] == [
        "0h%25201.mgf",
        "0h%202.mgf",
        "0h%203.mgf",
        "3d%201.mgf",
        "3d%202.mgf",
        "3d%203.mgf",
    ]
    assert all(row.matched_registry_file == row.data_file for row in exact)
    assert all(row.matched_registry_file == "" for row in conflicted)
    assert len(audit.conflicts) == 6
    assert {row.conflict_type for row in audit.conflicts} == {
        "unmatched_exact_filename"
    }


def test_real_sdrf_does_not_decode_or_rewrite_mgf_names() -> None:
    audit = _audit()

    conflicted_names = {
        row.data_file
        for row in audit.file_mappings
        if row.mapping_status == "conflict"
    }
    assert "0h201.mgf" not in conflicted_names
    assert "0h 201.mgf" not in conflicted_names
    assert "0h%25201.mgf" in conflicted_names
