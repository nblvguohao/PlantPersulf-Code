def test_fetch_script_accepts_metadata_only_command() -> None:
    from scripts.fetch_registered_data import build_parser

    arguments = build_parser().parse_args(["--metadata-only"])
    assert arguments.metadata_only is True


def test_cli_exposes_registry_audit_command() -> None:
    from plantpersulf.cli import build_parser

    arguments = build_parser().parse_args(["audit-registry"])
    assert arguments.command == "audit-registry"


def test_fetch_script_accepts_registered_download_command() -> None:
    from scripts.fetch_registered_data import build_parser

    arguments = build_parser().parse_args(
        [
            "--accession",
            "PXD006140",
            "--file-class",
            "metadata,results",
        ]
    )
    assert arguments.accession == "PXD006140"
    assert arguments.file_class == "metadata,results"


def test_cli_exposes_file_audit_command() -> None:
    from plantpersulf.cli import build_parser

    arguments = build_parser().parse_args(["audit-files", "--accession", "PXD006140"])
    assert arguments.command == "audit-files"
    assert arguments.accession == "PXD006140"
