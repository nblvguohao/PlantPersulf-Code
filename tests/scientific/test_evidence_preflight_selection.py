from pathlib import Path

import pytest
import yaml

from plantpersulf.download.registered import resolve_registered_selection

POLICY_PATH = Path("configs/evidence_preflight_v1.yaml")
SELECTION_PATH = Path("configs/download_selection.yaml")
FILES_REGISTRY = Path("data/registry/files.tsv")

EXPECTED = {
    "PXD024061": {"checksums": ("checksum.txt",)},
    "PXD035795": {
        "design": ("SDRF.txt",),
        "results": (
            "peptide.csv",
            "peptides_1_1_0.mzid.gz",
            "proteins.csv",
        ),
        "checksums": ("checksum.txt",),
    },
    "PXD039999": {"checksums": ("checksum.txt",)},
}


def test_evidence_preflight_policy_is_exact() -> None:
    from plantpersulf.evidence.policy import load_preflight_policy

    policy = load_preflight_policy(POLICY_PATH)

    assert policy.version == 1
    assert policy.studies == (
        "PXD006140",
        "PXD024061",
        "PXD035795",
        "PXD039999",
    )
    assert policy.evidence_classes == (
        "site_ms",
        "site_mutagenesis",
        "site_biochemical",
        "protein_level_only",
        "identification_only",
        "unresolved",
    )
    assert policy.large_file_threshold_bytes == 104_857_600
    assert policy.prohibited_label_terms == (
        "positive",
        "unlabeled",
        "negative",
    )


@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("version", 2, "version"),
        ("studies", ["PXD006140", "PXD006140"], "studies"),
        ("studies", ["NOT_A_PXD"], "studies"),
        ("evidence_classes", ["positive"], "evidence classes"),
        ("large_file_threshold_bytes", 100, "threshold"),
        ("prohibited_label_terms", ["negative"], "prohibited"),
    ),
)
def test_evidence_preflight_policy_rejects_changed_contract(
    tmp_path: Path,
    field: str,
    value: object,
    message: str,
) -> None:
    from plantpersulf.evidence.policy import load_preflight_policy

    policy = {
        "version": 1,
        "studies": [
            "PXD006140",
            "PXD024061",
            "PXD035795",
            "PXD039999",
        ],
        "evidence_classes": [
            "site_ms",
            "site_mutagenesis",
            "site_biochemical",
            "protein_level_only",
            "identification_only",
            "unresolved",
        ],
        "large_file_threshold_bytes": 104_857_600,
        "prohibited_label_terms": ["positive", "unlabeled", "negative"],
    }
    policy[field] = value
    path = tmp_path / "invalid_policy.yaml"
    path.write_text(yaml.safe_dump(policy), encoding="utf-8")

    with pytest.raises(RuntimeError, match=message):
        load_preflight_policy(path)


def test_stage_a_selection_resolves_exact_real_pride_files() -> None:
    observed: dict[str, dict[str, tuple[str, ...]]] = {}
    for accession, file_classes in EXPECTED.items():
        observed[accession] = {}
        for file_class in file_classes:
            selection = resolve_registered_selection(
                accession,
                (file_class,),
                SELECTION_PATH,
                FILES_REGISTRY,
            )
            observed[accession][file_class] = tuple(
                item.file_name for item in selection
            )
            for item in selection:
                assert item.record_type == "source_file"
                assert item.source_url.startswith(("https://", "ftp://"))
                assert item.size_bytes > 0
                assert item.remote_checksum_algorithm == "SHA1"
                assert len(item.remote_checksum) == 40
                assert item.file_name not in {
                    "txt_persulfproject.zip",
                    "proteinSeq.txt",
                    "arabidopsis_uniprot_072020_identified.fasta",
                }
    assert observed == EXPECTED
