import csv
import json
from pathlib import Path

import pytest

from plantpersulf.provenance.audit import assert_registered_input

REFERENCE_REGISTRY = Path("data/registry/reference_sequences.tsv")
REFERENCE_PATHS = (
    Path("data/registry/cache/uniprot/Q93VK9.fasta"),
    Path("data/registry/cache/uniprot/Q9ZW96.fasta"),
)
FIXTURE_MANIFEST = Path("tests/fixtures/real/PXD006140/source_manifest.json")


def _fixture_source(tmp_path: Path):  # type: ignore[no-untyped-def]
    from plantpersulf.proteomics.metadata import ProteomicsSource

    manifest = json.loads(FIXTURE_MANIFEST.read_text(encoding="utf-8"))
    fixture_path = FIXTURE_MANIFEST.parent / manifest["fixture_file"]
    registry_path = tmp_path / "fixture_registry.tsv"
    with registry_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("path", "sha256"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "path": fixture_path.resolve().as_posix(),
                "sha256": manifest["fixture_sha256"],
            }
        )
    return ProteomicsSource(
        study_accession="PXD006140",
        source_file=fixture_path,
        registry_file=registry_path,
        source_sha256=manifest["fixture_sha256"],
    )


@pytest.mark.parametrize("path", REFERENCE_PATHS)
def test_task4_coordinate_reference_is_registered(path: Path) -> None:
    assert_registered_input(path, REFERENCE_REGISTRY)
    accession = path.stem
    header = path.read_text(encoding="utf-8").splitlines()[0]
    assert header.startswith(f">sp|{accession}|")
    assert " OS=Arabidopsis thaliana " in header
    assert " SV=1" in header

    with REFERENCE_REGISTRY.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    matches = [
        row
        for row in rows
        if row["study_accession"] == "PXD006140"
        and row["repository"] == "UniProt"
        and row["protein_accession"] == accession
    ]
    assert len(matches) == 1
    assert matches[0]["sequence_version"] == "1"
    assert matches[0]["scientific_use"] == "coordinate_validation_only"
    assert matches[0]["source_url"] == (
        f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    )


def test_real_q93vk9_and_q9zw96_coordinates_match_registered_sequences(
    tmp_path: Path,
) -> None:
    from plantpersulf.proteomics.peptide_parser import parse_omssa
    from plantpersulf.proteomics.site_normalizer import (
        load_reference_sequence,
        normalize_sites,
    )

    parsed = parse_omssa(_fixture_source(tmp_path))
    references = {
        reference.accession: reference
        for reference in (
            load_reference_sequence(path, REFERENCE_REGISTRY)
            for path in REFERENCE_PATHS
        )
    }
    q93vk9 = next(psm for psm in parsed.psms if psm.spectrum_id == "29")
    q9zw96 = next(psm for psm in parsed.psms if psm.spectrum_id == "34")

    q93_sites = normalize_sites(q93vk9, references, "site_coordinates")
    q9_sites = normalize_sites(q9zw96, references, "site_coordinates")

    assert [
        (site.cys_position_in_peptide, site.cys_position_in_protein)
        for site in q93_sites.sites
    ] == [(10, 24)]
    assert [
        (site.cys_position_in_peptide, site.cys_position_in_protein)
        for site in q9_sites.sites
    ] == [(4, 5), (7, 8)]
    assert not q93_sites.conflicts
    assert not q9_sites.conflicts
    for site in (*q93_sites.sites, *q9_sites.sites):
        assert site.evidence_level == "psm_coordinate_only"
        assert site.source_sha256 == (
            "10fce0b8ac26f564c7c16879b0d6e32cfa14a8ee1116f0f7911cbc758650b43f"
        )
        assert all(
            "negative" not in str(value).lower()
            for value in vars(site).values()
        )


def test_real_record_without_registered_sequence_is_a_conflict(
    tmp_path: Path,
) -> None:
    from plantpersulf.proteomics.peptide_parser import parse_omssa
    from plantpersulf.proteomics.site_normalizer import (
        load_reference_sequence,
        normalize_sites,
    )

    parsed = parse_omssa(_fixture_source(tmp_path))
    q9zw96 = next(psm for psm in parsed.psms if psm.spectrum_id == "34")
    q93 = load_reference_sequence(REFERENCE_PATHS[0], REFERENCE_REGISTRY)

    normalized = normalize_sites(
        q9zw96,
        {q93.accession: q93},
        "site_coordinates",
    )

    assert not normalized.sites
    assert [issue.reason for issue in normalized.conflicts] == [
        "sequence_unavailable"
    ]


def test_protein_level_scope_never_creates_a_site(tmp_path: Path) -> None:
    from plantpersulf.proteomics.peptide_parser import parse_omssa
    from plantpersulf.proteomics.site_normalizer import normalize_sites

    q93vk9 = next(
        psm
        for psm in parse_omssa(_fixture_source(tmp_path)).psms
        if psm.spectrum_id == "29"
    )

    normalized = normalize_sites(q93vk9, {}, "protein_level_only")

    assert not normalized.sites
    assert not normalized.conflicts


def test_site_output_audit_detects_declared_hash_corruption(
    tmp_path: Path,
) -> None:
    from plantpersulf.proteomics.peptide_parser import (
        publish_parsed_proteomics,
    )
    from plantpersulf.proteomics.site_normalizer import audit_site_output

    output_root = tmp_path / "interim"
    output = output_root / "PXD006140/proteomics_parser_v1"
    publish_parsed_proteomics(
        (_fixture_source(tmp_path),),
        REFERENCE_PATHS,
        REFERENCE_REGISTRY,
        output,
    )
    clean = audit_site_output("PXD006140", output_root=output_root)
    assert clean.site_count == 3

    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["outputs"][0]["sha256"] = "0" * 64
    manifest_path.write_bytes(
        (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    )

    with pytest.raises(RuntimeError, match="output SHA256 mismatch"):
        audit_site_output("PXD006140", output_root=output_root)
