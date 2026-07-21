"""RED (Phase A2, Route S/H): PXD006140 Dataset S3 → residue-resolved site_ms.

First executable Phase A2 target from
``docs/superpowers/specs/2026-07-21-evidence-site-acquisition-design.md`` §5.5/§9.

Dataset S3's ``ident_peptides`` records the authors' own custom persulfidation
PTMs on cysteine — ``MOD:99998`` (Sulfide), ``MOD:99997``/``MOD:99996``
(CN-Biotin-Sulfide) — with a peptide-internal position (``MOD:99998 C26``) and
protein coordinates (``inferredCoords``). This test pins the extraction contract:

- ONLY a cysteine carrying an author-designated persulfidation PTM becomes a
  ``site_ms`` site; ``MOD:99999`` (MSBT control), ``MOD:00110`` (artifact) and
  unmodified cysteines never do;
- ``is_decoy`` rows are dropped entirely;
- protein position = ``inferredCoords`` span-start + peptide-position − 1, and the
  residue is verified to be C against the registered UniProt sequence;
- isoform multi-mapped spectra go to a conflict table, never silently first-picked;
- a multi-cysteine peptide yields a site ONLY at the modified cysteine.

Driven by a provenance-locked real micro-fixture cut verbatim from Dataset S3
(``tests/fixtures/real/PXD006140_dataset_s3/``; biological values unchanged).

Expected RED: ``plantpersulf.proteomics.persulfidation_sites`` does not exist yet,
so collection fails with ImportError. GREEN implements the minimal parser.
"""

from __future__ import annotations

import json
from pathlib import Path

from plantpersulf.proteomics.persulfidation_sites import (  # RED: module missing
    PERSULFIDATION_MODS,
    parse_dataset_s3_sites,
)

from plantpersulf.proteomics.metadata import ReferenceSequence
from plantpersulf.provenance.hashing import hash_file

FIXTURE = Path("tests/fixtures/real/PXD006140_dataset_s3")
SUBSET_TSV = FIXTURE / "ident_peptides_subset.tsv"
MANIFEST = json.loads((FIXTURE / "source_manifest.json").read_text(encoding="utf-8"))


def _references() -> dict[str, ReferenceSequence]:
    references: dict[str, ReferenceSequence] = {}
    for path in sorted(FIXTURE.glob("*.fasta")):
        lines = path.read_text(encoding="utf-8").splitlines()
        accession = lines[0].split("|")[1]
        references[accession] = ReferenceSequence(
            accession=accession,
            sequence="".join(line for line in lines[1:] if line),
            sequence_version="1",
            source_file=path,
            source_sha256=hash_file(path, "sha256"),
        )
    return references


def _normalized():  # type: ignore[no-untyped-def]
    return parse_dataset_s3_sites(
        source_tsv=SUBSET_TSV,
        references=_references(),
        source_sha256=MANIFEST["source_sha256"],
    )


def test_persulfidation_mod_codes_exclude_control_and_artifact() -> None:
    assert set(PERSULFIDATION_MODS) == {"99998", "99997", "99996"}
    assert "99999" not in PERSULFIDATION_MODS  # MSBT free-thiol control
    assert "00110" not in PERSULFIDATION_MODS  # L-cysteine methyl disulfide


def test_only_author_designated_persulfidation_cysteines_become_sites() -> None:
    result = _normalized()
    sites = {
        (s.protein_accession_raw, s.cys_position_in_protein) for s in result.sites
    }

    # O03042 C284 (MOD:99998) + O03042 C247 (MOD:99996) + Q9SLA0 C353 (MOD:99997)
    assert sites == {("O03042", 284), ("O03042", 247), ("Q9SLA0", 353)}
    assert all(s.evidence_level == "site_ms" for s in result.sites)
    assert all(s.reference_source_sha256 for s in result.sites)


def test_multicysteine_peptide_marks_only_the_modified_cysteine() -> None:
    result = _normalized()
    q9sla0 = [s for s in result.sites if s.protein_accession_raw == "Q9SLA0"]

    # GSY…C(21)…C(33)…; only the MOD:99997 C21 → protein 353 is a site, not C33/365.
    assert len(q9sla0) == 1
    site = q9sla0[0]
    assert site.cys_position_in_peptide == 21
    assert site.cys_position_in_protein == 353
    assert site.peptide_sequence[site.cys_position_in_peptide - 1] == "C"
    assert 365 not in {s.cys_position_in_protein for s in q9sla0}


def test_isoform_multimapped_spectrum_goes_to_conflicts_not_a_site() -> None:
    result = _normalized()

    assert not any(
        s.protein_accession_raw.startswith("P27140") for s in result.sites
    )
    reasons = {c.reason for c in result.conflicts}
    assert any("isoform" in r or "multi" in r for r in reasons)


def test_controls_artifacts_and_decoys_never_produce_sites() -> None:
    result = _normalized()
    site_proteins = {s.protein_accession_raw for s in result.sites}

    assert "Q944G9" not in site_proteins  # MOD:99999 MSBT control only
    assert "Q9ZNZ7" not in site_proteins  # MOD:00110 artifact only
    assert not any(
        s.protein_accession_raw.startswith(("rev_", "O65387")) for s in result.sites
    )
    # A decoy row must not leak a site anywhere.
    assert all("rev_" not in s.protein_accession_raw for s in result.sites)


def test_fixture_provenance_is_intact() -> None:
    assert MANIFEST["biological_values_modified"] is False
    assert MANIFEST["accession"] == "PXD006140"
    assert hash_file(SUBSET_TSV, "sha256") == MANIFEST["fixture_sha256"]
