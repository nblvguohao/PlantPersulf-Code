"""RED (Phase A2, study #2): PXD024061 MaxQuant Sites → residue-resolved site_ms.

Second Gate-1 site-level study from
``docs/superpowers/specs/2026-07-21-evidence-site-acquisition-design.md`` §5.6/§9.

PXD024061's deposited MaxQuant search output contains dedicated PTM-site tables
``Sulfide(C)Sites.txt`` and ``CianoBiotin(C)Sites.txt`` — the same lab's
tag-switch persulfidation chemistry as PXD006140. This test pins the extraction
contract for the MaxQuant ``…Sites.txt`` schema:

- keep only ``Amino acid == C``;
- exclude ``Reverse == "+"`` and ``Potential contaminant == "+"`` rows;
- apply the MaxQuant class-I threshold ``Localization prob >= 0.75`` — sites below
  it go to a flagged lower-confidence tier, never silently dropped;
- a razor site mapping to multiple proteins becomes a conflict, never first-picked;
- ``Position`` in the leading protein is verified to be C against the registered
  UniProt sequence, and the ``Sequence window`` must match that neighbourhood;
- emit ``evidence_level=site_ms`` for validated class-I sites.

Driven by provenance-locked real micro-fixtures cut verbatim from the two site
tables (range-extracted from the 2.42 GB archive; biological values unchanged).

Expected RED: ``parse_maxquant_persulfidation_sites`` / ``MAXQUANT_CLASS_I_THRESHOLD``
do not exist yet, so collection fails with ImportError.
"""

from __future__ import annotations

import json
from pathlib import Path

from plantpersulf.proteomics.metadata import ReferenceSequence
from plantpersulf.proteomics.persulfidation_sites import (  # RED: names missing
    MAXQUANT_CLASS_I_THRESHOLD,
    parse_maxquant_persulfidation_sites,
)
from plantpersulf.provenance.hashing import hash_file

FIXTURE = Path("tests/fixtures/real/PXD024061_maxquant_sites")
SULFIDE_TSV = FIXTURE / "sulfide_sites_subset.tsv"
CIANOBIOTIN_TSV = FIXTURE / "cianobiotin_sites_subset.tsv"
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


def _source_sha(modification_name: str) -> str:
    return next(
        s["member_sha256"]
        for s in MANIFEST["sources"]
        if s["modification_name"] == modification_name
    )


def _sulfide():  # type: ignore[no-untyped-def]
    return parse_maxquant_persulfidation_sites(
        source_tsv=SULFIDE_TSV,
        references=_references(),
        source_sha256=_source_sha("Sulfide"),
        modification_name="Sulfide",
    )


def _cianobiotin():  # type: ignore[no-untyped-def]
    return parse_maxquant_persulfidation_sites(
        source_tsv=CIANOBIOTIN_TSV,
        references=_references(),
        source_sha256=_source_sha("CN-Biotin-Sulfide"),
        modification_name="CN-Biotin-Sulfide",
    )


def _sites(result):  # type: ignore[no-untyped-def]
    return {(s.protein_accession_raw, s.cys_position_in_protein) for s in result.sites}


def test_class_i_threshold_is_maxquant_default() -> None:
    assert MAXQUANT_CLASS_I_THRESHOLD == 0.75


def test_only_class_i_single_protein_cysteines_become_sites() -> None:
    result = _sulfide()

    assert _sites(result) == {("A0A1P8B1I9", 2424), ("A0A1P8APX1", 401)}
    assert all(s.evidence_level == "site_ms" for s in result.sites)
    assert all(s.modification_name_raw == "Sulfide" for s in result.sites)
    assert all(s.reference_source_sha256 for s in result.sites)


def test_low_localization_sites_are_flagged_not_dropped_or_promoted() -> None:
    result = _sulfide()
    low = {
        (s.protein_accession_raw, s.cys_position_in_protein)
        for s in result.low_confidence
    }

    # prob 0.5 site is retained in the low-confidence tier, never a class-I site.
    assert ("A0A1I9LQH7", 1019) in low
    assert ("A0A1I9LQH7", 1019) not in _sites(result)
    # the just-above-threshold prob 0.8247 site IS class-I, not low-confidence.
    assert ("A0A1P8APX1", 401) not in low


def test_razor_multiprotein_site_is_a_conflict_not_a_site() -> None:
    result = _sulfide()

    assert not any(";" in s.protein_accession_raw for s in result.sites)
    assert not any(";" in s.protein_accession_raw for s in result.low_confidence)
    reasons = {c.reason for c in result.conflicts}
    assert any("razor" in r or "multi" in r for r in reasons)


def test_contaminants_are_excluded_and_never_become_sites() -> None:
    result = _sulfide()

    emitted = list(result.sites) + list(result.low_confidence)
    assert not any(s.protein_accession_raw.startswith("CON__") for s in emitted)
    assert any("contaminant" in issue.reason for issue in result.excluded)


def test_cianobiotin_sites_carry_the_labeled_persulfidation_modification() -> None:
    result = _cianobiotin()

    assert _sites(result) == {("Q8VZF1", 290)}
    assert all(s.modification_name_raw == "CN-Biotin-Sulfide" for s in result.sites)


def test_fixture_provenance_is_intact() -> None:
    assert MANIFEST["biological_values_modified"] is False
    assert MANIFEST["accession"] == "PXD024061"
    assert hash_file(SULFIDE_TSV, "sha256") == _source_subset_sha("Sulfide")
    assert hash_file(CIANOBIOTIN_TSV, "sha256") == _source_subset_sha(
        "CN-Biotin-Sulfide"
    )


def _source_subset_sha(modification_name: str) -> str:
    return next(
        s["subset_sha256"]
        for s in MANIFEST["sources"]
        if s["modification_name"] == modification_name
    )
