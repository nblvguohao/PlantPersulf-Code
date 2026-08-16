"""Streaming mzIdentML audit with exact reference resolution."""

from __future__ import annotations

import gzip
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.provenance.hashing import hash_file

EXPECTED_NAMESPACE = "http://psidev.info/psi/pi/mzIdentML/1.1"
COUNTED_ELEMENTS = (
    "DBSequence",
    "Peptide",
    "PeptideEvidence",
    "SpectrumIdentificationResult",
    "SpectrumIdentificationItem",
    "Modification",
)


@dataclass(frozen=True)
class _Modification:
    location: str
    residue: str
    monoisotopic_mass_delta: str
    cv_accession: str
    cv_name: str
    cv_value: str


@dataclass(frozen=True)
class _Peptide:
    sequence: str
    modifications: tuple[_Modification, ...]


@dataclass(frozen=True)
class _PeptideEvidence:
    peptide_ref: str
    db_sequence_ref: str
    start: str
    end: str
    is_decoy: str


@dataclass(frozen=True)
class MzidCandidate:
    spectrum_result_id: str
    spectrum_id: str
    spectrum_identification_item_id: str
    rank: str
    pass_threshold: str
    peptide_id: str
    peptide_sequence: str
    peptide_evidence_id: str
    db_sequence_id: str
    protein_accession: str
    modification_index: int
    location: str
    residue: str
    protein_position: str
    monoisotopic_mass_delta: str
    cv_accession: str
    cv_name: str
    cv_value: str
    evidence_class: str
    method_mapping_status: str
    conflict_status: str
    source_sha256: str
    source_locator: str


@dataclass(frozen=True)
class MzidConflict:
    conflict_type: str
    detail: str
    spectrum_result_id: str
    spectrum_identification_item_id: str
    peptide_id: str
    peptide_evidence_id: str
    db_sequence_id: str
    source_sha256: str
    source_locator: str


@dataclass(frozen=True)
class MzidAudit:
    namespace: str
    version: str
    structural_counts: dict[str, int]
    candidates: tuple[MzidCandidate, ...]
    conflicts: tuple[MzidConflict, ...]
    source_sha256: str


def _namespace(tag: str) -> str:
    if not tag.startswith("{") or "}" not in tag:
        return ""
    return tag[1 : tag.index("}")]


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def _required_attribute(element: ET.Element, name: str, context: str) -> str:
    value = element.get(name, "")
    if not value:
        raise RuntimeError(f"mzIdentML {context} lacks {name}")
    return value


def _integer(value: str, context: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise RuntimeError(f"mzIdentML has invalid {context}: {value}") from exc


def _append_conflict(
    conflicts: list[MzidConflict],
    conflict_type: str,
    detail: str,
    result_id: str,
    item_id: str,
    peptide_id: str,
    evidence_id: str,
    db_sequence_id: str,
    source_sha256: str,
) -> None:
    conflicts.append(
        MzidConflict(
            conflict_type=conflict_type,
            detail=detail,
            spectrum_result_id=result_id,
            spectrum_identification_item_id=item_id,
            peptide_id=peptide_id,
            peptide_evidence_id=evidence_id,
            db_sequence_id=db_sequence_id,
            source_sha256=source_sha256,
            source_locator=(
                f"peptides_1_1_0.mzid.gz:result={result_id}:item={item_id}"
            ),
        )
    )


def _peptide_from_element(element: ET.Element, qualified: str) -> _Peptide:
    sequence_element = element.find(f"{qualified}PeptideSequence")
    if sequence_element is None or not sequence_element.text:
        raise RuntimeError("mzIdentML Peptide lacks sequence")
    modifications: list[_Modification] = []
    for modification in element.findall(f"{qualified}Modification"):
        cv_param = modification.find(f"{qualified}cvParam")
        if cv_param is None:
            raise RuntimeError("mzIdentML Modification lacks cvParam")
        modifications.append(
            _Modification(
                location=_required_attribute(modification, "location", "Modification"),
                residue=modification.get("residues", ""),
                monoisotopic_mass_delta=_required_attribute(
                    modification,
                    "monoisotopicMassDelta",
                    "Modification",
                ),
                cv_accession=_required_attribute(cv_param, "accession", "cvParam"),
                cv_name=_required_attribute(cv_param, "name", "cvParam"),
                cv_value=cv_param.get("value", ""),
            )
        )
    return _Peptide(sequence_element.text, tuple(modifications))


def _validate_mapping(
    peptide: _Peptide,
    evidence: _PeptideEvidence,
    protein_sequence: str,
) -> tuple[str, ...]:
    problems: list[str] = []
    start = _integer(evidence.start, "PeptideEvidence start")
    end = _integer(evidence.end, "PeptideEvidence end")
    if start < 1 or end < start or end > len(protein_sequence):
        problems.append("peptide_evidence_range")
    elif protein_sequence[start - 1 : end] != peptide.sequence:
        problems.append("peptide_protein_sequence_mismatch")
    return tuple(problems)


def _validate_modification(
    modification: _Modification,
    peptide_sequence: str,
) -> tuple[str, ...]:
    problems: list[str] = []
    location = _integer(modification.location, "Modification location")
    if location < 0 or location > len(peptide_sequence) + 1:
        problems.append("modification_location_out_of_range")
        return tuple(problems)
    if 1 <= location <= len(peptide_sequence):
        observed = peptide_sequence[location - 1]
        if modification.residue and modification.residue != observed:
            problems.append("modification_residue_mismatch")
    elif modification.residue:
        problems.append("terminal_modification_has_residue")
    return tuple(problems)


def _process_result(
    element: ET.Element,
    qualified: str,
    peptides: dict[str, _Peptide],
    peptide_evidence: dict[str, _PeptideEvidence],
    db_sequences: dict[str, tuple[str, str]],
    candidates: list[MzidCandidate],
    conflicts: list[MzidConflict],
    source_sha256: str,
) -> None:
    result_id = _required_attribute(element, "id", "SpectrumIdentificationResult")
    spectrum_id = _required_attribute(
        element, "spectrumID", "SpectrumIdentificationResult"
    )
    for item in element.findall(f"{qualified}SpectrumIdentificationItem"):
        item_id = _required_attribute(item, "id", "SpectrumIdentificationItem")
        peptide_id = _required_attribute(
            item, "peptide_ref", "SpectrumIdentificationItem"
        )
        peptide = peptides.get(peptide_id)
        evidence_refs = [
            reference.get("peptideEvidence_ref", "")
            for reference in item.findall(f"{qualified}PeptideEvidenceRef")
        ]
        if peptide is None:
            _append_conflict(
                conflicts,
                "missing_peptide_reference",
                f"Peptide does not exist: {peptide_id}",
                result_id,
                item_id,
                peptide_id,
                "",
                "",
                source_sha256,
            )
            continue
        if not evidence_refs or any(not reference for reference in evidence_refs):
            _append_conflict(
                conflicts,
                "missing_peptide_evidence_reference",
                "SpectrumIdentificationItem has no complete PeptideEvidenceRef",
                result_id,
                item_id,
                peptide_id,
                "",
                "",
                source_sha256,
            )
            continue
        for evidence_id in evidence_refs:
            evidence = peptide_evidence.get(evidence_id)
            if evidence is None:
                _append_conflict(
                    conflicts,
                    "missing_peptide_evidence",
                    f"PeptideEvidence does not exist: {evidence_id}",
                    result_id,
                    item_id,
                    peptide_id,
                    evidence_id,
                    "",
                    source_sha256,
                )
                continue
            db_sequence_id = evidence.db_sequence_ref
            db_sequence = db_sequences.get(db_sequence_id)
            problems: list[str] = []
            if evidence.peptide_ref != peptide_id:
                problems.append("peptide_evidence_reference_mismatch")
            if db_sequence is None:
                problems.append("missing_db_sequence")
                protein_accession = ""
                protein_sequence = ""
            else:
                protein_accession, protein_sequence = db_sequence
                problems.extend(_validate_mapping(peptide, evidence, protein_sequence))
            for problem in problems:
                _append_conflict(
                    conflicts,
                    problem,
                    "exact mzIdentML reference or sequence mapping failed",
                    result_id,
                    item_id,
                    peptide_id,
                    evidence_id,
                    db_sequence_id,
                    source_sha256,
                )
            for index, modification in enumerate(peptide.modifications, start=1):
                modification_problems = list(
                    _validate_modification(modification, peptide.sequence)
                )
                for problem in modification_problems:
                    _append_conflict(
                        conflicts,
                        problem,
                        "localized modification does not match peptide sequence",
                        result_id,
                        item_id,
                        peptide_id,
                        evidence_id,
                        db_sequence_id,
                        source_sha256,
                    )
                location = _integer(modification.location, "Modification location")
                protein_position = ""
                if 1 <= location <= len(peptide.sequence):
                    protein_position = str(
                        _integer(evidence.start, "PeptideEvidence start") + location - 1
                    )
                locator = (
                    "peptides_1_1_0.mzid.gz:"
                    f"result={result_id}:item={item_id}:"
                    f"evidence={evidence_id}:modification={index}"
                )
                candidates.append(
                    MzidCandidate(
                        spectrum_result_id=result_id,
                        spectrum_id=spectrum_id,
                        spectrum_identification_item_id=item_id,
                        rank=_required_attribute(
                            item, "rank", "SpectrumIdentificationItem"
                        ),
                        pass_threshold=_required_attribute(
                            item,
                            "passThreshold",
                            "SpectrumIdentificationItem",
                        ),
                        peptide_id=peptide_id,
                        peptide_sequence=peptide.sequence,
                        peptide_evidence_id=evidence_id,
                        db_sequence_id=db_sequence_id,
                        protein_accession=protein_accession,
                        modification_index=index,
                        location=modification.location,
                        residue=modification.residue,
                        protein_position=protein_position,
                        monoisotopic_mass_delta=(modification.monoisotopic_mass_delta),
                        cv_accession=modification.cv_accession,
                        cv_name=modification.cv_name,
                        cv_value=modification.cv_value,
                        evidence_class="unresolved",
                        method_mapping_status="absent_v1",
                        conflict_status=(
                            "conflict" if problems or modification_problems else "clear"
                        ),
                        source_sha256=source_sha256,
                        source_locator=locator,
                    )
                )


def parse_mzidentml(path: Path, source_sha256: str) -> MzidAudit:
    """Stream registered mzIdentML and preserve exact identified modifications."""
    if hash_file(path, "sha256") != source_sha256:
        raise RuntimeError(f"mzIdentML source SHA256 mismatch: {path}")
    counts = {name: 0 for name in COUNTED_ELEMENTS}
    db_sequences: dict[str, tuple[str, str]] = {}
    peptides: dict[str, _Peptide] = {}
    peptide_evidence: dict[str, _PeptideEvidence] = {}
    candidates: list[MzidCandidate] = []
    conflicts: list[MzidConflict] = []
    namespace = ""
    version = ""
    with gzip.open(path, "rb") as handle:
        for event, element in ET.iterparse(handle, events=("start", "end")):
            if event == "start" and not namespace:
                namespace = _namespace(element.tag)
                version = element.get("version", "")
                if namespace != EXPECTED_NAMESPACE or version != "1.1.0":
                    raise RuntimeError("mzIdentML namespace or version mismatch")
                continue
            if event != "end":
                continue
            local_name = _local_name(element.tag)
            if local_name in counts:
                counts[local_name] += 1
            qualified = f"{{{namespace}}}"
            if local_name == "DBSequence":
                identifier = _required_attribute(element, "id", "DBSequence")
                sequence_element = element.find(f"{qualified}Seq")
                if sequence_element is None or not sequence_element.text:
                    raise RuntimeError("mzIdentML DBSequence lacks Seq")
                if identifier in db_sequences:
                    raise RuntimeError(f"duplicate DBSequence: {identifier}")
                db_sequences[identifier] = (
                    _required_attribute(element, "accession", "DBSequence"),
                    sequence_element.text,
                )
                element.clear()
            elif local_name == "Peptide":
                identifier = _required_attribute(element, "id", "Peptide")
                if identifier in peptides:
                    raise RuntimeError(f"duplicate Peptide: {identifier}")
                peptides[identifier] = _peptide_from_element(element, qualified)
                element.clear()
            elif local_name == "PeptideEvidence":
                identifier = _required_attribute(element, "id", "PeptideEvidence")
                if identifier in peptide_evidence:
                    raise RuntimeError(f"duplicate PeptideEvidence: {identifier}")
                peptide_evidence[identifier] = _PeptideEvidence(
                    peptide_ref=_required_attribute(
                        element, "peptide_ref", "PeptideEvidence"
                    ),
                    db_sequence_ref=_required_attribute(
                        element, "dBSequence_ref", "PeptideEvidence"
                    ),
                    start=_required_attribute(element, "start", "PeptideEvidence"),
                    end=_required_attribute(element, "end", "PeptideEvidence"),
                    is_decoy=_required_attribute(element, "isDecoy", "PeptideEvidence"),
                )
                element.clear()
            elif local_name == "SpectrumIdentificationResult":
                _process_result(
                    element,
                    qualified,
                    peptides,
                    peptide_evidence,
                    db_sequences,
                    candidates,
                    conflicts,
                    source_sha256,
                )
                element.clear()
    return MzidAudit(
        namespace=namespace,
        version=version,
        structural_counts=counts,
        candidates=tuple(candidates),
        conflicts=tuple(conflicts),
        source_sha256=source_sha256,
    )
