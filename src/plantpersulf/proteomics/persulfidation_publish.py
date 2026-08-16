"""Deterministic publish and audit of persulfidation site evidence.

Reads ``configs/persulfidation_sites_v1.yaml``, calls each study's registered
parser against hash-verified supplement sources and a SHA256-pinned UniProt
reference bundle, and writes atomically-replaceable TSV + manifest output under
``data/interim/<ACCESSION>/persulfidation_sites_v1/`` that the readiness v2 gate
can read. Re-running into an existing directory must produce byte-identical
output (or fail); ``audit_persulfidation_sites`` re-verifies the frozen output
against current source/reference hashes.
"""

from __future__ import annotations

import json
import re
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, cast

import yaml

from plantpersulf.proteomics.metadata import (
    ParseIssue,
    ReferenceSequence,
    SiteEvidence,
    SiteNormalization,
)
from plantpersulf.proteomics.persulfidation_sites import (
    MaxquantPersulfidationSites,
    parse_dataset_s3_sites,
    parse_maxquant_persulfidation_sites,
)
from plantpersulf.provenance.hashing import hash_file
from plantpersulf.provenance.supplementary import (
    SupplementarySource,
    audit_supplementary_sources,
)

SCHEMA_VERSION = 1
PARSER_VERSION = "plantpersulf/0.0.0"
PERSULF_SITE_FIELDS = (
    "study_accession",
    "sample_id",
    "source_file",
    "spectrum_id",
    "peptide_sequence",
    "modified_sequence",
    "protein_accession_raw",
    "protein_accession_canonical",
    "cys_position_in_peptide",
    "cys_position_in_protein",
    "modification_name_raw",
    "evidence_level",
    "quant_value",
    "quant_unit",
    "parser_version",
    "source_sha256",
    "reference_source_sha256",
)
ISSUE_FIELDS = (
    "source_file",
    "spectrum_id",
    "protein_accession_raw",
    "modified_sequence",
    "reason",
    "detail",
)
OUTPUT_NAMES = (
    "sites.tsv",
    "low_confidence.tsv",
    "conflicts.tsv",
    "excluded.tsv",
    "manifest.json",
)


@dataclass(frozen=True)
class PersulfidationPublishSummary:
    site_count: int = 0
    distinct_site_count: int = 0
    low_confidence_count: int = 0
    conflict_count: int = 0
    excluded_count: int = 0
    missing_sequence_count: int = 0


def _write_tsv(
    path: Path,
    fieldnames: tuple[str, ...],
    records: tuple[SiteEvidence | ParseIssue, ...],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        import csv

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for record in records:
            writer.writerow(asdict(record))


def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hash_file(path, "sha256")
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _load_config(config_path: Path) -> dict[str, Any]:
    try:
        loaded: object = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    except (OSError, yaml.YAMLError) as exc:
        raise RuntimeError(
            f"invalid persulfidation site config: {config_path}"
        ) from exc
    if not isinstance(loaded, dict) or loaded.get("version") != 1:
        raise RuntimeError("persulfidation site config requires version 1")
    return cast(dict[str, Any], loaded)


def _load_reference_bundle(bundle_info: dict[str, Any]) -> dict[str, ReferenceSequence]:
    path = Path(str(bundle_info["path"]))
    expected_sha = str(bundle_info["sha256"])
    actual = hash_file(path, "sha256")
    if actual != expected_sha:
        raise RuntimeError(f"reference bundle SHA256 mismatch: {path}")

    references: dict[str, ReferenceSequence] = {}
    cur_header = ""
    cur_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur_header:
                acc = cur_header.strip().split("|")[1]
                references[acc] = ReferenceSequence(
                    accession=acc,
                    sequence="".join(cur_lines),
                    sequence_version=_sequence_version(cur_header),
                    source_file=path,
                    source_sha256=actual,
                )
            cur_header = line
            cur_lines = []
        elif line:
            cur_lines.append(line)
    if cur_header:
        acc = cur_header.strip().split("|")[1]
        references[acc] = ReferenceSequence(
            accession=acc,
            sequence="".join(cur_lines),
            sequence_version=_sequence_version(cur_header),
            source_file=path,
            source_sha256=actual,
        )
    if len(references) != bundle_info.get("record_count", 0):
        raise RuntimeError("reference bundle record count mismatch")
    return references


def _sequence_version(header: str) -> str:
    match = re.search(r"SV=(\d+)", header)
    return match.group(1) if match else "1"


def _resolve_source(
    study: str,
    member: str,
    supplement_sources: tuple[SupplementarySource, ...],
) -> SupplementarySource:
    for source in supplement_sources:
        if source.study_accession == study and source.member == member:
            return source
    raise RuntimeError(f"supplementary source not found: {study} / {member}")


def _publish_006140(
    study_cfg: dict[str, Any],
    references: dict[str, ReferenceSequence],
    supplement_sources: tuple[SupplementarySource, ...],
) -> tuple[SiteNormalization, None]:
    member = str(study_cfg["source_member"])
    _resolve_source("PXD006140", member, supplement_sources)
    ident_path = Path(str(study_cfg["ident_peptides_path"]))
    expected = str(study_cfg["ident_peptides_sha256"])
    if hash_file(ident_path, "sha256") != expected:
        raise RuntimeError(f"ident_peptides SHA256 mismatch: {ident_path}")
    return parse_dataset_s3_sites(
        ident_path,
        references,
        expected,
    ), None


def _publish_024061(
    study_cfg: dict[str, Any],
    references: dict[str, ReferenceSequence],
    supplement_sources: tuple[SupplementarySource, ...],
) -> tuple[None, MaxquantPersulfidationSites]:
    all_sites: list[SiteEvidence] = []
    all_low: list[SiteEvidence] = []
    all_conflicts: list[ParseIssue] = []
    all_excluded: list[ParseIssue] = []
    for member_cfg in study_cfg["members"]:
        member = str(member_cfg["source_member"])
        mod_name = str(member_cfg["modification_name"])
        source = _resolve_source("PXD024061", member, supplement_sources)
        result = parse_maxquant_persulfidation_sites(
            source.local_path,
            references,
            source.sha256,
            mod_name,
        )
        all_sites.extend(result.sites)
        all_low.extend(result.low_confidence)
        all_conflicts.extend(result.conflicts)
        all_excluded.extend(result.excluded)
    return None, MaxquantPersulfidationSites(
        tuple(all_sites),
        tuple(all_low),
        tuple(all_conflicts),
        tuple(all_excluded),
    )


def publish_persulfidation_sites(
    study_accession: str,
    config_path: Path = Path("configs/persulfidation_sites_v1.yaml"),
    registry_dir: Path = Path("data/registry"),
    output_root: Path = Path("data/interim"),
) -> PersulfidationPublishSummary:
    config = _load_config(config_path)
    studies_cfg = config.get("studies")
    if not isinstance(studies_cfg, dict):
        raise RuntimeError("persulfidation site config requires a studies mapping")
    study_accession = study_accession.strip().upper()
    study_cfg = studies_cfg.get(study_accession)
    if not isinstance(study_cfg, dict):
        raise RuntimeError(f"study not registered: {study_accession}")

    bundle_info = config.get("reference_bundle")
    if not isinstance(bundle_info, dict):
        raise RuntimeError("persulfidation site config requires reference_bundle")
    references = _load_reference_bundle(bundle_info)

    supplement_sources = audit_supplementary_sources(
        registry_path=registry_dir / "supplementary_sources.tsv",
    )

    parser_label = str(study_cfg["parser"])
    if parser_label == "dataset_s3":
        s3_result, _ = _publish_006140(study_cfg, references, supplement_sources)
        sites = s3_result.sites
        conflicts = s3_result.conflicts
        low_confidence: tuple[SiteEvidence, ...] = ()
        excluded: tuple[ParseIssue, ...] = tuple(
            issue
            for issue in conflicts
            if "decoy" in issue.reason or "contaminant" in issue.reason.lower()
        )
        conflicts = tuple(issue for issue in conflicts if issue not in excluded)
    elif parser_label == "maxquant_sites":
        _, mq_result = _publish_024061(study_cfg, references, supplement_sources)
        sites = mq_result.sites
        low_confidence = mq_result.low_confidence
        conflicts = mq_result.conflicts
        excluded = mq_result.excluded
    else:
        raise RuntimeError(f"unknown persulfidation site parser: {parser_label}")

    summary = PersulfidationPublishSummary(
        site_count=len(sites),
        distinct_site_count=len(
            {(s.protein_accession_raw, s.cys_position_in_protein) for s in sites}
        ),
        low_confidence_count=len(low_confidence),
        conflict_count=len(conflicts),
        excluded_count=len(excluded),
        missing_sequence_count=sum(
            issue.reason in {"sequence_unavailable", "isoform_sequence_unavailable"}
            for issue in conflicts
        ),
    )

    output_directory = output_root / study_accession / "persulfidation_sites_v1"
    output_directory.parent.mkdir(parents=True, exist_ok=True)
    temporary_root = Path(
        tempfile.mkdtemp(
            prefix=f".{output_directory.name}-",
            dir=output_directory.parent,
        )
    )
    staged = temporary_root / output_directory.name
    staged.mkdir()
    try:
        _write_tsv(staged / "sites.tsv", PERSULF_SITE_FIELDS, sites)
        _write_tsv(staged / "low_confidence.tsv", PERSULF_SITE_FIELDS, low_confidence)
        _write_tsv(staged / "conflicts.tsv", ISSUE_FIELDS, conflicts)
        _write_tsv(staged / "excluded.tsv", ISSUE_FIELDS, excluded)

        output_hashes = [
            {"file": name, "sha256": hash_file(staged / name, "sha256")}
            for name in OUTPUT_NAMES
            if name != "manifest.json"
        ]
        manifest: dict[str, object] = {
            "schema_version": SCHEMA_VERSION,
            "study_accession": study_accession,
            "parser_version": PARSER_VERSION,
            "parser_label": parser_label,
            "reference_bundle_sha256": bundle_info["sha256"],
            "counts": asdict(summary),
            "outputs": output_hashes,
            "labels_created": False,
            "nondetection_labeled_negative": False,
            "biological_values_modified": False,
        }
        (staged / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        if output_directory.exists():
            current = _tree_hashes(output_directory)
            proposed = _tree_hashes(staged)
            if proposed != current:
                raise RuntimeError(
                    f"existing persulfidation site output differs: {output_directory}"
                )
            return summary
        staged.replace(output_directory)
        return summary
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def _manifest_mapping(value: object, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise RuntimeError(f"site manifest requires mapping: {context}")
    return cast(dict[str, Any], value)


def _manifest_items(value: object, context: str) -> list[Any]:
    if not isinstance(value, list):
        raise RuntimeError(f"site manifest requires list: {context}")
    return value


def _read_tsv(path: Path, fields: tuple[str, ...]) -> list[dict[str, str]]:
    import csv

    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != fields:
            raise RuntimeError(
                f"persulfidation site output has invalid columns: {path}"
            )
        rows = [dict(row) for row in reader]
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise RuntimeError(f"persulfidation site output has malformed row: {path}")
    return cast(list[dict[str, str]], rows)


def audit_persulfidation_sites(
    study_accession: str,
    output_directory: Path,
    reference_bundle_sha256: str | None = None,
) -> PersulfidationPublishSummary:
    """Re-verify frozen persulfidation site output."""
    observed = {path.name for path in output_directory.iterdir() if path.is_file()}
    if observed != set(OUTPUT_NAMES):
        raise RuntimeError("persulfidation site output file set mismatch")
    manifest_path = output_directory / "manifest.json"
    try:
        parsed: object = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(
            f"invalid persulfidation site manifest: {manifest_path}"
        ) from exc
    manifest = _manifest_mapping(parsed, "root")
    if (
        manifest.get("schema_version") != SCHEMA_VERSION
        or manifest.get("study_accession") != study_accession
    ):
        raise RuntimeError("persulfidation site manifest identity mismatch")
    if (
        manifest.get("labels_created") is not False
        or manifest.get("nondetection_labeled_negative") is not False
        or manifest.get("biological_values_modified") is not False
    ):
        raise RuntimeError(
            "persulfidation site manifest violates scientific integrity policy"
        )

    if reference_bundle_sha256 is not None and (
        manifest.get("reference_bundle_sha256") != reference_bundle_sha256
    ):
        raise RuntimeError("persulfidation site reference bundle sha256 mismatch")

    declared: set[str] = set()
    for raw_entry in _manifest_items(manifest.get("outputs"), "outputs"):
        entry = _manifest_mapping(raw_entry, "outputs[]")
        name = str(entry.get("file", ""))
        declared.add(name)
        path = output_directory / name
        if not path.is_file() or hash_file(path, "sha256") != entry.get("sha256"):
            raise RuntimeError(f"persulfidation site output SHA256 mismatch: {path}")
    if declared != set(OUTPUT_NAMES) - {"manifest.json"}:
        raise RuntimeError("persulfidation site output list mismatch")

    sites = _read_tsv(output_directory / "sites.tsv", PERSULF_SITE_FIELDS)
    low = _read_tsv(output_directory / "low_confidence.tsv", PERSULF_SITE_FIELDS)
    conflicts = _read_tsv(output_directory / "conflicts.tsv", ISSUE_FIELDS)
    excluded = _read_tsv(output_directory / "excluded.tsv", ISSUE_FIELDS)
    summary = PersulfidationPublishSummary(
        site_count=len(sites),
        distinct_site_count=len(
            {
                (row["protein_accession_raw"], int(row["cys_position_in_protein"]))
                for row in sites
            }
        ),
        low_confidence_count=len(low),
        conflict_count=len(conflicts),
        excluded_count=len(excluded),
        missing_sequence_count=sum(
            row["reason"] in {"sequence_unavailable", "isoform_sequence_unavailable"}
            for row in conflicts
        ),
    )
    if _manifest_mapping(manifest.get("counts"), "counts") != {
        key: value for key, value in asdict(summary).items()
    }:
        raise RuntimeError("persulfidation site manifest counts mismatch")
    if any(
        row["evidence_level"] != "site_ms"
        or "negative" in row["evidence_level"].lower()
        for row in (sites + low)
    ):
        raise RuntimeError(
            "persulfidation site output contains forbidden evidence level"
        )
    return summary
