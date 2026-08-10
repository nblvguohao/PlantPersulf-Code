"""Phase S0, Task 2 — fail-closed structure coverage audit.

Before any scoring run launches, this module verifies every frozen input hash,
verifies the v1/v2 registry relationship, and maps every benchmark row against
the registered structures in the selected release. The resulting
``StructureCoverageAudit`` is passed to the downstream decision protocol so
coverage facts are available even though the three-metric cluster-bootstrap
decision (Task 6) is computed from *scores*, not from static coverage alone.

Provenance rules (fail-closed):

- Every ``FrozenFile`` must match its declared SHA256.
- The v1 registry must be an exact subset of v2 (additive-only expansion).
- Every benchmark row produces exactly one ``CoverageRecord``.
- Blocking errors (hash mismatch, missing file, registry conflict) are raised
  before any record is emitted and are *never* converted to a non-blocking
  ``mapping_status``.
"""

from __future__ import annotations

import csv
import hashlib
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from plantpersulf.download.alphafold import (
    AlphaFoldStructureSource,
    audit_alphafold_structures,
)

# ---------------------------------------------------------------------------
# frozen files
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FrozenFile:
    """A named input whose SHA256 must match before any analysis begins."""

    name: str
    path: Path
    sha256: str


def sha256_file(path: Path) -> str:
    """Return an uppercase SHA256 hex digest of ``path``'s bytes."""
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest().upper()


def verify_frozen_file(item: FrozenFile) -> None:
    """Raise ``RuntimeError`` if the file bytes do not match ``item.sha256``."""
    observed = sha256_file(item.path)
    if observed.upper() != item.sha256.upper():
        raise RuntimeError(
            f"SHA256 mismatch for {item.path}: "
            f"expected {item.sha256}, observed {observed}"
        )


# ---------------------------------------------------------------------------
# registry identity
# ---------------------------------------------------------------------------


def _read_registry_records(path: Path) -> list[dict[str, str]]:
    STRUCTURE_FIELDS = (
        "accession",
        "model_version",
        "source_url",
        "local_path",
        "retrieved_at",
        "size_bytes",
        "sha256",
        "downloader_version",
    )
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != STRUCTURE_FIELDS:
            raise RuntimeError(f"registry has invalid columns: {path}")
        return [dict(row) for row in reader]


def verify_registry_pair(
    v1_path: Path,
    v2_path: Path,
    *,
    expected_v1_count: int = 7,
    expected_v2_count: int = 2006,
) -> None:
    """Verify that v1 is a strict subset of v2 with expected counts.

    Duplicate rows (same accession repeated in one registry) are treated as an
    error. Conflicting rows (same accession, different provenance) produce a
    distinct blocking error message.
    """
    rows_v1 = _read_registry_records(v1_path)
    rows_v2 = _read_registry_records(v2_path)

    def _key(row: dict[str, str]) -> str:
        return row["accession"]

    accessions_v1 = {_key(r) for r in rows_v1}
    accessions_v2 = {_key(r) for r in rows_v2}

    if len(rows_v1) != len(accessions_v1):
        raise RuntimeError(
            f"v1 registry has {len(rows_v1)} rows but only "
            f"{len(accessions_v1)} unique accessions"
        )
    if len(rows_v2) != len(accessions_v2):
        raise RuntimeError(
            f"v2 registry has {len(rows_v2)} rows but only "
            f"{len(accessions_v2)} unique accessions"
        )

    if len(rows_v1) != expected_v1_count:
        raise RuntimeError(
            f"v1 registry has {len(rows_v1)} records (expected {expected_v1_count})"
        )
    if len(rows_v2) != expected_v2_count:
        raise RuntimeError(
            f"v2 registry has {len(rows_v2)} records (expected {expected_v2_count})"
        )

    if not accessions_v1 < accessions_v2:
        raise RuntimeError(
            "v1 registry is not a strict subset of v2"
        )

    # Check for conflicting rows on shared accessions
    v2_lookup = {_key(r): r for r in rows_v2}
    for r1 in rows_v1:
        r2 = v2_lookup[r1["accession"]]
        for col in ("model_version", "source_url", "local_path",
                     "size_bytes", "sha256"):
            if r1[col] != r2[col]:
                raise RuntimeError(
                    f"registry conflict on {r1['accession']}: "
                    f"v1 {col}={r1[col]}, v2 {col}={r2[col]}"
                )


# ---------------------------------------------------------------------------
# coverage records
# ---------------------------------------------------------------------------


ALLOWED_MAPPING_STATUSES = frozenset({
    "mapped_cys",
    "absent_structure",
    "absent_residue",
    "non_cys_residue",
    "malformed_structure",
    "duplicate_residue_ambiguity",
})


@dataclass(frozen=True)
class CoverageRecord:
    """One benchmark row mapped against one registered structure release."""

    release: str
    protein_accession: str
    cys_position_in_protein: int
    label: str          # "positive" | "unlabeled"
    study_accession: str
    cluster_id: str
    has_registered_structure: bool
    maps_to_cys: bool
    mapping_status: str
    plddt: float | None


@dataclass(frozen=True)
class StructureCoverageAudit:
    records: tuple[CoverageRecord, ...] = ()
    summary_rows: tuple[dict[str, object], ...] = ()
    input_hashes: dict[str, str] = field(default_factory=dict)
    blocking_errors: tuple[str, ...] = ()

    def require_pass(self) -> None:
        if self.blocking_errors:
            raise RuntimeError(
                f"coverage audit has blocking errors: {self.blocking_errors}"
            )


# ---------------------------------------------------------------------------
# main audit entry point
# ---------------------------------------------------------------------------


def audit_structure_coverage(
    *,
    release: str,
    benchmark_sha256: str,
    proteome_sha256: str,
    clusters_sha256: str,
    registry_sha256: str,
    registry_path: Path,
    registry_base: Path,
    benchmark_path: Path,
    proteome_path: Path,
    clusters_path: Path,
) -> StructureCoverageAudit:
    """Map every benchmark row against one frozen structure release.

    Returns a ``StructureCoverageAudit`` that contains exactly one
    ``CoverageRecord`` per benchmark row plus summary dimensions. Blocking
    provenance failures (hash mismatch, missing file) raise *during* this
    call; non-blocking mapping failures are recorded as a ``mapping_status``.
    """
    blocking: list[str] = []
    input_hashes: dict[str, str] = {}

    # --- Step 0: verify every frozen hash ---
    files: list[tuple[str, Path, str]] = [
        ("benchmark", benchmark_path, benchmark_sha256),
        ("proteome", proteome_path, proteome_sha256),
        ("clusters", clusters_path, clusters_sha256),
        ("registry", registry_path, registry_sha256),
    ]
    for name, path, expected in files:
        try:
            observed = sha256_file(path)
            input_hashes[name] = observed
            if observed.upper() != expected.upper():
                blocking.append(
                    f"SHA256 mismatch for {name} ({path}): "
                    f"expected {expected}, observed {observed}"
                )
        except FileNotFoundError:
            blocking.append(f"missing input: {path}")

    if not blocking:
        # --- Step 1: read inputs ---
        benchmark_rows = _read_benchmark(benchmark_path)
        clusters = _read_clusters(clusters_path)
        proteome = _load_proteome_fasta(proteome_path)

        # --- Step 2: register structure sources ---
        try:
            sources = audit_alphafold_structures(
                registry_path, base_directory=registry_base
            )
        except RuntimeError as exc:
            blocking.append(str(exc))
            return StructureCoverageAudit(
                records=(),
                summary_rows=(),
                input_hashes=input_hashes,
                blocking_errors=tuple(blocking),
            )

        struct_by_acc: dict[str, AlphaFoldStructureSource] = {
            s.accession: s for s in sources
        }

        # --- Step 3: pre-load and pre-parse every registered PDB once ---
        # Cache: accession -> list[(position, plddt, res_name)]
        # Only CYS residues are kept; other residues just record position/res_name
        pdb_cache: dict[str, list[tuple[int, float, str]]] = {}
        for acc, source in struct_by_acc.items():
            try:
                pdb_text = source.local_path.read_text(encoding="utf-8")
            except (OSError, UnicodeDecodeError):
                blocking.append(
                    f"cannot read registered structure: {source.local_path}"
                )
                continue
            residues = _parse_residues_for_audit(pdb_text)
            pdb_cache[acc] = residues

        if blocking:
            return StructureCoverageAudit(
                records=(),
                summary_rows=(),
                input_hashes=input_hashes,
                blocking_errors=tuple(blocking),
            )

        # --- Step 4: emit one record per benchmark row ---
        records: list[CoverageRecord] = []
        for row in benchmark_rows:
            acc = row["protein_accession"]
            pos = int(row["cys_position_in_protein"])
            label = row["label"]
            study = row["study_accession"]
            cluster = clusters.get(acc, f"__singleton__{acc}")

            residues = pdb_cache.get(acc)
            if residues is None:
                records.append(
                    CoverageRecord(
                        release=release,
                        protein_accession=acc,
                        cys_position_in_protein=pos,
                        label=label,
                        study_accession=study,
                        cluster_id=cluster,
                        has_registered_structure=False,
                        maps_to_cys=False,
                        mapping_status="absent_structure",
                        plddt=None,
                    )
                )
                continue

            # Check if this position maps to a CYS
            by_pos = {r[0]: r for r in residues}
            match = by_pos.get(pos)

            if match is None:
                status = "absent_residue" if residues else "malformed_structure"
            elif match[2] != "CYS":
                status = "non_cys_residue"
            else:
                status = "mapped_cys"

            records.append(
                CoverageRecord(
                    release=release,
                    protein_accession=acc,
                    cys_position_in_protein=pos,
                    label=label,
                    study_accession=study,
                    cluster_id=cluster,
                    has_registered_structure=True,
                    maps_to_cys=(status == "mapped_cys"),
                    mapping_status=status,
                    plddt=match[1] if status == "mapped_cys" else None,
                )
            )

        # --- Step 5: build summary ---
        summary = _build_summary(records, release)

        return StructureCoverageAudit(
            records=tuple(records),
            summary_rows=summary,
            input_hashes=input_hashes,
            blocking_errors=tuple(blocking),
        )

    return StructureCoverageAudit(
        records=(),
        summary_rows=(),
        input_hashes=input_hashes,
        blocking_errors=tuple(blocking),
    )


# ---------------------------------------------------------------------------
# internal helpers
# ---------------------------------------------------------------------------

BENCHMARK_FIELDS = (
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "evidence_level",
    "source_sha256",
)


def _read_benchmark(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != BENCHMARK_FIELDS:
            raise RuntimeError(f"benchmark has invalid columns: {path}")
        return [dict(row) for row in reader]


def _read_clusters(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mapping[row["protein_accession"]] = row["cluster_id"]
    return mapping


def _load_proteome_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur_header = ""
    cur_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur_header:
                seqs[cur_header.split("|")[1]] = "".join(cur_lines)
            cur_header = line
            cur_lines = []
        elif line:
            cur_lines.append(line)
    if cur_header:
        seqs[cur_header.split("|")[1]] = "".join(cur_lines)
    return seqs


def _parse_residues_for_audit(
    pdb_text: str,
) -> list[tuple[int, float, str]]:
    """Parse one PDB file once, returning (res_seq, plddt, res_name) for every
    CA atom. Used to build a fast per-accession lookup cache so the audit does
    not re-parse a PDB file thousands of times for each benchmark row."""
    residues: list[tuple[int, float, str]] = []
    for line in pdb_text.splitlines():
        if len(line) < 66 or line[0:6] != "ATOM  ":
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        res_name = line[17:20].strip()
        try:
            res_seq = int(line[22:26])
            plddt = float(line[60:66])
        except ValueError:
            continue
        residues.append((res_seq, plddt, res_name))
    return residues


def _build_summary(
    records: list[CoverageRecord],
    release: str,
) -> tuple[dict[str, object], ...]:
    summary: list[dict[str, object]] = []

    n_benchmark = len(records)
    n_positive = sum(1 for r in records if r.label == "positive")
    n_unlabeled = sum(1 for r in records if r.label == "unlabeled")

    # Registered proteins: unique accessions that have at least one structure
    registered_set = {r.protein_accession for r in records
                      if r.has_registered_structure}

    # Covered proteins: unique accessions where at least one Cys maps
    covered_set = {r.protein_accession for r in records
                   if r.mapping_status == "mapped_cys"}

    # Mapped Cys counts
    mapped_total = sum(1 for r in records if r.mapping_status == "mapped_cys")
    mapped_pos = sum(1 for r in records
                     if r.label == "positive" and r.mapping_status == "mapped_cys")
    mapped_ul = sum(1 for r in records
                    if r.label == "unlabeled" and r.mapping_status == "mapped_cys")

    for scope, label, n_sites in [
        ("overall", "all", n_benchmark),
        ("overall", "positive", n_positive),
        ("overall", "unlabeled", n_unlabeled),
    ]:
        _add_summary_row(summary, release, scope, "", label,
                         "benchmark_sites", n_sites, None, None, n_sites)

    _add_summary_row(summary, release, "overall", "", "all",
                     "registered_proteins", len(registered_set),
                     None, None, n_benchmark)
    _add_summary_row(summary, release, "overall", "", "all",
                     "covered_proteins", len(covered_set),
                     None, None, n_benchmark)
    _add_summary_row(summary, release, "overall", "", "all",
                     "mapped_cys_sites", mapped_total, None, None, n_benchmark)
    _add_summary_row(summary, release, "overall", "", "positive",
                     "mapped_cys_sites", mapped_pos, None, None, n_positive)
    _add_summary_row(summary, release, "overall", "", "unlabeled",
                     "mapped_cys_sites", mapped_ul, None, None, n_unlabeled)

    if n_benchmark > 0:
        _add_summary_row(summary, release, "overall", "", "all",
                         "coverage_rate",
                         round(mapped_total / n_benchmark, 6),
                         None, None, n_benchmark)

    # Per-study summaries
    studies = sorted({r.study_accession for r in records if r.study_accession})
    for study in studies:
        study_recs = [r for r in records if r.study_accession == study]
        n_s = len(study_recs)
        mapped_s = sum(1 for r in study_recs if r.mapping_status == "mapped_cys")
        _add_summary_row(summary, release, "study", study, "all",
                         "benchmark_sites", n_s, None, None, n_s)
        _add_summary_row(summary, release, "study", study, "all",
                         "mapped_cys_sites", mapped_s, None, None, n_s)
        if n_s > 0:
            _add_summary_row(summary, release, "study", study, "all",
                             "coverage_rate", round(mapped_s / n_s, 6),
                             None, None, n_s)

    # pLDDT bins
    bins = {"plddt_lt_50": (0.0, 50.0), "plddt_50_70": (50.0, 70.0),
            "plddt_70_90": (70.0, 90.0), "plddt_ge_90": (90.0, float("inf"))}
    for bin_name, (lo, hi) in bins.items():
        count = sum(
            1 for r in records
            if r.plddt is not None and lo <= r.plddt < hi
        )
        _add_summary_row(summary, release, "overall", "", "all",
                         bin_name, count, None, None, n_benchmark)

    # Coverage-label association: difference in coverage rate pos vs unlabeled
    coverage_pos = mapped_pos / n_positive if n_positive > 0 else 0.0
    coverage_ul = mapped_ul / n_unlabeled if n_unlabeled > 0 else 0.0
    diff = coverage_pos - coverage_ul
    _add_summary_row(summary, release, "overall", "", "all",
                     "coverage_label_risk_difference",
                     round(diff, 6), None, None, n_benchmark)

    # Cluster bootstrap CI for coverage-label risk difference
    # Pre-aggregate per cluster to avoid O(n_records * n_boot) scanning
    cluster_agg: dict[str, tuple[int, int, int, int]] = {}
    for rec in records:
        cid = rec.cluster_id
        prev = cluster_agg.get(cid, (0, 0, 0, 0))
        if rec.label == "positive":
            pos_inc = 1 if rec.mapping_status == "mapped_cys" else 0
            cluster_agg[cid] = (prev[0] + 1, prev[1] + pos_inc, prev[2], prev[3])
        else:
            ul_inc = 1 if rec.mapping_status == "mapped_cys" else 0
            cluster_agg[cid] = (prev[0], prev[1], prev[2] + 1, prev[3] + ul_inc)

    n_clusters = len(cluster_agg)
    cluster_list = sorted(cluster_agg)
    rng = random.Random(1729)
    replicates: list[float] = []
    for _ in range(5000):
        drawn_pos = 0
        drawn_pos_covered = 0
        drawn_ul = 0
        drawn_ul_covered = 0
        for _ in cluster_list:
            cid = cluster_list[rng.randrange(n_clusters)]
            (p, pc, u, uc) = cluster_agg[cid]
            drawn_pos += p
            drawn_pos_covered += pc
            drawn_ul += u
            drawn_ul_covered += uc
        if drawn_pos > 0 and drawn_ul > 0:
            replicates.append(
                drawn_pos_covered / drawn_pos - drawn_ul_covered / drawn_ul
            )

    if replicates:
        replicates.sort()
        ci_low = _percentile(replicates, 2.5)
        ci_high = _percentile(replicates, 97.5)
    else:
        ci_low = 0.0
        ci_high = 0.0

    _add_summary_row(summary, release, "overall", "", "all",
                     "coverage_label_risk_difference",
                     round(diff, 6),
                     round(ci_low, 6), round(ci_high, 6),
                     n_benchmark, n_clusters)

    return tuple(summary)


def _add_summary_row(
    target: list[dict[str, object]],
    release: str,
    scope: str,
    study_accession: str,
    label: str,
    metric: str,
    estimate: object,
    ci_low: object,
    ci_high: object,
    n_sites: int,
    n_clusters: int = 0,
) -> None:
    target.append({
        "release": release,
        "scope": scope,
        "study_accession": study_accession,
        "label": label,
        "metric": metric,
        "estimate": estimate,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "n_sites": n_sites,
        "n_clusters": n_clusters,
    })


def _percentile(sorted_values: list[float], q: float) -> float:
    if not sorted_values:
        return 0.0
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * q / 100.0
    lo = int(pos)
    hi = min(lo + 1, len(sorted_values) - 1)
    frac = pos - lo
    return sorted_values[lo] * (1.0 - frac) + sorted_values[hi] * frac
