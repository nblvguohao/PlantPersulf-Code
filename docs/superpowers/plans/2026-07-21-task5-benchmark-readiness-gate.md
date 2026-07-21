# Task 5 Benchmark Readiness Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish a deterministic, fail-closed decision showing whether the registered real evidence is sufficient to start Task 5 positive–unlabeled benchmark construction, without creating any benchmark labels when no eligible persulfidation-site evidence exists.

**Architecture:** A frozen YAML policy defines only the evidence levels already permitted by the project specification and the required minimum of two distinct studies. A readiness builder audits the existing PXD006140 coordinate output and Cycle 2 content output, classifies their already-recorded evidence levels without reinterpretation, and writes study-level counts, explicit blockers, and a hashed manifest. A separate auditor rehashes every input and output and rejects any benchmark artifact or label column.

**Tech Stack:** Python 3.10, standard-library `csv` and `json`, PyYAML, existing Task 4 and Cycle 2 auditors, pytest, Ruff, mypy.

## Global Constraints

- Use only registered, checksum-verified real files and deterministic outputs derived from PXD006140, PXD024061, PXD035795, and PXD039999.
- Eligible site evidence levels are exactly `site_ms`, `site_mutagenesis`, and `site_biochemical`.
- `psm_coordinate_only`, `protein_level_only`, `identification_only`, and `unresolved` are never upgraded to eligible site evidence.
- Require at least two distinct studies with eligible site evidence before benchmark construction is permitted, as required by Gate 1 in `docs/PlantPersulf_Code_TDD_Codex.md`.
- Do not create `positive`, `unlabeled`, or `negative` rows; do not create `data/processed/benchmark_v1`; do not start Task 6, features, models, metrics, figures, candidates, or conclusions.
- PXD006140's three coordinate rows retain `evidence_level=psm_coordinate_only`; their methionine-oxidation annotations are not persulfidation calls.
- PXD035795 retains zero `site_ms`; its protein and peptide records remain non-site evidence and its candidate modifications remain unresolved.

---

### Task 1: Freeze the readiness policy and publish a fail-closed audit

**Files:**
- Create: `configs/benchmark_readiness_v1.yaml`
- Create: `src/plantpersulf/benchmark/__init__.py`
- Create: `src/plantpersulf/benchmark/readiness.py`
- Create: `tests/scientific/test_benchmark_readiness_gate.py`

**Interfaces:**
- Consumes: `data/interim/PXD006140/proteomics_parser_v1/{sites.tsv,manifest.json}`, `data/interim/evidence_preflight_v1/{study_inventory.tsv,evidence_records.tsv,candidate_modifications.tsv,manifest.json}`, and their existing auditors.
- Produces: `build_benchmark_readiness(policy_path: Path, parser_output_root: Path, content_output_directory: Path, output_directory: Path, registry_dir: Path) -> BenchmarkReadinessSummary` and `audit_benchmark_readiness(output_directory: Path) -> BenchmarkReadinessSummary`.
- Publishes: `study_readiness.tsv`, `blockers.tsv`, and `manifest.json` under `data/interim/benchmark_readiness_v1`.

- [ ] **Step 1: Write the failing real-evidence readiness tests**

```python
def test_current_real_evidence_forces_benchmark_stop(tmp_path: Path) -> None:
    parser_root, content_dir = _build_real_prerequisites(tmp_path)
    summary = build_benchmark_readiness(
        policy_path=Path("configs/benchmark_readiness_v1.yaml"),
        parser_output_root=parser_root,
        content_output_directory=content_dir,
        output_directory=tmp_path / "readiness",
        registry_dir=Path("data/registry"),
    )

    assert summary.decision == "STOP"
    assert summary.eligible_site_count == 0
    assert summary.eligible_study_count == 0
    assert summary.coordinate_only_count == 3
    assert summary.non_site_evidence_count == 10_787
    assert summary.unresolved_candidate_count == 25
    assert summary.required_study_count == 2


def test_readiness_output_contains_no_benchmark_or_label_artifact(
    tmp_path: Path,
) -> None:
    output = _build_readiness(tmp_path)
    assert {path.name for path in output.iterdir()} == {
        "study_readiness.tsv",
        "blockers.tsv",
        "manifest.json",
    }
    for tsv_path in output.glob("*.tsv"):
        header = tsv_path.read_text(encoding="utf-8").splitlines()[0]
        assert "label" not in header.lower()
    assert not (tmp_path / "processed" / "benchmark_v1").exists()
```

- [ ] **Step 2: Run RED and confirm the feature is missing**

Run: `python -m pytest tests/scientific/test_benchmark_readiness_gate.py -v`

Expected: FAIL because `plantpersulf.benchmark.readiness` and `configs/benchmark_readiness_v1.yaml` do not exist.

- [ ] **Step 3: Implement the minimal frozen policy and publisher**

Use this exact policy:

```yaml
version: 1
studies:
  - PXD006140
  - PXD024061
  - PXD035795
  - PXD039999
eligible_site_evidence_levels:
  - site_ms
  - site_mutagenesis
  - site_biochemical
coordinate_only_evidence_levels:
  - psm_coordinate_only
minimum_distinct_site_evidence_studies: 2
prohibited_record_labels:
  - positive
  - unlabeled
  - negative
```

Use immutable dataclasses and these exact TSV fields:

```python
STUDY_READINESS_FIELDS = (
    "study_accession",
    "eligible_site_count",
    "coordinate_only_count",
    "non_site_evidence_count",
    "unresolved_candidate_count",
    "content_audit_status",
    "readiness_status",
)

BLOCKER_FIELDS = (
    "scope",
    "study_accession",
    "blocker_code",
    "observed_value",
    "required_value",
    "source_locator",
)
```

Before reading rows, call `audit_site_output` and `audit_content_output`. Count only literal evidence values already present in the audited TSVs. Emit one study row for each frozen accession. Emit study blockers for coordinate-only or unsupported evidence and the global blocker `insufficient_distinct_site_evidence_studies`. Write with sorted rows and `\n` endings through a temporary directory, then atomically replace the output. The manifest must include all input and output SHA256 values plus:

```json
{
  "decision": "STOP",
  "benchmark_created": false,
  "labels_created": false,
  "nondetection_labeled_negative": false,
  "biological_values_modified": false
}
```

The decision may be `GO` only when the literal eligible-site rows span at least two frozen studies; the current registered evidence must produce `STOP`.

- [ ] **Step 4: Run GREEN and determinism checks**

Run:

```powershell
python -m pytest tests/scientific/test_benchmark_readiness_gate.py -v
```

Expected: readiness tests PASS and their scope assertion confirms that no benchmark fixture or label file exists.

- [ ] **Step 5: Commit the readiness publisher**

```powershell
git add configs/benchmark_readiness_v1.yaml src/plantpersulf/benchmark tests/scientific/test_benchmark_readiness_gate.py
git commit -m "feat: gate benchmark construction on eligible site evidence"
```

---

### Task 2: Add CLI entry points and enforce the STOP decision end to end

**Files:**
- Modify: `src/plantpersulf/cli.py`
- Modify: `tests/scientific/test_benchmark_readiness_gate.py`
- Create: `tests/release/test_benchmark_stop_has_no_labels.py`

**Interfaces:**
- Produces CLI commands `build-benchmark-readiness --version v1` and `audit-benchmark-readiness --version v1`.
- Both commands use `data/interim/PXD006140/proteomics_parser_v1`, `data/interim/evidence_preflight_v1`, and `data/interim/benchmark_readiness_v1` by default.

- [ ] **Step 1: Write failing CLI and release-gate tests**

```python
def test_cli_exposes_benchmark_readiness_without_building_benchmark() -> None:
    build = build_parser().parse_args(
        ["build-benchmark-readiness", "--version", "v1"]
    )
    audit = build_parser().parse_args(
        ["audit-benchmark-readiness", "--version", "v1"]
    )
    assert build.command == "build-benchmark-readiness"
    assert audit.command == "audit-benchmark-readiness"


def test_stop_manifest_is_label_free_and_fail_closed(tmp_path: Path) -> None:
    manifest = _build_and_load_manifest(tmp_path)
    assert manifest["decision"] == "STOP"
    assert manifest["benchmark_created"] is False
    assert manifest["labels_created"] is False
    assert manifest["nondetection_labeled_negative"] is False
```

- [ ] **Step 2: Run RED**

Run: `python -m pytest tests/scientific/test_benchmark_readiness_gate.py tests/release/test_benchmark_stop_has_no_labels.py -v`

Expected: FAIL because the CLI parsers and dispatch branches do not exist.

- [ ] **Step 3: Add minimal CLI parsing and dispatch**

`build-benchmark-readiness` first calls the existing real-data builders for PXD006140 and Cycle 2, then calls `build_benchmark_readiness`. `audit-benchmark-readiness` only audits the frozen readiness directory. Both print the summary dataclass as sorted JSON and return zero for a scientifically valid `STOP`; STOP is a result, not a software error.

- [ ] **Step 4: Run GREEN and the full completion gate**

```powershell
python -m plantpersulf.cli build-benchmark-readiness --version v1
python -m plantpersulf.cli audit-benchmark-readiness --version v1
python -m pytest tests/unit tests/scientific tests/release -q
python -m pytest -m network tests/integration -q
python -m ruff check .
python -m mypy src/plantpersulf
git ls-files data/raw data/interim data/processed results
```

Expected: both CLI commands report `decision=STOP`, every test and static check exits zero, and the final Git command prints nothing.

- [ ] **Step 5: Commit the CLI and release gate**

```powershell
git add src/plantpersulf/cli.py tests/scientific/test_benchmark_readiness_gate.py tests/release/test_benchmark_stop_has_no_labels.py
git commit -m "feat: publish fail-closed Task 5 readiness decision"
```

## Self-review

- Spec coverage: Task 5 admissible evidence levels, Gate 1's two-study minimum, no-negative rule, real-file hashes, deterministic outputs, and STOP behavior are each enforced by code and tests.
- Placeholder scan: the plan contains no invented biological values, placeholder result, fake label, fake metric, or unresolved implementation instruction.
- Type consistency: both tasks use the same `BenchmarkReadinessSummary`, output names, and default directories.
- Scope boundary: this cycle publishes only a readiness decision. It does not implement the original Task 5 benchmark files because the audited eligible-site count is zero.
