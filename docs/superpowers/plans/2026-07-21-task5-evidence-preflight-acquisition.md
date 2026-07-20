# Task 5 Evidence Preflight Acquisition Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Acquire and register the exact small Stage A PRIDE files, then publish a deterministic metadata-level evidence-preflight inventory without creating benchmark labels.

**Architecture:** Extend the existing exact-name download selection and reuse the fail-closed registered downloader. A new evidence-preflight module validates a versioned policy, joins the official PRIDE file registry to the central download registry, and atomically publishes metadata-level decisions. Content-level evidence parsing is a second cycle written only after the real Stage A bytes and schemas are available.

**Tech Stack:** Python 3.10, standard-library CSV/JSON/dataclasses/pathlib, PyYAML, pytest, existing PlantPersulf provenance and download modules.

## Global Constraints

- Execute only the acquisition and metadata-inventory cycle of the Task 5 preflight; do not build Task 5 labels or benchmark artifacts.
- Use only PXD006140, PXD024061, PXD035795, and PXD039999 official registered PRIDE metadata and exact selected files.
- Never assign `site_ms`, `site_mutagenesis`, `site_biochemical`, `positive`, `unlabeled`, or `negative` from filenames or metadata inventory.
- Preserve all biological values and identifiers; missing information is a reportable gap.
- Every downloaded input must resolve exactly once in `data/registry/files.tsv`, pass its repository SHA1 when present, and be registered with a local SHA256 in `data/registry/downloads.tsv`.
- Do not download raw spectra, peak files, PXD024061 `txt_persulfproject.zip`, search FASTA files, PXD039999 `proteinSeq.txt`, or publication supplements in this cycle.
- Apply strict RED, GREEN, REFACTOR order. Do not change implementation before the matching RED test has failed for the expected reason.

---

### Task 1: Freeze the Stage A policy and exact registered selection

**Files:**
- Create: `configs/evidence_preflight_v1.yaml`
- Modify: `configs/download_selection.yaml`
- Create: `src/plantpersulf/evidence/__init__.py`
- Create: `src/plantpersulf/evidence/policy.py`
- Create: `tests/scientific/test_evidence_preflight_selection.py`

**Interfaces:**
- Consumes: `resolve_registered_selection(accession, file_classes, selection_path, files_registry_path)` from `plantpersulf.download.registered`.
- Produces: `PreflightPolicy`, `load_preflight_policy(path: Path) -> PreflightPolicy`, and exact Stage A selection entries in `configs/download_selection.yaml`.

- [ ] **Step 1: Write the policy and selection RED test**

Create a scientific test that imports `load_preflight_policy`, requires policy version `1`, exact study order, threshold `104857600`, and the exact evidence classes from the design. Resolve the following exact selections against the real `files.tsv`:

```python
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
```

Assert every selected row is `record_type == "source_file"`, has a nonempty official URL, positive size, `remote_checksum_algorithm == "SHA1"`, and a 40-character repository checksum. Assert no selected filename has category `RAW` or `PEAK` and that `txt_persulfproject.zip`, `proteinSeq.txt`, and `arabidopsis_uniprot_072020_identified.fasta` are absent.

- [ ] **Step 2: Run the RED test**

Run:

```bash
python -m pytest tests/scientific/test_evidence_preflight_selection.py -v
```

Expected: collection fails because `plantpersulf.evidence.policy` and the new policy file do not exist. The failure must not be caused by missing real PRIDE registry rows.

- [ ] **Step 3: Implement the minimal policy loader and selection**

Implement immutable policy records:

```python
@dataclass(frozen=True)
class PreflightPolicy:
    version: int
    studies: tuple[str, ...]
    evidence_classes: tuple[str, ...]
    large_file_threshold_bytes: int
    prohibited_label_terms: tuple[str, ...]
```

`load_preflight_policy` must reject non-mappings, version other than `1`, duplicate/non-PXD studies, evidence classes other than the six approved design values, threshold other than `104857600`, and prohibited terms other than `positive`, `unlabeled`, and `negative`.

Create the YAML with study order PXD006140, PXD024061, PXD035795, PXD039999. Extend `download_selection.yaml` only with the exact `EXPECTED` mappings above; preserve the existing PXD006140 entries byte-for-byte except for YAML reserialization that is strictly necessary.

- [ ] **Step 4: Run GREEN and regression tests**

Run:

```bash
python -m pytest tests/scientific/test_evidence_preflight_selection.py tests/unit/test_registered_download_selection.py -v
python -m plantpersulf.cli audit-registry
```

Expected: all selected tests pass and registry counts remain 9 datasets, 180 files, 44 samples, and 10 publications.

- [ ] **Step 5: Commit the policy cycle**

```powershell
git add -- configs/evidence_preflight_v1.yaml configs/download_selection.yaml src/plantpersulf/evidence tests/scientific/test_evidence_preflight_selection.py
git commit -m "feat: freeze evidence preflight acquisition policy"
```

---

### Task 2: Download and register the real Stage A files

**Files:**
- Modify: `data/registry/downloads.tsv`
- Create ignored real inputs under: `data/raw/PXD024061/`, `data/raw/PXD035795/`, `data/raw/PXD039999/`
- Test: `tests/scientific/test_evidence_preflight_selection.py`

**Interfaces:**
- Consumes: existing `download_registered_files` and Task 1 selections.
- Produces: seven new provenance rows in `data/registry/downloads.tsv` and seven local files with exact registered bytes.

- [ ] **Step 1: Add the local-provenance RED test**

Extend the scientific test to call `audit_downloaded_files` for PXD024061, PXD035795, and PXD039999 using the real central registries. Require summaries of `(1, 1, 0)`, `(5, 5, 0)`, and `(1, 1, 0)` for selected, downloaded, and cached counts.

- [ ] **Step 2: Run the RED test**

Run:

```bash
python -m pytest tests/scientific/test_evidence_preflight_selection.py::test_stage_a_downloads_are_registered -v
```

Expected: FAIL because the seven selected files do not yet have rows in `downloads.tsv`; no test may fabricate local files.

- [ ] **Step 3: Download through the existing fail-closed path**

Run sequentially:

```bash
python scripts/fetch_registered_data.py --accession PXD024061 --file-class checksums
python scripts/fetch_registered_data.py --accession PXD035795 --file-class design,results,checksums
python scripts/fetch_registered_data.py --accession PXD039999 --file-class checksums
```

Each command must use the default official registry, central downloads registry, and raw directory. Record exact command output, file sizes, repository SHA1 values, computed SHA256 values, and any external failure. Do not retry by changing URLs or filenames outside the registered rows.

- [ ] **Step 4: Run GREEN audits**

Run:

```bash
python -m pytest tests/scientific/test_evidence_preflight_selection.py::test_stage_a_downloads_are_registered -v
python -m plantpersulf.cli audit-files --accession PXD024061
python -m plantpersulf.cli audit-files --accession PXD035795
python -m plantpersulf.cli audit-files --accession PXD039999
```

Expected: the scientific test passes; each CLI audit reports every selected source downloaded with zero metadata-cache selections.

- [ ] **Step 5: Commit the download provenance**

Commit only the registry row changes and test. Raw inputs stay ignored but must remain present locally.

```bash
git add data/registry/downloads.tsv tests/scientific/test_evidence_preflight_selection.py
git commit -m "data: register evidence preflight source files"
```

---

### Task 3: Publish a deterministic metadata-level inventory

**Files:**
- Create: `src/plantpersulf/evidence/preflight.py`
- Create: `tests/scientific/test_evidence_preflight_inventory.py`
- Modify: `src/plantpersulf/cli.py`
- Create ignored outputs under: `data/interim/evidence_preflight_v1/`

**Interfaces:**
- Consumes: `PreflightPolicy`, `files.tsv`, `downloads.tsv`, `datasets.tsv`, and `download_selection.yaml`.
- Produces: `build_metadata_preflight(...) -> PreflightSummary`, `audit_metadata_preflight(...) -> PreflightSummary`, CLI commands `build-evidence-preflight` and `audit-evidence-preflight`.

- [ ] **Step 1: Write deterministic inventory RED tests**

Tests must use the real registries and assert:

1. exactly four study rows;
2. every official PRIDE file receives one file-decision row;
3. the seven new Stage A files are `downloaded` with local SHA256;
4. the existing two PXD006140 OMSSA files are `downloaded_existing_scope`;
5. `txt_persulfproject.zip` is `deferred_large_not_stage_a` with size `2417963957`;
6. all RAW and PEAK rows are deferred and never selected;
7. `evidence_records.tsv` contains only its header in cycle 1;
8. `metadata_gaps.tsv` explicitly records missing PXD006140 sample mapping and `content_audit_pending` for Stage A result/design files;
9. two builds in different temporary roots have byte-identical output hashes;
10. no TSV schema contains a `label` field and no output claims that nondetection is negative; the manifest may contain only the required boolean integrity flags whose names mention labels or negatives.

Add a corruption test that changes a declared output hash in a temporary manifest and requires the audit to fail.

- [ ] **Step 2: Run the RED tests**

Run:

```bash
python -m pytest tests/scientific/test_evidence_preflight_inventory.py -v
```

Expected: FAIL because `plantpersulf.evidence.preflight` and the CLI commands do not exist.

- [ ] **Step 3: Implement fixed schemas and atomic publication**

Define immutable `FileDecision` and `PreflightSummary` records and fixed TSV columns. Join registry rows by `(dataset_accession, file_name)` only after enforcing uniqueness. Rehash every local selected file through its exact `downloads.tsv` row.

Decision order is fixed:

```text
registered metadata cache -> registered_metadata
selected and locally audited -> downloaded
existing PXD006140 approved result -> downloaded_existing_scope
size >= 104857600 -> deferred_large_not_stage_a
otherwise -> deferred_not_stage_a
```

Every file decision has `evidence_class=unresolved`. Output order is study order from the policy, then registry filename. Write all TSV files to a temporary sibling, compute hashes, write a sorted JSON manifest last, and atomically rename. Existing identical output is accepted; differing existing output is fatal.

The manifest includes source registry hashes, selection and policy hashes, output hashes, deterministic counts, `content_audit_complete=false`, `labels_created=false`, `nondetection_labeled_negative=false`, and `biological_values_modified=false`.

- [ ] **Step 4: Add CLI commands**

Add:

```bash
python -m plantpersulf.cli build-evidence-preflight --version v1
python -m plantpersulf.cli audit-evidence-preflight --version v1
```

Reject versions other than `v1`. Defaults point to the central registries and `data/interim/evidence_preflight_v1`.

- [ ] **Step 5: Run GREEN tests and real publication twice**

Run:

```bash
python -m pytest tests/scientific/test_evidence_preflight_inventory.py -v
python -m plantpersulf.cli build-evidence-preflight --version v1
python -m plantpersulf.cli audit-evidence-preflight --version v1
python -m plantpersulf.cli build-evidence-preflight --version v1
python -m plantpersulf.cli audit-evidence-preflight --version v1
```

Expected: tests pass and both real runs return identical counts and output hashes.

- [ ] **Step 6: Commit the metadata inventory cycle**

```powershell
git add -- src/plantpersulf/evidence src/plantpersulf/cli.py tests/scientific/test_evidence_preflight_inventory.py
git commit -m "feat: publish deterministic evidence preflight inventory"
```

---

### Task 4: Complete cycle-1 verification and schema handoff

**Files:**
- Modify: `docs/superpowers/specs/2026-07-21-task5-evidence-preflight-design.md` only if real schemas disprove an explicit design statement.
- Create: `docs/superpowers/specs/2026-07-21-task5-evidence-content-audit-design.md` after inspecting registered Stage A headers and format metadata.

**Interfaces:**
- Consumes: registered Stage A files and metadata preflight manifest.
- Produces: a cycle-1 completion report and a separate content-audit design grounded in observed real schemas; no content parser is implemented in this cycle.

- [ ] **Step 1: Inspect registered schemas without modifying values**

Record exact encodings, delimiters, headers, row counts, mzIdentML namespace/version, modification CV accessions/names, and SDRF columns from the seven registered files. Do not copy biological rows into documentation; cite file SHA256 plus structural summaries.

- [ ] **Step 2: Run all completion gates**

```bash
python -m pytest tests/unit tests/scientific tests/release -q
python -m pytest -m network tests/integration -q
python -m plantpersulf.cli audit-registry
python -m plantpersulf.cli audit-files --accession PXD006140
python -m plantpersulf.cli audit-files --accession PXD024061
python -m plantpersulf.cli audit-files --accession PXD035795
python -m plantpersulf.cli audit-files --accession PXD039999
python -m plantpersulf.cli audit-evidence-preflight --version v1
python -m ruff check .
python -m mypy src/plantpersulf scripts
git diff --check
```

Expected: all tests, provenance audits, output audit, lint, typing, and whitespace checks pass.

- [ ] **Step 3: Scientific-integrity self-review**

Confirm from staged paths and generated manifest that no benchmark, label, split, feature, model, metric, plot, candidate, or conclusion artifact exists; every local scientific input is SHA256-registered; and no nondetection or missing record is represented as negative.

- [ ] **Step 4: Commit the observed-schema content-audit design**

Write and self-review the separate design with exact observed schemas and conservative evidence rules. Commit it independently; stop before implementing its content parser unless that second cycle is explicitly continued under its own plan.

```bash
git add docs/superpowers/specs/2026-07-21-task5-evidence-content-audit-design.md
git commit -m "docs: design real-file evidence content audit"
```

- [ ] **Step 5: Report and stop**

Report every RED and GREEN command, downloaded files with official URLs and both repository checksum and SHA256, created/modified files, deterministic inventory counts, observed schema summary, integrity answers, commit SHAs, external blockers, and whether the Task 5 label-building gate is satisfied. Do not start benchmark implementation.
