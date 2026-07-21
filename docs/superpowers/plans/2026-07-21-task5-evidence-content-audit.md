# Task 5 Real-File Evidence Content Audit Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Audit the registered PXD035795 Stage A files at content level, preserve their exact biological strings and reference conflicts, and publish deterministic non-label evidence tables without starting benchmark construction.

**Architecture:** A method-source registry gates scientific interpretation independently from file parsing. Focused positional SDRF, exact-schema CSV, and streaming mzIdentML readers return immutable records; a content publisher joins only exact identifiers, emits every conflict, and defaults every localized modification to `unresolved` because version 1 contains no reviewed mass-to-persulfidation mapping. The existing metadata preflight remains the base inventory and its integrity flags remain label-free.

**Tech Stack:** Python 3.10, standard-library `csv`, `gzip`, `json`, `xml.etree.ElementTree`, existing atomic downloader/hash registry, PyYAML, pytest, Ruff, mypy.

## Global Constraints

- Use only registered real files from PXD006140, PXD024061, PXD035795, and PXD039999; do not synthesize biological values or labels.
- Do not treat nondetection, zero, blank quantification, or an unmatched file as a negative.
- Do not infer a site from CSV PTM text, protein-level records, approximate mass agreement, URL decoding, or filename rewriting.
- `site_ms` is forbidden in this cycle unless a versioned reviewed method mapping exists; `configs/evidence_method_mappings_v1.yaml` intentionally contains no mappings.
- Preserve the original source SHA256 and a reproducible locator for every emitted biological record.
- Do not create benchmark files, labels, splits, features, models, metrics, figures, candidates, or conclusions.

---

### Task 1: Register and audit the official PXD035795 method source

**Files:**
- Create: `configs/evidence_method_sources_v1.yaml`
- Create: `configs/evidence_method_mappings_v1.yaml`
- Create: `data/registry/evidence_methods.tsv`
- Create: `src/plantpersulf/evidence/methods.py`
- Modify: `src/plantpersulf/cli.py`
- Test: `tests/scientific/test_evidence_method_registry.py`

**Interfaces:**
- Consumes: `download_verified_file(DownloadRequest)`, DOI `10.1111/nph.18838`, University of Seville repository record `ce5a8e5e-6bf0-4564-a0d5-c6dcb4eeaa6f`, and its public PDF bitstream URL.
- Produces: `acquire_method_source(accession: str, config_path: Path, registry_path: Path, raw_root: Path) -> MethodSource`; `audit_method_sources(config_path: Path, registry_path: Path) -> tuple[MethodSource, ...]`; CLI commands `acquire-evidence-method` and `audit-evidence-methods`.

- [ ] **Step 1: Write the failing registry and no-auto-upgrade tests**

```python
def test_pxd035795_method_source_must_be_registered_and_hashed() -> None:
    from plantpersulf.evidence.methods import audit_method_sources

    sources = audit_method_sources(
        Path("configs/evidence_method_sources_v1.yaml"),
        Path("data/registry/evidence_methods.tsv"),
    )
    source = next(row for row in sources if row.study_accession == "PXD035795")
    assert source.identifier == "10.1111/nph.18838"
    assert len(source.sha256) == 64
    assert source.local_path.is_file()


def test_version_one_has_no_reviewed_site_upgrade_mapping() -> None:
    mapping = yaml.safe_load(
        Path("configs/evidence_method_mappings_v1.yaml").read_text(encoding="utf-8")
    )
    assert mapping == {"version": 1, "mappings": []}
```

- [ ] **Step 2: Run RED and verify missing module/config/registry is the reason**

Run: `python -m pytest tests/scientific/test_evidence_method_registry.py -v`

Expected: FAIL importing `plantpersulf.evidence.methods` or finding the reviewed config/registry, before any scientific record can be upgraded.

- [ ] **Step 3: Implement the minimal fail-closed source registry**

Use this exact TSV schema:

```python
METHOD_FIELDS = (
    "study_accession", "identifier_type", "identifier", "repository",
    "repository_record_id", "official_url", "local_path", "retrieved_at",
    "size_bytes", "sha256", "method_scope", "data_use_status",
    "downloader_version",
)
```

The YAML allow-list fixes the accession, DOI, official URL, destination name, observed HTTP size `11657462`, method scope `dimedone_switch_lc_ms_ms`, and use status `method_audit_only`. Acquisition must use a temporary sibling through `download_verified_file`, require the exact byte size, compute SHA256, append exactly one sorted registry row, and refuse identity drift. Audit must rehash the local PDF and reject missing, duplicate, malformed, or mismatching rows.

- [ ] **Step 4: Acquire the real PDF and run GREEN**

Run:

```powershell
python -m plantpersulf.cli acquire-evidence-method --accession PXD035795
python -m plantpersulf.cli audit-evidence-methods
python -m pytest tests/scientific/test_evidence_method_registry.py -v
```

Expected: one registered, size-checked method PDF; all method-registry tests PASS; mappings remain empty.

- [ ] **Step 5: Commit**

```powershell
git add configs/evidence_method_sources_v1.yaml configs/evidence_method_mappings_v1.yaml data/registry/evidence_methods.tsv src/plantpersulf/evidence/methods.py src/plantpersulf/cli.py tests/scientific/test_evidence_method_registry.py
git commit -m "data: register audited PXD035795 method source"
```

---

### Task 2: Parse the SDRF positionally and expose exact mapping conflicts

**Files:**
- Create: `src/plantpersulf/evidence/sdrf.py`
- Test: `tests/scientific/test_pxd035795_sdrf_content.py`

**Interfaces:**
- Consumes: registered `data/raw/PXD035795/SDRF.txt` and `data/registry/files.tsv`.
- Produces: `parse_sdrf(path: Path, source_sha256: str, files_registry_path: Path) -> SdrfAudit`, containing `rows`, `file_mappings`, `conflicts`, `duplicate_header_positions`, and `schema`.

- [ ] **Step 1: Write the failing positional parsing test against the real bytes**

```python
def test_real_sdrf_preserves_duplicate_headers_and_exact_file_mapping() -> None:
    audit = parse_sdrf(SDRF_PATH, SDRF_SHA256, FILES_PATH)
    assert len(audit.header) == 26
    assert len(audit.rows) == 6
    assert audit.duplicate_header_positions[
        "comment[modification parameters]"
    ] == (19, 20, 21, 22, 23)
    assert len(audit.file_mappings) == 12
    assert sum(row.mapping_status == "exact" for row in audit.file_mappings) == 6
    assert sum(row.mapping_status == "conflict" for row in audit.file_mappings) == 6
```

Also assert all six RAW values match registry names exactly, all encoded MGF-like values remain unchanged and conflicted, and five modification strings are preserved in position order.

- [ ] **Step 2: Run RED**

Run: `python -m pytest tests/scientific/test_pxd035795_sdrf_content.py -v`

Expected: FAIL because `plantpersulf.evidence.sdrf` does not exist.

- [ ] **Step 3: Implement the minimal positional parser**

Read UTF-8 with optional BOM using `csv.reader(delimiter="\t")`. Require the exact 26-field registered header and six rows. Split `comment[data file]` only on commas; do not URL-decode, strip embedded percent sequences, replace spaces, or perform fuzzy matching. Emit a mapping row for each component and an explicit `unmatched_exact_filename` conflict when the `(PXD035795, file_name)` registry key is absent.

- [ ] **Step 4: Run GREEN and regressions**

Run:

```powershell
python -m pytest tests/scientific/test_pxd035795_sdrf_content.py -v
python -m pytest tests/scientific/test_evidence_preflight_inventory.py -q
```

Expected: all tests PASS.

- [ ] **Step 5: Commit**

```powershell
git add src/plantpersulf/evidence/sdrf.py tests/scientific/test_pxd035795_sdrf_content.py
git commit -m "feat: audit PXD035795 SDRF mappings positionally"
```

---

### Task 3: Preserve exact peptide and protein CSV content without inferring sites

**Files:**
- Create: `src/plantpersulf/evidence/csv_content.py`
- Test: `tests/scientific/test_pxd035795_csv_content.py`

**Interfaces:**
- Consumes: registered `peptide.csv` and `proteins.csv` plus their exact SHA256 values.
- Produces: `parse_peptide_csv(path: Path, source_sha256: str) -> CsvAudit` and `parse_protein_csv(path: Path, source_sha256: str) -> CsvAudit`; each audit contains exact `header`, `rows`, `declared_modifications`, `schema`, and a record type.

- [ ] **Step 1: Write the failing exact-schema tests**

```python
def test_real_peptide_csv_preserves_all_rows_and_has_no_inferred_site() -> None:
    audit = parse_peptide_csv(PEPTIDE_PATH, PEPTIDE_SHA256)
    assert len(audit.header) == 19
    assert len(audit.rows) == 9326
    assert all(record.site_localization_status == "not_available" for record in audit.rows)
    assert all(record.evidence_class == "identification_only" for record in audit.rows)


def test_real_protein_csv_is_never_promoted_to_site_level() -> None:
    audit = parse_protein_csv(PROTEIN_PATH, PROTEIN_SHA256)
    assert len(audit.header) == 28
    assert len(audit.rows) == 1461
    assert all(record.site_localization_status == "not_applicable" for record in audit.rows)
    assert not any(record.evidence_class == "site_ms" for record in audit.rows)
```

Tests additionally compare representative raw row tuples to `csv.reader` output and require blank/zero/area strings to remain byte-derived strings, not imputed values.

- [ ] **Step 2: Run RED**

Run: `python -m pytest tests/scientific/test_pxd035795_csv_content.py -v`

Expected: FAIL because the exact CSV reader module is missing.

- [ ] **Step 3: Implement minimal exact readers**

Use `csv.reader` with UTF-8-SIG and require the literal headers shown in the design. Require 9326 and 1461 rows respectively. Preserve every field as a string tuple and serialize it later as compact JSON. Enumerate PTM names only by splitting the raw PTM field on semicolons; never interpret a residue coordinate. Peptide records remain `identification_only`; protein records may become `protein_level_only` only in the publisher after a method-source audit succeeds.

- [ ] **Step 4: Run GREEN**

Run: `python -m pytest tests/scientific/test_pxd035795_csv_content.py -v`

Expected: all tests PASS.

- [ ] **Step 5: Commit**

```powershell
git add src/plantpersulf/evidence/csv_content.py tests/scientific/test_pxd035795_csv_content.py
git commit -m "feat: preserve exact PXD035795 CSV evidence content"
```

---

### Task 4: Stream mzIdentML references and localized candidate modifications

**Files:**
- Create: `src/plantpersulf/evidence/mzidentml.py`
- Test: `tests/scientific/test_pxd035795_mzidentml_content.py`

**Interfaces:**
- Consumes: registered `peptides_1_1_0.mzid.gz` with SHA256 `62105dfde42d9675bc5dc7c83969d971eec57c9e02983659cd2df96ec4d1b5fe`.
- Produces: `parse_mzidentml(path: Path, source_sha256: str) -> MzidAudit`, containing exact structural counts, `candidates`, `reference_conflicts`, and namespace/version.

- [ ] **Step 1: Write the failing real-file structural and policy tests**

```python
def test_real_mzid_preserves_registered_structure_and_candidates() -> None:
    audit = parse_mzidentml(MZID_PATH, MZID_SHA256)
    assert audit.namespace == "http://psidev.info/psi/pi/mzIdentML/1.1"
    assert audit.version == "1.1.0"
    assert audit.structural_counts == {
        "DBSequence": 5,
        "Peptide": 11409,
        "PeptideEvidence": 95,
        "SpectrumIdentificationResult": 661,
        "SpectrumIdentificationItem": 699,
        "Modification": 124,
    }
    assert audit.candidates
    assert {row.evidence_class for row in audit.candidates} == {"unresolved"}
    assert not any(row.evidence_class == "site_ms" for row in audit.candidates)
```

Also require exact `Peptide_ref`, `PeptideEvidence_ref`, and `DBSequence_ref` resolution; pass-threshold and rank preservation; 1-based modification range checking; sequence/residue disagreement and missing references in `reference_conflicts`; and source locators containing stable mzIdentML IDs.

- [ ] **Step 2: Run RED**

Run: `python -m pytest tests/scientific/test_pxd035795_mzidentml_content.py -v`

Expected: FAIL because the streaming mzIdentML reader is missing.

- [ ] **Step 3: Implement the minimal streaming parser**

Use `gzip.open(..., "rb")` with `ElementTree.iterparse`. Index DBSequence, Peptide, and PeptideEvidence identifiers exactly, then resolve each SpectrumIdentificationItem. Emit one candidate per identification-item × peptide-evidence × localized modification mapping, preserving blank residues, numeric mass strings, CV accession/name, rank, threshold, and all raw identifiers. Never compare masses with a tolerance in this parser and always set `evidence_class="unresolved"` and `method_mapping_status="absent_v1"`.

- [ ] **Step 4: Run GREEN and check memory-safe behavior**

Run: `python -m pytest tests/scientific/test_pxd035795_mzidentml_content.py -v`

Expected: all tests PASS and the exact registered structural counts are retained.

- [ ] **Step 5: Commit**

```powershell
git add src/plantpersulf/evidence/mzidentml.py tests/scientific/test_pxd035795_mzidentml_content.py
git commit -m "feat: audit PXD035795 mzIdentML candidate modifications"
```

---

### Task 5: Publish and audit the deterministic Cycle 2 inventory

**Files:**
- Create: `src/plantpersulf/evidence/content.py`
- Modify: `src/plantpersulf/evidence/preflight.py`
- Modify: `src/plantpersulf/cli.py`
- Create: `tests/scientific/test_evidence_content_audit.py`
- Modify: `tests/release/test_no_synthetic_scientific_artifacts.py`

**Interfaces:**
- Consumes: audited method registry, empty reviewed mapping config, four registered PXD035795 design/result files, Cycle 1 policy/registries, and the readers from Tasks 2–4.
- Produces: `build_content_audit(...) -> ContentAuditSummary`, `audit_content_output(...) -> ContentAuditSummary`, and CLI commands `build-evidence-content-audit --version v1` and `audit-evidence-content --version v1`.

- [ ] **Step 1: Write failing publication, integrity, and determinism tests**

Require the final directory to contain the six Cycle 1 files plus:

```python
CONTENT_OUTPUT_NAMES = {
    "candidate_modifications.tsv",
    "sample_file_mapping.tsv",
    "mapping_conflicts.tsv",
    "content_schema.tsv",
}
```

Use this evidence schema:

```python
EVIDENCE_FIELDS = (
    "dataset_accession", "record_id", "record_type", "source_file",
    "source_sha256", "source_locator", "evidence_class", "method",
    "method_source_sha256", "protein_accession", "peptide_sequence",
    "modification_raw", "raw_values_json", "sample_mapping_status",
    "quantification_status", "site_localization_status",
)
```

Assert two independent builds have identical tree hashes; all source/method/output hashes audit; PXD035795 has no pending content gap; PXD024061/PXD039999 remain explicitly unsupported in Stage A; `content_audit_complete=true`; `labels_created=false`; `nondetection_labeled_negative=false`; `biological_values_modified=false`; no TSV column is named `label`; and no row has `evidence_class=site_ms` while mappings are empty.

- [ ] **Step 2: Run RED**

Run: `python -m pytest tests/scientific/test_evidence_content_audit.py -v`

Expected: FAIL because the content publisher and CLI do not exist.

- [ ] **Step 3: Implement the minimal atomic publisher**

Build the Cycle 1 inventory into a temporary base directory, replace the four PXD035795 `content_audit_pending` gaps with schema/mapping conflict results, update study status, and write sorted TSVs with `\n` endings. Peptide CSV rows become `identification_only`; protein CSV rows become `protein_level_only` only after the registered PDF hash passes; mzIdentML candidates remain `unresolved` in the candidate table and do not become evidence positives. Quantification status remains `unresolved_no_registered_column_mapping`, and all raw CSV fields are JSON strings without numeric conversion. Set `content_audit_complete=true` only after all four registered files have successful schema rows and all integrity audits pass.

- [ ] **Step 4: Run GREEN, repeatability, and CLI audits**

Run:

```powershell
python -m pytest tests/scientific/test_evidence_content_audit.py -v
python -m plantpersulf.cli build-evidence-content-audit --version v1
python -m plantpersulf.cli audit-evidence-content --version v1
python -m plantpersulf.cli build-evidence-content-audit --version v1
python -m plantpersulf.cli audit-evidence-content --version v1
```

Expected: tests PASS, both builds are byte-identical, and both audits print the same summary.

- [ ] **Step 5: Run the full completion gate**

```powershell
python -m pytest tests/unit tests/scientific tests/release -q
python -m pytest -m network tests/integration -q
python -m plantpersulf.cli audit-registry
python -m plantpersulf.cli audit-files --accession PXD006140
python -m plantpersulf.cli audit-files --accession PXD024061
python -m plantpersulf.cli audit-files --accession PXD035795
python -m plantpersulf.cli audit-files --accession PXD039999
python -m plantpersulf.cli audit-evidence-methods
python -m plantpersulf.cli audit-evidence-content --version v1
python -m plantpersulf.cli parse-proteomics --accession PXD006140
python -m plantpersulf.cli audit-sites --accession PXD006140
python -m ruff check .
python -m mypy src/plantpersulf
```

Expected: every command exits zero; no labels, splits, model artifacts, or unregistered inputs exist.

- [ ] **Step 6: Commit**

```powershell
git add src/plantpersulf/evidence/content.py src/plantpersulf/evidence/preflight.py src/plantpersulf/cli.py tests/scientific/test_evidence_content_audit.py tests/release/test_no_synthetic_scientific_artifacts.py
git commit -m "feat: publish deterministic real-file content audit"
```

## Self-review

- Spec coverage: method provenance, all four selected PXD035795 design/result files, duplicate SDRF headers, exact conflicts, raw CSV preservation, mzIdentML references, empty method mapping, deterministic publication, and scientific-integrity flags each have a task and test.
- Placeholder scan: no incomplete implementation instruction, guessed biological value, or placeholder result is present.
- Type consistency: Tasks 2–4 return immutable audit objects consumed only by Task 5; method acquisition/audit interfaces and all CLI names are consistent throughout.
- Scope boundary: no benchmark label or model interface is introduced. A zero-`site_ms` result is an allowed and reportable scientific outcome.
