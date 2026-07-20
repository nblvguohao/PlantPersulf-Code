# Task 2 Checksummed Downloads Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Download an explicitly approved set of real repository files, fail closed on transport or integrity errors, record local SHA256 provenance, and audit every downloaded byte.

**Architecture:** A repository-neutral atomic downloader streams into a temporary file and verifies size plus repository checksum before publication. A registered-download layer resolves only exact allowlisted file names against `data/registry/files.tsv`, writes a separate `data/registry/downloads.tsv`, and powers the existing fetch script and a new `audit-files` CLI command.

**Tech Stack:** Python 3.10 standard library, PyYAML, pytest, Ruff, mypy, official PRIDE FTP/HTTPS endpoints.

## Global Constraints

- Execute Task 2 only; do not parse biological values or begin Task 3.
- Use only real public files already tied to an accession and official source URL.
- Never leave a final file after HTTP, size, SHA1, or SHA256 failure.
- Store downloaded data under ignored `data/raw/`; commit provenance, code, and tests only.
- Record accession, repository, exact source file, retrieval time, URL, size, remote checksum, local SHA256, usage statement, and downloader version.
- PXD006140 `results` means exactly the two OMSSA text files named in `configs/download_selection.yaml`; no category or extension inference is permitted.

---

### Task 1: Atomic fail-closed downloader

**Files:**
- Create: `src/plantpersulf/download/base.py`
- Create: `src/plantpersulf/provenance/hashing.py`
- Create: `tests/unit/test_download_failure_is_fatal.py`

**Interfaces:**
- Produces: `DownloadRequest`, `DownloadResult`, and `download_verified_file(request)`.
- Produces: `hash_file(path, algorithm)` for SHA1 and SHA256 streaming hashes.

- [ ] **Step 1: Write the HTTP failure test**

Use a local `ThreadingHTTPServer` handler returning HTTP 503. Call `download_verified_file` and assert `RuntimeError`, an absent destination, and no temporary files.

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/unit/test_download_failure_is_fatal.py::test_http_failure_is_fatal_and_leaves_no_file -v`

Expected: `ModuleNotFoundError: No module named 'plantpersulf.download.base'`.

- [ ] **Step 3: Implement minimal atomic transport**

Define immutable request/result dataclasses. Open the official URL with a project User-Agent, stream bytes into a temporary sibling file, and delete that temporary file on every exception before raising `RuntimeError`.

- [ ] **Step 4: Verify GREEN**

Run the Step 2 command and require one passing test.

- [ ] **Step 5: Write checksum and size failure tests**

Serve non-biological marker bytes from a local HTTP server. Supply an incorrect expected size or checksum and assert neither case publishes the destination.

- [ ] **Step 6: Verify RED**

Run: `python -m pytest tests/unit/test_download_failure_is_fatal.py -v`

Expected: the new integrity assertions fail because transport-only code publishes the file.

- [ ] **Step 7: Implement integrity verification**

Stream SHA1 and SHA256 while downloading, compare byte count and expected remote checksum, then use `Path.replace` only after every check succeeds.

- [ ] **Step 8: Verify GREEN**

Run the Step 6 command and require all tests to pass.

---

### Task 2: Official small-file integration test

**Files:**
- Create: `tests/integration/test_real_file_download.py`

**Interfaces:**
- Consumes: `download_verified_file` and the PXD035795 SDRF row from `data/registry/files.tsv`.

- [ ] **Step 1: Write the real-network test**

Locate registered file `SDRF.txt` for PXD035795, download its official 4,671-byte source URL, verify the registered SHA1, and assert the returned SHA256 equals a fresh local SHA256 calculation.

- [ ] **Step 2: Verify RED**

Run: `python -m pytest -m network tests/integration/test_real_file_download.py -v`.

Expected: fail on the unimplemented or incomplete verified-download behavior.

- [ ] **Step 3: Add only the compatibility required by the official endpoint**

Keep source URL, expected size, and checksum sourced from the committed registry; do not substitute a local fixture.

- [ ] **Step 4: Verify GREEN**

Run the Step 2 command and require one passing real-network test.

---

### Task 3: Exact registered selection, manifest, and audit CLI

**Files:**
- Create: `configs/download_selection.yaml`
- Create: `src/plantpersulf/download/registered.py`
- Create: `tests/unit/test_registered_download_selection.py`
- Modify: `scripts/fetch_registered_data.py`
- Modify: `src/plantpersulf/cli.py`

**Interfaces:**
- Produces: `download_registered_files(accession, file_classes, ...)`.
- Produces: `audit_downloaded_files(accession, ...)`.
- Produces: TSV columns `dataset_accession`, `repository`, `file_name`, `file_class`, `source_url`, `retrieved_at`, `size_bytes`, `remote_checksum_algorithm`, `remote_checksum`, `path`, `sha256`, `license_or_usage`, `downloader_version`.

- [ ] **Step 1: Write selection and CLI tests**

Assert that PXD006140 `results` resolves exactly the two approved OMSSA names, unknown classes fail, the fetch parser accepts `--accession` plus `--file-class`, and the CLI accepts `audit-files --accession`.

- [ ] **Step 2: Verify RED**

Run: `python -m pytest tests/unit/test_registered_download_selection.py tests/unit/test_task1_cli.py -v`.

Expected: missing registered-download module and unsupported arguments/command.

- [ ] **Step 3: Implement exact selection and manifest updates**

Load the allowlist and source registry, require exact accession/name/URL identity, download only selected remote rows, and atomically rewrite `downloads.tsv`. Reuse an existing local file only after its registered SHA256, size, and source identity pass audit.

- [ ] **Step 4: Implement script and CLI dispatch**

Keep `--metadata-only` behavior intact. Add the authorized download mode and `audit-files`, both returning JSON summaries and nonzero exit on any mismatch.

- [ ] **Step 5: Verify GREEN**

Run the Step 2 command and require all tests to pass.

---

### Task 4: Real download, complete verification, and independent commit

**Files:**
- Create during execution: `data/registry/downloads.tsv`
- Create locally but keep ignored: `data/raw/PXD006140/omssa.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt`
- Create locally but keep ignored: `data/raw/PXD006140/omssa.ne.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt`

- [ ] **Step 1: Run the authorized acquisition**

Run: `python scripts/fetch_registered_data.py --accession PXD006140 --file-class metadata,results`.

Expected: exactly two result downloads plus validation of the existing metadata cache; no other PXD006140 file is fetched.

- [ ] **Step 2: Audit downloaded files**

Run: `python -m plantpersulf.cli audit-files --accession PXD006140`.

Expected: both local result files match official size/SHA1 and registered local SHA256.

- [ ] **Step 3: Run all completion gates**

Run Task 2 unit and network tests, all fast/scientific/release tests, `ruff check .`, `mypy src/plantpersulf scripts`, `python -m plantpersulf.cli audit-registry`, and `git diff --check` excluding intentional terminal tabs in TSV rows.

- [ ] **Step 4: Inspect provenance and scope**

Confirm exactly two downloaded result rows, two ignored raw files, no Task 3 paths, and no synthetic biology, labels, metrics, splits, or model artifacts.

- [ ] **Step 5: Commit independently**

Commit message: `feat: add checksummed real-data download pipeline`.
