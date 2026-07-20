# Task 3 Provenance-Locked Real Fixtures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate deterministic, tests-only fixtures from three already registered real public sources without modifying biological values or performing new downloads.

**Architecture:** `scripts/build_real_fixtures.py` owns a fixed mapping from the three approved accessions to their audited source files and extraction methods. It verifies every source through the existing provenance registry, writes byte-preserving fixture output and a deterministic JSON manifest through temporary siblings, and exposes small functions that scientific tests can call directly.

**Tech Stack:** Python 3.10 standard library, existing PlantPersulf provenance hashing/audit code, pytest, Ruff, mypy.

## Global Constraints

- Execute Task 3 only; do not download data, parse scientific fields, or begin Task 4.
- Use PXD006140 source SHA256 `f8d052626f9f792495c785b7b42d657e2effb3a982c2f8ba561c7c320f37bdcb` from `data/registry/downloads.tsv`.
- Use PXD051570 source SHA256 `84a5d83e099910a780721a3984c753a35f1e25f781675f5a2e5097174a336168` and GSE163745 source SHA256 `20f52bd62689003b8d7185167c6e18ec70b9e1178da3342134a447c138c5f228` from `data/registry/files.tsv`.
- Preserve PXD006140 header plus the next 32 physical lines as original bytes; copy both metadata caches byte-for-byte.
- Store fixtures only under `tests/fixtures/real`; mark them `tests_only` and never route them into training or scientific outputs.
- Every failure is fatal. Never substitute another file, accession, generated record, or placeholder.
- Complete Task 3 as one implementation commit after the already approved design commit.

---

### Task 1: Write the complete RED scientific contract

**Files:**
- Create: `tests/scientific/test_fixture_provenance.py`

**Interfaces:**
- Consumes later interface: `FixtureSource`, `build_fixture`, and `build_real_fixtures` from `scripts.build_real_fixtures`.
- Establishes the exact output and rejection behavior before production code exists.

- [ ] **Step 1: Write the manifest completeness test**

Add a test that finds `tests/fixtures/real/**/source_manifest.json`, requires exactly the accessions `PXD006140`, `PXD051570`, and `GSE163745`, and checks:

```python
assert manifest["schema_version"] == 1
assert manifest["accession"]
assert len(manifest["source_sha256"]) == 64
assert len(manifest["fixture_sha256"]) == 64
assert manifest["extraction_command"] == (
    "python scripts/build_real_fixtures.py "
    "--accessions PXD006140,PXD051570,GSE163745"
)
assert manifest["biological_values_modified"] is False
assert manifest["use"] == "tests_only"
```

- [ ] **Step 2: Write byte-fidelity tests**

For PXD006140, open the registered source in binary mode, read exactly 33 physical lines, and compare those bytes with `omssa_head_32.csv`. For PXD051570 and GSE163745, compare the fixture bytes with their complete registered cache bytes. Recompute every source and fixture SHA256 and compare with its manifest.

- [ ] **Step 3: Write deterministic regeneration test**

Import `build_real_fixtures`, build all three accessions into a pytest temporary output root twice, and compare the complete relative-path-to-SHA256 mapping after each run:

```python
build_real_fixtures(APPROVED_ACCESSIONS, output_root=output_root)
first = tree_hashes(output_root)
build_real_fixtures(APPROVED_ACCESSIONS, output_root=output_root)
second = tree_hashes(output_root)
assert first == second
```

- [ ] **Step 4: Write fail-closed source test**

Create only non-biological marker bytes and a registry with the wrong SHA256 in `tmp_path`. Construct a `FixtureSource` pointing to it, call `build_fixture`, assert `RuntimeError`, and assert that neither the fixture nor its manifest exists.

- [ ] **Step 5: Verify RED**

Run:

```bash
python -m pytest tests/scientific/test_fixture_provenance.py -v
```

Expected: failures because `scripts.build_real_fixtures` and `tests/fixtures/real` do not exist. Confirm failures are missing Task 3 behavior rather than syntax or test-data errors.

---

### Task 2: Implement the audited deterministic fixture builder

**Files:**
- Create: `scripts/build_real_fixtures.py`
- Create during command execution: `tests/fixtures/real/PXD006140/omssa_head_32.csv`
- Create during command execution: `tests/fixtures/real/PXD006140/source_manifest.json`
- Create during command execution: `tests/fixtures/real/PXD051570/PXD051570.metadata.json`
- Create during command execution: `tests/fixtures/real/PXD051570/source_manifest.json`
- Create during command execution: `tests/fixtures/real/GSE163745/GSE163745.metadata.json`
- Create during command execution: `tests/fixtures/real/GSE163745/source_manifest.json`

**Interfaces:**
- Produces immutable `FixtureSource` with accession, repository, registry path, source path, expected source SHA256, fixture name, extraction method, and optional physical-line count.
- Produces `build_fixture(source, output_root=Path("tests/fixtures/real"), extraction_command=EXTRACTION_COMMAND) -> Path`.
- Produces `build_real_fixtures(accessions, output_root=Path("tests/fixtures/real")) -> FixtureBuildSummary`.
- Produces CLI `python scripts/build_real_fixtures.py --accessions <comma-separated>` with JSON output.

- [ ] **Step 1: Define fixed source specifications**

Create `FixtureSource` and an immutable mapping equivalent to:

```python
SOURCES = {
    "PXD006140": FixtureSource(
        accession="PXD006140",
        source_repository="PRIDE",
        source_registry=Path("data/registry/downloads.tsv"),
        source_file=Path(
            "data/raw/PXD006140/"
            "omssa.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt"
        ),
        expected_source_sha256=(
            "f8d052626f9f792495c785b7b42d657e2effb3a982c2f8ba561c7c320f37bdcb"
        ),
        fixture_name="omssa_head_32.csv",
        extraction_method="binary_header_plus_32_lines",
        physical_line_count=33,
    ),
}
```

Add the two metadata specifications with `byte_for_byte_copy` and no line limit. Keep all paths repository-relative and do not discover alternatives dynamically.

- [ ] **Step 2: Implement source verification and extraction**

`build_fixture` must first call:

```python
assert_registered_input(source.source_file, source.source_registry)
actual_source_sha256 = hash_file(source.source_file, "sha256")
if actual_source_sha256 != source.expected_source_sha256:
    raise RuntimeError("fixture source SHA256 differs from approved Task 3 source")
```

For `binary_header_plus_32_lines`, read 33 lines with `readline()` from an `rb` handle and fail if any line is empty before the limit. For `byte_for_byte_copy`, stream chunks from source to temporary fixture without JSON parsing or serialization.

- [ ] **Step 3: Implement deterministic manifest and atomic publication**

Write `source_manifest.json` with exactly the design fields, using `json.dumps(..., indent=2, sort_keys=True)` plus one newline. Calculate `fixture_sha256` from the completed temporary fixture, write a temporary manifest, replace the fixture, replace the manifest, and delete remaining temporary files in `finally`.

Before replacing an existing manifest, parse it and require `schema_version == 1`; malformed or unsupported manifests are fatal. Do not add build timestamps.

- [ ] **Step 4: Implement accession validation and CLI**

Parse comma-separated accessions, strip and uppercase each value, reject empty/duplicate/unsupported inputs, preserve requested order, and build each requested fixture. Return:

```python
@dataclass(frozen=True)
class FixtureBuildSummary:
    accession_count: int
    fixture_count: int
```

Print `json.dumps(asdict(summary), sort_keys=True)` and exit nonzero on any exception.

- [ ] **Step 5: Generate the approved fixtures**

Run:

```bash
python scripts/build_real_fixtures.py \
  --accessions PXD006140,PXD051570,GSE163745
```

Expected JSON:

```json
{"accession_count": 3, "fixture_count": 3}
```

- [ ] **Step 6: Verify GREEN**

Run:

```bash
python -m pytest tests/scientific/test_fixture_provenance.py -v
```

Expected: all Task 3 provenance, fidelity, determinism, and rejection tests pass.

---

### Task 3: Refactor and complete scientific-integrity verification

**Files:**
- Modify only if tests expose duplication: `scripts/build_real_fixtures.py`
- Modify only if assertions need clearer diagnostics: `tests/scientific/test_fixture_provenance.py`

**Interfaces:**
- Preserve all Task 2 public names and manifest fields.

- [ ] **Step 1: Refactor only while GREEN**

Extract shared temporary-file cleanup or binary-copy helpers if the implementation contains duplication. Do not add new source types, extraction modes, accessions, or configuration formats.

- [ ] **Step 2: Re-run deterministic build**

Run the builder twice, then run the deterministic scientific test that compares the complete fixture tree across consecutive builds:

```bash
python scripts/build_real_fixtures.py --accessions PXD006140,PXD051570,GSE163745
python scripts/build_real_fixtures.py --accessions PXD006140,PXD051570,GSE163745
python -m pytest tests/scientific/test_fixture_provenance.py -v
```

- [ ] **Step 3: Run all completion gates**

```bash
python -m pytest tests/scientific/test_fixture_provenance.py -v
python -m pytest tests/unit tests/scientific tests/release -q
python -m pytest -m network tests/integration -q
python -m plantpersulf.cli audit-registry
python -m plantpersulf.cli audit-files --accession PXD006140
ruff check .
mypy src/plantpersulf scripts
git diff --check
```

Require zero failures, zero Ruff errors, zero mypy errors, and successful Task 1/2 audits.

- [ ] **Step 4: Inspect scientific scope**

Confirm exactly three fixture data files and three manifests; verify no files under model, benchmark, split, processed-data, results, or Task 4 parser paths; verify temporary policy-marker bytes are absent from `tests/fixtures/real`.

- [ ] **Step 5: Commit Task 3 implementation**

```bash
git add scripts/build_real_fixtures.py \
  tests/fixtures/real \
  tests/scientific/test_fixture_provenance.py \
  docs/superpowers/plans/2026-07-20-task3-real-fixtures.md
git commit -m "test: add provenance-locked real-data fixtures"
```
