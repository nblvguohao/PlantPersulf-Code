# Task 3 Provenance-Locked Real Fixtures Design

## Objective

Build small, deterministic test fixtures exclusively from already downloaded or cached real public data. Every fixture must retain a complete link to its registered source and must be limited to software testing; no fixture may enter training, scientific analysis, figures, or conclusions.

## Approved approach

Use a hybrid source strategy:

- PXD006140 uses an exact byte-preserving prefix of a downloaded and checksummed OMSSA result file.
- PXD051570 uses a byte-for-byte copy of its registered official iProX metadata cache.
- GSE163745 uses a byte-for-byte copy of its registered official GEO metadata cache.

This approach satisfies the documented three-accession acceptance command without acquiring new files during Task 3. The two metadata fixtures are metadata-only and must never be described as raw experimental data.

## Alternatives considered

1. Copy metadata caches for all three accessions. This is smallest but provides no real peptide-result fixture for the next parser task.
2. Download PXD051570 and GSE163745 experimental source files first. This could provide richer fixtures but would reopen data acquisition outside Task 3.
3. Use the approved hybrid approach. It stays inside the Task 3 boundary and provides one genuine result fixture plus two genuine metadata fixtures.

## Source files

### PXD006140

- Source registry: `data/registry/downloads.tsv`
- Source file: `data/raw/PXD006140/omssa.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt`
- Registered SHA256: `f8d052626f9f792495c785b7b42d657e2effb3a982c2f8ba561c7c320f37bdcb`
- Extraction: copy the CSV header and the next 32 newline-terminated physical lines as raw bytes, preserving delimiters, quoting, biological text, numeric values, and line endings.
- Fixture name: `omssa_head_32.csv`

The builder must not interpret target/decoy status, peptide sequence, modification names, scores, or accessions.

### PXD051570

- Source registry: `data/registry/files.tsv`
- Source file: `data/registry/cache/iprox/PXD051570.json`
- Registered SHA256: `84a5d83e099910a780721a3984c753a35f1e25f781675f5a2e5097174a336168`
- Extraction: byte-for-byte copy.
- Fixture name: `PXD051570.metadata.json`

This fixture is official repository metadata, not proteomic raw data or a processed result table.

### GSE163745

- Source registry: `data/registry/files.tsv`
- Source file: `data/registry/cache/geo/GSE163745.json`
- Registered SHA256: `20f52bd62689003b8d7185167c6e18ec70b9e1178da3342134a447c138c5f228`
- Extraction: byte-for-byte copy.
- Fixture name: `GSE163745.metadata.json`

This fixture is official GEO Series/Sample metadata, not an expression matrix or sequencing reads.

## Output layout

```text
tests/fixtures/real/
├── PXD006140/
│   ├── omssa_head_32.csv
│   └── source_manifest.json
├── PXD051570/
│   ├── PXD051570.metadata.json
│   └── source_manifest.json
└── GSE163745/
    ├── GSE163745.metadata.json
    └── source_manifest.json
```

The builder is `scripts/build_real_fixtures.py` and accepts:

```bash
python scripts/build_real_fixtures.py \
  --accessions PXD006140,PXD051570,GSE163745
```

Unknown, duplicate, empty, or unsupported accessions are fatal.

## Manifest contract

Each `source_manifest.json` contains:

```text
schema_version
accession
source_repository
source_registry
source_file
source_sha256
fixture_file
fixture_sha256
extraction_method
extraction_command
biological_values_modified
use
```

Required values and rules:

- `schema_version` is `1`.
- `source_registry` and `source_file` are repository-relative paths.
- `source_sha256` must equal the current registered digest before extraction.
- `fixture_sha256` is calculated from the completed temporary fixture before publication and must still match after publication.
- `extraction_method` is `binary_header_plus_32_lines` for PXD006140 and `byte_for_byte_copy` for the metadata fixtures.
- `extraction_command` is the exact documented builder command above.
- `biological_values_modified` is the JSON boolean `false`.
- `use` is `tests_only`.

Timestamps are excluded so repeated builds from unchanged sources are byte-for-byte deterministic.

## Data flow

1. Parse and normalize the requested accession list without reordering it.
2. Resolve each accession through the fixed Task 3 source specification.
3. Call the existing registered-input audit using the appropriate source registry.
4. Recompute and compare the source SHA256 with the fixed expected Task 3 source digest.
5. Extract to a temporary sibling file and prepare a temporary manifest.
6. Compute the fixture SHA256 and place it in the temporary manifest.
7. Atomically replace the fixture and then atomically replace its manifest. Any interruption between replacements is detectable because the retained manifest cannot validate the new fixture.
8. Return a machine-readable summary containing fixture count and accession count.

No network request is permitted in the fixture builder.

## Failure behavior

The command fails closed when:

- a requested accession is unsupported;
- a source path or registry is missing;
- the source is not registered;
- source size or SHA256 differs from provenance;
- fewer than 32 PXD006140 data records are available;
- a temporary write, hash, or atomic replacement fails, in which case the command is fatal and the next audit must reject any fixture/manifest mismatch;
- an existing manifest uses an unsupported schema.

Failures before replacement publish neither temporary output. An interruption between the fixture and manifest replacements may leave a detectable mismatch, which remains unusable because provenance tests and audits must reject it.

An already valid fixture may be replaced only by deterministic regeneration from the same audited source. The builder must not fall back to another file, another accession, or generated content.

## Testing strategy

The first RED test is `tests/scientific/test_fixture_provenance.py` and requires at least one real manifest plus complete accession, source SHA256, fixture SHA256, extraction command, `biological_values_modified: false`, and `use: tests_only` fields.

Additional tests verify:

- all three approved accessions have exactly one manifest and one fixture;
- current source and fixture hashes match each manifest;
- PXD006140 contains exactly one header line plus 32 original physical data lines;
- both metadata fixtures are byte-for-byte equal to their registered sources;
- a temporary unregistered source is rejected and produces no fixture;
- two consecutive builds produce identical fixture and manifest hashes.

Tests that exercise rejection may use non-biological marker bytes only in pytest temporary directories. Those bytes must never be written under `tests/fixtures/real`.

## Non-goals

Task 3 does not:

- download new scientific files;
- parse or normalize peptides, proteins, modifications, expression values, or sample conditions;
- remove decoys or contaminants;
- create labels, benchmark tables, splits, metrics, plots, or conclusions;
- use fixtures in training or publication outputs;
- start Task 4.

## Acceptance

```bash
python scripts/build_real_fixtures.py \
  --accessions PXD006140,PXD051570,GSE163745
pytest tests/scientific/test_fixture_provenance.py -v
```

Completion additionally requires all fast tests, scientific-integrity tests, release tests, Ruff, mypy, Task 1 registry audit, and Task 2 downloaded-file audit to pass.
