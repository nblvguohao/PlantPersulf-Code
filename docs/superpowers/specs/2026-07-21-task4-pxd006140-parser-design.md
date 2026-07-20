# Task 4 PXD006140 Proteomics Parser Design

## Status and scope

This design covers Task 4 only: parsing the two already registered PXD006140
OMSSA result files into traceable peptide-spectrum-match records and producing
only sequence-verified cysteine coordinates. It does not build the PU benchmark,
freeze splits, extract model features, train a model, or create candidate lists.

Task 4 uses the real PXD006140 inputs already registered by Task 2 and the
tests-only PXD006140 fixture frozen by Task 3. It adds two narrowly scoped
official UniProt FASTA records, Q93VK9 and Q9ZW96, solely to prove the coordinate
verification path against real sequences. Downloading the complete Arabidopsis
reference proteome remains Task 7.

## Observed source facts

- `omssa.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt` is registered in
  `data/registry/downloads.tsv` with SHA256
  `f8d052626f9f792495c785b7b42d657e2effb3a982c2f8ba561c7c320f37bdcb`.
- `omssa.ne.20150821_01_AAroca_TMT6plex.cmpd.mgf.txt` is registered in
  `data/registry/downloads.tsv` with SHA256
  `9bbf2f88f51ab9a80b62315dd2b61f86f0d4290d8457591d28912a252af971ea`.
- The files are comma-delimited OMSSA exports whose headers contain spaces
  after delimiters. Modified residues may be represented with lower-case
  letters in the `Peptide` field, while `Mods` retains the raw modification
  description.
- OMSSA `Start` and `Stop` values are one-based and inclusive for the two
  target records checked against current official UniProt sequences:
  Q93VK9/`ETmASLGLICEK` and Q9ZW96/`RmSCNGCRVLR`.
- The Task 3 fixture also contains a reversed Q9FKH0 record. It is a real decoy
  record and must be excluded with an audit reason, not treated as biology.
- PXD006140 currently has no sample rows in `data/registry/samples.tsv`.
  Therefore Task 4 must report missing sample metadata and must not manufacture
  a sample identifier or condition.
- The OMSSA exports contain search-result PSMs. Their presence alone is not
  sufficient to label a cysteine as experimentally confirmed persulfidation.

## Approaches considered

### Selected: limited, versioned UniProt evidence bundle

Retrieve Q93VK9 and Q9ZW96 FASTA bytes from official UniProt REST endpoints,
store them under `data/registry/cache/uniprot`, and register each file with its
source URL, retrieval time, byte size, sequence version from the FASTA header,
and SHA256 in `data/registry/reference_sequences.tsv`. Scientific tests then
run offline against registered bytes.

This is reproducible and satisfies the Task 4 coordinate test without advancing
the complete proteome acquisition assigned to Task 7.

### Rejected: live UniProt requests during parsing or tests

This would make scientific tests depend on network availability and current
database state. It would also prevent an exact source-byte audit.

### Rejected: complete Arabidopsis reference proteome acquisition

This would give broader coordinate coverage, but it is Task 7 scope and would
substantially expand storage, versioning, and mapping decisions before Task 4
has passed review.

## Registered inputs

Task 4 consumes only inputs that pass `assert_registered_input`:

1. Both PXD006140 OMSSA result files through `data/registry/downloads.tsv`.
2. Q93VK9 and Q9ZW96 FASTA files through
   `data/registry/reference_sequences.tsv`.
3. The PXD006140 tests-only fixture through its Task 3 source manifest.

The dedicated reference registry has exact columns `study_accession`,
`protein_accession`, `repository`, `source_url`, `retrieved_at`,
`sequence_version`, `path`, `size_bytes`, `sha256`, and `scientific_use`. The two
rows use `study_accession=PXD006140`, `repository=UniProt`, and
`scientific_use=coordinate_validation_only`. This describes their role as
supporting reference inputs without registering them as independent proteomics
studies or modifying the Task 1 metadata registry. Their SHA256 values are
calculated from the downloaded response bytes; a file is unusable until the
corresponding registry row exists and its SHA256 audit passes.

No unregistered fallback sequence, copied sequence literal, inferred residue,
or alternate accession may be substituted after a retrieval or mapping failure.

## Components and interfaces

### `proteomics/metadata.py`

Defines immutable records for:

- `ProteomicsSource`: study accession, source path, source SHA256, parser name,
  and evidence scope;
- `PeptideSpectrumMatch`: one source row and one exact protein mapping;
- `SiteEvidence`: one sequence-verified cysteine coordinate;
- `ParseConflict`: the original source identity, reason code, and unmodified raw
  values needed for audit;
- `ParseSummary`: deterministic counts for parsed, verified, excluded,
  conflicted, and metadata-missing records.

These types contain strings, integers, optional decimal strings, and provenance
fields only. They do not assign benchmark labels.

### `proteomics/peptide_parser.py`

Provides a streaming OMSSA CSV parser. It:

- opens registered files with `utf-8-sig`, `newline=""`, and
  `skipinitialspace=True`;
- validates the exact required source columns before yielding any records;
- preserves the original `Peptide`, `Mods`, defline, spectrum number, start,
  stop, and source filename;
- derives an unmodified peptide sequence only by upper-casing the real peptide
  characters; it never adds, removes, or guesses residues;
- retains repeated spectrum/peptide mappings as separate rows;
- rejects malformed numeric coordinates into the conflict stream rather than
  coercing them;
- emits no output when the header is invalid or the source SHA256 audit fails.

The parser supports OMSSA in Task 4. Other search-engine formats fail closed.
Raw modification names are preserved so later versioned normalization can map
different engines without silently declaring equivalence.

### `proteomics/site_normalizer.py`

Provides conservative protein and site normalization:

- parses `sp|ACCESSION|...` and `tr|ACCESSION|...` target deflines;
- classifies `rev_`, reversed, and configured contaminant prefixes before any
  biological record is emitted;
- preserves an exact UniProt accession when present;
- does not collapse an isoform suffix to its canonical parent; isoform mappings
  without an exact registered sequence become conflicts;
- checks `stop - start + 1 == len(peptide_sequence)`;
- computes one-based peptide cysteine position as `offset + 1`;
- computes one-based protein cysteine position as `start + offset`;
- verifies the complete peptide slice and cysteine residue against the exact
  registered FASTA sequence before emitting `SiteEvidence`;
- reports `sequence_unavailable`, `peptide_length_conflict`,
  `peptide_sequence_conflict`, `unparseable_accession`, and
  `isoform_sequence_unavailable` without attempting correction.

If a source is declared `protein_level_only`, the parser may emit its PSM or
protein mapping but `site_normalizer` emits no site evidence. The test for this
policy reuses a real PXD006140 row and changes only the software evidence-scope
setting; it does not create a biological record.

## Scientific semantics

Task 4 uses `evidence_level=psm_coordinate_only` for verified coordinates. This
means that a real search-result peptide contains a sequence-verified cysteine at
the reported position. It does not mean the cysteine is a persulfidation-positive
site.

No output contains `negative`. Records without verified site evidence remain
unverified or protein-level evidence; they are not converted into negative
examples. Task 5 must apply a separately reviewed experimental-evidence policy
before assigning any positive or unlabeled benchmark status.

Quantitative fields remain blank because the selected OMSSA exports do not
provide an audited abundance value or unit. Missing values are not filled with
zero. `sample_id` remains blank because no registered PXD006140 sample mapping
exists; the omission is repeated in `missing_metadata.tsv`.

## Outputs

`python -m plantpersulf.cli parse-proteomics --accession PXD006140` writes an
atomic, deterministic tree under the ignored directory
`data/interim/PXD006140/proteomics_parser_v1/`:

- `psms.tsv`: every valid target PSM/protein mapping, including records that
  cannot yet be sequence verified;
- `sites.tsv`: one row per cysteine that passes complete sequence verification;
- `conflicts.tsv`: coordinate, accession, isoform, sequence, and schema
  conflicts with original values and reason codes;
- `excluded.tsv`: decoy and contaminant rows with explicit reasons;
- `missing_metadata.tsv`: missing sample condition and other metadata gaps;
- `manifest.json`: parser version, ordered input paths and SHA256 values,
  reference-sequence SHA256 values, output SHA256 values, schema version, and
  deterministic row counts.

The site table contains the required Task 4 fields:

```text
study_accession
sample_id
source_file
spectrum_id
peptide_sequence
modified_sequence
protein_accession_raw
protein_accession_canonical
cys_position_in_peptide
cys_position_in_protein
modification_name_raw
evidence_level
quant_value
quant_unit
parser_version
source_sha256
```

`psms.tsv` additionally retains raw start, stop, defline, search-engine
accession, and mapping status. Conflict and exclusion tables retain enough raw
identity to trace every decision to the source row.

Files are written to temporary siblings first. A failed source, schema,
registry, parse, or hash check publishes no partial output. Re-running the same
parser version with identical registered inputs produces byte-identical files
and manifest; no run timestamp is embedded in deterministic outputs.

## CLI audit behavior

`python -m plantpersulf.cli audit-sites --accession PXD006140`:

- verifies the parser manifest schema and accession;
- rehashes every source, reference, and generated output;
- verifies required columns and parser version;
- rejects any site lacking a registered source SHA256;
- rejects any emitted site whose peptide or cysteine does not match its
  registered reference sequence;
- rejects forbidden `negative` evidence;
- reports counts for PSMs, verified sites, conflicts, excluded rows, missing
  sequences, and missing sample metadata.

Conflicts and missing sequences are reportable limitations, not silent success.
The audit exits nonzero for corrupt provenance or an invalid emitted site, but
does not fail merely because an explicitly reported record lacks a sequence.

## TDD sequence

1. RED: real fixture parsing preserves peptide and modification bytes, preserves
   repeated mappings, excludes the real decoy with a reason, and reports missing
   sample metadata.
2. GREEN: implement only the OMSSA parser and immutable record types needed for
   those tests.
3. RED: Q93VK9 and Q9ZW96 coordinates must match registered official sequences;
   a real record without an available registered sequence must enter conflicts;
   a real row under `protein_level_only` must not create a site.
4. GREEN: implement minimal FASTA loading and conservative site normalization.
5. RED: CLI parsing and auditing must be deterministic, fail closed on a changed
   source SHA256, and publish no partial tree.
6. GREEN: implement the two CLI commands, deterministic writers, manifest, and
   audit.
7. REFACTOR: remove duplication only after all Task 4 tests are green.

Tests may create temporary directories, malformed headers, and policy-only
configuration states. Biological rows used by tests must come unchanged from
the registered PXD006140 fixture. No test may invent a peptide, protein,
coordinate, modification, measurement, or label.

## Completion gates

Task 4 is complete only when all of the following pass:

```text
python -m pytest tests/unit/test_peptide_parser_real_fixture.py \
  tests/scientific/test_site_coordinates_match_sequence.py -v
python -m plantpersulf.cli parse-proteomics --accession PXD006140
python -m plantpersulf.cli audit-sites --accession PXD006140
python -m pytest tests/unit tests/scientific tests/release -q
python -m pytest -m network tests/integration -q
python -m plantpersulf.cli audit-registry
python -m plantpersulf.cli audit-files --accession PXD006140
ruff check .
mypy src/plantpersulf scripts
git diff --check
```

The generated interim tree is inspected for deterministic hashes and scientific
scope but is not committed. The implementation commit contains source code,
tests, registered UniProt cache bytes, the dedicated reference-sequence
registry, and the implementation plan. Task 5 is not started.

## Known limitations and Task 5 gate

- Sequence verification in Task 4 is deliberately limited to Q93VK9 and
  Q9ZW96. All other accessions remain explicitly `sequence_unavailable` until a
  later approved reference-sequence acquisition task.
- PXD006140 lacks a registered sample-to-file mapping, so sample conditions
  remain missing.
- The OMSSA exports do not by themselves establish site-specific
  persulfidation evidence or audited quantification.
- Task 4 produces no benchmark labels, model inputs, performance values, plots,
  candidates, or biological conclusions.

Therefore a successful Task 4 implementation does not automatically satisfy the
Task 5 scientific start gate. Task 5 requires separate review of the Task 4
commit and confirmation that available experimental evidence can support
positive-versus-unlabeled labels without converting nondetection into a
negative.
