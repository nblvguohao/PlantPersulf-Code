# Task 5 Real-File Evidence Content Audit Design

## Scope

This is cycle 2 of the Task 5 evidence preflight. It interprets the registered
Stage A files structurally and conservatively. It does not build benchmark
labels, splits, features, models, metrics, figures, candidates, or biological
conclusions.

The cycle may emit audited `identification_only`, `protein_level_only`, or
`unresolved` records. It may emit `site_ms` only if an exact site record and a
registered experimental-method source jointly satisfy the reviewed evidence
rule. Nondetection and missing evidence never become negative.

## Registered real-file observations

All observations below were made from the exact local files in
`data/registry/downloads.tsv` after repository SHA1 and local SHA256 audit.

### PXD024061 checksum inventory

- File: `checksum.txt`
- SHA256: `38bb488f0c30ae11b5f8f50d7eb35350ffa19649a4c30825047b9494167073d3`
- Encoding/delimiter: UTF-8 with optional BOM, tab-delimited path and SHA1
- Structure: one `#Checksum File` header and 26 checksum rows
- Limitation: the Stage A bytes contain no processed biological result; the
  2.42 GB `txt_persulfproject.zip` remains deferred.

### PXD035795 SDRF

- File: `SDRF.txt`
- SHA256: `335934a0c9b307c99cbdd1c6be05d9da3813d4ee0d9645006f3cbf528a9c5fd1`
- Encoding/delimiter: UTF-8 with optional BOM, tab-delimited
- Structure: 26 columns and 6 assay rows
- Design: six source/assay names, three registered biological-replicate values,
  label-free acquisition, and two treatment values representing active and
  non-photorespiratory conditions
- File mapping: each assay names one RAW and one MGF-like value; RAW filenames
  resolve exactly in the official registry, while encoded MGF-like values do
  not resolve exactly and must enter mapping conflicts
- Repeated headers: five `comment[modification parameters]` columns occur at
  distinct positions and must be parsed positionally rather than collapsed by
  a dictionary reader
- Named search modifications include NBF_N, NBF_K, NBF_C, DCP, and dcp-ac.
  These names describe search parameters and are not by themselves proof of a
  persulfidated cysteine.

### PXD035795 peptide and protein CSV

- `peptide.csv` SHA256:
  `63f2df50a2748d5614a8cfad16cf19c4d7dfac4825d67136382b1fc1fb62ab47`
- `proteins.csv` SHA256:
  `cd816699fd3427634c347bb0eb52ee6392c47218bc2bcee5cd55ef367e2115ba`
- Encoding/delimiter: UTF-8 with optional BOM, comma-delimited
- `peptide.csv`: 19 columns and 9326 data rows, including raw peptide,
  accession, PTM, six replicate-area fields, and two group-level fields
- `proteins.csv`: 28 columns and 1461 data rows, including protein accession,
  peptide counts, PTM summary, replicate areas, and group-level areas
- PTM fields are free-text semicolon-separated summaries with many distinct
  values. A CSV PTM value has no independently verified residue location and
  therefore cannot alone produce `site_ms`.
- Protein rows cannot be promoted to site-level evidence.

### PXD035795 mzIdentML

- File: `peptides_1_1_0.mzid.gz`
- Compressed size/SHA256: 184107 bytes,
  `62105dfde42d9675bc5dc7c83969d971eec57c9e02983659cd2df96ec4d1b5fe`
- Decompressed size: 1967088 bytes, exactly matching the PRIDE metadata size
- Namespace/version: `http://psidev.info/psi/pi/mzIdentML/1.1`, version `1.1.0`
- Structural counts: 5 DBSequence, 11409 Peptide, 95 PeptideEvidence, 661
  SpectrumIdentificationResult, 699 SpectrumIdentificationItem, and 124
  Modification elements
- All identified modifications have a location. Observed modification residues
  are blank, C, or K.
- The only modification CV term is `MS:1001460` / `unknown modification`.
- Search masses include 163.0012 on blank/C/K, 168.0786 on C, 196.08 on C,
  and 394.1557 on C. Identified modification elements contain the first three
  masses but not 394.1557.
- Approximate mass agreement with an SDRF name is not sufficient evidence. The
  mapping must be supported by a registered method or search-parameter source.

### PXD035795 and PXD039999 checksums

- PXD035795 checksum SHA256:
  `81033d0e37a129ed6d59601d8662e964af94a4f66825571fa0a0ed4145fd137c`;
  22 tab-delimited SHA1 rows
- PXD039999 checksum SHA256:
  `429e956508a6166ca8c51fdbe886aba8104f5c2550e9a99d88f8f268ec53eac0`;
  13 tab-delimited SHA1 rows
- PXD039999 Stage A contains no processed result or sample-design file.

## Architecture

### Method-source registry

Before assigning any site-level evidence, register the PXD035795 publication
method and any official supplement used to define the chemical tags. A new
`data/registry/evidence_methods.tsv` records study accession, DOI or repository
identifier, official URL, local path, retrieval time, size, SHA256, method
scope, and data-use status. Public supplementary files are downloaded through
the same temporary-file and hash-first policy as repository files.

If no registered method source explicitly supports a DCP/NBF mass-to-evidence
mapping, all modified peptide candidates remain `unresolved` or
`identification_only`.

### SDRF reader

Read the header and rows as positional lists. Require exactly 26 columns and six
rows for the registered version. Preserve duplicate modification headers with
their column indices. Emit exact source, assay, replicate, treatment, label,
instrument, and data-file strings.

Split the two comma-separated data-file values without URL decoding or filename
rewriting. Resolve each component exactly against `files.tsv`; RAW values map to
assays, while unmatched MGF-like values enter `mapping_conflicts.tsv`. A later
reviewed mapping table may resolve them, but the parser cannot silently remove
percent characters or spaces.

### CSV readers

Require the exact registered headers and row counts. Preserve peptide,
accession, PTM, quality/significance, replicate areas, and group areas as raw
strings. Split PTM summaries only to enumerate declared names; do not infer
site positions from them. Protein rows are reported as `protein_level_only`
only if the registered experimental method supports that interpretation;
otherwise they remain `identification_only`.

### mzIdentML reader

Stream the gzip XML under the exact mzIdentML 1.1 namespace. Preserve spectrum
identity, peptide identity, DBSequence/PeptideEvidence references, modification
location, residue, mass delta, CV accession/name, rank, and pass-threshold
status. Resolve references exactly and retain every peptide-to-protein mapping.
Malformed references, absent sequences, or out-of-range locations are explicit
conflicts.

The reader emits `candidate_modification` records with `evidence_class` set to
`unresolved`. A separate rule joins a candidate to a registered method mapping.
Only an exact rule for the same study, residue, mass tolerance, chemical tag,
and experimental scope may upgrade the record to `site_ms`.

### Quantification and sample mapping

Replicate-area values remain raw strings until each CSV column maps to one SDRF
assay/treatment and a unit/normalization description is registered. Missing or
ambiguous mappings leave quantification status unresolved; they are never
filled with zero. Group fields cannot replace replicate-level mapping.

## Evidence policy

For PXD035795, `site_ms` requires all of:

1. a registered mzIdentML spectrum identification that passes its reported
   threshold;
2. an exact peptide, protein mapping, modification location, cysteine residue,
   and source SHA256;
3. a registered PXD035795 method source that explicitly maps the observed
   chemical tag/mass to persulfidation-site evidence;
4. no unresolved peptide/protein/reference conflict;
5. an evidence locator sufficient to reproduce the record.

Failure of any item results in `unresolved` or `identification_only`, never a
negative. `protein.csv` cannot satisfy these site-level conditions.

PXD006140 retains `psm_coordinate_only`; PXD024061 and PXD039999 remain
`unresolved` in this cycle because Stage A contains no inspected processed
result capable of establishing a site.

## Outputs

Cycle 2 updates the deterministic `data/interim/evidence_preflight_v1/` tree and
adds:

- populated `evidence_records.tsv` for content-audited records;
- `candidate_modifications.tsv` for unresolved localized modifications;
- `sample_file_mapping.tsv` and `mapping_conflicts.tsv`;
- `content_schema.tsv` with input SHA256, encoding, delimiter/namespace,
  headers, row/element counts, and parser version.

The manifest changes to `content_audit_complete=true` only when every selected
Stage A design/result file has a supported parser or an explicit unsupported
reason. Integrity flags remain `labels_created=false`,
`nondetection_labeled_negative=false`, and
`biological_values_modified=false`.

## TDD sequence

1. RED: registered method evidence is absent, so no candidate modification may
   become `site_ms`.
2. GREEN: register and audit the official method source without yet changing
   evidence classes.
3. RED/GREEN: positional SDRF parsing preserves all duplicate headers and
   exposes exact mapping conflicts.
4. RED/GREEN: exact CSV schemas and row counts are preserved with no inferred
   sites.
5. RED/GREEN: mzIdentML reference and modification coordinates are parsed from
   the registered gzip bytes; all candidates initially remain unresolved.
6. RED/GREEN: only reviewed method mappings may upgrade an exact candidate to
   `site_ms`; tests require failure when the method source or mapping is absent.
7. GREEN: deterministic publication, source/output hash audit, full regression,
   and scientific-integrity scan.

## Completion gate

The content audit is complete when all selected PXD035795 design/result files
are parsed or explicitly rejected, every emitted record has registered source
and method hashes, mapping conflicts are visible, and repeated builds are
byte-identical.

Task 5 label construction remains blocked if no record meets the complete
`site_ms`, `site_mutagenesis`, or `site_biochemical` contract. Even if
PXD035795 yields site-level records, a second independent site-level study is
still the recommended practical prerequisite for downstream leave-study-out
evaluation.

