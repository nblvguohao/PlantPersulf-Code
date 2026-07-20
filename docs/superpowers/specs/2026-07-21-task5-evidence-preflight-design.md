# Task 5 Evidence Preflight Design

## Status and scope

This is a prerequisite data-completeness task for Task 5. It does not create
benchmark labels, positive or negative examples, train/validation/test splits,
features, models, metrics, plots, candidates, or biological conclusions.

The scope is limited to the four registered Arabidopsis persulfidation studies:
PXD006140, PXD024061, PXD035795, and PXD039999. Every scientific input must be
an already registered local file or an official public file that is downloaded,
hashed with SHA256, and registered before inspection.

## Scientific question

Determine, without inventing evidence, whether the currently public files can
support a site-level positive-versus-unlabeled benchmark. In particular, the
preflight must distinguish:

- `site_ms`: an official result explicitly localizes experimental evidence to
  a cysteine site;
- `site_mutagenesis`: a published or formally supplied site-mutagenesis result;
- `site_biochemical`: a published or formally supplied site-specific
  biochemical result;
- `protein_level_only`: persulfidation evidence is assigned only to a protein;
- `identification_only`: a peptide/protein identification lacks explicit
  persulfidation-site evidence;
- `unresolved`: required files or metadata are absent or cannot be mapped.

`psm_coordinate_only` from Task 4 is not promoted to any positive evidence
class. Nondetection, missing sequence, missing quantification, and conflicts are
never converted into experimental negatives.

## Approaches considered

### Selected: evidence-first inventory across all four studies

Refresh official project metadata, classify public files by likely evidence
role, download only small high-value design/result/checksum files, and publish a
deterministic preflight report. Large archives remain pending until their need
and expected evidence value are documented.

This approach preserves cross-study visibility while minimizing bandwidth and
prevents a large download from being mistaken for scientific progress.

### Rejected: PXD006140-only benchmark construction

PXD006140 currently lacks registered sample mapping and its OMSSA exports do not
by themselves establish site-specific persulfidation. Treating the three
sequence-verified cysteine coordinates as positives would violate the evidence
policy and would not support downstream leave-study-out validation.

### Rejected: bulk download of all raw and search files

The registered inventory includes multi-gigabyte raw files and a 2.42 GB search
archive. Downloading them before evidence triage would consume substantial time
and storage without proving that a site-level result exists.

## Inputs and staged acquisition

The immutable starting inputs are `data/registry/datasets.tsv`, `files.tsv`,
`publications.tsv`, their registered official metadata caches, and the two
already downloaded PXD006140 OMSSA result files.

Stage A acquires only these official PRIDE files:

- PXD035795: `SDRF.txt`, `peptide.csv`, `proteins.csv`,
  `peptides_1_1_0.mzid.gz`, and `checksum.txt`;
- PXD024061: `checksum.txt`;
- PXD039999: `checksum.txt`;
- PXD006140: no additional file, because the two approved OMSSA result files
  and official metadata cache are already downloaded and registered.

Stage A must not acquire PXD024061 `txt_persulfproject.zip`, raw mass spectra,
peak files, search FASTA files, or PXD039999 `proteinSeq.txt`. A later
publication-supplement acquisition is a separate subproject and is not silently
folded into this PRIDE-file preflight. Each excluded or large file receives an
explicit decision row containing its size, reason deferred, expected evidence
value, and the evidence needed to authorize its later download.

## Components and outputs

### Selection policy

`configs/evidence_preflight_v1.yaml` pins the four accessions, accepted evidence
classes, a large-file threshold of 104857600 bytes, review-priority rules, and
prohibited label terms. `configs/download_selection.yaml` remains the exact
download authority and is extended only with the Stage A file names and classes.

### Inventory builder

`src/plantpersulf/evidence/preflight.py` reads the official registry and emits a
deterministic inventory. It never infers an evidence class from a filename
alone. Filename/category rules may assign `review_priority`, but
`evidence_class` remains `unresolved` until content-level checks succeed.

### Registered downloader

The existing fail-closed downloader and central `data/registry/downloads.tsv`
are reused. A selected remote file must resolve exactly once in `files.tsv`;
download to a temporary sibling; verify a repository checksum when one is
available; compute SHA256; then add its exact path, source URL, retrieval time,
size, and SHA256 to `downloads.tsv`. Failed or partial downloads are not
registered.

The download registry preserves both `registry_size_bytes` from the official
PRIDE metadata response and the actual transferred `size_bytes`. When an
official repository checksum is present, the checksum is the byte-identity
gate and a size difference is retained rather than silently corrected. When no
repository checksum exists, the registry size remains mandatory. This handles
compressed PRIDE objects whose API size describes a different representation
without weakening provenance.

### Content evidence audit

Format-specific readers may inspect SDRF, CSV/TSV, mzIdentML, and checksum text.
They preserve raw identifiers and values. Every
classification row records the source file SHA256, source row/sheet locator,
method, condition mapping status, quantification status, site-localization
status, and a reason code.

No parser may assign `site_ms` merely because an identified peptide contains a
cysteine or because its parent study concerns persulfidation.

### Deterministic report

The ignored directory `data/interim/evidence_preflight_v1/` contains:

- `study_inventory.tsv`: one row per study with registered metadata status;
- `file_decisions.tsv`: one row per official file, including selected,
  downloaded, deferred, or rejected status and reason;
- `evidence_records.tsv`: only content-audited evidence records;
- `metadata_gaps.tsv`: missing sample, condition, quantification, mapping, or
  site-localization facts;
- `large_file_queue.tsv`: explicit decisions for files outside Stage A;
- `manifest.json`: config hash, registry hashes, input/output hashes, parser
  versions, deterministic counts, and integrity flags.

The manifest must state `labels_created=false`,
`nondetection_labeled_negative=false`, and
`biological_values_modified=false`.

## TDD and failure behavior

Implementation is split into two evidence-preserving cycles. Cycle 1 freezes
the exact Stage A selection, downloads and registers the real files, and emits
the metadata-level inventory. Cycle 2 begins only after those registered bytes
are available; its RED tests are written against the actual file schemas and
then add the content evidence audit. This ordering prevents tests or parsers
from assuming biological columns that have not been observed.

RED tests use only registered real metadata caches and provenance-locked real
fixtures. They must first demonstrate missing selection, missing inventory
behavior, and rejection of unregistered or hash-mismatched inputs.

GREEN implementation is the minimum code needed to:

1. classify the existing official inventory without labels;
2. select only exact Stage A allow-listed files;
3. download and register selected files atomically;
4. audit supported content conservatively;
5. publish byte-identical outputs for identical registered inputs.

Any ambiguous mapping, unsupported file format, missing sample design, missing
quantification, or absent site-localization evidence becomes a report row. It
does not become a guessed value and does not abort unrelated study inventory.
Registry, checksum, schema, or manifest corruption is fatal.

## Completion and Task 5 gate

This preflight is complete when all four studies have deterministic inventory
rows, every Stage A selection is downloaded or has an explicit external
blocker, every inspected file is registered by SHA256, and the report states
which studies provide site-level, protein-level-only, identification-only, or
unresolved evidence.

Task 5 remains blocked unless at least one traceable site-level positive source
exists. For a scientifically useful downstream leave-study-out design, the
recommended practical gate is two independent studies with site-level positive
evidence. Known mechanism cards remain excluded from training and cannot be
counted as an independent discovery study.
