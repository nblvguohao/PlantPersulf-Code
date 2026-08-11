# Blind candidate-validation protocol

This is a policy template, not a candidate list or experimental result. A real
release requires a signed, human-readable statistical analysis plan (SAP) and
the matching machine-readable frozen manifest.

## Blinding and ownership

Only the key custodian receives `private_key.tsv`; assay operators receive
`blind_table.tsv`, which exposes the blind ID, site key, protein accession and
matched stratum but never rank, score, tier or arm. The key is unblinded only
after the full raw-file package and the pre-specified primary analysis are
returned.

## Arms and assay rules

The SAP fixes all sizes and matching variables before release: tomato-model
targets (T), optional admitted cross-crop targets (X), matched-unlabeled
controls, and process controls. All arms use identical sample handling,
assays, exclusions, batches and operator procedures. The cross-crop arm may
not exceed one quarter of non-control candidates.

## Data return and analysis

The laboratory returns every attempted result, raw files, batch/operator/
instrument metadata, exclusions and failure reasons. The SAP defines hit
criteria, the primary enrichment comparison, secondary analyses and failure
handling. No post-result filtering, candidate replacement or reranking is
permitted; unsuccessful candidates remain in the returned record.

## Release boundary

This repository policy creates no real K, power calculation, candidate table,
or biological claim. A real release is a separate approved operation after
model, candidate table, SAP and collaborator protocol are all frozen.
