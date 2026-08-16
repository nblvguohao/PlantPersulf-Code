# Tomato model-input provenance v1

## Frozen sequence universe

`data/raw/references/tomato_ref_proteome_v1.fasta` is the frozen
36,988-entry UniProtKB *Solanum lycopersicum* taxonomy snapshot queried as
`taxonomy_id:4081`, not the narrower current reference-proteome query. Its
SHA256 is
`1a46efe6461b239ec80b915a4c23ea39bbe39e37d43f9cc7c38b648523efca0f`.

The local file creation time, 2026-07-21 19:30:36 Asia/Shanghai, is retained
as the available retrieval record. The original UniProt release label was not
captured and is not reconstructed by guesswork. On 2026-08-11, the current
official query returned 37,015 entries: every accession and amino-acid
sequence in this frozen file matched exactly, with 27 current-only entries.
The current response is therefore not substituted for the frozen snapshot.

This resolves source identity and sequence integrity while preserving the
coordinate universe used for KIAE271 parsing. Any future analysis that adopts
the current snapshot must use a new input ID, recompute coordinate validation
and clustering, and never overwrite this v1 input.

## Homology clusters

`data/processed/clusters/tomato_proteome_clusters_v1.tsv` has SHA256
`48750f12bb1f7a09c5f192dd069f568f10b4c1407a13e46c356bb7e9ba996e2f`.
It is reproduced from the preserved MMseqs2 log
`tmp/tomato_local_v1/mmseqs.log` using MMseqs2
`6f45232ac8daca14e354ae320a4359056ec524c2`:

```text
easy-cluster data/raw/references/tomato_ref_proteome_v1.fasta tmp/tomato_local_v1/tomato_clust tmp/tomato_local_v1/mmseqs_tmp --min-seq-id 0.3 -c 0.5 --cov-mode 0 --threads 8
```

The `createtsv` output was converted from representative/member to
member/representative, given the header `protein_accession\tcluster_id`, and
written with LF line endings. Rebuilding this transformation reproduces the
registered SHA256 exactly.

No biological labels, scores, or candidate ranks are created by this record.
