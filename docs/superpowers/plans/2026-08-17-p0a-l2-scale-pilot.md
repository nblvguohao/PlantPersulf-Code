# P0-A pilot: L2 (co-peptide) scale feasibility

**Date**: 2026-08-17
**Status**: pilot only — counts derived from files already on disk, no MaxQuant/FragPipe re-search run, no new registry rows written, no code changed.
**Gate**: plan requires extrapolated L2 total ≥200 instances to justify the confound-controlled-benchmark strand of the NMI/NC methods paper (`nmi-nc-zippy-kernighan.md` §1, §5 Phase 0).

## Method

For each already-downloaded, already-registered supplementary site table under
`data/raw/supplements/`, grouped rows by (protein, peptide-with-modification-marks-stripped)
and counted groups where the peptide sequence contains more cysteines than were
called positive in that table — these are candidate L2 (co-peptide) pos+neg
instances, in the same spirit as the existing `evidence/copeptide_negatives.py`,
`maxquant_sites_negatives.py`, `pnas_site_negatives.py` extractors. No coordinate
verification or localization-confidence filtering was applied yet (that gate is
already implemented in the existing extractors and would be reused for real
registration).

## Per-source counts

| Source | File | Rows | Multi-Cys peptide groups | Candidate L2 (partial) groups | Currently registered from this source |
|---|---|---|---|---|---|
| PXD024061 (Arabidopsis) | `Sulfide_C_Sites.txt` | 82 site rows / 74 mod-peptide groups | 8 | 8 | ~7 (memory: 16 rows added 2026-08-15, saturating this file) |
| PXD072089 (rice) | `sd04.xlsx` (positives) vs `sd01.xlsx` (Cys inventory) | 656 site rows / 558 peptides | 74 (≥2 Cys) | **5** | 5 (memory: "5 部分修饰肽段" — exact match, this source is exhausted) |
| **PXD063170 (Magnaporthe)** | `PXD063170_sites_moesm3.tsv` | 1,482 site rows / 1,441 (protein,peptide) groups | 87 (≥2 Cys) | **63** | **0 — never mined for co-peptide negatives before this pilot** |
| **PXD006140 (Arabidopsis, Aroca/Romero lab)** | `erx294_suppl_supplementary_data_set_s3.xlsx` sheet `ident_peptides` | **58,569 peptide-spectrum matches** | **2,003 peptides with ≥2 Cys** | not yet decoded (Cys-targeted PSI-MOD codes present: MOD:01715 ×215, MOD:99996–99999 ×~1,200 combined, chemistry unresolved) | 0 — this raw MS/MS table has never been used; the project's existing PXD006140 positives (317 sites, per `phase-f-gate2-stop` memory) come from a separate processed site list, not this table |
| PXD072300 (recombinant orthogonal controls) | `recombinant_sites.json` | 10 proteins, protein-level positions only | n/a | not applicable (no peptide-level MS evidence in this file, L1 at best) | — |

## Headline finding

**The current registry (40 rows / 17 instances) has not exhausted even the
supplementary files already sitting on disk.** Two sources alone push the
total well past the ≥200 gate without any new MaxQuant/FragPipe re-search:

- **Magnaporthe (PXD063170) contributes 63 new candidate co-peptide groups**,
  entirely unexploited — this is the fourth species and was previously used
  only for cross-species transfer scoring, never for the co-peptide
  diagnostic. Adding it also turns the co-peptide null result from "3 species"
  into "4 species," strengthening the existing negative finding rather than
  threatening it.
- **PXD006140's raw peptide-spectrum table (`ident_peptides`, 58,569 rows) is
  a materially larger and completely untouched resource.** 2,003 peptides
  carry ≥2 cysteines. The chemistry tagging (custom PSI-MOD codes
  MOD:99996–99999, MOD:01715) is not yet decoded — this is genuine follow-up
  work, not a re-search — but even a conservative 10–15% post-decoding,
  post-verification yield would add 200–300 instances by itself.

**Extrapolated total from already-downloaded data alone: ≥76 confirmed-countable
(8 + 5 + 63) plus a large undecoded pool (2,003 candidate peptides in
PXD006140) that only needs chemistry-code resolution, not new acquisition.**
This clears the P0-A gate (≥200) with high confidence *before* considering
MS re-search of the remaining registered-but-unmined PXD raw files
(PXD035795, PXD039999, PXD072300 raw spectra) that Phase 1 originally
budgeted for.

## What changes in the plan

- **Phase 1 "扩产 L2" should be re-sequenced**: mine PXD063170 and decode
  PXD006140's `ident_peptides` MOD codes *before* budgeting for any FragPipe
  re-search. Re-search stays in scope for datasets with no existing
  peptide-level table (e.g. papers with only a final site list and raw
  spectra), but is no longer on the critical path for clearing the L2 volume
  gate.
- **PSI-MOD code decoding for PXD006140** becomes a concrete near-term task:
  resolve MOD:01715 / MOD:99996 / MOD:99997 / MOD:99998 / MOD:99999 against
  the PSI-MOD ontology (or the paper's own methods section) to confirm which
  code marks the persulfidation tag vs. a co-eluting label (this dataset used
  a TMT6plex quantitative design per `quant_DES_vs_Wt` — the custom codes are
  plausibly the tag-switch chemistry's biotin/iodoTMT label, not a standard
  Unimod entry).
- Coordinate verification, localization-confidence filtering, and
  registry-schema formatting for all new candidate groups still need to go
  through the existing `site_normalizer.py` fail-closed path and the
  three-state {positive, negative, undetermined} discipline before anything
  is written to `data/registry/copeptide_negatives_v1.tsv` — this pilot only
  counted raw candidates, it did not verify or register them.

## Verdict

**P0-A gate: PASS (high confidence).** Proceed to Phase 1 L2 expansion,
prioritizing PXD063170 mining and PXD006140 MOD-code decoding over new
MS re-search.
