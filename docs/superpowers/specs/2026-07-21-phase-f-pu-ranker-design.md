# Phase F — Structure-aware PU ranker + strict external validation (design)

> Codex Task 9→10. Implements the roadmap Phase F. Scope: build the ranker and
> the external-validation machinery, then judge **Gate 2** honestly. The point
> is not to make metrics look good — it is to produce reproducible evidence for
> a GO/STOP decision and make dishonest wording structurally impossible.

## Decisive constraint

Only two studies carry site-level positives: `PXD006140` (317) and
`PXD024061` (73). **Both are from the same laboratory (Romero/Gotor, U.
Sevilla), the same tag-switch chemistry, the same species (Arabidopsis).** They
are not *independent* held-out studies. Existing baselines are near-random
across studies (leave-study-out `test_ap` ≈ 0.02). Therefore Gate 2 condition 1
("≥2 independent held-out studies beat the no-learning baseline") cannot be
satisfied by current public data, and the expected honest outcome is
`GATE2_STOP` → roadmap **Phase Z** (data-resource / evidence-audit deliverable).

## Components

- **Ranker** (`models/structure_ranker.py`, `models/calibration.py`): minimal
  gated-fusion of sequence / frozen-ESM / structure / study-context branches;
  Elkan-Noto two-step PU training (reuses `models/pu_risk`); MC-dropout
  uncertainty; Platt calibration. Integrity is structural: missing structure is
  masked at the encoder output (never mean-imputed; structure scaler fit on
  present rows only), and an unseen held-out study is gated to zero.
- **Runner wiring** (`scripts/run_experiment.py`): `structure_ranker` model +
  multi-branch feature assembly + the full Codex ablation matrix
  (`configs/experiments/pu_ranker_v1.yaml`; ESM-free variant
  `pu_ranker_seqstruct_v1.yaml`). ESM extraction is skipped when no ablation
  needs it.
- **External validation** (`evaluation/{bootstrap,permutation,external_validation}.py`,
  `scripts/validate_external.py`): cluster-level bootstrap CI, permutation test,
  canonical no-leakage leave-study-out partition, known-control recovery with
  leakage detection + independent-unit counting + failed/unmappable reporting.
- **Gate 2** (`evaluation/conclusion_gate.py`, `configs/gate2_v1.yaml`):
  mechanical five-condition evaluation → `GATE2_GO`/`GATE2_STOP`; honesty
  verifiers that forbid predictive wording under STOP and require every claim to
  cite an existing table. `studies_are_independent: false` is frozen until
  genuinely independent studies enter the benchmark.

## Honesty rules (enforced by tests)

- `test_predictive_claim_requires_gate2_go`: STOP + predictive wording ⇒ fail.
- `test_claims_have_supporting_tables`: a claim citing a missing table ⇒ fail.
- `test_missing_structure_mask_is_respected`: masked structure never leaks.
- Unproven conditions (no effect-CI, no structure gain, no cluster-robustness)
  fail rather than being assumed — `validate_external.py` records gaps, never
  fills them.

## Two-track evaluation policy: leave-study-out (primary) + Split A cluster
## split (supplementary, literature-comparable)

**Question this section answers**: published cysteine-PTM predictors (e.g.
Sul-BertGRU, Bioinformatics 2025, doi:10.1093/bioinformatics/btaf078; pCysMod,
Front Cell Dev Biol 2021, doi:10.3389/fcell.2021.617366) evaluate with a
within-integrated-dataset split, not a cross-study split. Is our
leave-study-out requirement self-imposed rigor with no basis, or a
recognised concern?

**It is a recognised concern, not self-imposed.** Evidence:

- Whalen, Schreiber & Noble, *Nat Rev Genet* 2021 (doi:10.1038/s41576-021-00434-9),
  "Navigating the pitfalls of applying machine learning in genomics" — the
  field's own landmark review of how data structure (incl. homology) inflates
  reported ML performance in genomics/proteomics.
- Mahmood et al., *Hum Genomics* 2017 (doi:10.1186/s40246-017-0104-8) —
  concrete empirical demonstration: variant-effect predictor AUCs collapse
  from optimistic published numbers to 0.52–0.75 when re-evaluated on truly
  independent, non-circular functional datasets.
- iSNO-PseAAC, *PLoS ONE* 2013 (doi:10.1371/journal.pone.0055844) — even
  cysteine-PTM prediction specifically has required an explicit
  sequence-identity cutoff to control homology bias since at least 2013; this
  is not a foreign constraint, it is the field's own accepted baseline
  practice — we are applying a *stronger* version of it (cross-study/lab/
  chemistry, not just cross-sequence-identity) because our benchmark's
  specific weakness (only 2 studies, same lab) sits exactly at the level
  sequence-identity dedup does not control for.
- MD-HIT, *npj Comput Mater* 2024 (doi:10.1038/s41524-024-01468-1) — states
  the general principle plainly: redundancy-controlled evaluation numbers
  are lower, but "better reflect models' true prediction capability." A lower
  number under a stricter split is not a worse model — it is a more honest one.

**Resulting policy**: strictness should match the strength of the claim, not
be maximised unconditionally.

| Track | Split | Purpose | Admissible for Gate 2? |
|---|---|---|---|
| Primary | leave-study-out (`pu_ranker_v1.yaml`) | Judge cross-study predictive value | Yes — the only admissible source |
| Supplementary | Split A cluster split (`pu_ranker_cluster_v1.yaml`) | Literature-comparable "is there any learnable signal" number, reported side-by-side with the primary track | **No, structurally never** |

`scripts/validate_external.py::_parse_fold_study` only recognises model names
of the shape `leave_<study>_out|...`; a cluster-split experiment's rows never
match this pattern and are silently excluded from `_collect_fold_metrics` —
this is locked in by
`tests/release/test_gate2_ignores_within_dataset_split_metrics.py`, so Gate 2
cannot ingest Split A numbers even if `--model-release` is pointed at the
wrong experiment by mistake.

**Reporting rule**: any manuscript/collaboration document may state the Split
A number for literature comparison, but must accompany it with the
`pu_ranker_cluster_v1.yaml` `limitation` text verbatim, and must never cite it
as evidence for cross-study, cross-lab, or cross-species generalisation.

## Environment note

`fair-esm` + ESM-2 weights are required for any ablation with `use_esm: true`
and are absent from the current dev machine; the full `pu_ranker_v1` run
therefore executes on a machine with those extras. The ESM-free path, the
runner wiring, the statistical primitives, the orchestrator, and the Gate-2
decision are all verified here on real data (sequence + 7 real AlphaFold
structures + registered controls).
