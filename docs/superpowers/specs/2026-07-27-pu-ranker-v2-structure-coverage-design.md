# Phase S0 — Frozen Structure-Coverage v2 Validation Design

**Date:** 2026-07-27  
**Status:** approved for specification; implementation requires a separate
reviewed plan  
**Parent objective:** develop a credible state-of-the-art persulfidation
predictor without promoting an exploratory run or a biased benchmark result
into a scientific claim

## 1. Decision

The first SOTA-program task is a controlled validation of the expanded
AlphaFold structure coverage. It does not introduce the proposed
PersulfFormer-SAR architecture and does not compare against external
predictors yet.

The experiment answers one question:

> When the structure registry expands from the seven structures used by the
> frozen `pu_ranker_v1` release to the 2,006 structures currently registered,
> does real structure information provide a stable, leakage-controlled gain
> over sequence-only prediction in both leave-study-out folds?

The informal rerun described in `docs/phase_f_gate2_decision.md` is excluded
from all evidence. It used a mutable default registry, was not designed as an
experiment, and produced large seed instability. No number from that run may
be copied into the new result record.

## 2. Scientific boundary

This task can establish that expanded structure coverage is or is not worth
carrying into later SOTA model development. It cannot:

- change `configs/gate2_v1.yaml`;
- change `studies_are_independent: false`;
- turn the same-laboratory Arabidopsis folds into independent validation;
- use the rice or fungal transfer tracks for model or threshold selection;
- overwrite `pu_ranker_v1` inputs, scores, decisions, or documentation;
- claim SOTA, general prediction, cross-laboratory transfer, or biological
  mechanism;
- add tomato context, a graph neural network, SAR-PU, external predictor data,
  or new scientific data acquisition.

Gate 2 remains `GATE2_STOP` regardless of the numerical result because its
independent-study condition is structurally locked. The binding downgrade
statement remains:

> 当前公开数据不足以证明跨研究预测能力，模型仅用于候选组织与假设生成。  
> Current public data are insufficient to demonstrate cross-study predictive
> ability; the model is used only for candidate organisation and hypothesis
> generation.

## 3. Frozen inputs

The implementation must create immutable release snapshots rather than rely
on the mutable default path
`data/registry/alphafold_structures.tsv`.

| Role | Source | Frozen property |
|---|---|---|
| Benchmark | `data/processed/benchmark_v1/sites.tsv` | SHA256 `386178330DF17897606650FBD8C356BA66132CD1C3DB10F8E33AD53F46D5904C` |
| Reference proteome | `data/raw/references/arabidopsis_ref_proteome_v1.fasta` | SHA256 `51559016416634D52E2C92F6C30BB4297149841B01A60F683394CACF9D1033BF` |
| Homology clusters | `data/processed/clusters/protein_clusters_v2.tsv` | SHA256 `E13D16AEE81E0B98C9687BC98F8565122FBC027249D9E020A336A4F509671978` |
| Legacy experiment | `configs/experiments/pu_ranker_v1.yaml` | SHA256 `079D17A8BD91D77E2CC5C138FA8D414500B843DD5E17E16177D1561BE3097AE0` |
| Structure registry v1 | the seven-row registry in the parent of commit `c827277` | exact seven records; every local structure file must pass its recorded SHA256 |
| Structure registry v2 | current `data/registry/alphafold_structures.tsv` | 2,006 records; SHA256 `BEE2D30ED28D754A7162283BB3C6080928DBC6A4CFA563958C48BD0146188075` |

The release snapshots will be stored as:

- `data/registry/releases/alphafold_structures_release_v1.tsv`;
- `data/registry/releases/alphafold_structures_release_v2.tsv`.

The v1 snapshot is provenance-preserving extraction from Git history. The v2
snapshot is a byte-for-byte copy of the audited current registry. Biological
values and source metadata must not be edited while creating either snapshot.

## 4. Pre-model coverage audit

No training may start until both structure snapshots pass a fail-closed audit.
The audit must produce one machine-readable table and one manifest containing
the input hashes and code commit.

For each snapshot and each study, report:

- registered protein count;
- registry rows with duplicate or conflicting accessions;
- local files present, missing, or hash-mismatched;
- benchmark proteins covered;
- benchmark cysteine rows that map to a CYS residue;
- positive and unlabeled site coverage separately;
- mapping failure counts: absent structure, absent residue, non-CYS residue,
  malformed structure, and duplicate-residue ambiguity;
- pLDDT distribution and the fraction below 50;
- coverage by protein cluster.

The audit must explicitly test whether structure availability is associated
with the observed PU label, overall and within each discovery study. This is a
coverage-bias diagnostic, not evidence that availability causes the label.
Effect sizes and cluster-aware confidence intervals must be reported; a
significant difference does not abort the experiment, but it forces the
coverage-only control described below and must appear in the result
limitations.

Training is blocked if:

- an input hash differs from the frozen design;
- a registered local structure file is absent or fails SHA256;
- an accession has conflicting registry records;
- the v1 snapshot is not an exact seven-record subset of v2;
- a benchmark, split, seed list, or model hyperparameter differs between the
  two coverage arms.

## 5. Experiment design

### 5.1 Controlled variable

The only scientific input allowed to differ between the paired release arms is
the structure registry snapshot:

- `structcover_v1`: seven registered structures;
- `structcover_v2`: 2,006 registered structures.

Both arms use:

- the same benchmark rows;
- the same leave-study-out folds (`PXD006140`, `PXD024061`);
- the same real MMseqs2 cluster assignments;
- the same train/validation partition;
- the same 1:20 unlabeled-to-positive subsampling;
- seeds `[0, 1, 2, 3, 4]`;
- identical sequence features, scaler fitting, model parameters, epochs, and
  evaluation code.

The structure registry path must become an explicit config value propagated
through feature assembly. Production scoring must never call the default
registry path implicitly.

### 5.2 Targeted arms

The new experiment is ESM-free and study-context-free so it isolates structure
coverage without the known high-dimensional ESM instability or an unseen-study
branch.

| Arm | Sequence | Coverage indicator | Contact proxy | pLDDT | Purpose |
|---|---:|---:|---:|---:|---|
| `sequence_only` | yes | no | no | no | common reference |
| `sequence_coverage_only` | yes | yes | no | no | detects benefit from structure availability alone |
| `sequence_contact` | yes | implicit | yes | no | contact contribution |
| `sequence_plddt` | yes | implicit | no | yes | confidence contribution |
| `sequence_contact_plddt` | yes | implicit | yes | yes | complete current structure branch |

`sequence_coverage_only` uses the real `has_structure` flag as a single
feature. It must not be described as a biological structure model. It tests
whether registry membership alone carries study, annotation, abundance, or
accession-selection bias.

No hyperparameter search is allowed in this task. The current ranker
hyperparameters are copied unchanged from the frozen v1 config. If numerical
failure occurs, the run fails and is reported; parameters are not altered
after seeing test performance.

## 6. Evaluation

### 6.1 Primary estimand

For each held-out study, estimate the paired difference in Average Precision:

```text
AP(sequence_contact_plddt, structcover_v2)
-
AP(sequence_only)
```

The primary uncertainty interval is a paired bootstrap over real MMseqs2
protein clusters. Point estimates, all five seeds, the seed mean, and the 95%
cluster-bootstrap confidence interval are mandatory.

### 6.2 Secondary estimands

Report, without selecting the model from them:

- v2 complete structure versus v1 complete structure;
- contact-only and pLDDT-only gains versus sequence-only;
- complete structure versus coverage-only;
- Recall@K, MRR, and enrichment over the fold-specific PU base rate;
- results within the structure-covered subset;
- results stratified by pLDDT `<50`, `50–70`, `70–90`, and `>=90`;
- top-cluster removal sensitivity;
- fixed-seed permutation test;
- score and rank stability across seeds;
- risk-coverage curves using existing model uncertainty.

The covered-subset comparison must use aligned sites across the compared arms.
It cannot compare v1's small covered subset against v2's different covered
subset and call the difference a model gain.

### 6.3 Stability decision

The structure signal is marked `STRUCTURE_SIGNAL_STABLE` only when all
conditions hold:

1. `sequence_contact_plddt` exceeds `sequence_only` in both held-out studies;
2. the paired cluster-bootstrap 95% CI for the AP difference excludes zero in
   both studies;
3. the direction is positive for every one of the five frozen seeds in both
   studies;
4. complete structure exceeds `sequence_coverage_only` in both studies;
5. the gain remains positive after removing the single highest-contributing
   protein cluster;
6. no input, split, or registry audit fails.

Otherwise the decision is `STRUCTURE_SIGNAL_UNSTABLE`.

This decision is a development gate:

- `STABLE`: structure may enter the later PersulfFormer-SAR design;
- `UNSTABLE`: the SOTA program proceeds with external-predictor reproduction
  and selection-aware PU development, but no larger 3D graph branch is
  justified from current evidence.

Neither outcome changes Gate 2 or authorizes a SOTA claim.

## 7. Software boundaries

The implementation will introduce focused units:

- `src/plantpersulf/evaluation/structure_coverage_audit.py`  
  Audits frozen registries, local-file hashes, residue mapping, label coverage,
  and coverage-selection bias. It does not train a model.

- `configs/experiments/pu_ranker_structcover_v2.yaml`  
  Freezes both registry snapshots, exact input hashes, folds, seeds, model
  parameters, targeted arms, and the output directory.

- `scripts/score_structure_coverage.py`  
  Runs the paired v1/v2 arms from the frozen config and emits aligned per-site
  scores. It refuses an existing output directory unless an explicit
  verification-only mode is used.

- `scripts/validate_structure_coverage.py`  
  Computes the predefined comparisons and emits
  `structure_coverage_decision.json`.

- `docs/phase_s0_structure_coverage_decision.md`  
  Records the complete result, failed conditions, limitations, exact commands,
  data hashes, and commit. It is created only after the experiment runs.

Existing feature extraction and scoring code may be refactored only enough to
accept an explicit registry path. `pu_ranker_v1` behavior and defaults must
remain unchanged and be protected by regression tests.

## 8. Output contract

All generated scientific outputs live under:

```text
results/experiments/pu_ranker_structcover_v2/
├── coverage_audit.tsv
├── coverage_audit_manifest.json
├── scored/
│   ├── structcover_v1.tsv
│   └── structcover_v2.tsv
├── metrics.tsv
├── paired_effects.tsv
├── seed_stability.tsv
├── cluster_sensitivity.tsv
├── structure_coverage_decision.json
└── manifest.json
```

Every manifest includes:

- config SHA256;
- benchmark, proteome, cluster, and both registry SHA256 values;
- exact registry paths;
- code commit;
- package versions;
- fixed seeds;
- run timestamp;
- output-file SHA256 values.

No result path may point into
`results/experiments/pu_ranker_v1` or
`results/external_validation/pu_ranker_v1`.

## 9. TDD and integrity tests

Implementation follows strict RED → GREEN → REFACTOR cycles. Required tests
include:

1. explicit registry injection changes only structure lookup and never reads
   the mutable default registry;
2. v1 snapshot has exactly seven records and is an exact subset of v2;
3. checksum mismatch, missing local file, and conflicting accession each fail
   the coverage audit;
4. structure mapping retains every benchmark row and assigns an explicit
   missingness reason;
5. positive and unlabeled coverage are reported separately by study;
6. both coverage arms use identical rows, splits, seeds, subsampling, and
   hyperparameters;
7. the coverage-only arm cannot consume pLDDT or contact values;
8. the decision fails closed when any of its six conditions is absent;
9. the new run cannot overwrite v1 or an existing v2 output;
10. identical input, config, and seed reproduce identical per-site scores;
11. `configs/gate2_v1.yaml` remains unchanged and Gate 2 remains STOP;
12. every result-table claim cites an existing generated file.

Pure software-policy tests may use temporary paths and non-biological marker
text. Tests of structure coverage, labels, mappings, or scientific outputs
must use registered real records and must not fabricate biological values.

## 10. Failure handling

- Missing or corrupt registered structure: fail before training; no substitute
  file is generated.
- Mapping failure for an otherwise valid structure: retain the site with an
  explicit reason and `has_structure=false`.
- Non-finite score or metric: fail the affected run and the release decision.
- Seed-specific training failure: retain the failure record; do not drop the
  seed or replace it.
- Output already exists: fail rather than overwrite.
- Result below expectation: record `STRUCTURE_SIGNAL_UNSTABLE`; do not change
  thresholds, splits, seeds, or model parameters.

## 11. Completion gate

The task is complete only after:

- target RED and GREEN commands are recorded;
- all unit, fast scientific-integrity, and release tests pass;
- Ruff and strict mypy pass;
- both registry snapshots and all scientific inputs pass SHA256 audit;
- the paired experiment is reproducible from one documented command;
- the decision document reports all seeds, both studies, all failed cases,
  coverage bias, and limitations;
- the completion report follows the template in
  `docs/PlantPersulf_Code_TDD_Codex.md`;
- the implementation has its own reviewed commit SHA.

The next SOTA-program task is external Sul-BertGRU reproduction and the frozen
legacy/trustworthy dual leaderboard. It cannot begin until this structure
coverage task is reviewed.
