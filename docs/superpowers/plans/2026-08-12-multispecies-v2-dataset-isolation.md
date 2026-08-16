# Multispecies v2 Dataset Isolation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the Task 9.3 v2 data materialization boundary so all positive and unlabeled Cys rows inherit their frozen global-protein split, with within-partition PU panels and no v1 mutation.

**Architecture:** Keep `multispecies_v1` and its PANTHER grouping untouched. Add a v2-specific site-row materializer which reads registered positive evidence and reference proteins, resolves global MMseqs2 membership, assigns the immutable frozen split, then samples unlabeled rows inside each already-isolated partition. The experiment runner will expose development folds only and will continue to lock frozen-test labels and metrics by default.

**Tech Stack:** Python 3, PyYAML, pytest, registered FASTA/TSV inputs, immutable MMseqs2 cluster and split tables.

## Global Constraints

- Do not modify v1 configurations, data, splits, or historical results.
- A site key is `(species, protein_accession, cys_position)` and duplicate study evidence is merged into `study_accessions`.
- The primary split unit is global MMseqs2 `cluster_id`; PANTHER may only remain an unlabelled evolutionary feature.
- Unlabelled Cys are split first and sampled independently by partition at 20 per positive.
- Standardization, imputation, feature selection, PU prior fitting, and tuning must receive development-training rows only.
- Report Arabidopsis, rice, and tomato together; Magnaporthe is a separate pressure-test section.
- Frozen tests remain unreadable and unscorable unless a matching `test_unlock.json` exists.

---

### Task 1: Materialize split-bound v2 site rows

**Files:**
- Create: `src/plantpersulf/proteomics/multispecies_v2_dataset.py`
- Create: `tests/scientific/test_multispecies_v2_dataset.py`

**Interfaces:**
- Consumes: registered positive records, reference-proteome Cys positions, `GlobalClusterRow`, and `FrozenMultispeciesSplit`.
- Produces: `MultispeciesV2SiteRow` and `build_split_bound_v2_rows(...)`.

- [ ] **Step 1: Write failing tests** for duplicate-study merging, complete Cys inheritance to one split, and missing cluster membership failure.
- [ ] **Step 2: Run tests** and confirm import/API failure.
- [ ] **Step 3: Implement** the minimal immutable row type and mapping/materialization behavior.
- [ ] **Step 4: Run tests** and confirm they pass.

### Task 2: Sample PU panels only after split assignment

**Files:**
- Modify: `src/plantpersulf/proteomics/multispecies_v2_dataset.py`
- Modify: `tests/scientific/test_multispecies_v2_dataset.py`

**Interfaces:**
- Consumes: split-bound v2 rows.
- Produces: `sample_v2_unlabeled_panels(rows, per_positive=20, seed=...)`.

- [ ] **Step 1: Write failing tests** that assert no unlabeled key crosses partitions, every positive is retained, and sampling is deterministic.
- [ ] **Step 2: Run tests** and confirm missing behavior.
- [ ] **Step 3: Implement** independent within-partition sampling.
- [ ] **Step 4: Run tests** and confirm they pass.

### Task 3: Expose a development-only preparation API

**Files:**
- Modify: `src/plantpersulf/workflows/multispecies_v2.py`
- Modify: `tests/scientific/test_multispecies_v2_evaluation.py`

**Interfaces:**
- Consumes: frozen split and split-bound rows.
- Produces: development train/validation data for one fixed fold; never returns frozen-test rows by default.

- [ ] **Step 1: Write failing tests** proving train/validation contain no test rows and train-only fitting input cannot include validation rows.
- [ ] **Step 2: Run tests** and confirm failure.
- [ ] **Step 3: Implement** a narrow development-fold preparation function.
- [ ] **Step 4: Run tests** and confirm they pass.

### Task 4: Lock v2 contract and verify integration

**Files:**
- Modify: `configs/experiments/multispecies_v2_global_clusters_v2.yaml`
- Modify: `tests/unit/test_multispecies_v2_config.py`
- Modify: `scripts/train_multispecies_v2.py`

**Interfaces:**
- Consumes: Task 1–3 APIs and v2 configuration.
- Produces: default development-only preflight with an explicit data-isolation declaration.

- [ ] **Step 1: Write failing configuration and CLI tests** for v2-only inputs, `per_positive: 20`, and disabled test scoring by default.
- [ ] **Step 2: Run tests** and confirm failure.
- [ ] **Step 3: Implement** only the configuration/CLI wiring necessary to validate development-only preparation.
- [ ] **Step 4: Run focused, unit, scientific, release, registry, file-audit, leakage-audit, Ruff, and mypy checks; record unavailable verifier commands explicitly rather than claiming success.**
