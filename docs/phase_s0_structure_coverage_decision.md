# Phase S0 -- Structure-Coverage v2 Validation Decision Record

**Date**: 2026-08-10
**Decision**: **`STRUCTURE_SIGNAL_UNSTABLE`** (5/6 conditions passed)
**Gate 2**: Remains `GATE2_STOP` (no gate override performed or possible)
**Evidence**: `results/experiments/pu_ranker_structcover_v2/`

## What was tested

The PU ranker was scored under two frozen AlphaFold structure-registry
releases -- v1 (7 accessions, the coverage that existed at Gate 2 freeze)
and v2 (2,006 accessions, the expanded registry from cross-species work) --
on identical benchmark rows, leave-study-out folds, seeds, subsamples, and
hyperparameters. Five ESM-free ablation arms were compared:

| Arm | Description |
|---|---|
| `sequence_only` | [hydrophobicity, cys_density] only |
| `sequence_coverage_only` | sequence + binary structure-registered mask |
| `sequence_contact` | sequence + contact_number_proxy from structure |
| `sequence_plddt` | sequence + pLDDT from structure |
| `sequence_contact_plddt` | sequence + contact_proxy + pLDDT (full arm) |

All arms disable ESM and study-context branches. The primary comparison is
`sequence_contact_plddt` (v2) minus `sequence_only` (v2) within each
held-out study, paired on the same rows and real protein clusters.

## Coverage audit

| Release | Registered proteins | Benchmark sites mapped to Cys |
|---|---|---|
| structcover_v1 | 7 | 81 / 395,878 |
| structcover_v2 | 2,006 | 4,525 / 395,878 |

The expansion increases structure-annotated benchmark sites by 55x
(81 to 4,525). Both audits pass: every registered PDB file verifies
against its SHA256; no blocking errors.

## Six-condition decision protocol

| # | Condition | Result | Detail |
|---|---|---|---|
| 1 | `full_exceeds_sequence_both_studies` | **PASS** | PXD006140: full 0.145 > seq 0.057. PXD024061: full 0.160 > seq 0.048 |
| 2 | `cluster_ci_excludes_zero_both_studies` | **PASS** | PXD006140: delta +0.088, 95% CI [+0.084, +0.133]. PXD024061: delta +0.112, 95% CI [+0.064, +0.177] |
| 3 | `all_seed_directions_positive_both_studies` | **PASS** | All 10 seeds (2 studies x 5 seeds) show a positive full-vs-sequence-only AP delta |
| 4 | `full_exceeds_coverage_only_both_studies` | **FAIL** | PXD006140: full 0.145 < coverage_only 0.278. PXD024061: full 0.160 > coverage_only 0.041. The full arm does not outperform the binary "structure registered" mask in the PXD006140 fold |
| 5 | `top_cluster_removed_gain_positive_both_studies` | **PASS** | PXD006140: gain without Q3E937 (825 rows) = +0.092. PXD024061: gain without Q3E937 (160 rows) = +0.114 |
| 6 | `all_audits_pass` | **PASS** | Both release audits produce zero blocking errors |

**Result**: 5/6 conditions pass. One condition fails:
`full_exceeds_coverage_only_both_studies`. Decision is
`STRUCTURE_SIGNAL_UNSTABLE`.

## Why condition 4 fails

In the PXD006140 hold-out fold, the binary "has a registered AlphaFold
structure" signal alone (`sequence_coverage_only`) achieves AP 0.278 --
nearly double the full arm's AP of 0.145. The contact_number_proxy and
pLDDT features not only fail to add marginal value over that binary mask
in this fold, they actively degrade performance.

The pattern reverses for PXD024061 (full 0.160 > coverage_only 0.041),
so the structure-content features DO carry signal -- but that signal is
not stable across folds.

In essence, with 2,006 proteins covered, the simple act of *having* a
structure is a strong covariate that correlates with which cysteines are
persulfidation-positive in one study but not the other, and the
structure-content features (contact/pLDDT) do not reliably improve on
that binary baseline.

## Per-study AP table (structcover_v2)

| Study | sequence_only | coverage_only | contact+pLDDT |
|---|---|---|---|
| PXD006140 | 0.057 | **0.278** | 0.145 |
| PXD024061 | 0.048 | 0.041 | **0.160** |

## Seed-level delta_vs_sequence_only

| Study | seed=0 | seed=1 | seed=2 | seed=3 | seed=4 |
|---|---|---|---|---|---|
| PXD006140 | +0.013 | +0.135 | +0.258 | +0.467 | +0.433 |
| PXD024061 | +0.350 | +0.202 | +0.094 | +0.012 | +0.187 |

All seed deltas are positive -- the structure gain is consistent in
direction across seeds for the full-vs-sequence comparison.

## Seed-level delta_vs_coverage_only

| Study | seed=0 | seed=1 | seed=2 | seed=3 | seed=4 |
|---|---|---|---|---|---|
| PXD006140 | **-0.330** | +0.145 | +0.160 | +0.206 | **-0.044** |
| PXD024061 | +0.343 | +0.215 | +0.093 | +0.009 | +0.202 |

Two of five seeds are negative for PXD006140 -- the instability of the
full arm relative to coverage-only is seed-level, not just fold-level.

## Cluster sensitivity

Both studies survive removal of the top protein cluster (Q3E937):

| Study | Full AP | AP without Q3E937 | Retention ratio | Gain without Q3E937 |
|---|---|---|---|---|
| PXD006140 | 0.145 | 0.150 | 103.1% | +0.092 |
| PXD024061 | 0.160 | 0.162 | 101.3% | +0.114 |

Removing the largest cluster *increases* AP slightly in both folds,
confirming the structure gain is not driven by one dominant homology
cluster.

## Frozen input hashes (verified)

| Input | SHA256 |
|---|---|
| `benchmark_v1/sites.tsv` | `38617833...` |
| `arabidopsis_ref_proteome_v1.fasta` | `51559016...` |
| `protein_clusters_v2.tsv` | `E13D16AE...` |
| `pu_ranker_v1.yaml` | `079D17A8...` |
| `alphafold_structures_release_v1.tsv` | `E09D18D3...` |
| `alphafold_structures_release_v2.tsv` | `BEE2D30E...` |

(Full SHA256 values recorded in
`configs/experiments/pu_ranker_structcover_v2.yaml` and verified by
`scripts/score_structure_coverage.py --verify-only`.)

## Binding limitations

1. Both training studies are from the same laboratory (Romero/Gotor,
   U. Sevilla), the same tag-switch persulfidation chemistry, and the
   same species (Arabidopsis thaliana). Structure-coverage gain is
   demonstrated only within this specific experimental framework.

2. The `STRUCTURE_SIGNAL_UNSTABLE` decision is a development-level
   diagnostic, not a gate override. Gate 2 remains `GATE2_STOP`.

## Next-task consequence

Structure features may enter later architecture only when the structure
signal is demonstrably stable (six-condition pass). Until then, no larger
3D branch (graph neural networks, attention over residue graphs,
geometry-aware embeddings) is justified by the evidence collected here.

The binary coverage signal (`sequence_coverage_only`) is strong for
PXD006140 (AP 0.278) alone -- this warrants a separate investigation of
why structure availability correlates so strongly with persulfidation
positives in that study, but that investigation is outside the scope of
this frozen validation.

## Reproducibility commands

```powershell
# Preflight (hash-only)
$env:PYTHONPATH = "src;."
python scripts/score_structure_coverage.py --verify-only

# Full scoring (deterministic CPU, ~5-10 min)
$env:PYTHONPATH = "src;."
python scripts/score_structure_coverage.py --config configs/experiments/pu_ranker_structcover_v2.yaml

# Validation + decision
$env:PYTHONPATH = "src;."
python scripts/validate_structure_coverage.py --config configs/experiments/pu_ranker_structcover_v2.yaml --results results/experiments/pu_ranker_structcover_v2
```
