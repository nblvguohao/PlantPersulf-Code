# Task 9.6 Unified Statistical Analysis and Claim Gate Design

## Goal

Give the multispecies v2 dual-track evaluation (Task 9.1-9.5) one honest,
testable claim gate: compute the frozen-test statistical report the strict
homology-cluster track will need once a model is frozen, evaluate whether the
literature-comparable random-protein track (Task 9.4/9.5, already produces
per-seed macro AP) actually leads its baselines, and only then select one of
three mandated claim classes. No Gate 2 evidence, no external/prospective
wording, and no frozen-test scoring happen as a side effect of this task.

## Scope and constraints

- Applies only to Task 9.6 statistical-report and claim-gate infrastructure.
  It does not select or freeze a model, does not unlock or score the frozen
  20% test set, does not start any wet-lab or blind-cohort work, and does not
  change any frozen split, benchmark, or Task 9.1-9.5 artifact.
- The strict-track statistical policy (protein-cluster bootstrap over the
  model-vs-best-baseline paired AP delta) is frozen at exactly 10,000
  replicates, seed `20260811` — the same seed already used to freeze the
  strict 20% split, so the whole strict-track story reproduces from one
  documented constant. No caller may silently substitute a different
  `n_boot`/seed for a claim-bearing result.
- The literature-comparable track's "leads" signal is a plain paired mean
  comparison across the same 10 seeds already produced by Task 9.4/9.5
  (`three_crop_macro_average_precision` per seed/model in
  `literature_random_protein/summary.json`); it is intentionally not a CI —
  the spec only mandates a CI for the strict track.
- `claim_class_for_multispecies_result` gains a required `literature_leads`
  argument. Previously the function silently fell back to the
  "literature-comparable" claim whenever the strict-track criteria failed,
  even with zero evidence the literature track ever won anything. That is a
  real defect against the project's fail-closed discipline (red line 15: "达
  不到门槛时主动降级结论，不得包装升级") and is fixed by adding a third,
  strictly weaker claim class, `no_supported_multispecies_claim`, used when
  neither track's evidence clears its bar.
- Forbidden outward phrases ("外部泛化", "跨作物泛化", "前瞻验证" and their
  English equivalents) are unconditionally forbidden at this project stage:
  no tomato blind cohort has returned data yet, so there is no state in which
  they would be permitted. The guard therefore takes no "is validation
  complete" flag — that decision belongs to a future task once a blind
  cohort actually exists.
- Gate 2 isolation is already structurally guaranteed by
  `scripts/validate_external.py::_parse_fold_study`, which only recognises
  `leave_<study>_out`-tagged model names as leave-study-out folds (see
  `tests/release/test_gate2_ignores_within_dataset_split_metrics.py`). Task
  9.6 extends that existing lock-in test with the two multispecies v2 track
  tags (`strict_cluster_holdout`, `literature_random_protein`) rather than
  inventing a parallel isolation mechanism.

## Architecture

`plantpersulf.evaluation.multispecies_reporting` gains five additions on top
of the existing `SpeciesMetric` / `summarize_species_metrics` /
`claim_class_for_multispecies_result`:

1. `compute_species_metric(scored, species)` turns raw `(score, label)` pairs
   into a `SpeciesMetric` using the existing PU-appropriate
   `metrics.average_precision` / `recall_at_k` / `mean_reciprocal_rank`. Fails
   closed (`ValueError`/`RuntimeError`) on empty input or no positives, same
   as the underlying metric functions.
2. `pooled_average_precision(scored_by_species, species)` pools raw rows
   across species into one ranking and reports AP on the pooled list — the
   "池化结果" the spec asks for alongside the macro average that
   `summarize_species_metrics` already reports.
3. `literature_track_leads(candidate_macro_aps, baseline_macro_aps)` — a pure
   paired-mean comparison; requires equal non-empty length (same 10 seeds).
4. `strict_track_bootstrap_delta(model_scored, baseline_scored)` freezes the
   10,000-replicate/seed-`20260811` cluster-bootstrap policy on top of the
   existing `effect_size.paired_cluster_bootstrap_delta_ci`.
5. `find_forbidden_external_claims(text)` scans outward text for the
   unconditionally-forbidden phrase list.

`claim_class_for_multispecies_result` gains `literature_leads: bool` and a
third return value `NO_SUPPORTED_CLAIM`; `claim_statement_for_class(cls)` maps
each of the three claim classes to its exact mandated wording (bilingual,
matching the spec text verbatim for the two positive classes and a downgrade
statement for the negative one, in the same spirit as
`conclusion_gate.DOWNGRADE_STATEMENT`).

`build_multispecies_statistical_report(...)` is the single orchestrator a
future frozen-test run will call: given raw scored rows per species, primary
species tuple, pressure species tuple, the already-computed strict-track
bootstrap inputs, and the two literature-track macro-AP series, it returns
one dict containing the full species report, the pooled AP, the literature
lead flag, the selected claim class and its mandated statement, and an
explicit `gate2_eligible: False` marker.

## Data flow

```text
raw (score, label, cluster) per species   -> compute_species_metric (x N)
                                           -> summarize_species_metrics
                                           -> pooled_average_precision
literature_random_protein/summary.json    -> literature_track_leads
strict-track model/baseline scored rows   -> strict_track_bootstrap_delta
                                           -> claim_class_for_multispecies_result
                                           -> claim_statement_for_class
                                           -> build_multispecies_statistical_report
```

## Failure behavior

- Empty scored input, no positives, or missing species all fail closed
  (existing `metrics.py` / `summarize_species_metrics` behavior, reused
  as-is).
- Mismatched-length or empty literature macro-AP series raise `ValueError`
  before any lead/lag conclusion is drawn.
- An unknown claim class passed to `claim_statement_for_class` raises
  `ValueError` — there is no silent default statement.
- `no_supported_multispecies_claim` is the fail-closed default whenever
  neither track's evidence clears its bar; it is never upgraded implicitly.

## Tests and acceptance

Tests use only temporary policy markers and pre-computed synthetic
(score, label[, cluster]) tuples — the same pattern already used throughout
`tests/unit/test_multispecies_reporting.py` and
`tests/unit/test_multispecies_v2_completion.py`. They must prove: correct
per-species metric computation; pooled AP differs from macro AP when species
base rates differ; the three-way claim decision (homology-robust,
literature-comparable, no-claim) including the previously-missing no-claim
branch; the frozen bootstrap policy constants; the forbidden-phrase scanner;
and that the two multispecies v2 track tags remain invisible to
`scripts/validate_external.py::_parse_fold_study`.

Completion requires the specified unit, scientific, release, registry,
file-audit, v4 leakage-audit, Ruff and mypy commands. No claim, metric, or
result is generated against real frozen-test data by this task.

## Self-review

The design confines changes to Task 9.6 claim-gate/statistics
infrastructure, fixes a real fail-open defect in the existing claim function,
freezes the mandated bootstrap policy as named constants instead of ad hoc
call-site values, and reuses the existing Gate 2 isolation mechanism instead
of duplicating it. It does not select a model, unlock the frozen test, or
introduce any placeholder or unapproved wording.
