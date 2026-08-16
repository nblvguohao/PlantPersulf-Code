# Task 9.6 Unified Statistical Analysis and Claim Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add frozen-test per-species/pooled statistics, a frozen strict-track
bootstrap policy, a literature-track lead comparison, a three-way claim gate
(fixing the current fail-open "no evidence still claims improvement" defect),
and a forbidden-external-claim phrase guard to
`plantpersulf.evaluation.multispecies_reporting`; extend the existing Gate 2
isolation lock-in test to cover the two multispecies v2 track tags.

**Architecture:** Pure functions only, composed by one orchestrator
(`build_multispecies_statistical_report`). No I/O, no model selection, no
frozen-test unlock.

**Tech Stack:** Python 3.10, stdlib `dataclasses`, existing
`plantpersulf.evaluation.metrics` / `bootstrap` / `effect_size`, pytest,
Ruff, mypy.

## Global Constraints

- Task scope is Task 9.6 claim-gate/statistics infrastructure only; do not
  select or freeze a model, unlock the frozen test, or touch any frozen
  split/benchmark file.
- The strict-track bootstrap policy is frozen at `n_boot=10_000`,
  `seed=20260811` as named constants, never an ad hoc call-site value.
- `literature_track_leads` compares paired per-seed means over the same 10
  seeds; mismatched or empty series fail closed.
- `claim_class_for_multispecies_result` must return
  `no_supported_multispecies_claim` whenever neither track's evidence clears
  its bar — this is the primary defect fix and needs its own explicit test.
- Tests use only temporary/synthetic (score, label[, cluster]) tuples, never
  real biological records.

---

### Task 1: Per-species and pooled frozen-test metrics

**Files:**
- Modify: `src/plantpersulf/evaluation/multispecies_reporting.py`
- Modify: `tests/unit/test_multispecies_reporting.py`

**Interfaces:**
- `compute_species_metric(scored: list[tuple[float, str]], species: str) -> SpeciesMetric`
- `pooled_average_precision(scored_by_species: dict[str, list[tuple[float, str]]], species: tuple[str, ...]) -> float`

- [ ] **Step 1: Write the failing per-species metric test**

```python
def test_compute_species_metric_matches_underlying_pu_metrics() -> None:
    scored = [(0.9, "positive"), (0.8, "unlabeled"), (0.7, "positive"), (0.1, "unlabeled")]
    metric = compute_species_metric(scored, "arabidopsis")
    assert metric.species == "arabidopsis"
    assert metric.base_rate == 0.5
    assert metric.recall_at_50 == 1.0
    assert metric.mean_reciprocal_rank == 1.0
```

- [ ] **Step 2: Verify RED** — `compute_species_metric` does not exist yet.

- [ ] **Step 3: Implement `compute_species_metric`** using
  `metrics.average_precision` / `recall_at_k(50)` / `recall_at_k(200)` /
  `mean_reciprocal_rank`; raise on empty input or no positives.

- [ ] **Step 4: Write the failing pooled-AP test**

```python
def test_pooled_average_precision_differs_from_macro_average() -> None:
    scored_by_species = {
        "arabidopsis": [(0.9, "positive"), (0.1, "unlabeled")],
        "rice": [(0.2, "positive")] + [(0.8, "unlabeled")] * 9,
    }
    pooled = pooled_average_precision(scored_by_species, ("arabidopsis", "rice"))
    assert 0.0 < pooled < 1.0
```

- [ ] **Step 5: Verify RED, implement, verify GREEN.**

- [ ] **Step 6: Commit**

```powershell
git add src/plantpersulf/evaluation/multispecies_reporting.py tests/unit/test_multispecies_reporting.py
git commit -m "feat: add per-species and pooled frozen-test metrics for v2"
```

### Task 2: Literature-track lead comparison and frozen strict-track bootstrap policy

**Files:**
- Modify: `src/plantpersulf/evaluation/multispecies_reporting.py`
- Modify: `tests/unit/test_multispecies_reporting.py`

**Interfaces:**
- `literature_track_leads(candidate_macro_aps: list[float], baseline_macro_aps: list[float]) -> bool`
- `STRICT_BOOTSTRAP_N_BOOT = 10_000`, `STRICT_BOOTSTRAP_SEED = 20260811`
- `strict_track_bootstrap_delta(model_scored: ClusteredScored, baseline_scored: ClusteredScored) -> BootstrapResult`

- [ ] **Step 1: Write the failing lead-comparison test**

```python
def test_literature_track_leads_compares_paired_seed_means() -> None:
    assert literature_track_leads([0.5, 0.6], [0.4, 0.4]) is True
    assert literature_track_leads([0.3, 0.3], [0.5, 0.5]) is False
    with pytest.raises(ValueError, match="equal length"):
        literature_track_leads([0.5], [0.4, 0.4])
```

- [ ] **Step 2: Verify RED, implement, verify GREEN.**

- [ ] **Step 3: Write the failing frozen-bootstrap-policy test**

```python
def test_strict_track_bootstrap_delta_uses_frozen_policy() -> None:
    model = [(0.9, "positive", "C1"), (0.2, "unlabeled", "C2")]
    baseline = [(0.6, "positive", "C1"), (0.5, "unlabeled", "C2")]
    result = strict_track_bootstrap_delta(model, baseline)
    assert result.n_boot == STRICT_BOOTSTRAP_N_BOOT == 10_000
    assert STRICT_BOOTSTRAP_SEED == 20260811
```

- [ ] **Step 4: Verify RED, implement (thin wrapper over
  `effect_size.paired_cluster_bootstrap_delta_ci`), verify GREEN.**

- [ ] **Step 5: Commit**

```powershell
git add src/plantpersulf/evaluation/multispecies_reporting.py tests/unit/test_multispecies_reporting.py
git commit -m "feat: freeze v2 strict-track bootstrap policy and literature lead check"
```

### Task 3: Three-way claim gate and forbidden-phrase guard

**Files:**
- Modify: `src/plantpersulf/evaluation/multispecies_reporting.py`
- Modify: `tests/unit/test_multispecies_reporting.py`

**Interfaces:**
- `claim_class_for_multispecies_result(*, strict_ci_lower, species_deltas, same_frozen_inputs, literature_leads) -> str`
- `NO_SUPPORTED_CLAIM`, `HOMOLOGY_ROBUST_CLAIM`, `LITERATURE_COMPARABLE_CLAIM`
- `claim_statement_for_class(claim_class: str) -> str`
- `find_forbidden_external_claims(text: str) -> list[str]`

- [ ] **Step 1: Write the failing no-evidence test (the actual defect fix)**

```python
def test_claim_gate_returns_no_supported_claim_without_either_track_win() -> None:
    """Previously this silently fell back to the literature-improvement claim
    even though the literature track never won anything — a fail-open bug."""
    assert claim_class_for_multispecies_result(
        strict_ci_lower=-0.02,
        species_deltas={"arabidopsis": -0.01, "rice": 0.0, "tomato": -0.02},
        same_frozen_inputs=True,
        literature_leads=False,
    ) == NO_SUPPORTED_CLAIM
```

- [ ] **Step 2: Verify RED** — the current two-way function has no
  `literature_leads` parameter and always falls back to the improvement
  claim.

- [ ] **Step 3: Implement the three-way gate**, add `literature_leads`
  as a required keyword, add `NO_SUPPORTED_CLAIM`, update existing callers
  (only the test file) for the new required argument.

- [ ] **Step 4: Verify GREEN**, including the two pre-existing
  `test_claim_requires_strict_track_ci_and_all_three_positive_deltas`
  assertions updated with an explicit `literature_leads=True`.

- [ ] **Step 5: Write the failing claim-statement test**

```python
def test_claim_statement_matches_mandated_wording() -> None:
    assert "文献同口径" in claim_statement_for_class(LITERATURE_COMPARABLE_CLAIM)
    assert "对未见同源家族稳健" in claim_statement_for_class(HOMOLOGY_ROBUST_CLAIM)
    assert "两条轨道均未证明" in claim_statement_for_class(NO_SUPPORTED_CLAIM)
    with pytest.raises(ValueError, match="unknown"):
        claim_statement_for_class("made_up_class")
```

- [ ] **Step 6: Verify RED, implement, verify GREEN.**

- [ ] **Step 7: Write the failing forbidden-phrase test**

```python
def test_forbidden_external_claim_phrases_are_flagged() -> None:
    assert find_forbidden_external_claims("这是外部泛化结果") == ["外部泛化"]
    assert find_forbidden_external_claims("within-dataset improvement only") == []
```

- [ ] **Step 8: Verify RED, implement, verify GREEN.**

- [ ] **Step 9: Commit**

```powershell
git add src/plantpersulf/evaluation/multispecies_reporting.py tests/unit/test_multispecies_reporting.py
git commit -m "fix: close the no-evidence claim-gate gap and add forbidden-phrase guard"
```

### Task 4: Statistical report orchestrator and Gate 2 isolation lock-in

**Files:**
- Modify: `src/plantpersulf/evaluation/multispecies_reporting.py`
- Modify: `tests/unit/test_multispecies_reporting.py`
- Modify: `tests/release/test_gate2_ignores_within_dataset_split_metrics.py`

**Interfaces:**
- `build_multispecies_statistical_report(*, scored_by_species, primary_species, pressure_species, strict_ci_lower, species_deltas, same_frozen_inputs, candidate_literature_macro_aps, baseline_literature_macro_aps) -> dict[str, object]`

- [ ] **Step 1: Write the failing orchestrator test**

```python
def test_statistical_report_composes_species_pooled_and_claim(...) -> None:
    report = build_multispecies_statistical_report(
        scored_by_species={...},
        primary_species=("arabidopsis", "rice", "tomato"),
        pressure_species=("magnaporthe",),
        strict_ci_lower=0.01,
        species_deltas={"arabidopsis": 0.01, "rice": 0.02, "tomato": 0.03},
        same_frozen_inputs=True,
        candidate_literature_macro_aps=[0.5, 0.6],
        baseline_literature_macro_aps=[0.4, 0.4],
    )
    assert report["claim_class"] == HOMOLOGY_ROBUST_CLAIM
    assert report["gate2_eligible"] is False
    assert "primary_pooled_average_precision" in report
```

- [ ] **Step 2: Verify RED, implement, verify GREEN.**

- [ ] **Step 3: Extend the Gate 2 isolation lock-in test**

```python
def test_multispecies_v2_track_tags_are_not_recognised_as_a_fold() -> None:
    assert _parse_fold_study("strict_cluster_holdout|structure_ranker:full") is None
    assert _parse_fold_study("literature_random_protein|structure_ranker:full") is None
```

- [ ] **Step 4: Verify RED is actually GREEN already** (the existing prefix
  check already excludes these tags) — this step locks the guarantee in as a
  regression test rather than fixing a defect; note this explicitly in the
  commit body.

- [ ] **Step 5: Run the Task 9.6 completion gate**

```powershell
pytest tests/unit -q
pytest tests/scientific -q
pytest tests/release -q
python -m plantpersulf.cli audit-registry
python -m plantpersulf.cli audit-files
python -m plantpersulf.cli audit-leakage --split-version multispecies_strict_v4
ruff check .
mypy src/plantpersulf
```

- [ ] **Step 6: Commit**

```powershell
git add src/plantpersulf/evaluation/multispecies_reporting.py tests/unit/test_multispecies_reporting.py tests/release/test_gate2_ignores_within_dataset_split_metrics.py
git commit -m "feat: add v2 statistical report orchestrator and Gate 2 isolation lock-in"
```

## Self-review

The plan adds only pure statistics/claim-gate functions on top of already
frozen split/benchmark infrastructure, fixes one concrete fail-open defect
(no-evidence results silently claiming a literature-comparable improvement),
freezes the mandated bootstrap policy as named constants, and extends an
existing Gate 2 isolation test rather than inventing a parallel mechanism. It
does not select a model, unlock the frozen test, or introduce any external/
prospective wording.
