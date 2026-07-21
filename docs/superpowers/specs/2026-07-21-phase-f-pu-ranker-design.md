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

## Environment note

`fair-esm` + ESM-2 weights are required for any ablation with `use_esm: true`
and are absent from the current dev machine; the full `pu_ranker_v1` run
therefore executes on a machine with those extras. The ESM-free path, the
runner wiring, the statistical primitives, the orchestrator, and the Gate-2
decision are all verified here on real data (sequence + 7 real AlphaFold
structures + registered controls).
