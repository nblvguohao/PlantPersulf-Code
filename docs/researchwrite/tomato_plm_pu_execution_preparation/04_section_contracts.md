# Section contracts

## Section: Phase A — v2 acceptance

- Purpose: define artifact validity and the only allowed interpretation of the ongoing run.
- Inputs: frozen config, implementation commit, registered tomato inputs, result manifest.
- Allowed claims: run completed reproducibly; candidate admitted or fallback selected.
- Forbidden claims: external validation, final candidate release, general prediction.
- Required evidence: file hashes, OOF completeness, admission intervals, all-run metrics.
- Validation: no overwrite; no best-seed reporting; panel is primary; failures retained.

## Section: Phase B — v3 PLM residual

- Purpose: define a bounded comparison that tests added representation value.
- Inputs: accepted v2 report, registered PLM weights, immutable split generator.
- Allowed claims: stable paired improvement if the existing admission rule passes.
- Forbidden claims: PLM attention as mechanism; tomato test labels used for tuning.
- Required evidence: fold-local preprocessing, ablations, repeated OOF scores, resource manifest.
- Validation: all model-selection dependencies are train-fold-only; structure defaults off; hard negatives absent.

## Section: Phase C — blind wet-lab interface

- Purpose: define prerequisites and immutable handoff fields without inventing sample sizes.
- Inputs: one selected final model, power analysis, collaborator-frozen protocol.
- Allowed claims: pre-registered candidate prioritization and later blind enrichment.
- Forbidden claims: choosing K from favorable ranks; reranking after wet-lab results.
- Required evidence: blind IDs, hidden key, matched-unlabeled controls, process controls, endpoint and analysis plan.
- Validation: all results retained; one site per protein by default; post-hit mechanism work cannot alter the primary release.

## Section: Stopping rules

- Purpose: make negative and technically invalid outcomes actionable.
- Inputs: manifests, admission decision, power calculation, complete wet-lab return.
- Allowed claims: downgrade, fallback, or stop according to a frozen branch.
- Forbidden claims: silent retries, sample deletion, threshold relaxation.
- Required evidence: reason code and immutable failed artifact where applicable.
- Validation: every branch ends in a defined artifact and permitted claim class.

