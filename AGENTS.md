# PlantPersulf-Code execution contract

These instructions apply to the entire repository.

## Task sequencing

- Read `docs/PlantPersulf_Code_TDD_Codex.md` before making changes.
- Execute only the task explicitly approved by the project reviewer.
- Do not begin a later task, data acquisition, scientific analysis, or model
  development before the current task is reviewed.
- Use strict test-driven development for every behavior change: write the
  smallest failing test, confirm the failure is caused by missing behavior,
  implement only enough to pass, then refactor while all tests remain green.

## Scientific integrity

- Scientific inputs must be real, public, traceable data or formally supplied
  collaboration data with a sample sheet, experimental design, and raw files.
- Do not create synthetic or random biological values, simulated labels,
  fabricated samples, placeholder results, invented metrics, or substitute
  files after a download or parse failure.
- Do not treat an unobserved persulfidation site as an experimental negative.
  The primary task is positive-unlabeled unless explicit negative evidence is
  registered.
- Do not use a scientific input unless its accession, source file, and SHA256
  are registered. Do not guess missing fields or silently resolve ID conflicts.
- Do not tune on a test set or overwrite a frozen data split.
- Do not use the same known mechanism for training and independent validation.
- Do not rerank a frozen prospective candidate list after wet-lab results.
- Pure software-policy tests may create temporary paths, malformed files, and
  non-biological marker text. They must never enter scientific outputs.

## Completion gate

For each task, run the task tests, all fast tests, scientific integrity tests,
release tests, Ruff, and mypy. Report the RED and GREEN commands and outputs,
data provenance, integrity checks, limitations, and commit SHA using the
template in the TDD document.
