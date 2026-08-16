# Task 9.5 Runtime Integrity Design

## Goal

Make each multispecies v2 development or literature-comparison task resumable
only when its scientific inputs, configuration and code are byte-identical,
while recording sufficient execution provenance for a release audit.

## Scope and constraints

- Applies only to Task 9.5 runtime infrastructure. It does not change a
  frozen split, create test metrics, alter model architecture, or begin Task
  9.6 statistical analysis.
- A task is a `(track, fold, seed, model)` tuple. Independent tuples may be
  assigned to different configured GPUs; this is task-level parallelism, not
  DDP.
- The only permitted OOM mitigation is reducing `score_batch_size`. The event
  must be logged with the original and replacement batch size; model,
  optimizer, data panel and split remain unchanged.
- No checkpoint may be reused unless its recorded input, config and code
  SHA256 values equal the current task fingerprint exactly.
- Manifests must include command, environment/dependency versions, device/GPU,
  code revision, dirty state, configured GPU map, and SHA256 values for all
  inputs, configuration, split, checkpoints and outputs.

## Architecture

`multispecies_v2` will expose a small runtime-contract module used by the v2
entrypoint. It has four responsibilities.

1. `TaskFingerprint` deterministically represents the task identity plus the
   input/config/code hashes and is serialized into every checkpoint and
   manifest record.
2. `append_training_event` validates all required JSONL fields, emits one
   complete newline-delimited JSON object with an atomic append boundary, and
   rejects partial or non-finite metric values.
3. `write_task_checkpoint` writes a task payload and fingerprint to a
   temporary sibling then atomically renames it. `load_resumable_checkpoint`
   rejects absent, malformed, incomplete, mismatched, or output-hash-invalid
   checkpoints; no best-effort recovery exists.
4. `write_run_manifest` expands the existing manifest to a stable schema and
   validates that every referenced existing artifact hashes to the declared
   value before atomically writing the manifest.

The CLI constructs one runtime contract from the approved v2 config and gives
it to each development-fold task. The default execution remains serial when
`max_parallel_tasks` is one. GPU dispatch is deterministic round-robin over
the configured `gpu_map`; an empty mapping leaves device selection unchanged.

## Data flow

```text
approved inputs + config + frozen split + git revision
                    -> TaskFingerprint
                    -> JSONL event / atomic checkpoint
                    -> resume gate (exact fingerprint required)
                    -> hash-verified run manifest
```

Checkpoint contents contain training state and hashes, never frozen test
labels. The strict track may only receive development rows under the existing
test-unlock boundary.

## Failure behavior

- Missing event fields, invalid device mapping, non-finite values, malformed
  JSON, missing artifacts, or any hash mismatch fail closed with `RuntimeError`
  or `ValueError` before a task resumes.
- A dirty git worktree is not hidden: its boolean state and stable dirty-file
  list are recorded. The code revision remains the exact `HEAD` identifier.
- An OOM without a positive smaller configured batch size fails. A successful
  retry records both attempts in JSONL and manifest deviations.
- Existing non-identical checkpoint or manifest paths are never overwritten.

## Tests and acceptance

Tests use only temporary policy markers and do not generate biological values.
They must prove: deterministic fingerprints; required JSONL schema; atomic
checkpoint creation; rejection of config/input/code mismatch; rejection of a
tampered checkpoint/output hash; deterministic GPU assignment; OOM batch-size
deviation recording; and release failure when manifest/log required fields are
absent. Existing v2 tests continue to prove that test labels are unavailable
without unlock.

Completion requires the specified unit, scientific, release, registry,
file-audit, v4 leakage-audit, Ruff and mypy commands. The Task 9.5 manifest
will describe runtime provenance only; no claims or new performance metrics
are introduced.

## Self-review

The design confines changes to Task 9.5, names every runtime boundary, defines
the only recovery condition, and prohibits a fallback that could change a
scientific result. It contains no placeholder data, score, or unapproved model
behavior.
