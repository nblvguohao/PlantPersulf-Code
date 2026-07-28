# Frozen Structure-Coverage v2 Validation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Determine, under frozen inputs and a fail-closed protocol, whether
expanding the registered AlphaFold coverage from 7 to 2,006 proteins produces
a stable structure-specific gain over sequence-only and structure-availability
controls in both leave-study-out folds.

**Architecture:** Freeze two immutable structure-registry releases, inject the
selected release explicitly into the existing feature path, audit provenance
and coverage before training, run five targeted ESM-free arms with identical
rows/splits/seeds/hyperparameters, then make a predefined six-condition
development decision. The experiment is additive: it cannot mutate
`pu_ranker_v1`, Gate 2, external-validation outputs, or the mutable registry.

**Tech Stack:** Python 3.11, PyTorch, NumPy, PyYAML, existing
`plantpersulf` feature/model/evaluation modules, pytest, Ruff, mypy, Git.

## Global Constraints

- Read `docs/PlantPersulf_Code_TDD_Codex.md` before every implementation
  session and follow its completion-report template.
- Work in strict RED -> GREEN -> REFACTOR cycles. Run and record the smallest
  failing test before production code.
- Stop after every numbered task and request project-reviewer approval before
  beginning the next task.
- Use only the registered public scientific inputs named in the approved
  design. Do not synthesize biological rows, labels, structures, metrics, or
  replacement files.
- Preserve positive-unlabeled semantics. `unlabeled` must never be renamed or
  interpreted as an experimental negative.
- Never tune on held-out studies or change frozen seeds, arms, thresholds, or
  hyperparameters after observing results.
- Never overwrite an existing output directory. Do not delete or alter the
  user-owned untracked `tmp/` directory.
- Keep `configs/gate2_v1.yaml` byte-identical and retain `GATE2_STOP`.
- The informal 2,006-structure rerun is not an input and none of its numbers
  may enter this experiment.

## Frozen Scientific Contract

| Input | Required SHA256 |
|---|---|
| `data/processed/benchmark_v1/sites.tsv` | `386178330DF17897606650FBD8C356BA66132CD1C3DB10F8E33AD53F46D5904C` |
| `data/raw/references/arabidopsis_ref_proteome_v1.fasta` | `51559016416634D52E2C92F6C30BB4297149841B01A60F683394CACF9D1033BF` |
| `data/processed/clusters/protein_clusters_v2.tsv` | `E13D16AEE81E0B98C9687BC98F8565122FBC027249D9E020A336A4F509671978` |
| `configs/experiments/pu_ranker_v1.yaml` | `079D17A8BD91D77E2CC5C138FA8D414500B843DD5E17E16177D1561BE3097AE0` |
| release-v1 registry | `E09D18D334A3CF41233E59EB6CF66C9611EF9B13CD407C30770F608674128326` |
| release-v2 registry | `BEE2D30ED28D754A7162283BB3C6080928DBC6A4CFA563958C48BD0146188075` |

The release-v1 bytes must be extracted from
`c827277^:data/registry/alphafold_structures.tsv`; release-v2 must be a
byte-for-byte copy of the currently audited registry. Do not reconstruct either
file by reserializing parsed rows.

---

## Task 1: Freeze the Registry Releases and Inject Registry Paths Explicitly

**Files:**

- Create:
  `data/registry/releases/alphafold_structures_release_v1.tsv`
- Create:
  `data/registry/releases/alphafold_structures_release_v2.tsv`
- Create:
  `tests/scientific/test_structure_registry_release_snapshots.py`
- Create:
  `tests/scientific/test_structure_registry_injection.py`
- Modify: `scripts/run_experiment.py`

- [ ] **Step 1: Create the release snapshots without editing scientific
  values**

Use Git plumbing for v1 and a byte-preserving filesystem copy for v2. Run the
v1 extraction end-to-end in `cmd` so PowerShell does not decode and re-encode
the blob. Before adding the files, calculate SHA256 and row counts:

```powershell
New-Item -ItemType Directory -Force data/registry/releases
cmd /d /c "git show c827277^^:data/registry/alphafold_structures.tsv > data\registry\releases\alphafold_structures_release_v1.tsv"
Copy-Item -LiteralPath data/registry/alphafold_structures.tsv `
  -Destination data/registry/releases/alphafold_structures_release_v2.tsv
Get-FileHash -Algorithm SHA256 data/registry/releases/*.tsv
```

Do not continue unless the two hashes exactly match the frozen contract.

- [ ] **Step 2: Write the snapshot-integrity RED test**

The test must verify exact hashes, 7 and 2,006 unique accessions, and v1 as an
exact record subset of v2:

```python
def test_frozen_structure_registries_have_expected_identity() -> None:
    v1 = RELEASES / "alphafold_structures_release_v1.tsv"
    v2 = RELEASES / "alphafold_structures_release_v2.tsv"
    assert sha256(v1) == V1_SHA256
    assert sha256(v2) == V2_SHA256
    rows_v1 = read_registry_records(v1)
    rows_v2 = read_registry_records(v2)
    assert len(rows_v1) == len({row["protein_accession"] for row in rows_v1}) == 7
    assert len(rows_v2) == len({row["protein_accession"] for row in rows_v2}) == 2006
    assert {canonical_record(row) for row in rows_v1} <= {
        canonical_record(row) for row in rows_v2
    }
```

Run:

```powershell
python -m pytest tests/scientific/test_structure_registry_release_snapshots.py -q
```

Expected RED: import failure because the release-registry reader/constants do
not yet exist, or missing release files if the snapshots were intentionally
held back until after the test was authored.

- [ ] **Step 3: Implement the smallest snapshot reader used by the test**

Place private parsing helpers in the test if production code does not yet need
them. Do not add a second scientific registry parser merely to satisfy the
test; production auditing in Task 2 will use
`plantpersulf.download.alphafold.audit_alphafold_structures`.

- [ ] **Step 4: Write the explicit-injection RED test**

The test must monkeypatch only the registry auditor, not biological feature
values. It should record the received path and return an empty audited
registry so no fabricated structure is created:

```python
def test_feature_assembly_uses_only_explicit_registry(monkeypatch, tmp_path) -> None:
    requested: list[Path] = []

    def record_registry(path: Path):
        requested.append(path)
        return ()

    monkeypatch.setattr(run_experiment, "audit_alphafold_structures", record_registry)
    explicit = tmp_path / "software-policy-registry.tsv"
    run_experiment._structure_feature_vectors_with_mask(
        REAL_BENCHMARK_ROWS[:1],
        REAL_PROTEOME,
        tmp_path,
        "injection",
        structure_registry_path=explicit,
    )
    assert requested == [explicit]
```

Run:

```powershell
python -m pytest tests/scientific/test_structure_registry_injection.py -q
```

Expected RED: `_structure_feature_vectors_with_mask()` rejects
`structure_registry_path`.

- [ ] **Step 5: Propagate the explicit path through existing feature assembly**

Change only these interfaces in `scripts/run_experiment.py`:

```python
def _structure_feature_vectors_with_mask(
    rows: list[BenchmarkRow],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
    *,
    structure_registry_path: Path = Path(
        "data/registry/alphafold_structures.tsv"
    ),
) -> tuple[np.ndarray, np.ndarray]:
    sources = audit_alphafold_structures(structure_registry_path)
    ...


def _structure_feature_vectors(
    rows: list[BenchmarkRow],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
    *,
    structure_registry_path: Path = Path(
        "data/registry/alphafold_structures.tsv"
    ),
) -> np.ndarray:
    ...


def _build_branch_features(
    rows: list[BenchmarkRow],
    proteome_path: Path,
    scratch_dir: Path,
    tag: str,
    *,
    need_esm: bool = True,
    structure_registry_path: Path = Path(
        "data/registry/alphafold_structures.tsv"
    ),
) -> BranchFeatures:
    ...
```

All new experiment callers must pass the path explicitly. Existing v1 callers
retain their default behavior, protected by the current release regression
tests.

- [ ] **Step 6: Run GREEN and regression tests**

```powershell
python -m pytest tests/scientific/test_structure_registry_release_snapshots.py tests/scientific/test_structure_registry_injection.py tests/scientific/test_missing_structure_mask_is_respected.py tests/scientific/test_model_release_is_reproducible.py -q
```

- [ ] **Step 7: Commit and stop for review**

```powershell
git add -- data/registry/releases scripts/run_experiment.py tests/scientific/test_structure_registry_release_snapshots.py tests/scientific/test_structure_registry_injection.py
git commit -m "feat: freeze and inject structure registry releases"
```

Record RED/GREEN output, hashes, limitations, and commit SHA using the TDD
template. Do not begin Task 2 until approved.

---

## Task 2: Implement the Fail-Closed Pre-Model Coverage Audit

**Files:**

- Create:
  `src/plantpersulf/evaluation/structure_coverage_audit.py`
- Create:
  `tests/unit/test_structure_coverage_audit_policy.py`
- Create:
  `tests/scientific/test_structure_coverage_audit.py`

- [ ] **Step 1: Write RED tests for the result schema and fail-closed rules**

Define the public API in the tests:

```python
@dataclass(frozen=True)
class FrozenFile:
    path: Path
    sha256: str


@dataclass(frozen=True)
class CoverageRecord:
    release: str
    protein_accession: str
    cys_position_in_protein: int
    label: str
    study_accession: str
    cluster_id: str
    has_registered_structure: bool
    maps_to_cys: bool
    mapping_status: str
    plddt: float | None


@dataclass(frozen=True)
class StructureCoverageAudit:
    records: tuple[CoverageRecord, ...]
    summary_rows: tuple[dict[str, object], ...]
    input_hashes: dict[str, str]
    blocking_errors: tuple[str, ...]

    def require_pass(self) -> None: ...
```

Allowed `mapping_status` values are exactly:

```python
{
    "mapped_cys",
    "absent_structure",
    "absent_residue",
    "non_cys_residue",
    "malformed_structure",
    "duplicate_residue_ambiguity",
}
```

Software-policy tests may copy a real registered row into a temporary file and
then alter its checksum field, remove its local file, or duplicate it with a
conflicting path. They must not invent a biological accession, residue, label,
or measurement.

Run:

```powershell
python -m pytest tests/unit/test_structure_coverage_audit_policy.py tests/scientific/test_structure_coverage_audit.py -q
```

Expected RED: module import failure.

- [ ] **Step 2: Implement frozen-file verification and registry identity
  checks**

Add:

```python
def sha256_file(path: Path) -> str: ...

def verify_frozen_file(item: FrozenFile) -> None:
    observed = sha256_file(item.path)
    if observed.upper() != item.sha256.upper():
        raise StructureCoverageAuditError(
            f"SHA256 mismatch for {item.path}: expected {item.sha256}, observed {observed}"
        )

def verify_registry_pair(
    v1_path: Path,
    v2_path: Path,
    *,
    expected_v1_count: int = 7,
    expected_v2_count: int = 2006,
) -> None: ...
```

Treat duplicate identical rows as an error because registry identity must be
one accession to one source. Treat duplicate conflicting rows as a distinct,
explicit blocking error. Call the existing AlphaFold auditor so every local
file is verified against its registered size and SHA256 before mapping.

- [ ] **Step 3: Implement retained site-level mapping records**

Add:

```python
def audit_structure_coverage(
    *,
    release: str,
    benchmark: FrozenFile,
    proteome: FrozenFile,
    clusters: FrozenFile,
    registry: FrozenFile,
) -> StructureCoverageAudit: ...
```

Required behavior:

- emit exactly one `CoverageRecord` for every benchmark row;
- join clusters by `protein_accession` and fail on any missing assignment;
- keep positives and unlabeled rows distinct;
- for positive rows, keep the real `study_accession`;
- represent unlabeled study context as the literal empty source field already
  present in the benchmark, not a guessed study;
- inspect registered structures with the existing residue parser;
- catch only expected mapping/parser exceptions and convert them to an allowed
  mapping status; never convert checksum, missing-file, or registry-conflict
  errors into non-blocking mapping failures.

- [ ] **Step 4: Implement summary and coverage-bias rows**

Summary dimensions must be machine-readable:

```text
release
scope                  # overall | study
study_accession
label                  # positive | unlabeled | all
metric                 # registered_proteins, covered_proteins, benchmark_sites,
                       # mapped_cys_sites, each mapping_status, plddt_lt_50,
                       # coverage_rate, coverage_label_risk_difference
estimate
ci_low
ci_high
n_sites
n_clusters
```

For coverage-label association, calculate the difference in registered
structure availability between observed positives and unlabeled sites. Obtain
the 95% interval by resampling real protein clusters with frozen seed `1729`
and `5000` replicates. This is a diagnostic only and must not set the final
stable/unstable decision.

- [ ] **Step 5: Run GREEN tests**

```powershell
python -m pytest tests/unit/test_structure_coverage_audit_policy.py tests/scientific/test_structure_coverage_audit.py -q
```

Assert on the real release snapshots that all registered files pass, all
benchmark rows are retained, every status is allowed, and positive/unlabeled
coverage is separately present.

- [ ] **Step 6: Commit and stop for review**

```powershell
git add -- src/plantpersulf/evaluation/structure_coverage_audit.py tests/unit/test_structure_coverage_audit_policy.py tests/scientific/test_structure_coverage_audit.py
git commit -m "feat: audit frozen structure coverage"
```

---

## Task 3: Freeze the Paired Experiment Config and Validate Controlled Variables

**Files:**

- Create:
  `configs/experiments/pu_ranker_structcover_v2.yaml`
- Create:
  `src/plantpersulf/evaluation/structure_coverage_config.py`
- Create:
  `tests/unit/test_structure_coverage_config.py`
- Create:
  `tests/release/test_structure_coverage_keeps_gate2_stop.py`

- [ ] **Step 1: Write the config-validation RED tests**

Tests must require:

- both exact release paths, hashes, and counts;
- the four frozen non-registry input hashes;
- studies `PXD006140` and `PXD024061`;
- seeds `[0, 1, 2, 3, 4]`;
- unlabeled-to-positive ratio `20` and subsample seed `12345`;
- copied ranker parameters `hidden=16`, `dropout=0.2`, `epochs=200`,
  `lr=0.05`, `n_mc_dropout=16`;
- exactly five allowed arms;
- `use_esm=false` and `use_study_context=false` for every arm;
- a distinct new output directory;
- rejection if any controlled variable differs across release arms.

Also hash `configs/gate2_v1.yaml` in the test and assert the current Gate 2
decision remains `GATE2_STOP`.

Run:

```powershell
python -m pytest tests/unit/test_structure_coverage_config.py tests/release/test_structure_coverage_keeps_gate2_stop.py -q
```

Expected RED: missing config and validation module.

- [ ] **Step 2: Add the frozen YAML**

Use this structure, with no hyperparameter search section:

```yaml
version: 2
experiment:
  name: pu_ranker_structcover_v2
  benchmark_labels:
    path: data/processed/benchmark_v1/sites.tsv
    sha256: 386178330DF17897606650FBD8C356BA66132CD1C3DB10F8E33AD53F46D5904C
  reference_proteome:
    path: data/raw/references/arabidopsis_ref_proteome_v1.fasta
    sha256: 51559016416634D52E2C92F6C30BB4297149841B01A60F683394CACF9D1033BF
  protein_clusters:
    path: data/processed/clusters/protein_clusters_v2.tsv
    sha256: E13D16AEE81E0B98C9687BC98F8565122FBC027249D9E020A336A4F509671978
  legacy_config:
    path: configs/experiments/pu_ranker_v1.yaml
    sha256: 079D17A8BD91D77E2CC5C138FA8D414500B843DD5E17E16177D1561BE3097AE0
registries:
  structcover_v1:
    path: data/registry/releases/alphafold_structures_release_v1.tsv
    sha256: E09D18D334A3CF41233E59EB6CF66C9611EF9B13CD407C30770F608674128326
    records: 7
  structcover_v2:
    path: data/registry/releases/alphafold_structures_release_v2.tsv
    sha256: BEE2D30ED28D754A7162283BB3C6080928DBC6A4CFA563958C48BD0146188075
    records: 2006
splits:
  mode: leave_study_out
  studies: [PXD006140, PXD024061]
evaluation:
  seeds: [0, 1, 2, 3, 4]
  cluster_bootstrap_seed: 1729
  cluster_bootstrap_replicates: 5000
  permutation_seed: 2718
  permutation_replicates: 5000
subsample:
  unlabeled_per_positive: 20
  seed: 12345
ranker:
  hidden: 16
  dropout: 0.2
  epochs: 200
  lr: 0.05
  n_mc_dropout: 16
arms:
  sequence_only:
    use_esm: false
    use_structure: false
    use_plddt: false
    use_accessibility: false
    use_study_context: false
  sequence_coverage_only:
    use_esm: false
    use_structure: true
    use_plddt: false
    use_accessibility: false
    use_study_context: false
  sequence_contact:
    use_esm: false
    use_structure: true
    use_plddt: false
    use_accessibility: true
    use_study_context: false
  sequence_plddt:
    use_esm: false
    use_structure: true
    use_plddt: true
    use_accessibility: false
    use_study_context: false
  sequence_contact_plddt:
    use_esm: false
    use_structure: true
    use_plddt: true
    use_accessibility: true
    use_study_context: false
output:
  directory: results/experiments/pu_ranker_structcover_v2
```

- [ ] **Step 3: Implement strict config loading**

Add:

```python
@dataclass(frozen=True)
class StructureCoverageExperimentConfig:
    ...

def load_structure_coverage_config(path: Path) -> StructureCoverageExperimentConfig:
    ...

def validate_controlled_variables(config: StructureCoverageExperimentConfig) -> None:
    ...
```

Reject unknown arms, missing fields, duplicated seeds, non-positive replicate
counts, an output inside v1/external-validation directories, or any frozen
hash mismatch. Validation must run before output-directory creation.

- [ ] **Step 4: Run GREEN tests**

```powershell
python -m pytest tests/unit/test_structure_coverage_config.py tests/release/test_structure_coverage_keeps_gate2_stop.py -q
```

- [ ] **Step 5: Commit and stop for review**

```powershell
git add -- configs/experiments/pu_ranker_structcover_v2.yaml src/plantpersulf/evaluation/structure_coverage_config.py tests/unit/test_structure_coverage_config.py tests/release/test_structure_coverage_keeps_gate2_stop.py
git commit -m "feat: freeze structure coverage experiment config"
```

---

## Task 4: Isolate the Five Structure Ablation Arms

**Files:**

- Modify: `src/plantpersulf/models/structure_ranker.py`
- Create:
  `tests/unit/test_structure_coverage_arms.py`
- Modify only if necessary:
  `tests/scientific/test_missing_structure_mask_is_respected.py`

- [ ] **Step 1: Write RED tests for feature isolation**

Use deterministic non-biological marker arrays in this pure software-policy
test. Verify output invariance:

```python
def test_coverage_only_cannot_consume_contact_or_plddt() -> None:
    first = project_structure_inputs(
        structure=np.array([[11.0, 22.0]]),
        mask=np.array([[1.0]]),
        arm="sequence_coverage_only",
    )
    second = project_structure_inputs(
        structure=np.array([[111.0, 222.0]]),
        mask=np.array([[1.0]]),
        arm="sequence_coverage_only",
    )
    np.testing.assert_array_equal(first.values, second.values)
    np.testing.assert_array_equal(first.values, np.array([[1.0]]))
```

Also assert:

- `sequence_only` exposes no structure tensor;
- `sequence_contact` changes with contact but not pLDDT;
- `sequence_plddt` changes with pLDDT but not contact;
- `sequence_contact_plddt` changes with both;
- a missing structure always produces an all-zero projected value and retains
  mask `0`.

Run:

```powershell
python -m pytest tests/unit/test_structure_coverage_arms.py -q
```

Expected RED: `project_structure_inputs` does not exist.

- [ ] **Step 2: Implement an explicit projection boundary**

Add:

```python
@dataclass(frozen=True)
class ProjectedStructureInputs:
    values: np.ndarray
    mask: np.ndarray

def project_structure_inputs(
    *,
    structure: np.ndarray,
    mask: np.ndarray,
    arm: str,
) -> ProjectedStructureInputs:
    ...
```

The coverage-only value is the real binary `structure_mask`, not a learned
proxy and not either scientific structure column. Do not rely solely on
zeroing inputs before a biased linear layer; this boundary must make prohibited
values inaccessible to the arm.

If the existing `StructureRanker` needs a one-dimensional structure encoder
for coverage-only, make input width explicit while preserving the current
two-dimensional default and its serialized v1 behavior.

- [ ] **Step 3: Run GREEN and v1 regression tests**

```powershell
python -m pytest tests/unit/test_structure_coverage_arms.py tests/scientific/test_missing_structure_mask_is_respected.py tests/scientific/test_model_release_is_reproducible.py -q
```

- [ ] **Step 4: Commit and stop for review**

```powershell
git add -- src/plantpersulf/models/structure_ranker.py tests/unit/test_structure_coverage_arms.py tests/scientific/test_missing_structure_mask_is_respected.py
git commit -m "feat: isolate structure coverage ablation arms"
```

---

## Task 5: Build the Non-Overwriting Paired Scoring Runner

**Files:**

- Create: `scripts/score_structure_coverage.py`
- Create: `tests/unit/test_score_structure_coverage.py`
- Create:
  `tests/scientific/test_structure_coverage_score_reproducibility.py`

- [ ] **Step 1: Write RED tests for fail-closed orchestration**

The unit tests must prove:

- config and both audits run before training;
- an existing output path raises `FileExistsError`;
- v1 and external-validation paths are rejected;
- each release receives its own explicit registry path;
- fold rows, subsample rows, seeds, model parameters, and arms are identical
  across releases;
- score rows are aligned by
  `(held_out_study, seed, arm, protein_accession, cys_position_in_protein)`;
- any missing/non-finite score fails the run and is not silently dropped.

Run:

```powershell
python -m pytest tests/unit/test_score_structure_coverage.py -q
```

Expected RED: scoring module import failure.

- [ ] **Step 2: Implement the runner interfaces**

Add:

```python
SCORE_FIELDS = (
    "coverage_release",
    "held_out_study",
    "seed",
    "arm",
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "cluster_id",
    "has_registered_structure",
    "mapping_status",
    "plddt_bin",
    "score",
    "uncertainty",
)

def prepare_output_directory(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    path.mkdir(parents=True)

def run_structure_coverage_scoring(
    config_path: Path,
    *,
    verify_only: bool = False,
) -> Path:
    ...
```

Required order:

1. load and validate frozen config;
2. verify all frozen hashes;
3. verify the v1/v2 registry relationship;
4. audit both releases and call `require_pass()`;
5. in `verify_only`, emit no directory and stop successfully;
6. ensure output is absent and create it once;
7. build frozen folds/subsamples once;
8. loop release -> fold -> seed -> arm, passing the release registry path
   explicitly into `_build_branch_features`;
9. write sorted per-release TSVs and the audit tables;
10. hash every output and write the manifest last.

Write to a temporary sibling directory and atomically rename it only after all
runs pass. On failure, retain a text failure record outside the scientific
result directory; never publish a partial result tree as complete.

- [ ] **Step 3: Add the CLI**

```python
def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("configs/experiments/pu_ranker_structcover_v2.yaml"),
    )
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    run_structure_coverage_scoring(args.config, verify_only=args.verify_only)
    return 0
```

- [ ] **Step 4: Run GREEN unit tests and the real preflight**

```powershell
python -m pytest tests/unit/test_score_structure_coverage.py -q
python scripts/score_structure_coverage.py --verify-only
```

The preflight must print both release hashes, registered file counts, benchmark
site counts, blocking-error count zero, and `training_started=false`.

- [ ] **Step 5: Add and run the real reproducibility test**

The test must use registered benchmark/structure data and fixed seed. Compare
two score-table byte streams produced in separate temporary directories for a
single frozen fold/seed/arm execution. It may reduce the number of epochs only
if that reduced value is a test-only parameter that cannot enter the release
config or scientific outputs.

```powershell
python -m pytest tests/scientific/test_structure_coverage_score_reproducibility.py -q
```

- [ ] **Step 6: Commit and stop for review**

```powershell
git add -- scripts/score_structure_coverage.py tests/unit/test_score_structure_coverage.py tests/scientific/test_structure_coverage_score_reproducibility.py
git commit -m "feat: score paired structure coverage releases"
```

Do not launch the full paired scientific run until this runner and its real
preflight are reviewed.

---

## Task 6: Implement Predefined Metrics and the Six-Condition Decision

**Files:**

- Create:
  `src/plantpersulf/evaluation/structure_coverage_decision.py`
- Create: `scripts/validate_structure_coverage.py`
- Create:
  `tests/unit/test_structure_coverage_decision.py`
- Create:
  `tests/scientific/test_structure_coverage_validation.py`

- [ ] **Step 1: Write the fail-closed decision RED tests**

Use policy booleans, not invented biological measurements:

```python
def test_decision_requires_every_condition() -> None:
    passed = {
        "full_exceeds_sequence_both_studies": True,
        "cluster_ci_excludes_zero_both_studies": True,
        "all_seed_directions_positive_both_studies": True,
        "full_exceeds_coverage_only_both_studies": True,
        "top_cluster_removed_gain_positive_both_studies": True,
        "all_audits_pass": True,
    }
    assert decide_structure_signal(passed).status == "STRUCTURE_SIGNAL_STABLE"
    for condition in passed:
        failed = passed | {condition: False}
        decision = decide_structure_signal(failed)
        assert decision.status == "STRUCTURE_SIGNAL_UNSTABLE"
        assert condition in decision.failed_conditions
```

Reject missing/unknown conditions. Run:

```powershell
python -m pytest tests/unit/test_structure_coverage_decision.py -q
```

Expected RED: decision module import failure.

- [ ] **Step 2: Implement exact decision types**

```python
REQUIRED_CONDITIONS = (
    "full_exceeds_sequence_both_studies",
    "cluster_ci_excludes_zero_both_studies",
    "all_seed_directions_positive_both_studies",
    "full_exceeds_coverage_only_both_studies",
    "top_cluster_removed_gain_positive_both_studies",
    "all_audits_pass",
)

@dataclass(frozen=True)
class StructureSignalDecision:
    status: Literal[
        "STRUCTURE_SIGNAL_STABLE",
        "STRUCTURE_SIGNAL_UNSTABLE",
    ]
    conditions: dict[str, bool]
    failed_conditions: tuple[str, ...]
    gate2_status: Literal["GATE2_STOP"] = "GATE2_STOP"

def decide_structure_signal(
    conditions: Mapping[str, bool],
) -> StructureSignalDecision: ...
```

- [ ] **Step 3: Write scientific validation RED tests**

Tests should read real score rows from a checked-in minimal release fixture
only if such a fixture is produced from registered inputs and records its
source hashes. Otherwise, generate the fixture during the test using the
fixed-seed real-data scoring helper; do not hand-author scores.

Assert calculations for:

- AP, Recall@K, MRR, and enrichment;
- paired AP delta by real protein cluster;
- seed-specific directions;
- complete structure versus coverage-only;
- top-contributing-cluster removal;
- aligned covered-subset comparisons;
- pLDDT bins `<50`, `50-<70`, `70-<90`, `>=90`;
- fixed-seed permutation;
- seed score/rank stability;
- uncertainty risk-coverage;
- failure on unaligned or non-finite rows.

- [ ] **Step 4: Implement validation and table output**

Add:

```python
def validate_structure_coverage_results(
    *,
    config_path: Path,
    result_directory: Path,
) -> StructureSignalDecision:
    ...
```

Write exactly:

```text
metrics.tsv
paired_effects.tsv
seed_stability.tsv
cluster_sensitivity.tsv
structure_coverage_decision.json
manifest.json
```

Use existing `average_precision`, `recall_at_k`,
`mean_reciprocal_rank`, `paired_cluster_bootstrap_delta_ci`, and
`top_cluster_dominance` implementations. Extend them only when an exact
predefined statistic is absent. The primary comparison is:

```text
structcover_v2 / sequence_contact_plddt
minus
sequence_only
```

within each held-out study, paired on the same rows and real cluster IDs.

The CLI must refuse to validate a partial result, an unexpected hash, or an
already finalized decision file.

- [ ] **Step 5: Run GREEN tests**

```powershell
python -m pytest tests/unit/test_structure_coverage_decision.py tests/scientific/test_structure_coverage_validation.py -q
```

- [ ] **Step 6: Commit and stop for review**

```powershell
git add -- src/plantpersulf/evaluation/structure_coverage_decision.py scripts/validate_structure_coverage.py tests/unit/test_structure_coverage_decision.py tests/scientific/test_structure_coverage_validation.py
git commit -m "feat: validate structure coverage stability"
```

---

## Task 7: Run the Frozen Experiment and Write the Phase S0 Decision

**Files:**

- Generate only:
  `results/experiments/pu_ranker_structcover_v2/`
- Create after successful validation:
  `docs/phase_s0_structure_coverage_decision.md`

- [ ] **Step 1: Record a clean pre-run integrity snapshot**

```powershell
git status --short
Get-FileHash -Algorithm SHA256 data/processed/benchmark_v1/sites.tsv
Get-FileHash -Algorithm SHA256 data/raw/references/arabidopsis_ref_proteome_v1.fasta
Get-FileHash -Algorithm SHA256 data/processed/clusters/protein_clusters_v2.tsv
Get-FileHash -Algorithm SHA256 configs/experiments/pu_ranker_v1.yaml
Get-FileHash -Algorithm SHA256 data/registry/releases/alphafold_structures_release_v1.tsv
Get-FileHash -Algorithm SHA256 data/registry/releases/alphafold_structures_release_v2.tsv
```

Abort on any mismatch. Confirm the only unrelated untracked item remains
`tmp/`, and do not stage it.

- [ ] **Step 2: Run the immutable scorer once**

```powershell
python scripts/score_structure_coverage.py --config configs/experiments/pu_ranker_structcover_v2.yaml
```

Do not rerun into the same path. A computational failure is a failed frozen run
to report, not authorization to tune parameters.

- [ ] **Step 3: Run validation once**

```powershell
python scripts/validate_structure_coverage.py --config configs/experiments/pu_ranker_structcover_v2.yaml --results results/experiments/pu_ranker_structcover_v2
```

- [ ] **Step 4: Verify the result manifest**

Recalculate every output SHA256, confirm row alignment, confirm all five seeds
and both studies are present, and independently recompute the six decision
booleans from the generated tables.

- [ ] **Step 5: Write the decision document from generated evidence**

`docs/phase_s0_structure_coverage_decision.md` must contain:

1. exact decision (`STRUCTURE_SIGNAL_STABLE` or
   `STRUCTURE_SIGNAL_UNSTABLE`);
2. all six conditions with pass/fail evidence paths;
3. both study-specific primary AP deltas and cluster-bootstrap intervals;
4. every seed direction;
5. complete-versus-coverage-only results;
6. top-cluster-removal results;
7. coverage audit and coverage-label-bias diagnostic;
8. all predefined secondary analyses;
9. frozen input/output hashes and exact commands;
10. explicit statement that Gate 2 remains `GATE2_STOP`;
11. the binding limitation that both studies are from the same laboratory,
    species, and chemistry;
12. the next-task consequence:
    structure may enter later architecture only if stable; otherwise no larger
    3D branch is justified.

Every numerical sentence must cite a generated TSV/JSON path. Do not round a
confidence bound across zero.

- [ ] **Step 6: Run the full completion gate**

Use the exact commands mandated in the repository TDD document. At minimum:

```powershell
python -m pytest tests/unit -q
python -m pytest tests/scientific -q
python -m pytest tests/release -q
python -m pytest -m "not slow" -q
python -m ruff check .
python -m mypy --strict src scripts
```

If repository configuration defines narrower canonical commands, run those as
well and report both. Any failure blocks completion.

- [ ] **Step 7: Commit scientific outputs and decision, then stop**

```powershell
git add -- results/experiments/pu_ranker_structcover_v2 docs/phase_s0_structure_coverage_decision.md
git commit -m "results: validate expanded structure coverage"
```

The completion report must include RED and GREEN commands/output, provenance,
integrity checks, limitations, decision, and commit SHA. Do not begin
Sul-BertGRU reproduction or PersulfFormer-SAR development until this task is
reviewed.

---

## Final Review Checklist

- [ ] Release-v1 is exactly 7 real registered records and an exact subset of
  release-v2.
- [ ] Release-v2 is exactly 2,006 real registered records.
- [ ] Every scientific input and registered local structure passes SHA256.
- [ ] Registry selection is explicit in every new scoring path.
- [ ] All five arms are ESM-free and study-context-free.
- [ ] Coverage-only cannot consume contact or pLDDT values.
- [ ] The paired releases use identical rows, folds, seeds, subsampling, and
  ranker parameters.
- [ ] No output was overwritten and no v1 result changed.
- [ ] Every benchmark row has one retained mapping record and reason.
- [ ] Both held-out studies and all five seeds appear in every required
  comparison.
- [ ] The six-condition decision fails closed.
- [ ] Gate 2 remains `GATE2_STOP`.
- [ ] No SOTA, independent-validation, or mechanism claim is made.
- [ ] Full pytest, Ruff, and strict mypy gates pass.
- [ ] The reviewer approved each task before the next one began.
