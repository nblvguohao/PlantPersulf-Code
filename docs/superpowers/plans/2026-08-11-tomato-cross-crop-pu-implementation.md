# Tomato Blind and Cross-Crop PU Ranker Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a leakage-safe tomato PU site ranker, a gated target-label-free cross-crop arm, and a hash-frozen blind candidate release workflow without weakening the existing Gate 2 or scientific-integrity controls.

**Architecture:** A small additive PU ranker combines within-protein pairwise ranking with non-negative PU risk and train-fold-only scaling. Structure is a separately admitted residual, observation propensity is fail-closed until a complete discovery table exists, and cross-crop transfer uses a source-only dependency graph that rejects every tomato label. Candidate release is a later hard-gated task that consumes only an admitted frozen model and a pre-registered analysis plan.

**Tech Stack:** Python 3.10–3.11, PyTorch 2.x CPU for differentiable PU loss, scikit-learn 1.5+ for frozen baselines, PyYAML 6.x, pytest 8+, Ruff, mypy; no new dependency is permitted by this plan.

## Global Constraints

- Follow `docs/PlantPersulf_Code_TDD_Codex.md` and the repository `AGENTS.md` exactly.
- Execute only one task after explicit project-reviewer approval; stop after its commit and report RED/GREEN evidence before requesting approval for the next task.
- Do not start this plan while the current unreviewed task or dirty overlapping edits remain unresolved.
- Scientific inputs must be real, public or formally supplied, registered by accession/source file/SHA256, and parsed fail-closed.
- Never create synthetic biological sequences, samples, site labels, metrics, or replacement files. Pure software tests may use numeric arrays and non-biological marker IDs only; scientific tests must use registered real micro-fixtures.
- Labels are exactly `positive` or `unlabeled`; no implementation may introduce a hard biological `negative` label.
- Every scaler, feature selector, projection, class-prior sensitivity decision, and hyperparameter fit is train-fold-only.
- Test labels are evaluated exactly once and never used for model or feature selection.
- The future tomato lockbox and blind results never enter training, tuning, model admission, or reranking.
- `S_bio` cannot consume study, laboratory, species one-hot, accession, structure availability, pLDDT, observation status, PPI, AF3, Boltz-2, PLIP, docking, or MD scores.
- `q_obs` remains disabled until a complete registered discovery-observation table exists; this plan does not invent or approximate that table.
- The cross-crop arm may fit and tune on registered non-tomato source data only. Existing tomato knowledge means the claim is target-label-free fitting/tuning, not a target-naive study design.
- Random-protein, development-CV, cross-crop, and blind-candidate artifacts cannot change the existing Gate 2 decision unless a separately approved Gate 2 task satisfies its frozen rules.
- Structure residual defaults to disabled (`delta=0`) and can be enabled only by the admission task defined below.
- MSA conservation, frozen protein-language embeddings, fitted SAR propensity, second-crop wet lab, and mechanism simulation are outside this first implementation plan. Each requires a new reviewer-approved plan after its input gate exists.
- Blind-release policy Task 10 cannot start until the model, K, power analysis, blind protocol, and collaborator conditions are separately approved and frozen.
- Repository audit on 2026-08-11 found KIAE271 in `data/registry/supplementary_sources.tsv`, but did not find `tomato_ref_proteome_v1.fasta` or `tomato_proteome_clusters_v1.tsv` in a path/SHA256 registry. Before Task 7, a separately approved provenance task must verify the reference-proteome accession/source and the exact MMseqs2 generation command, then register both files. Neither this plan nor later code may infer that provenance.

## File Map

| Responsibility | File |
|---|---|
| Biological feature whitelist and six-dimensional core vector | `src/plantpersulf/features/site_biology.py` |
| Pairwise and nnPU losses | `src/plantpersulf/models/ranking_loss.py` |
| Deterministic additive PU model | `src/plantpersulf/models/additive_pu_ranker.py` |
| Observation-propensity fail-closed policy | `src/plantpersulf/models/observation_propensity.py` |
| Optional structure-only residual | `src/plantpersulf/models/structure_residual.py` |
| Structure residual admission | `src/plantpersulf/evaluation/structure_admission.py` |
| Rank intervals and applicability envelope | `src/plantpersulf/evaluation/rank_ensemble.py` |
| Candidate-vs-baseline admission | `src/plantpersulf/evaluation/model_admission.py` |
| Novel tomato candidate registry | `src/plantpersulf/proteomics/tomato_candidate_registry.py` |
| Tomato v2 orchestration | `src/plantpersulf/workflows/tomato_ranker_v2.py` |
| Tomato v2 CLI/config | `scripts/run_tomato_ranker_v2.py`, `configs/experiments/tomato_ranker_v2.yaml` |
| Target-label-free dependency firewall/workflow | `src/plantpersulf/workflows/cross_crop_target_label_free.py` |
| Cross-crop CLI/config | `scripts/run_cross_crop_target_label_free_v1.py`, `configs/experiments/cross_crop_target_label_free_v1.yaml` |
| Blind release policy and writer | `src/plantpersulf/reporting/candidate_release.py` |
| Blind release CLI/protocol | `scripts/build_candidate_release_v1.py`, `docs/wetlab_validation_matrix.md` |

## Completion Gate Used by Every Task

Run these commands after the task-specific GREEN test and before each commit:

```powershell
pytest tests/unit -q
pytest tests/scientific -q
pytest tests/release -q
python -m plantpersulf.cli audit-registry
python -m plantpersulf.cli audit-files
python -m plantpersulf.cli audit-benchmark --version v1
python -m plantpersulf.cli audit-leakage --split-version v1
ruff check .
mypy src/plantpersulf
```

Expected: every pytest suite passes; every audit exits `0`; Ruff and mypy report no errors. Network integration tests remain a separately approved operation because they access external repositories.

---

### Task 1: Biological Feature Contract

**Reviewer gate:** Approve only the feature schema and real-input extraction test. This task does not train a model.

**Files:**
- Create: `src/plantpersulf/features/site_biology.py`
- Create: `tests/unit/test_site_biology_policy.py`
- Create: `tests/scientific/test_site_biology_real_input.py`

**Interfaces:**
- Consumes: `plantpersulf.features.sequence.SequenceFeatureRow`
- Produces: `BIOLOGY_FEATURE_NAMES: tuple[str, ...]`, `SiteBiologyVector`, `build_site_biology_vector(row: SequenceFeatureRow) -> SiteBiologyVector`, `build_site_biology_vector_from_sequence(accession, position, sequence, radius=10) -> SiteBiologyVector`

- [ ] **Step 1: Write the policy RED test**

```python
from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES


def test_biology_features_exclude_shortcut_fields() -> None:
    assert BIOLOGY_FEATURE_NAMES == (
        "hydrophobicity",
        "protein_cys_density",
        "local_positive_charge_density",
        "local_negative_charge_density",
        "local_cys_density",
        "local_sequence_entropy",
    )
    forbidden = {"protein_length", "plddt", "has_structure", "study_id"}
    assert forbidden.isdisjoint(BIOLOGY_FEATURE_NAMES)
```

- [ ] **Step 2: Run the policy test and verify RED**

Run: `pytest tests/unit/test_site_biology_policy.py -v`

Expected: FAIL during import because `plantpersulf.features.site_biology` does not exist.

- [ ] **Step 3: Write a scientific RED test driven by a registered real sequence**

```python
from pathlib import Path

from plantpersulf.features.sequence import SequenceFeatureRow
from plantpersulf.features.site_biology import build_site_biology_vector
from plantpersulf.provenance.hashing import hash_file


def test_registered_q9zw96_builds_six_finite_core_features() -> None:
    fasta = Path("data/registry/cache/uniprot/Q9ZW96.fasta")
    assert hash_file(fasta, "sha256") == (
        "f4a9ba53775e0433e3efe92599d023a5f54f8cf8995806cba42db27a3bb78b19"
    )
    sequence = "".join(fasta.read_text(encoding="utf-8").splitlines()[1:])
    position = sequence.index("C") + 1
    radius = 10
    start = max(0, position - 1 - radius)
    end = min(len(sequence), position + radius)
    window = "X" * max(0, radius - position + 1) + sequence[start:end]
    window += "X" * (2 * radius + 1 - len(window))
    row = SequenceFeatureRow(
        protein_accession="Q9ZW96",
        cys_position=position,
        label="unlabeled",
        flanking_window=window,
        hydrophobicity=0.0,
        cys_density=sequence.count("C") / len(sequence),
        protein_length=len(sequence),
        local_positive_charge_density=sum(aa in "KR" for aa in window) / len(window),
    )
    vector = build_site_biology_vector(row)
    assert len(vector.values) == 6
    assert all(value == value for value in vector.values)
```

- [ ] **Step 4: Implement the minimal feature contract**

```python
from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass

from plantpersulf.features.sequence import KYTE_DOOLITTLE, SequenceFeatureRow

BIOLOGY_FEATURE_NAMES = (
    "hydrophobicity",
    "protein_cys_density",
    "local_positive_charge_density",
    "local_negative_charge_density",
    "local_cys_density",
    "local_sequence_entropy",
)


@dataclass(frozen=True)
class SiteBiologyVector:
    protein_accession: str
    cys_position: int
    values: tuple[float, ...]


def _density(window: str, residues: frozenset[str]) -> float:
    return sum(aa in residues for aa in window) / len(window) if window else 0.0


def _entropy(window: str) -> float:
    residues = [aa for aa in window if aa != "X"]
    if not residues:
        return 0.0
    counts = Counter(residues)
    return -sum((n / len(residues)) * math.log2(n / len(residues)) for n in counts.values())


def build_site_biology_vector(row: SequenceFeatureRow) -> SiteBiologyVector:
    values = (
        row.hydrophobicity,
        row.cys_density,
        row.local_positive_charge_density,
        _density(row.flanking_window, frozenset({"D", "E"})),
        _density(row.flanking_window, frozenset({"C"})),
        _entropy(row.flanking_window),
    )
    return SiteBiologyVector(row.protein_accession, row.cys_position, values)


def build_site_biology_vector_from_sequence(
    accession: str,
    position: int,
    sequence: str,
    radius: int = 10,
) -> SiteBiologyVector:
    if position < 1 or position > len(sequence) or sequence[position - 1] != "C":
        raise ValueError("candidate coordinate is not a cysteine")
    start = max(0, position - 1 - radius)
    end = min(len(sequence), position + radius)
    left = "X" * max(0, radius - position + 1)
    window = left + sequence[start:end]
    window += "X" * (2 * radius + 1 - len(window))
    hydrophobicity = sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in window) / len(window)
    row = SequenceFeatureRow(
        protein_accession=accession,
        cys_position=position,
        label="unlabeled",
        flanking_window=window,
        hydrophobicity=hydrophobicity,
        cys_density=sequence.count("C") / len(sequence),
        protein_length=len(sequence),
        local_positive_charge_density=_density(window, frozenset({"K", "R"})),
    )
    return build_site_biology_vector(row)
```

- [ ] **Step 5: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_site_biology_policy.py tests/scientific/test_site_biology_real_input.py -v`

Expected: PASS; the real input hash remains pinned and exactly six finite features are produced.

- [ ] **Step 6: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass without changing any scientific artifact.

- [ ] **Step 7: Commit Task 1 only**

```powershell
git add src/plantpersulf/features/site_biology.py tests/unit/test_site_biology_policy.py tests/scientific/test_site_biology_real_input.py
git commit -m "feat: define shortcut-free site biology features"
```

Stop and submit the TDD evidence template for reviewer approval.

---

### Task 2: Pairwise and Non-Negative PU Loss

**Reviewer gate:** Task 1 must be accepted. This task implements pure numerical loss functions and no biological output.

**Files:**
- Create: `src/plantpersulf/models/ranking_loss.py`
- Create: `tests/unit/test_ranking_loss.py`
- Create: `tests/scientific/test_unlabeled_never_becomes_negative.py`

**Interfaces:**
- Produces: `within_protein_pairs(labels, protein_ids)`, `nnpu_logistic_risk(logits, labels, class_prior)`, and `combined_pu_ranking_loss(logits, labels, protein_ids, class_prior, pairwise_weight)`
- Labels accepted: exact strings `positive` and `unlabeled`

- [ ] **Step 1: Write RED tests for PU semantics and pairing**

```python
import pytest
import torch

from plantpersulf.models.ranking_loss import (
    combined_pu_ranking_loss,
    within_protein_pairs,
)


def test_pairs_are_within_protein_and_positive_to_unlabeled() -> None:
    pairs = within_protein_pairs(
        ["positive", "unlabeled", "unlabeled", "positive"],
        ["group-a", "group-a", "group-b", "group-b"],
    )
    assert pairs == ((0, 1), (3, 2))


def test_negative_label_is_rejected() -> None:
    with pytest.raises(ValueError, match="forbidden"):
        combined_pu_ranking_loss(
            torch.tensor([0.0, 1.0]),
            ["positive", "negative"],
            ["group-a", "group-a"],
            class_prior=0.5,
            pairwise_weight=1.0,
        )
```

- [ ] **Step 2: Run tests and verify RED**

Run: `pytest tests/unit/test_ranking_loss.py -v`

Expected: FAIL during import because `ranking_loss.py` does not exist.

- [ ] **Step 3: Implement exact nnPU and pairwise loss**

```python
from __future__ import annotations

import torch

ALLOWED_LABELS = frozenset({"positive", "unlabeled"})


def _validate_labels(labels: list[str]) -> None:
    forbidden = set(labels) - ALLOWED_LABELS
    if forbidden:
        raise ValueError(f"forbidden label(s): {sorted(forbidden)}")


def within_protein_pairs(
    labels: list[str], protein_ids: list[str]
) -> tuple[tuple[int, int], ...]:
    _validate_labels(labels)
    if len(labels) != len(protein_ids):
        raise ValueError("labels and protein_ids length mismatch")
    return tuple(
        (i, j)
        for i, label_i in enumerate(labels)
        for j, label_j in enumerate(labels)
        if label_i == "positive"
        and label_j == "unlabeled"
        and protein_ids[i] == protein_ids[j]
    )


def nnpu_logistic_risk(
    logits: torch.Tensor, labels: list[str], class_prior: float
) -> torch.Tensor:
    _validate_labels(labels)
    if not 0.0 < class_prior < 1.0:
        raise ValueError("class_prior must be in (0, 1)")
    positive = torch.tensor([x == "positive" for x in labels], dtype=torch.bool)
    unlabeled = ~positive
    if not positive.any() or not unlabeled.any():
        raise ValueError("nnPU requires positive and unlabeled rows")
    pos_positive = torch.nn.functional.softplus(-logits[positive]).mean()
    pos_negative = torch.nn.functional.softplus(logits[positive]).mean()
    unl_negative = torch.nn.functional.softplus(logits[unlabeled]).mean()
    negative_risk = unl_negative - class_prior * pos_negative
    return class_prior * pos_positive + torch.clamp(negative_risk, min=0.0)


def combined_pu_ranking_loss(
    logits: torch.Tensor,
    labels: list[str],
    protein_ids: list[str],
    class_prior: float,
    pairwise_weight: float,
) -> torch.Tensor:
    risk = nnpu_logistic_risk(logits, labels, class_prior)
    pairs = within_protein_pairs(labels, protein_ids)
    if not pairs:
        return risk
    pair_loss = torch.stack(
        [torch.nn.functional.softplus(-(logits[i] - logits[j])) for i, j in pairs]
    ).mean()
    return risk + pairwise_weight * pair_loss
```

- [ ] **Step 4: Add the scientific policy test**

```python
from pathlib import Path


def test_new_model_source_has_no_hard_negative_label_literal() -> None:
    source = Path("src/plantpersulf/models/ranking_loss.py").read_text(encoding="utf-8")
    assert '"negative"' not in source
    assert "'negative'" not in source
```

- [ ] **Step 5: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_ranking_loss.py tests/scientific/test_unlabeled_never_becomes_negative.py -v`

Expected: PASS; invalid biological negative labels fail closed.

- [ ] **Step 6: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 7: Commit Task 2 only**

```powershell
git add src/plantpersulf/models/ranking_loss.py tests/unit/test_ranking_loss.py tests/scientific/test_unlabeled_never_becomes_negative.py
git commit -m "feat: add pairwise non-negative PU loss"
```

Stop and submit the TDD evidence template for reviewer approval.

---

### Task 3: Deterministic Additive PU Ranker

**Reviewer gate:** Task 2 must be accepted. No structure or cross-crop logic is allowed in this task.

**Files:**
- Create: `src/plantpersulf/models/additive_pu_ranker.py`
- Create: `tests/unit/test_additive_pu_ranker.py`
- Create: `tests/scientific/test_additive_ranker_train_only_scaling.py`

**Interfaces:**
- Consumes: numeric core feature rows, PU labels, protein IDs, `BIOLOGY_FEATURE_NAMES`
- Produces: `AdditivePuConfig`, `AdditivePuModel`, `fit_additive_pu_ranker(train_features, train_labels, train_protein_ids, feature_names, config) -> AdditivePuModel`, `AdditivePuModel.score(features) -> tuple[float, ...]`, and JSON-safe `to_dict`/`from_dict`

- [ ] **Step 1: Write the deterministic RED test**

```python
from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)


def test_same_seed_produces_identical_additive_scores() -> None:
    features = [[0.0, 1.0], [1.0, 0.0], [0.2, 0.8], [0.8, 0.2]]
    labels = ["positive", "unlabeled", "positive", "unlabeled"]
    proteins = ["group-a", "group-a", "group-b", "group-b"]
    config = AdditivePuConfig(class_prior=0.5, seed=7, epochs=50)
    first = fit_additive_pu_ranker(features, labels, proteins, ("f1", "f2"), config)
    second = fit_additive_pu_ranker(features, labels, proteins, ("f1", "f2"), config)
    assert first.score(features) == second.score(features)
    restored = type(first).from_dict(first.to_dict())
    assert restored.score(features) == first.score(features)
```

- [ ] **Step 2: Run tests and verify RED**

Run: `pytest tests/unit/test_additive_pu_ranker.py -v`

Expected: FAIL during import because `additive_pu_ranker.py` does not exist.

- [ ] **Step 3: Implement the train-only scaler and model contract**

```python
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import torch

from plantpersulf.models.ranking_loss import combined_pu_ranking_loss
from plantpersulf.models.traditional import TrainOnlyScaler


@dataclass(frozen=True)
class AdditivePuConfig:
    class_prior: float
    seed: int
    epochs: int = 200
    learning_rate: float = 0.02
    pairwise_weight: float = 1.0
    l1: float = 0.001
    l2: float = 0.001


@dataclass(frozen=True)
class AdditivePuModel:
    feature_names: tuple[str, ...]
    scaler: TrainOnlyScaler
    weights: tuple[float, ...]
    intercept: float

    def score(self, features: list[list[float]]) -> tuple[float, ...]:
        scaled = self.scaler.transform(features)
        return tuple(
            self.intercept + sum(w * x for w, x in zip(self.weights, row, strict=True))
            for row in scaled
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "feature_names": list(self.feature_names),
            "scaler_mean": list(self.scaler.mean),
            "scaler_std": list(self.scaler.std),
            "weights": list(self.weights),
            "intercept": self.intercept,
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "AdditivePuModel":
        return cls(
            feature_names=tuple(str(x) for x in value["feature_names"]),
            scaler=TrainOnlyScaler(
                mean=tuple(float(x) for x in value["scaler_mean"]),
                std=tuple(float(x) for x in value["scaler_std"]),
            ),
            weights=tuple(float(x) for x in value["weights"]),
            intercept=float(value["intercept"]),
        )


def fit_additive_pu_ranker(
    train_features: list[list[float]],
    train_labels: list[str],
    train_protein_ids: list[str],
    feature_names: tuple[str, ...],
    config: AdditivePuConfig,
) -> AdditivePuModel:
    if not train_features:
        raise ValueError("train_features must not be empty")
    if not (
        len(train_features) == len(train_labels) == len(train_protein_ids)
    ):
        raise ValueError("training row count mismatch")
    if len(feature_names) != len(train_features[0]) or any(
        len(row) != len(feature_names) for row in train_features
    ):
        raise ValueError("feature_names width mismatch")
    torch.manual_seed(config.seed)
    torch.use_deterministic_algorithms(True)
    scaler = TrainOnlyScaler.fit(train_features)
    x = torch.tensor(scaler.transform(train_features), dtype=torch.float32)
    weights = torch.zeros(x.shape[1], requires_grad=True)
    intercept = torch.zeros(1, requires_grad=True)
    optimizer = torch.optim.Adam([weights, intercept], lr=config.learning_rate)
    for _ in range(config.epochs):
        optimizer.zero_grad()
        logits = x @ weights + intercept
        loss = combined_pu_ranking_loss(
            logits,
            train_labels,
            train_protein_ids,
            config.class_prior,
            config.pairwise_weight,
        )
        loss = loss + config.l1 * weights.abs().sum() + config.l2 * weights.square().sum()
        loss.backward()
        optimizer.step()
    return AdditivePuModel(
        feature_names=feature_names,
        scaler=scaler,
        weights=tuple(float(v) for v in weights.detach()),
        intercept=float(intercept.detach().item()),
    )
```

- [ ] **Step 4: Write the train-only scaling scientific test**

```python
from plantpersulf.models.additive_pu_ranker import (
    AdditivePuConfig,
    fit_additive_pu_ranker,
)


def test_prediction_distribution_cannot_change_fitted_scaler() -> None:
    model = fit_additive_pu_ranker(
        [[0.0], [1.0], [0.2], [0.8]],
        ["positive", "unlabeled", "positive", "unlabeled"],
        ["a", "a", "b", "b"],
        ("marker",),
        AdditivePuConfig(class_prior=0.5, seed=3, epochs=20),
    )
    frozen = (model.scaler.mean, model.scaler.std)
    model.score([[1_000_000.0], [-1_000_000.0]])
    assert (model.scaler.mean, model.scaler.std) == frozen
```

- [ ] **Step 5: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_additive_pu_ranker.py tests/scientific/test_additive_ranker_train_only_scaling.py -v`

Expected: PASS with bit-for-bit identical scores for a fixed seed.

- [ ] **Step 6: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 7: Commit Task 3 only**

```powershell
git add src/plantpersulf/models/additive_pu_ranker.py tests/unit/test_additive_pu_ranker.py tests/scientific/test_additive_ranker_train_only_scaling.py
git commit -m "feat: add deterministic additive PU ranker"
```

Stop and submit the TDD evidence template for reviewer approval.

---

### Task 4: Observation-Propensity Fail-Closed Policy

**Reviewer gate:** Task 3 must be accepted. This task does not fit SAR-PU because no complete registered discovery-observation table is currently available.

**Files:**
- Create: `src/plantpersulf/models/observation_propensity.py`
- Create: `tests/unit/test_observation_propensity_policy.py`
- Create: `tests/scientific/test_observation_propensity_requires_registered_table.py`

**Interfaces:**
- Produces: `ObservationPropensityStatus`, `resolve_observation_propensity(table_path, registry_path) -> ObservationPropensityStatus`
- Current valid release state: `enabled=False`, reason `complete_discovery_table_not_supplied`

- [ ] **Step 1: Write RED tests for absent and unregistered tables**

```python
from pathlib import Path

import pytest

from plantpersulf.models.observation_propensity import resolve_observation_propensity


def test_missing_table_disables_propensity_without_substitution() -> None:
    status = resolve_observation_propensity(None, None)
    assert status.enabled is False
    assert status.reason == "complete_discovery_table_not_supplied"


def test_unregistered_table_is_rejected(tmp_path: Path) -> None:
    marker = tmp_path / "policy-marker.tsv"
    marker.write_text("policy-marker\n", encoding="utf-8")
    with pytest.raises(RuntimeError, match="registry"):
        resolve_observation_propensity(marker, None)
```

- [ ] **Step 2: Run tests and verify RED**

Run: `pytest tests/unit/test_observation_propensity_policy.py -v`

Expected: FAIL during import because the policy module does not exist.

- [ ] **Step 3: Implement the fail-closed status resolver**

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from plantpersulf.provenance.audit import assert_registered_input


@dataclass(frozen=True)
class ObservationPropensityStatus:
    enabled: bool
    reason: str
    source_sha256: str | None


def resolve_observation_propensity(
    table_path: Path | None,
    registry_path: Path | None,
) -> ObservationPropensityStatus:
    if table_path is None:
        return ObservationPropensityStatus(
            enabled=False,
            reason="complete_discovery_table_not_supplied",
            source_sha256=None,
        )
    if registry_path is None:
        raise RuntimeError("observation table requires a path/SHA256 registry")
    assert_registered_input(table_path, registry_path)
    raise RuntimeError(
        "registered table supplied but fitted SAR-PU is outside the approved task"
    )
```

- [ ] **Step 4: Add source scanning for forbidden fallback behavior**

```python
from pathlib import Path


def test_propensity_module_has_no_imputation_or_random_fallback() -> None:
    source = Path("src/plantpersulf/models/observation_propensity.py").read_text(
        encoding="utf-8"
    )
    assert "random" not in source.lower()
    assert "imput" not in source.lower()
```

- [ ] **Step 5: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_observation_propensity_policy.py tests/scientific/test_observation_propensity_requires_registered_table.py -v`

Expected: PASS; no table yields an explicit disabled status and supplied-but-unapproved data fail closed.

- [ ] **Step 6: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 7: Commit Task 4 only**

```powershell
git add src/plantpersulf/models/observation_propensity.py tests/unit/test_observation_propensity_policy.py tests/scientific/test_observation_propensity_requires_registered_table.py
git commit -m "feat: enforce fail-closed observation propensity policy"
```

Stop and submit the TDD evidence template for reviewer approval.

---

### Task 5: Optional Structure Residual and Admission

**Reviewer gate:** Tasks 1–4 must be accepted. This task may implement a residual arm but cannot enable it in a release.

**Files:**
- Create: `src/plantpersulf/models/structure_residual.py`
- Create: `src/plantpersulf/evaluation/structure_admission.py`
- Create: `tests/unit/test_structure_residual.py`
- Create: `tests/scientific/test_structure_residual_cannot_learn_coverage.py`

**Interfaces:**
- Consumes: frozen sequence scores, up to three structure scalar features, quality mask, PU labels, protein IDs
- Produces: `score_structure_residual(weights, intercept, structure_features, quality_mask)`, `StructureAdmissionDecision`, `admit_structure_residual(paired_structure_deltas, paired_coverage_only_deltas)`

- [ ] **Step 1: Write RED tests for missing-structure neutrality**

```python
from plantpersulf.models.structure_residual import score_structure_residual


def test_missing_structure_gets_exact_zero_residual() -> None:
    residual = score_structure_residual(
        weights=(1.0, -1.0),
        intercept=0.5,
        structure_features=[[2.0, 1.0], [0.0, 0.0]],
        quality_mask=[True, False],
    )
    assert residual[1] == 0.0
```

- [ ] **Step 2: Write the coverage-only rejection RED test**

```python
from plantpersulf.evaluation.structure_admission import admit_structure_residual


def test_coverage_only_gain_blocks_structure_admission() -> None:
    decision = admit_structure_residual(
        paired_structure_deltas=(0.03, 0.02, 0.04, 0.01, 0.03),
        paired_coverage_only_deltas=(0.03, 0.02, 0.04, 0.01, 0.03),
    )
    assert decision.admitted is False
    assert decision.delta == 0
```

- [ ] **Step 3: Run tests and verify RED**

Run: `pytest tests/unit/test_structure_residual.py tests/scientific/test_structure_residual_cannot_learn_coverage.py -v`

Expected: FAIL during import because both modules are absent.

- [ ] **Step 4: Implement missing-neutral residual scoring**

```python
from __future__ import annotations


def score_structure_residual(
    weights: tuple[float, ...],
    intercept: float,
    structure_features: list[list[float]],
    quality_mask: list[bool],
) -> tuple[float, ...]:
    if len(structure_features) != len(quality_mask):
        raise ValueError("structure features and quality mask length mismatch")
    if len(weights) > 3:
        raise ValueError("structure residual is limited to three features")
    return tuple(
        intercept + sum(w * x for w, x in zip(weights, row, strict=True))
        if available
        else 0.0
        for row, available in zip(structure_features, quality_mask, strict=True)
    )
```

- [ ] **Step 5: Implement the explicit admission rule**

```python
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class StructureAdmissionDecision:
    admitted: bool
    delta: int
    reason: str


def admit_structure_residual(
    paired_structure_deltas: tuple[float, ...],
    paired_coverage_only_deltas: tuple[float, ...],
) -> StructureAdmissionDecision:
    if len(paired_structure_deltas) < 5:
        return StructureAdmissionDecision(False, 0, "insufficient_repeated_folds")
    structure_stable = min(paired_structure_deltas) > 0.0
    coverage_has_gain = max(paired_coverage_only_deltas) > 0.0
    if structure_stable and not coverage_has_gain:
        return StructureAdmissionDecision(True, 1, "matched_structure_gain")
    return StructureAdmissionDecision(False, 0, "coverage_or_instability_failure")
```

- [ ] **Step 6: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_structure_residual.py tests/scientific/test_structure_residual_cannot_learn_coverage.py -v`

Expected: PASS; every missing-structure row gets zero residual and coverage-only gain blocks admission.

- [ ] **Step 7: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 8: Commit Task 5 only**

```powershell
git add src/plantpersulf/models/structure_residual.py src/plantpersulf/evaluation/structure_admission.py tests/unit/test_structure_residual.py tests/scientific/test_structure_residual_cannot_learn_coverage.py
git commit -m "feat: add gated structure residual policy"
```

Stop and submit the TDD evidence template for reviewer approval.

---

### Task 6: Rank Uncertainty, Applicability, and Model Admission

**Reviewer gate:** Task 5 must be accepted. The admission layer consumes only out-of-fold predictions.

**Files:**
- Create: `src/plantpersulf/evaluation/rank_ensemble.py`
- Create: `src/plantpersulf/evaluation/model_admission.py`
- Create: `tests/unit/test_rank_ensemble.py`
- Create: `tests/scientific/test_model_admission_uses_paired_oof_runs.py`

**Interfaces:**
- Produces: `RankSummary`, `summarize_rank_runs(site_keys, run_scores)`, `ApplicabilityEnvelope.fit(train_features)`, `AdmissionDecision`, `admit_candidate_model(candidate_runs, baseline_runs)`

- [ ] **Step 1: Write the RED tests**

```python
from plantpersulf.evaluation.model_admission import admit_candidate_model
from plantpersulf.evaluation.rank_ensemble import summarize_rank_runs


def test_rank_summary_uses_run_distribution() -> None:
    summary = summarize_rank_runs(
        ("site-a", "site-b"),
        ((0.9, 0.1), (0.8, 0.2), (0.7, 0.3)),
    )
    assert summary[0].site_key == "site-a"
    assert summary[0].median_rank == 1.0


def test_candidate_must_improve_every_required_metric() -> None:
    candidate = {
        f"r{i}f0": (0.3, 1.5, 0.1 if i == 0 else 0.3)
        for i in range(5)
    }
    baseline = {f"r{i}f0": (0.1, 1.2, 0.2) for i in range(5)}
    blocks = {f"r{i}f0": f"r{i}" for i in range(5)}
    decision = admit_candidate_model(candidate, baseline, blocks, n_boot=500)
    assert decision.admitted is False
    assert "mrr" in decision.reason
```

- [ ] **Step 2: Run tests and verify RED**

Run: `pytest tests/unit/test_rank_ensemble.py tests/scientific/test_model_admission_uses_paired_oof_runs.py -v`

Expected: FAIL during import because the modules do not exist.

- [ ] **Step 3: Implement deterministic rank summaries**

```python
from __future__ import annotations

from dataclasses import dataclass
from statistics import median


@dataclass(frozen=True)
class RankSummary:
    site_key: str
    median_score: float
    median_rank: float
    best_rank: int
    worst_rank: int


def summarize_rank_runs(
    site_keys: tuple[str, ...], run_scores: tuple[tuple[float, ...], ...]
) -> tuple[RankSummary, ...]:
    ranks_by_run = []
    for scores in run_scores:
        order = sorted(range(len(scores)), key=lambda i: (-scores[i], site_keys[i]))
        rank = {idx: place for place, idx in enumerate(order, start=1)}
        ranks_by_run.append(rank)
    return tuple(
        RankSummary(
            site_key=key,
            median_score=median(run[i] for run in run_scores),
            median_rank=median(ranks[i] for ranks in ranks_by_run),
            best_rank=min(ranks[i] for ranks in ranks_by_run),
            worst_rank=max(ranks[i] for ranks in ranks_by_run),
        )
        for i, key in enumerate(site_keys)
    )
```

- [ ] **Step 4: Implement paired admission without test-set access**

```python
from __future__ import annotations

import random
from dataclasses import dataclass


@dataclass(frozen=True)
class AdmissionDecision:
    admitted: bool
    reason: str
    paired_intervals: dict[str, tuple[float, float]]


def _percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    return ordered[int((len(ordered) - 1) * fraction)]


def _paired_block_interval(
    deltas_by_block: list[float], n_boot: int, seed: int
) -> tuple[float, float]:
    if len(deltas_by_block) < 5:
        raise ValueError("at least five repeated-CV blocks are required")
    if n_boot < 500:
        raise ValueError("n_boot must be at least 500")
    rng = random.Random(seed)
    means = [
        sum(rng.choice(deltas_by_block) for _ in deltas_by_block)
        / len(deltas_by_block)
        for _ in range(n_boot)
    ]
    return _percentile(means, 0.025), _percentile(means, 0.975)


def admit_candidate_model(
    candidate_runs: dict[str, tuple[float, float, float]],
    baseline_runs: dict[str, tuple[float, float, float]],
    run_blocks: dict[str, str],
    n_boot: int = 2_000,
    seed: int = 20_260_811,
) -> AdmissionDecision:
    run_keys = set(candidate_runs)
    if run_keys != set(baseline_runs) or run_keys != set(run_blocks):
        raise ValueError("candidate, baseline, and block keys must match")
    metric_names = ("recall_at_k", "enrichment", "mrr")
    intervals: dict[str, tuple[float, float]] = {}
    for metric_index, metric_name in enumerate(metric_names):
        block_deltas = [
            sum(
                candidate_runs[key][metric_index]
                - baseline_runs[key][metric_index]
                for key in sorted(run_keys)
                if run_blocks[key] == block
            )
            / sum(run_blocks[key] == block for key in run_keys)
            for block in sorted(set(run_blocks.values()))
        ]
        interval = _paired_block_interval(block_deltas, n_boot, seed + metric_index)
        intervals[metric_name] = interval
        leave_one_block_out = [
            sum(value for j, value in enumerate(block_deltas) if j != i)
            / (len(block_deltas) - 1)
            for i in range(len(block_deltas))
        ]
        if interval[0] <= 0.0 or min(leave_one_block_out) <= 0.0:
            return AdmissionDecision(
                False,
                f"{metric_name}_conservative_interval_or_stability_failed",
                intervals,
            )
    return AdmissionDecision(
        True,
        "all_required_paired_intervals_exclude_zero",
        intervals,
    )
```

- [ ] **Step 5: Add `ApplicabilityEnvelope` using the train-fold feature range only**

```python
@dataclass(frozen=True)
class ApplicabilityEnvelope:
    lower: tuple[float, ...]
    upper: tuple[float, ...]

    @classmethod
    def fit(cls, train_features: list[list[float]]) -> "ApplicabilityEnvelope":
        columns = list(zip(*train_features, strict=True))
        return cls(
            lower=tuple(min(column) for column in columns),
            upper=tuple(max(column) for column in columns),
        )

    def contains(self, row: list[float]) -> bool:
        return all(lo <= value <= hi for value, lo, hi in zip(row, self.lower, self.upper, strict=True))
```

- [ ] **Step 6: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_rank_ensemble.py tests/scientific/test_model_admission_uses_paired_oof_runs.py -v`

Expected: PASS; admission uses matched OOF run IDs and fails if any required metric is unstable.

- [ ] **Step 7: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 8: Commit Task 6 only**

```powershell
git add src/plantpersulf/evaluation/rank_ensemble.py src/plantpersulf/evaluation/model_admission.py tests/unit/test_rank_ensemble.py tests/scientific/test_model_admission_uses_paired_oof_runs.py
git commit -m "feat: add rank uncertainty and admission gates"
```

Stop and submit the TDD evidence template for reviewer approval.

---

### Task 7: Novel Tomato Candidate Registry

**Reviewer gate:** Tasks 1–6 must be accepted, and the separate tomato-input provenance task described in Global Constraints must have been reviewed. This task enumerates candidates only; it does not acquire, score, or release them.

**Files:**
- Create: `src/plantpersulf/proteomics/tomato_candidate_registry.py`
- Create: `tests/unit/test_tomato_candidate_registry_contract.py`
- Create: `tests/scientific/test_tomato_candidate_registry_real_inputs.py`

**Interfaces:**
- Consumes: registered KIAE271 supplement and tomato reference proteome
- Produces: `TomatoCandidate`, `TomatoCandidateRegistry`, `build_tomato_candidate_registry(kiae271_xlsx, proteome_path, supplementary_registry, model_input_registry)`
- Excludes: all coordinate-verified known positives and every ambiguous/unverified KIAE271-touched Cys

- [ ] **Step 1: Write the registry contract RED test**

```python
from dataclasses import fields

from plantpersulf.proteomics.tomato_candidate_registry import TomatoCandidate


def test_candidate_contract_has_no_label_field() -> None:
    assert tuple(field.name for field in fields(TomatoCandidate)) == (
        "protein_accession",
        "cys_position",
        "site_key",
    )
```

- [ ] **Step 2: Run the unit test and verify RED**

Run: `pytest tests/unit/test_tomato_candidate_registry_contract.py -v`

Expected: FAIL during import because `tomato_candidate_registry.py` does not exist.

- [ ] **Step 3: Implement deterministic enumeration and provenance**

```python
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from plantpersulf.features.sequence import _load_proteome
from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites
from plantpersulf.proteomics.tomato_local_dataset import kiae271_excluded_keys
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file


@dataclass(frozen=True)
class TomatoCandidate:
    protein_accession: str
    cys_position: int
    site_key: str


@dataclass(frozen=True)
class TomatoCandidateRegistry:
    candidates: tuple[TomatoCandidate, ...]
    kiae271_sha256: str
    proteome_sha256: str
    excluded_known_positive_count: int
    excluded_ambiguous_count: int


def build_tomato_candidate_registry(
    kiae271_xlsx: Path,
    proteome_path: Path,
    supplementary_registry: Path,
    model_input_registry: Path,
) -> TomatoCandidateRegistry:
    assert_registered_input(kiae271_xlsx, supplementary_registry)
    assert_registered_input(proteome_path, model_input_registry)
    proteome = _load_proteome(proteome_path)
    positives = {
        (site.protein_accession, site.cys_position)
        for site in parse_kiae271_sites(kiae271_xlsx, proteome).sites
    }
    ambiguous = set(kiae271_excluded_keys(kiae271_xlsx, proteome))
    candidates = tuple(
        TomatoCandidate(accession, position, f"{accession}:C{position}")
        for accession, sequence in sorted(proteome.items())
        for position, residue in enumerate(sequence, start=1)
        if residue == "C"
        and (accession, position) not in positives
        and (accession, position) not in ambiguous
    )
    return TomatoCandidateRegistry(
        candidates=candidates,
        kiae271_sha256=hash_file(kiae271_xlsx, "sha256"),
        proteome_sha256=hash_file(proteome_path, "sha256"),
        excluded_known_positive_count=len(positives),
        excluded_ambiguous_count=len(ambiguous),
    )
```

- [ ] **Step 4: Write the registered-real-input scientific test**

```python
from pathlib import Path

from plantpersulf.features.sequence import _load_proteome
from plantpersulf.proteomics.kiae271_sites import parse_kiae271_sites
from plantpersulf.proteomics.tomato_candidate_registry import (
    build_tomato_candidate_registry,
)
from plantpersulf.provenance.audit import assert_registered_input


def test_real_registry_excludes_every_known_positive_and_keeps_only_cys() -> None:
    xlsx = Path("data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx")
    fasta = Path("data/raw/references/tomato_ref_proteome_v1.fasta")
    supplementary_registry = Path("data/registry/supplementary_sources.tsv")
    model_input_registry = Path("data/registry/model_inputs.tsv")
    assert_registered_input(xlsx, supplementary_registry)
    assert_registered_input(fasta, model_input_registry)
    proteome = _load_proteome(fasta)
    positive_keys = {
        (site.protein_accession, site.cys_position)
        for site in parse_kiae271_sites(xlsx, proteome).sites
    }
    registry = build_tomato_candidate_registry(
        xlsx, fasta, supplementary_registry, model_input_registry
    )
    keys = {(row.protein_accession, row.cys_position) for row in registry.candidates}
    assert keys.isdisjoint(positive_keys)
    assert all(proteome[accession][position - 1] == "C" for accession, position in keys)
```

- [ ] **Step 5: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_tomato_candidate_registry_contract.py tests/scientific/test_tomato_candidate_registry_real_inputs.py -v`

Expected: PASS; every candidate is a real tomato Cys and no known/ambiguous KIAE271 site enters the pool.

- [ ] **Step 6: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 7: Commit Task 7 only**

```powershell
git add src/plantpersulf/proteomics/tomato_candidate_registry.py tests/unit/test_tomato_candidate_registry_contract.py tests/scientific/test_tomato_candidate_registry_real_inputs.py
git commit -m "feat: add provenance-locked tomato candidate registry"
```

Stop and submit the TDD evidence template for reviewer approval.

---

### Task 8: Tomato v2 Repeated-Cluster Workflow and Frozen Score Release

**Reviewer gate:** Tasks 1–7 must be accepted. The reviewer must separately approve running the real experiment after approving its workflow implementation.

**Files:**
- Create: `src/plantpersulf/workflows/__init__.py`
- Create: `src/plantpersulf/workflows/tomato_ranker_v2.py`
- Create: `scripts/run_tomato_ranker_v2.py`
- Create: `configs/experiments/tomato_ranker_v2.yaml`
- Create: `tests/unit/test_tomato_ranker_v2_config.py`
- Create: `tests/scientific/test_tomato_ranker_v2_split_and_prior_policy.py`
- Create: `tests/release/test_tomato_ranker_v2_cannot_change_gate2.py`

**Interfaces:**
- Consumes: `build_tomato_pu_rows`, `assign_grouped_folds`, Task 1 features, Task 3 model, Task 4 status, Task 5/6 gates
- Produces: `select_release_model(candidate_admitted)`, `run_tomato_ranker_v2(config_path: Path, output_dir: Path) -> dict[str, object]`, per-site OOF scores, `summary.json`, `manifest.json`

- [ ] **Step 1: Write config-policy RED tests**

```python
from pathlib import Path

import yaml

from plantpersulf.workflows.tomato_ranker_v2 import select_release_model


def test_tomato_v2_config_freezes_primary_arena_and_structure_default() -> None:
    cfg = yaml.safe_load(Path("configs/experiments/tomato_ranker_v2.yaml").read_text())
    assert cfg["evaluation"]["primary_arena"] == "panel"
    assert cfg["evaluation"]["split"] == "repeated_homology_cluster_cv"
    assert cfg["evaluation"]["folds"] == 5
    assert cfg["evaluation"]["repetitions"] == 5
    assert cfg["structure"]["delta_default"] == 0
    assert cfg["observation_propensity"]["enabled"] is False
    assert cfg["models"]["fallback_on_admission_failure"] is True
    assert cfg["claim_class"] == "development_candidate_ranking_not_gate2"
    assert select_release_model(False) == "pu_logistic"
    assert select_release_model(True) == "additive_pu"
```

- [ ] **Step 2: Run tests and verify RED**

Run: `pytest tests/unit/test_tomato_ranker_v2_config.py -v`

Expected: FAIL because the config does not exist.

- [ ] **Step 3: Create the frozen configuration**

```yaml
experiment:
  name: tomato_ranker_v2
claim_class: development_candidate_ranking_not_gate2
data:
  kiae271_xlsx: data/raw/supplements/KIAE271_SUPPL/kiae271_DSs.xlsx
  reference_proteome: data/raw/references/tomato_ref_proteome_v1.fasta
  cluster_file: data/processed/clusters/tomato_proteome_clusters_v1.tsv
  supplementary_registry: data/registry/supplementary_sources.tsv
  model_input_registry: data/registry/model_inputs.tsv
  arena_ratio: 20
  subsample_seed: 12345
features:
  core:
    - hydrophobicity
    - protein_cys_density
    - local_positive_charge_density
    - local_negative_charge_density
    - local_cys_density
    - local_sequence_entropy
evaluation:
  primary_arena: panel
  audit_arena: proteome
  split: repeated_homology_cluster_cv
  folds: 5
  repetitions: 5
  model_seeds: [0, 1, 2, 3, 4]
  development_k: 50
pu:
  prior_multipliers: [1.0, 1.5, 2.0]
  prior_cap: 0.5
models:
  candidate: additive_pu
  baseline: pu_logistic
  fallback_on_admission_failure: true
structure:
  delta_default: 0
observation_propensity:
  enabled: false
output:
  include_per_site_oof_scores: true
  include_input_sha256: true
```

- [ ] **Step 4: Write the split/prior scientific RED test**

```python
from plantpersulf.workflows.tomato_ranker_v2 import sensitivity_priors


def test_prior_grid_starts_at_observed_positive_fraction() -> None:
    assert sensitivity_priors(10, 100, (1.0, 1.5, 2.0), 0.5) == (0.1, 0.15, 0.2)
```

- [ ] **Step 5: Implement deterministic config validation, input loading, and the prior grid**

```python
from __future__ import annotations

import csv
import json
import statistics
from pathlib import Path
from typing import Any

import yaml

from plantpersulf.evaluation.metrics import average_precision, mean_reciprocal_rank, recall_at_k
from plantpersulf.evaluation.model_admission import admit_candidate_model
from plantpersulf.evaluation.rank_ensemble import ApplicabilityEnvelope, summarize_rank_runs
from plantpersulf.features.sequence import _load_proteome
from plantpersulf.features.site_biology import (
    BIOLOGY_FEATURE_NAMES,
    build_site_biology_vector_from_sequence,
)
from plantpersulf.models.additive_pu_ranker import AdditivePuConfig, fit_additive_pu_ranker
from plantpersulf.models.traditional import pu_logistic_regression_scores
from plantpersulf.proteomics.tomato_local_dataset import (
    ARENA_PANEL,
    ARENA_PROTEOME,
    assign_grouped_folds,
    build_tomato_pu_rows,
)
from plantpersulf.proteomics.tomato_candidate_registry import build_tomato_candidate_registry
from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file


def sensitivity_priors(
    n_positive: int,
    n_total: int,
    multipliers: tuple[float, ...],
    cap: float,
) -> tuple[float, ...]:
    observed = n_positive / n_total
    return tuple(sorted({min(cap, observed * multiplier) for multiplier in multipliers}))


def _cluster_map(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return {
            row["protein_accession"]: row["cluster_id"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def _matrix(rows: list[Any], proteome: dict[str, str]) -> list[list[float]]:
    return [
        list(
            build_site_biology_vector_from_sequence(
                row.protein_accession,
                row.cys_position,
                proteome[row.protein_accession],
            ).values
        )
        for row in rows
    ]


def _run_metrics(scores: list[float], labels: list[str], k: int) -> tuple[float, float, float]:
    pairs = list(zip(scores, labels, strict=True))
    ap = average_precision(pairs)
    recall = recall_at_k(pairs, k)
    mrr = mean_reciprocal_rank(pairs)
    if ap is None or recall is None or mrr is None:
        raise RuntimeError("held-out fold lacks a positive")
    observed = labels.count("positive") / len(labels)
    return recall, ap / observed, mrr


def _lower_quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    index = int((len(ordered) - 1) * fraction)
    return ordered[index]


def _percentile_column(scores: tuple[float, ...]) -> tuple[float, ...]:
    order = sorted(range(len(scores)), key=lambda i: (-scores[i], i))
    denominator = max(1, len(scores) - 1)
    percentile = [0.0] * len(scores)
    for rank, index in enumerate(order):
        percentile[index] = 1.0 - rank / denominator
    return tuple(percentile)


def select_release_model(candidate_admitted: bool) -> str:
    return "additive_pu" if candidate_admitted else "pu_logistic"


def run_tomato_ranker_v2(config_path: Path, output_dir: Path) -> dict[str, object]:
    cfg: dict[str, Any] = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    if cfg["claim_class"] != "development_candidate_ranking_not_gate2":
        raise RuntimeError("tomato v2 cannot emit a Gate 2 claim")
    if (
        cfg["evaluation"]["primary_arena"] != "panel"
        or cfg["evaluation"]["split"] != "repeated_homology_cluster_cv"
        or int(cfg["evaluation"]["folds"]) != 5
        or int(cfg["evaluation"]["repetitions"]) != 5
    ):
        raise RuntimeError("primary arena or split differs from the frozen policy")
    if cfg["structure"]["delta_default"] != 0:
        raise RuntimeError("structure residual must default to disabled")
    if cfg["observation_propensity"]["enabled"] is not False:
        raise RuntimeError("observation propensity lacks an admissible input table")
    if (
        cfg["models"]["candidate"] != "additive_pu"
        or cfg["models"]["baseline"] != "pu_logistic"
        or cfg["models"]["fallback_on_admission_failure"] is not True
    ):
        raise RuntimeError("model comparison or fallback policy differs from frozen config")
    if output_dir.exists():
        raise FileExistsError(output_dir)
    xlsx_path = Path(cfg["data"]["kiae271_xlsx"])
    proteome_path = Path(cfg["data"]["reference_proteome"])
    cluster_path = Path(cfg["data"]["cluster_file"])
    supplementary_registry = Path(cfg["data"]["supplementary_registry"])
    model_input_registry = Path(cfg["data"]["model_input_registry"])
    for path in (
        xlsx_path,
        proteome_path,
        cluster_path,
        supplementary_registry,
        model_input_registry,
    ):
        if not path.is_file():
            raise RuntimeError(f"required registered input missing: {path}")
    assert_registered_input(xlsx_path, supplementary_registry)
    assert_registered_input(proteome_path, model_input_registry)
    assert_registered_input(cluster_path, model_input_registry)
    proteome = _load_proteome(proteome_path)
    clusters = _cluster_map(cluster_path)
    missing_clusters = sorted(set(proteome) - set(clusters))
    if missing_clusters:
        raise RuntimeError(
            f"registered cluster table misses {len(missing_clusters)} proteins"
        )
    feature_names = tuple(cfg["features"]["core"])
    if feature_names != BIOLOGY_FEATURE_NAMES:
        raise RuntimeError("core feature order differs from the frozen biological contract")
    candidate_runs: dict[str, tuple[float, float, float]] = {}
    baseline_runs: dict[str, tuple[float, float, float]] = {}
    oof_rows: list[dict[str, object]] = []
    for arena in (ARENA_PANEL, ARENA_PROTEOME):
        rows = build_tomato_pu_rows(
            xlsx_path,
            proteome,
            arena=arena,
            ratio=int(cfg["data"]["arena_ratio"]),
            seed=int(cfg["data"]["subsample_seed"]),
        )
        features = _matrix(rows, proteome)
        for repetition in range(int(cfg["evaluation"]["repetitions"])):
            folded = assign_grouped_folds(
                rows,
                clusters,
                n_folds=int(cfg["evaluation"]["folds"]),
                seed=int(cfg["data"]["subsample_seed"]) + repetition,
            )
            for fold in range(int(cfg["evaluation"]["folds"])):
                train_idx = [i for i, row in enumerate(folded) if row.fold != fold]
                test_idx = [i for i, row in enumerate(folded) if row.fold == fold]
                train_x = [features[i] for i in train_idx]
                test_x = [features[i] for i in test_idx]
                train_y = [folded[i].label for i in train_idx]
                test_y = [folded[i].label for i in test_idx]
                train_proteins = [folded[i].protein_accession for i in train_idx]
                envelope = ApplicabilityEnvelope.fit(train_x)
                priors = sensitivity_priors(
                    train_y.count("positive"),
                    len(train_y),
                    tuple(float(x) for x in cfg["pu"]["prior_multipliers"]),
                    float(cfg["pu"]["prior_cap"]),
                )
                candidate_columns: list[tuple[float, ...]] = []
                baseline_columns: list[list[float]] = []
                for seed in cfg["evaluation"]["model_seeds"]:
                    baseline_columns.append(
                        pu_logistic_regression_scores(train_x, train_y, test_x, int(seed))
                    )
                    for prior in priors:
                        model = fit_additive_pu_ranker(
                            train_x,
                            train_y,
                            train_proteins,
                            feature_names,
                            AdditivePuConfig(class_prior=prior, seed=int(seed)),
                        )
                        candidate_columns.append(model.score(test_x))
                candidate_scores = [
                    statistics.median(column[i] for column in candidate_columns)
                    for i in range(len(test_idx))
                ]
                baseline_scores = [
                    statistics.median(column[i] for column in baseline_columns)
                    for i in range(len(test_idx))
                ]
                run_key = f"{arena}_r{repetition}f{fold}"
                candidate_runs[run_key] = _run_metrics(
                    candidate_scores, test_y, int(cfg["evaluation"]["development_k"])
                )
                baseline_runs[run_key] = _run_metrics(
                    baseline_scores, test_y, int(cfg["evaluation"]["development_k"])
                )
                for local_index, row_index in enumerate(test_idx):
                    row = folded[row_index]
                    oof_rows.append(
                        {
                            "arena": arena,
                            "run": run_key,
                            "protein_accession": row.protein_accession,
                            "cys_position": row.cys_position,
                            "label": row.label,
                            "candidate_score": candidate_scores[local_index],
                            "baseline_score": baseline_scores[local_index],
                            "in_domain": envelope.contains(test_x[local_index]),
                        }
                    )
    panel_candidate = {key: value for key, value in candidate_runs.items() if key.startswith("panel_")}
    panel_baseline = {key: value for key, value in baseline_runs.items() if key.startswith("panel_")}
    panel_blocks = {
        key: key.rsplit("f", maxsplit=1)[0]
        for key in panel_candidate
    }
    decision = admit_candidate_model(
        panel_candidate,
        panel_baseline,
        panel_blocks,
        seed=int(cfg["data"]["subsample_seed"]),
    )
    rank_summaries: dict[str, list[dict[str, object]]] = {}
    for arena in (ARENA_PANEL, ARENA_PROTEOME):
        site_keys = tuple(
            sorted(
                {
                    f"{row['protein_accession']}:C{row['cys_position']}"
                    for row in oof_rows
                    if row["arena"] == arena
                }
            )
        )
        run_scores = []
        for repetition in range(int(cfg["evaluation"]["repetitions"])):
            by_site = {
                f"{row['protein_accession']}:C{row['cys_position']}": float(row["candidate_score"])
                for row in oof_rows
                if row["arena"] == arena
                and str(row["run"]).startswith(f"{arena}_r{repetition}f")
            }
            run_scores.append(tuple(by_site[key] for key in site_keys))
        rank_summaries[arena] = [
            {
                "site_key": item.site_key,
                "median_score": item.median_score,
                "median_rank": item.median_rank,
                "best_rank": item.best_rank,
                "worst_rank": item.worst_rank,
            }
            for item in summarize_rank_runs(site_keys, tuple(run_scores))
        ]
    selected_model = select_release_model(decision.admitted)
    release_models: list[dict[str, object]] = []
    candidate_score_rows: list[dict[str, object]] = []
    full_rows = build_tomato_pu_rows(
        xlsx_path,
        proteome,
        arena=ARENA_PANEL,
        ratio=int(cfg["data"]["arena_ratio"]),
        seed=int(cfg["data"]["subsample_seed"]),
    )
    full_x = _matrix(full_rows, proteome)
    full_y = [row.label for row in full_rows]
    full_proteins = [row.protein_accession for row in full_rows]
    full_priors = sensitivity_priors(
        full_y.count("positive"),
        len(full_y),
        tuple(float(x) for x in cfg["pu"]["prior_multipliers"]),
        float(cfg["pu"]["prior_cap"]),
    )
    registry = build_tomato_candidate_registry(
        xlsx_path,
        proteome_path,
        supplementary_registry,
        model_input_registry,
    )
    candidate_x = [
        list(
            build_site_biology_vector_from_sequence(
                row.protein_accession,
                row.cys_position,
                proteome[row.protein_accession],
            ).values
        )
        for row in registry.candidates
    ]
    percentile_columns: list[tuple[float, ...]] = []
    if selected_model == "additive_pu":
        for seed in cfg["evaluation"]["model_seeds"]:
            for prior in full_priors:
                model = fit_additive_pu_ranker(
                    full_x,
                    full_y,
                    full_proteins,
                    feature_names,
                    AdditivePuConfig(class_prior=prior, seed=int(seed)),
                )
                release_models.append(
                    {"seed": int(seed), "class_prior": prior, "model": model.to_dict()}
                )
                percentile_columns.append(_percentile_column(model.score(candidate_x)))
    else:
        for seed in cfg["evaluation"]["model_seeds"]:
            baseline_scores = pu_logistic_regression_scores(
                full_x, full_y, candidate_x, int(seed)
            )
            release_models.append(
                {
                    "seed": int(seed),
                    "model_type": "pu_logistic",
                    "refit_recipe": "traditional.pu_logistic_regression_scores",
                }
            )
            percentile_columns.append(_percentile_column(tuple(baseline_scores)))
    release_envelope = ApplicabilityEnvelope.fit(full_x)
    for index, candidate in enumerate(registry.candidates):
        per_model = [column[index] for column in percentile_columns]
        candidate_score_rows.append(
            {
                "site_key": candidate.site_key,
                "protein_accession": candidate.protein_accession,
                "cys_position": candidate.cys_position,
                "median_percentile": statistics.median(per_model),
                "conservative_percentile": _lower_quantile(per_model, 0.2),
                "in_domain": release_envelope.contains(candidate_x[index]),
            }
        )
    model_release_created = bool(percentile_columns and candidate_score_rows)
    if not model_release_created:
        raise RuntimeError("selected model produced no frozen candidate scores")
    output_dir.mkdir(parents=True, exist_ok=False)
    summary: dict[str, object] = {
        "experiment": cfg["experiment"]["name"],
        "claim_class": cfg["claim_class"],
        "candidate_model_admitted": decision.admitted,
        "admission_reason": decision.reason,
        "admission_paired_intervals": decision.paired_intervals,
        "candidate_runs": candidate_runs,
        "baseline_runs": baseline_runs,
        "structure_delta": 0,
        "selected_model": selected_model,
        "model_release_created": model_release_created,
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "oof_scores.json").write_text(
        json.dumps(oof_rows, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "rank_summaries.json").write_text(
        json.dumps(rank_summaries, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if model_release_created:
        (output_dir / "frozen_models.json").write_text(
            json.dumps(release_models, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (output_dir / "candidate_scores.json").write_text(
            json.dumps(candidate_score_rows, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    manifest = {
        "complete": True,
        "config_sha256": hash_file(config_path, "sha256"),
        "kiae271_sha256": hash_file(xlsx_path, "sha256"),
        "proteome_sha256": hash_file(proteome_path, "sha256"),
        "clusters_sha256": hash_file(cluster_path, "sha256"),
        "feature_names": feature_names,
        "candidate_model_admitted": decision.admitted,
        "selected_model": selected_model,
        "model_release_created": model_release_created,
        "frozen_models_sha256": hash_file(
            output_dir / "frozen_models.json", "sha256"
        ),
        "candidate_scores_sha256": hash_file(
            output_dir / "candidate_scores.json", "sha256"
        ),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary
```

- [ ] **Step 6: Add the thin CLI**

```python
from pathlib import Path

from plantpersulf.workflows.tomato_ranker_v2 import run_tomato_ranker_v2


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run_tomato_ranker_v2(args.config, args.output_dir)


if __name__ == "__main__":
    main()
```

- [ ] **Step 7: Add the Gate 2 firewall test**

```python
from plantpersulf.evaluation.external_validation import _parse_fold_study


def test_tomato_v2_fold_names_are_not_gate2_evidence() -> None:
    assert _parse_fold_study("tomato_v2_r0f0|additive_pu") is None
```

- [ ] **Step 8: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_tomato_ranker_v2_config.py tests/scientific/test_tomato_ranker_v2_split_and_prior_policy.py tests/release/test_tomato_ranker_v2_cannot_change_gate2.py -v`

Expected: PASS; no scientific model run is required by these policy tests.

- [ ] **Step 9: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 10: Commit Task 8 implementation only**

```powershell
git add src/plantpersulf/workflows scripts/run_tomato_ranker_v2.py configs/experiments/tomato_ranker_v2.yaml tests/unit/test_tomato_ranker_v2_config.py tests/scientific/test_tomato_ranker_v2_split_and_prior_policy.py tests/release/test_tomato_ranker_v2_cannot_change_gate2.py
git commit -m "feat: add gated tomato PU ranker workflow"
```

Stop. Running the real configuration is a distinct reviewer approval and must produce its own registered result manifest and TDD report.

---

### Task 9: Cross-Crop Target-Label-Free Arm

**Reviewer gate:** Task 8 implementation and its real run must be reviewed. This task cannot read any tomato positive-label table, including KIAE271.

**Approved graded-evidence amendment (2026-08-11; supersedes any conflicting
Task 9 threshold, split, holdout, or claim text below):**

- Task 9 is an auxiliary innovation arm. The final generalization evidence is
  the preregistered, blinded tomato candidate validation, not internal source
  cross-validation.
- Use three labeled source species—Arabidopsis thaliana, Oryza sativa, and
  Magnaporthe oryzae—and unlabeled Solanum lycopersicum as the fourth, target
  species. Magnaporthe is a source-domain distance stress test and must not be
  described as a crop.
- Require at least one registered study per source species. Two or more studies
  per species is the ideal evidence tier, not an absolute admission barrier.
  Emit `within_species_replication_supported=false` globally whenever any
  source species is below the ideal tier, plus a per-species status map.
- Within every source study, run repeated 10-fold positive-unlabeled
  cross-validation grouped by homology cluster when registered cluster IDs are
  available, otherwise by protein. All sites from a protein or cluster stay in
  one fold. Five repetitions are frozen for v1. Random site-level splitting is
  forbidden.
- Missing-value imputation, standardization, class-prior sensitivity, and any
  model selection must be fitted or derived from training rows only. Frozen
  hyperparameters may be used instead of inner-fold tuning; the outer held-out
  fold must never select them.
- When registered batch, condition, or independent-experiment IDs exist, add
  complete leave-one-batch/condition-out evaluations.
- Cross-source evaluation uses leave-one-source-domain-out, where one domain is
  exactly one `species-study` batch. It supports a cross-source domain-transfer
  statement only. With singleton-study species it cannot identify a pure
  species effect and cannot support a universal cross-crop generalization
  claim.
- The manifest must preserve fold-level within-source metrics, all source-domain
  metrics, study counts, preprocessing scope, and the three negative/conditional
  claim flags. Tomato labels remain structurally unavailable during fit and
  tuning.
- PLM-v3 remains a later, independent method-development enhancement and is not
  required for this Task 9 evidence tier.

**Files:**
- Create: `src/plantpersulf/workflows/cross_crop_target_label_free.py`
- Create: `scripts/run_cross_crop_target_label_free_v1.py`
- Create: `configs/experiments/cross_crop_target_label_free_v1.yaml`
- Create: `tests/unit/test_target_label_free_dependency_graph.py`
- Create: `tests/scientific/test_cross_crop_training_excludes_tomato_labels.py`
- Create: `tests/release/test_cross_crop_result_is_not_gate2.py`

**Interfaces:**
- Produces: `SourceBatch`, `TargetCandidateBatch`, `validate_target_label_free_sources(sources, target_species)`, `run_cross_crop_target_label_free(config_path, source_batches, target_candidates, output_dir)`, `scores.json`, and a last-written `manifest.json`
- Source species: registered Arabidopsis, rice, and Magnaporthe batches; the
  latter is a distance stress-test domain rather than a crop-transfer replicate
- Target input: every tomato reference-proteome Cys and its six core features, materialized from the registered proteome only, without using Task 7 or any KIAE271-derived exclusion during scoring. The later frozen release may apply the pre-registered “novel site” eligibility filter only after X scores are immutable.

- [ ] **Step 1: Write the dependency-firewall RED test**

```python
import pytest

from plantpersulf.workflows.cross_crop_target_label_free import (
    SourceBatch,
    validate_target_label_free_sources,
)
from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES


def test_tomato_source_batch_is_rejected() -> None:
    source = SourceBatch(
        species="Solanum lycopersicum",
        study_accession="policy-marker-study",
        feature_names=BIOLOGY_FEATURE_NAMES,
        features=((0.0,) * 6, (1.0,) * 6),
        labels=("positive", "unlabeled"),
        protein_ids=("group-a", "group-a"),
        source_sha256="marker-sha",
    )
    with pytest.raises(RuntimeError, match="target label leakage"):
        validate_target_label_free_sources((source,), target_species="Solanum lycopersicum")
```

- [ ] **Step 2: Run tests and verify RED**

Run: `pytest tests/unit/test_target_label_free_dependency_graph.py -v`

Expected: FAIL during import because the workflow module does not exist.

- [ ] **Step 3: Implement the source-batch firewall**

```python
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SourceBatch:
    species: str
    study_accession: str
    feature_names: tuple[str, ...]
    features: tuple[tuple[float, ...], ...]
    labels: tuple[str, ...]
    protein_ids: tuple[str, ...]
    source_sha256: str


@dataclass(frozen=True)
class TargetCandidateBatch:
    site_keys: tuple[str, ...]
    feature_names: tuple[str, ...]
    features: tuple[tuple[float, ...], ...]
    source_sha256: str


def validate_target_label_free_sources(
    sources: tuple[SourceBatch, ...], target_species: str
) -> None:
    from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES

    if not sources:
        raise RuntimeError("at least one registered source batch is required")
    source_keys = [(source.species, source.study_accession) for source in sources]
    if len(source_keys) != len(set(source_keys)):
        raise RuntimeError("duplicate source species-study batch")
    for source in sources:
        if source.species == target_species:
            raise RuntimeError("target label leakage: target species in source batches")
        if not source.source_sha256:
            raise RuntimeError("source batch lacks registered SHA256")
        if not source.study_accession:
            raise RuntimeError("source batch lacks study accession")
        if source.feature_names != BIOLOGY_FEATURE_NAMES:
            raise RuntimeError("source feature order differs from biological contract")
        if set(source.labels) - {"positive", "unlabeled"}:
            raise RuntimeError("source batch contains a forbidden label")
        if not (len(source.features) == len(source.labels) == len(source.protein_ids)):
            raise RuntimeError("source batch row count mismatch")
        if not source.features or any(
            len(row) != len(source.feature_names) for row in source.features
        ):
            raise RuntimeError("source feature width mismatch")


@dataclass(frozen=True)
class SourceTransferDecision:
    admitted: bool
    reason: str


def source_transfer_gate(sources: tuple[SourceBatch, ...]) -> SourceTransferDecision:
    import statistics

    from plantpersulf.evaluation.metrics import average_precision
    from plantpersulf.models.additive_pu_ranker import (
        AdditivePuConfig,
        fit_additive_pu_ranker,
    )
    from plantpersulf.models.traditional import pu_logistic_regression_scores
    from plantpersulf.workflows.tomato_ranker_v2 import sensitivity_priors

    if len({source.species for source in sources}) < 2:
        return SourceTransferDecision(False, "fewer_than_two_source_species")
    for holdout_name, training, held_out in leave_one_source_domain_out(sources):
        train_x = [list(row) for batch in training for row in batch.features]
        train_y = [label for batch in training for label in batch.labels]
        train_proteins = [
            f"{batch.species}|{batch.study_accession}|{protein}"
            for batch in training
            for protein in batch.protein_ids
        ]
        test_x = [list(row) for batch in held_out for row in batch.features]
        test_y = [label for batch in held_out for label in batch.labels]
        priors = sensitivity_priors(
            train_y.count("positive"), len(train_y), (1.0, 1.5, 2.0), 0.5
        )
        feature_names = training[0].feature_names
        candidate_ap = []
        baseline_ap = []
        for seed in (0, 1, 2, 3, 4):
            baseline_scores = pu_logistic_regression_scores(
                train_x, train_y, test_x, seed
            )
            baseline_value = average_precision(
                list(zip(baseline_scores, test_y, strict=True))
            )
            if baseline_value is None:
                raise RuntimeError("held-out source batch has no positive")
            baseline_ap.append(baseline_value)
            for prior in priors:
                model = fit_additive_pu_ranker(
                    train_x,
                    train_y,
                    train_proteins,
                    feature_names,
                    AdditivePuConfig(class_prior=prior, seed=seed),
                )
                value = average_precision(
                    list(zip(model.score(test_x), test_y, strict=True))
                )
                if value is None:
                    raise RuntimeError("held-out source batch has no positive")
                candidate_ap.append(value)
        if min(candidate_ap) <= max(baseline_ap):
            return SourceTransferDecision(
                False,
                f"unstable_source_transfer:{holdout_name}:"
                f"candidate_median={statistics.median(candidate_ap):.6f}",
            )
    return SourceTransferDecision(True, "all_source_holdouts_stably_better")
```

- [ ] **Step 4: Create a source-only config**

```yaml
experiment:
  name: cross_crop_target_label_free_v1
claim_class: target_label_free_transfer_not_gate2
target_species: Solanum lycopersicum
source_species:
  - Arabidopsis thaliana
  - Oryza sativa
  - Magnaporthe oryzae
minimum_source_species: 3
minimum_total_source_studies: 4
minimum_source_studies_per_species: 1
ideal_source_studies_per_species: 2
forbidden_target_label_sources:
  - KIAE271_SUPPL
  - tomato_local_v1
  - tomato_ranker_v2
within_source_validation:
  split: repeated_grouped_protein_or_homology_cluster_cv
  folds: 10
  repetitions: 5
  grouping_preference: homology_cluster_then_protein
  preprocessing_fit_scope: training_fold_only
  missing_value_policy: training_fold_median
  leave_batch_condition_out_when_available: true
evaluation:
  source_selection: leave_one_source_domain_out
  source_domain_unit: species_study
  target_labels_visible_during_fit: false
  target_labels_visible_during_tuning: false
claims:
  pure_species_effect_supported: false
  universal_cross_crop_generalization_supported: false
  final_generalization_requires_tomato_blind_validation: true
applicability:
  min_target_in_domain_fraction: 0.80
wetlab:
  enabled: false
  max_noncontrol_fraction: 0.25
```

- [ ] **Step 5: Implement source-only fit and target scoring over registered materialized batches**

Scientific source-table materialization is a separate reviewer-approved data task because it writes new registered artifacts. This workflow accepts only already materialized batches with verified SHA256 values, fits `AdditivePuModel`, and scores a target batch with no labels. Its signature has no target-label parameter:

```python
def run_cross_crop_target_label_free(
    config_path: Path,
    source_batches: tuple[SourceBatch, ...],
    target_candidates: TargetCandidateBatch,
    output_dir: Path,
) -> dict[str, object]:
    import json
    import statistics
    from pathlib import Path

    import yaml

    from plantpersulf.evaluation.rank_ensemble import ApplicabilityEnvelope
    from plantpersulf.features.site_biology import BIOLOGY_FEATURE_NAMES
    from plantpersulf.models.additive_pu_ranker import (
        AdditivePuConfig,
        fit_additive_pu_ranker,
    )
    from plantpersulf.provenance.hashing import hash_file
    from plantpersulf.workflows.tomato_ranker_v2 import (
        _lower_quantile,
        _percentile_column,
        sensitivity_priors,
    )

    cfg = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    target_species = str(cfg["target_species"])
    if cfg["claim_class"] != "target_label_free_transfer_not_gate2":
        raise RuntimeError("cross-crop workflow cannot emit a Gate 2 claim")
    if (
        cfg["evaluation"]["target_labels_visible_during_fit"] is not False
        or cfg["evaluation"]["target_labels_visible_during_tuning"] is not False
        or cfg["wetlab"]["enabled"] is not False
        or float(cfg["wetlab"]["max_noncontrol_fraction"]) > 0.25
    ):
        raise RuntimeError("target-label or wet-lab policy differs from frozen config")
    validate_target_label_free_sources(source_batches, target_species)
    expected_sources = set(str(x) for x in cfg["source_species"])
    actual_sources = {batch.species for batch in source_batches}
    if actual_sources != expected_sources:
        raise RuntimeError("source species differ from frozen config")
    minimum_studies = int(cfg["minimum_source_studies_per_species"])
    for species in expected_sources:
        study_count = len(
            {
                batch.study_accession
                for batch in source_batches
                if batch.species == species
            }
        )
        if study_count < minimum_studies:
            raise RuntimeError("source species lacks independent study support")
    if output_dir.exists():
        raise FileExistsError(output_dir)
    if not target_candidates.source_sha256:
        raise RuntimeError("target candidate batch lacks registered SHA256")
    if len(target_candidates.site_keys) != len(target_candidates.features):
        raise RuntimeError("target candidate row count mismatch")
    if target_candidates.feature_names != BIOLOGY_FEATURE_NAMES:
        raise RuntimeError("target feature order differs from biological contract")
    if len(set(target_candidates.site_keys)) != len(target_candidates.site_keys):
        raise RuntimeError("duplicate target site key")
    if not target_candidates.features or any(
        len(row) != len(target_candidates.feature_names)
        for row in target_candidates.features
    ):
        raise RuntimeError("target feature width mismatch")
    source_decision = source_transfer_gate(source_batches)
    train_x = [list(row) for batch in source_batches for row in batch.features]
    train_y = [label for batch in source_batches for label in batch.labels]
    train_proteins = [
        f"{batch.species}|{batch.study_accession}|{protein}"
        for batch in source_batches
        for protein in batch.protein_ids
    ]
    priors = sensitivity_priors(
        train_y.count("positive"), len(train_y), (1.0, 1.5, 2.0), 0.5
    )
    score_columns: list[tuple[float, ...]] = []
    target_x = [list(row) for row in target_candidates.features]
    feature_names = source_batches[0].feature_names
    for seed in (0, 1, 2, 3, 4):
        for prior in priors:
            model = fit_additive_pu_ranker(
                train_x,
                train_y,
                train_proteins,
                feature_names,
                AdditivePuConfig(class_prior=prior, seed=seed),
            )
            score_columns.append(_percentile_column(model.score(target_x)))
    median_scores = [
        statistics.median(column[i] for column in score_columns)
        for i in range(len(target_x))
    ]
    conservative_scores = [
        _lower_quantile([column[i] for column in score_columns], 0.2)
        for i in range(len(target_x))
    ]
    envelope = ApplicabilityEnvelope.fit(train_x)
    in_domain = [envelope.contains(row) for row in target_x]
    target_in_domain_fraction = sum(in_domain) / len(in_domain)
    applicability_passed = target_in_domain_fraction >= float(
        cfg["applicability"]["min_target_in_domain_fraction"]
    )
    output_dir.mkdir(parents=True, exist_ok=False)
    score_payload: dict[str, object] = {
        "scores": {
            site_key: {
                "median_percentile": median_scores[i],
                "conservative_percentile": conservative_scores[i],
                "in_domain": in_domain[i],
            }
            for i, site_key in enumerate(target_candidates.site_keys)
        },
    }
    (output_dir / "scores.json").write_text(
        json.dumps(score_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    result: dict[str, object] = {
        "complete": True,
        "claim_class": cfg["claim_class"],
        "target_species": target_species,
        "sources": [
            {
                "species": batch.species,
                "study_accession": batch.study_accession,
                "sha256": batch.source_sha256,
            }
            for batch in source_batches
        ],
        "target_candidate_sha256": target_candidates.source_sha256,
        "target_label_hash": None,
        "source_gate_admitted": source_decision.admitted,
        "source_gate_reason": source_decision.reason,
        "target_in_domain_fraction": target_in_domain_fraction,
        "applicability_passed": applicability_passed,
        "wetlab_eligible": source_decision.admitted and applicability_passed,
        "scores_sha256": hash_file(output_dir / "scores.json", "sha256"),
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result
```

The registered materialized source batches must be generated later from the existing public Arabidopsis and rice parsers without tomato labels. The target batch must be enumerated from the registered tomato proteome alone and must not consume KIAE271 or the Task 7 registry. Their creation, registry insertion, and SHA256 audit are one distinct scientific-data task and cannot be folded silently into this model task.

- [ ] **Step 6: Add a thin CLI that accepts one hash-verified materialized request**

The request JSON contains `sources` and `target_candidates` fields matching the dataclasses above. It may contain target features and site keys but has no target-label field. A real run must provide a registry containing the request path and its SHA256; the CLI rejects an unregistered request before parsing it.

```python
from __future__ import annotations

import argparse
import json
from pathlib import Path

from plantpersulf.provenance.audit import assert_registered_input
from plantpersulf.provenance.hashing import hash_file
from plantpersulf.workflows.cross_crop_target_label_free import (
    SourceBatch,
    TargetCandidateBatch,
    run_cross_crop_target_label_free,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--request-registry", type=Path, required=True)
    parser.add_argument("--request-sha256", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    assert_registered_input(args.request, args.request_registry)
    if hash_file(args.request, "sha256") != args.request_sha256:
        raise RuntimeError("materialized request SHA256 mismatch")
    raw = json.loads(args.request.read_text(encoding="utf-8"))
    if set(raw) != {"sources", "target_candidates"}:
        raise RuntimeError("materialized request schema mismatch")
    expected_source_fields = {
        "species",
        "study_accession",
        "feature_names",
        "features",
        "labels",
        "protein_ids",
        "source_sha256",
    }
    if any(set(item) != expected_source_fields for item in raw["sources"]):
        raise RuntimeError("source batch schema mismatch")
    sources = tuple(
        SourceBatch(
            species=str(item["species"]),
            study_accession=str(item["study_accession"]),
            feature_names=tuple(str(x) for x in item["feature_names"]),
            features=tuple(tuple(float(x) for x in row) for row in item["features"]),
            labels=tuple(str(x) for x in item["labels"]),
            protein_ids=tuple(str(x) for x in item["protein_ids"]),
            source_sha256=str(item["source_sha256"]),
        )
        for item in raw["sources"]
    )
    target = raw["target_candidates"]
    if set(target) != {
        "site_keys",
        "feature_names",
        "features",
        "source_sha256",
    }:
        raise RuntimeError("target label leakage or schema mismatch")
    target_batch = TargetCandidateBatch(
        site_keys=tuple(str(x) for x in target["site_keys"]),
        feature_names=tuple(str(x) for x in target["feature_names"]),
        features=tuple(tuple(float(x) for x in row) for row in target["features"]),
        source_sha256=str(target["source_sha256"]),
    )
    run_cross_crop_target_label_free(
        args.config, sources, target_batch, args.output_dir
    )


if __name__ == "__main__":
    main()
```

- [ ] **Step 7: Add the scientific no-target-label test**

```python
import inspect

from plantpersulf.workflows.cross_crop_target_label_free import (
    run_cross_crop_target_label_free,
)


def test_cross_crop_runner_has_no_target_label_parameter() -> None:
    names = set(inspect.signature(run_cross_crop_target_label_free).parameters)
    assert not any("target_label" in name for name in names)
```

- [ ] **Step 8: Add the Gate 2 firewall test**

```python
from plantpersulf.evaluation.external_validation import _parse_fold_study


def test_cross_crop_fold_names_are_not_gate2_evidence() -> None:
    assert _parse_fold_study("cross_crop_source_holdout_rice|additive_pu") is None
```

- [ ] **Step 9: Run task tests and verify GREEN**

Run: `pytest tests/unit/test_target_label_free_dependency_graph.py tests/scientific/test_cross_crop_training_excludes_tomato_labels.py tests/release/test_cross_crop_result_is_not_gate2.py -v`

Expected: PASS; the fit/tune API has no target-label channel.

- [ ] **Step 10: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 11: Commit Task 9 only**

```powershell
git add src/plantpersulf/workflows/cross_crop_target_label_free.py scripts/run_cross_crop_target_label_free_v1.py configs/experiments/cross_crop_target_label_free_v1.yaml tests/unit/test_target_label_free_dependency_graph.py tests/scientific/test_cross_crop_training_excludes_tomato_labels.py tests/release/test_cross_crop_result_is_not_gate2.py
git commit -m "feat: add target-label-free cross-crop arm"
```

Stop. A real run and any decision to allocate blind candidates require separate reviewer approval.

---

### Task 10: Hash-Frozen Blind Release Policy Library

**Reviewer gate:** Tasks 8 and 9 must be accepted. This task implements and tests the immutable release policy with non-biological marker files only. Building a real candidate release remains outside this plan until K, power, matching variables, collaborator protocol, and an admitted model are frozen.

**Files:**
- Create: `src/plantpersulf/reporting/__init__.py`
- Create: `src/plantpersulf/reporting/candidate_release.py`
- Create: `scripts/build_candidate_release_v1.py`
- Create: `docs/wetlab_validation_matrix.md`
- Create: `tests/release/test_candidate_release_traceability.py`
- Create: `tests/release/test_candidate_release_cannot_be_reranked.py`
- Create: `tests/release/test_cross_crop_quota_is_bounded.py`

**Interfaces:**
- Consumes: policy-selected frozen T score release (A if admitted, otherwise B), optional admitted X transfer manifest, pre-registered K/power/SAP document, and a preselected frozen candidate table
- Produces: `ReleaseInputs`, `CandidateRelease`, `release_id_from_inputs`, `blind_id`, and `build_candidate_release` for a preselected, frozen candidate table

- [ ] **Step 1: Write the traceability RED test**

```python
from pathlib import Path

import pytest

from plantpersulf.reporting.candidate_release import ReleaseInputs, build_candidate_release


def test_release_requires_frozen_analysis_plan(tmp_path: Path) -> None:
    inputs = ReleaseInputs(
        tomato_model_manifest=tmp_path / "model.json",
        cross_crop_model_manifest=None,
        candidate_table=tmp_path / "scores.tsv",
        analysis_plan=None,
        total_noncontrol_k=20,
        cross_crop_k=0,
        matched_unlabeled_k=20,
        process_control_k=2,
    )
    with pytest.raises(RuntimeError, match="analysis plan"):
        build_candidate_release(inputs, tmp_path / "release")
```

- [ ] **Step 2: Write the quota and no-overwrite RED tests**

```python
import json
from pathlib import Path

import pytest

from plantpersulf.reporting.candidate_release import ReleaseInputs, build_candidate_release


def test_cross_crop_quota_cannot_exceed_one_quarter() -> None:
    with pytest.raises(RuntimeError, match="25%"):
        ReleaseInputs(
            tomato_model_manifest=Path("model-marker"),
            cross_crop_model_manifest=Path("cross-crop-marker"),
            candidate_table=Path("score-marker"),
            analysis_plan=Path("plan-marker"),
            total_noncontrol_k=20,
            cross_crop_k=6,
            matched_unlabeled_k=20,
            process_control_k=2,
        ).validate()


def test_existing_release_directory_is_never_overwritten(tmp_path: Path) -> None:
    output = tmp_path / "release"
    output.mkdir()
    with pytest.raises(FileExistsError):
        build_candidate_release(ReleaseInputs(
            tomato_model_manifest=tmp_path / "model-marker",
            cross_crop_model_manifest=tmp_path / "cross-crop-marker",
            candidate_table=tmp_path / "score-marker",
            analysis_plan=tmp_path / "plan-marker",
            total_noncontrol_k=4,
            cross_crop_k=1,
            matched_unlabeled_k=4,
            process_control_k=1,
        ), output)


def test_cross_crop_candidates_require_an_admitted_transfer_manifest(
    tmp_path: Path,
) -> None:
    tomato = tmp_path / "tomato.json"
    table = tmp_path / "candidates.tsv"
    plan = tmp_path / "plan.md"
    tomato.write_text(
        '{"complete": true, "selected_model": "pu_logistic", '
        '"model_release_created": true}\n',
        encoding="utf-8",
    )
    table.write_text("policy-marker\n", encoding="utf-8")
    plan.write_text(
        json.dumps(
            {
                "status": "frozen",
                "power_passed": True,
                "collaborator_protocol_frozen": True,
                "cross_crop_power_preserved": True,
                "total_noncontrol_k": 4,
                "cross_crop_k": 1,
                "matched_unlabeled_k": 4,
                "process_control_k": 1,
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(RuntimeError, match="transfer manifest"):
        ReleaseInputs(
            tomato_model_manifest=tomato,
            cross_crop_model_manifest=None,
            candidate_table=table,
            analysis_plan=plan,
            total_noncontrol_k=4,
            cross_crop_k=1,
            matched_unlabeled_k=4,
            process_control_k=1,
        ).validate()
```

- [ ] **Step 3: Run tests and verify RED**

Run: `pytest tests/release/test_candidate_release_traceability.py tests/release/test_cross_crop_quota_is_bounded.py -v`

Expected: FAIL during import because the reporting module does not exist.

- [ ] **Step 4: Implement the immutable input contract**

```python
from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.provenance.hashing import hash_file


@dataclass(frozen=True)
class ReleaseInputs:
    tomato_model_manifest: Path
    cross_crop_model_manifest: Path | None
    candidate_table: Path
    analysis_plan: Path | None
    total_noncontrol_k: int
    cross_crop_k: int
    matched_unlabeled_k: int
    process_control_k: int

    def validate(self) -> None:
        if self.analysis_plan is None:
            raise RuntimeError("frozen analysis plan is required")
        if self.total_noncontrol_k <= 0:
            raise RuntimeError("total_noncontrol_k must be positive")
        if self.cross_crop_k < 0:
            raise RuntimeError("cross_crop_k cannot be negative")
        if self.matched_unlabeled_k < 0 or self.process_control_k < 0:
            raise RuntimeError("control arm sizes cannot be negative")
        if self.cross_crop_k / self.total_noncontrol_k > 0.25:
            raise RuntimeError("cross-crop arm exceeds the 25% cap")
        for path in (self.tomato_model_manifest, self.candidate_table, self.analysis_plan):
            if path is None or not path.exists():
                raise RuntimeError(f"release input missing: {path}")
        assert self.analysis_plan is not None
        analysis_plan = json.loads(self.analysis_plan.read_text(encoding="utf-8"))
        frozen_sizes = {
            "total_noncontrol_k": self.total_noncontrol_k,
            "cross_crop_k": self.cross_crop_k,
            "matched_unlabeled_k": self.matched_unlabeled_k,
            "process_control_k": self.process_control_k,
        }
        if (
            analysis_plan.get("status") != "frozen"
            or analysis_plan.get("power_passed") is not True
            or analysis_plan.get("collaborator_protocol_frozen") is not True
            or any(analysis_plan.get(key) != value for key, value in frozen_sizes.items())
        ):
            raise RuntimeError("release arguments differ from the frozen analysis plan")
        if (
            self.cross_crop_k > 0
            and analysis_plan.get("cross_crop_power_preserved") is not True
        ):
            raise RuntimeError("analysis plan does not preserve cross-crop power")
        tomato_manifest = json.loads(
            self.tomato_model_manifest.read_text(encoding="utf-8")
        )
        if (
            tomato_manifest.get("complete") is not True
            or tomato_manifest.get("selected_model")
            not in {"additive_pu", "pu_logistic"}
            or tomato_manifest.get("model_release_created") is not True
        ):
            raise RuntimeError("tomato model manifest has no selected release model")
        if self.cross_crop_k > 0:
            if self.cross_crop_model_manifest is None:
                raise RuntimeError("cross-crop candidates require a transfer manifest")
            if not self.cross_crop_model_manifest.exists():
                raise RuntimeError("cross-crop transfer manifest is missing")
            cross_crop_manifest = json.loads(
                self.cross_crop_model_manifest.read_text(encoding="utf-8")
            )
            if (
                cross_crop_manifest.get("complete") is not True
                or cross_crop_manifest.get("claim_class")
                != "target_label_free_transfer_not_gate2"
                or cross_crop_manifest.get("source_gate_admitted") is not True
                or cross_crop_manifest.get("applicability_passed") is not True
                or cross_crop_manifest.get("wetlab_eligible") is not True
                or "target_label_hash" not in cross_crop_manifest
                or cross_crop_manifest.get("target_label_hash") is not None
            ):
                raise RuntimeError("cross-crop transfer manifest is not admitted")
        elif self.cross_crop_model_manifest is not None:
            raise RuntimeError("cross-crop manifest supplied while cross_crop_k is zero")


@dataclass(frozen=True)
class CandidateRelease:
    release_id: str
    output_dir: Path


def release_id_from_inputs(inputs: ReleaseInputs) -> str:
    inputs.validate()
    cross_crop_manifest = inputs.cross_crop_model_manifest
    material = "|".join(
        (
            hash_file(inputs.tomato_model_manifest, "sha256"),
            hash_file(cross_crop_manifest, "sha256")
            if cross_crop_manifest is not None else "",
            hash_file(inputs.candidate_table, "sha256"),
            hash_file(inputs.analysis_plan, "sha256") if inputs.analysis_plan else "",
            str(inputs.total_noncontrol_k),
            str(inputs.cross_crop_k),
            str(inputs.matched_unlabeled_k),
            str(inputs.process_control_k),
        )
    )
    return "candidate-release-" + hashlib.sha256(material.encode("utf-8")).hexdigest()[:16]
```

- [ ] **Step 5: Implement deterministic blind IDs and release hashing**

Use `hashlib.sha256` over the frozen release ID plus stable site key; never use Python's randomized `hash()`:

```python
def blind_id(release_id: str, site_key: str) -> str:
    digest = hashlib.sha256(f"{release_id}|{site_key}".encode("utf-8")).hexdigest()
    return f"BLIND-{digest[:12].upper()}"
```

The preselected candidate table has the exact columns `site_key`, `protein_accession`, `arm`, `selection_order`, and `matched_stratum`. Allowed arms are `T`, `X`, `matched_unlabeled`, and `process_control`. Implement the writer exactly as follows; candidate selection and matching are intentionally not performed here because their variables must come from the future frozen SAP:

```python
RELEASE_COLUMNS = (
    "site_key",
    "protein_accession",
    "arm",
    "selection_order",
    "matched_stratum",
)


def build_candidate_release(inputs: ReleaseInputs, output_dir: Path) -> CandidateRelease:
    if output_dir.exists():
        raise FileExistsError(output_dir)
    release_id = release_id_from_inputs(inputs)
    assert inputs.analysis_plan is not None
    with inputs.candidate_table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != RELEASE_COLUMNS:
            raise RuntimeError("candidate table schema mismatch")
        rows = [dict(row) for row in reader]
    allowed = {"T", "X", "matched_unlabeled", "process_control"}
    if {row["arm"] for row in rows} - allowed:
        raise RuntimeError("candidate table contains an unknown arm")
    t_rows = [row for row in rows if row["arm"] == "T"]
    x_rows = [row for row in rows if row["arm"] == "X"]
    matched_rows = [row for row in rows if row["arm"] == "matched_unlabeled"]
    process_control_rows = [row for row in rows if row["arm"] == "process_control"]
    if len(t_rows) + len(x_rows) != inputs.total_noncontrol_k:
        raise RuntimeError("non-control candidate count differs from frozen K")
    if len(x_rows) != inputs.cross_crop_k:
        raise RuntimeError("cross-crop candidate count differs from frozen K")
    if len(matched_rows) != inputs.matched_unlabeled_k:
        raise RuntimeError("matched-unlabeled count differs from frozen K")
    if len(process_control_rows) != inputs.process_control_k:
        raise RuntimeError("process-control count differs from frozen K")
    proteins = [row["protein_accession"] for row in (*t_rows, *x_rows)]
    if len(proteins) != len(set(proteins)):
        raise RuntimeError("more than one non-control site selected per protein")
    site_keys = [row["site_key"] for row in rows]
    if len(site_keys) != len(set(site_keys)):
        raise RuntimeError("candidate table contains a duplicate site key")
    output_dir.mkdir(parents=True, exist_ok=False)
    private_rows = [
        {**row, "blind_id": blind_id(release_id, row["site_key"])}
        for row in sorted(rows, key=lambda row: (int(row["selection_order"]), row["site_key"]))
    ]
    with (output_dir / "private_key.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("blind_id", *RELEASE_COLUMNS), delimiter="\t")
        writer.writeheader()
        writer.writerows(private_rows)
    with (output_dir / "blind_table.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "blind_id",
                "site_key",
                "protein_accession",
                "matched_stratum",
            ),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(
            {
                "blind_id": row["blind_id"],
                "site_key": row["site_key"],
                "protein_accession": row["protein_accession"],
                "matched_stratum": row["matched_stratum"],
            }
            for row in private_rows
        )
    cross_crop_manifest = inputs.cross_crop_model_manifest
    manifest = {
        "complete": True,
        "release_id": release_id,
        "tomato_model_manifest_sha256": hash_file(
            inputs.tomato_model_manifest, "sha256"
        ),
        "cross_crop_model_manifest_sha256": (
            hash_file(cross_crop_manifest, "sha256")
            if cross_crop_manifest is not None else None
        ),
        "candidate_table_sha256": hash_file(inputs.candidate_table, "sha256"),
        "analysis_plan_sha256": hash_file(inputs.analysis_plan, "sha256"),
        "total_noncontrol_k": inputs.total_noncontrol_k,
        "cross_crop_k": inputs.cross_crop_k,
        "matched_unlabeled_k": inputs.matched_unlabeled_k,
        "process_control_k": inputs.process_control_k,
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    hash_lines = [
        f"{hash_file(output_dir / name, 'sha256')}  {name}"
        for name in ("private_key.tsv", "blind_table.tsv", "manifest.json")
    ]
    (output_dir / "SHA256SUMS").write_text("\n".join(hash_lines) + "\n", encoding="utf-8")
    return CandidateRelease(release_id=release_id, output_dir=output_dir)
```

- [ ] **Step 6: Add the content-addressed zero-rerank test**

```python
import json
from pathlib import Path

from plantpersulf.reporting.candidate_release import (
    ReleaseInputs,
    release_id_from_inputs,
)


def test_changed_candidate_order_gets_a_new_release_id(tmp_path: Path) -> None:
    model = tmp_path / "model.json"
    plan = tmp_path / "plan.md"
    table = tmp_path / "candidates.tsv"
    model.write_text(
        '{"complete": true, "selected_model": "pu_logistic", '
        '"model_release_created": true}\n',
        encoding="utf-8",
    )
    plan.write_text(
        json.dumps(
            {
                "status": "frozen",
                "power_passed": True,
                "collaborator_protocol_frozen": True,
                "cross_crop_power_preserved": False,
                "total_noncontrol_k": 1,
                "cross_crop_k": 0,
                "matched_unlabeled_k": 0,
                "process_control_k": 0,
            }
        ),
        encoding="utf-8",
    )
    header = "site_key\tprotein_accession\tarm\tselection_order\tmatched_stratum\n"
    table.write_text(header + "site-a\tprotein-a\tT\t1\tstratum-a\n", encoding="utf-8")
    first = ReleaseInputs(
        tomato_model_manifest=model,
        cross_crop_model_manifest=None,
        candidate_table=table,
        analysis_plan=plan,
        total_noncontrol_k=1,
        cross_crop_k=0,
        matched_unlabeled_k=0,
        process_control_k=0,
    )
    first_id = release_id_from_inputs(first)
    table.write_text(header + "site-a\tprotein-a\tT\t2\tstratum-a\n", encoding="utf-8")
    second = ReleaseInputs(
        tomato_model_manifest=model,
        cross_crop_model_manifest=None,
        candidate_table=table,
        analysis_plan=plan,
        total_noncontrol_k=1,
        cross_crop_k=0,
        matched_unlabeled_k=0,
        process_control_k=0,
    )
    assert release_id_from_inputs(second) != first_id
```

Wet-lab results require a separate future `prospective_validation_v1` plan referencing `candidate_release_v1`; no API in this task may update an existing release.

- [ ] **Step 7: Write the complete blind protocol document**

Write explicit sections for visible blind fields, hidden-key owner, T/X/matched-unlabeled/process-control arms, hit definition, identical assay/exclusion rules, full-result return, batch/operator/raw-file capture, unblinding condition, primary and secondary enrichment tests, failure handling, and prohibition on post-result reranking. Require a signed human-readable SAP plus the machine-readable frozen manifest consumed by `ReleaseInputs`. Every section states that real values come only from that later SAP; this policy document supplies no biological values.

- [ ] **Step 8: Add the explicit-argument CLI with no overwrite mode**

```python
from __future__ import annotations

import argparse
from pathlib import Path

from plantpersulf.reporting.candidate_release import ReleaseInputs, build_candidate_release


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tomato-model-manifest", type=Path, required=True)
    parser.add_argument("--cross-crop-model-manifest", type=Path)
    parser.add_argument("--candidate-table", type=Path, required=True)
    parser.add_argument("--analysis-plan", type=Path, required=True)
    parser.add_argument("--total-noncontrol-k", type=int, required=True)
    parser.add_argument("--cross-crop-k", type=int, required=True)
    parser.add_argument("--matched-unlabeled-k", type=int, required=True)
    parser.add_argument("--process-control-k", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    inputs = ReleaseInputs(
        tomato_model_manifest=args.tomato_model_manifest,
        cross_crop_model_manifest=args.cross_crop_model_manifest,
        candidate_table=args.candidate_table,
        analysis_plan=args.analysis_plan,
        total_noncontrol_k=args.total_noncontrol_k,
        cross_crop_k=args.cross_crop_k,
        matched_unlabeled_k=args.matched_unlabeled_k,
        process_control_k=args.process_control_k,
    )
    build_candidate_release(inputs, args.output_dir)


if __name__ == "__main__":
    main()
```

- [ ] **Step 9: Run task tests and verify GREEN**

Run: `pytest tests/release/test_candidate_release_traceability.py tests/release/test_candidate_release_cannot_be_reranked.py tests/release/test_cross_crop_quota_is_bounded.py -v`

Expected: PASS; an existing or altered release cannot be overwritten or silently treated as the original.

- [ ] **Step 10: Run the global completion gate**

Run every command under “Completion Gate Used by Every Task”.

Expected: all commands pass.

- [ ] **Step 11: Commit the policy library without building or tagging a real release**

```powershell
git add src/plantpersulf/reporting scripts/build_candidate_release_v1.py docs/wetlab_validation_matrix.md tests/release/test_candidate_release_traceability.py tests/release/test_candidate_release_cannot_be_reranked.py tests/release/test_cross_crop_quota_is_bounded.py
git commit -m "feat: add immutable blind release policy"
```

Stop and submit the policy-task TDD evidence. A real release requires a new plan after the prerequisites in this task's reviewer gate exist.

---

## Plan Self-Review Coverage Map

| Approved design requirement | Implemented by |
|---|---|
| Shortcut-free biological feature core | Task 1 |
| Pairwise + nnPU objective; no hard negatives | Task 2 |
| Small deterministic elastic-net-like additive model | Task 3 |
| q_obs disabled without complete registered table | Task 4 |
| Missing structure neutral; structure coverage cannot admit residual | Task 5 |
| Repetition-blocked paired intervals, rank distribution, applicability | Task 6 |
| Unregistered tomato proteome/cluster fail closed | Pre-Task-7 provenance gate and Tasks 7–8 |
| Novel candidate registry excludes known and ambiguous KIAE271 sites | Task 7 |
| Repeated homology-cluster tomato development; A-to-B fallback; Gate 2 firewall | Task 8 |
| Source-only target-label-free cross-crop arm with leave-species/leave-study and applicability gates | Task 9 |
| X arm bounded by power/cap and isolated from T | Tasks 9–10 |
| Machine-readable frozen SAP, all arm sizes, immutable blind release, zero reranking | Task 10 |
| MSA, PLM, fitted SAR-PU | Explicitly not admitted without new data/reviewer task |
| PPLM/AF3/docking/MD/QM-MM | Explicitly post-hit and outside this software plan |

## Execution Boundary

This plan is documentation only. It does not authorize Task 1 automatically. The reviewer must approve Task 1 explicitly, then each later task in order. Real model runs, target-label-free runs, and the candidate release each require a distinct approval beyond approval of their software task.
