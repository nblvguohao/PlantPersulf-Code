"""Task 9 — minimal structure-aware PU ranker with gated multi-branch fusion.

Architecture (the smallest thing the Codex permits — no heterogeneous graphs,
no end-to-end 650M fine-tuning, every branch ablatable):

    sequence branch  ┐
    frozen-ESM branch├─ per-branch encoder ─ gated fusion (missingness-masked)
    structure branch ┘        │                        └→ PU ranking head
    study-context branch ─────┘                             (+ MC-dropout σ)

Two integrity properties are baked into the *structure* of the model rather
than left to convention:

1. **Missing structure is masked, never mean-imputed.** A row whose
   ``structure_mask`` is ``False`` has its structure-branch *output* forced to
   zero (``e_str = enc(struct) * mask``). The raw structure values therefore
   cannot influence the score or the fitted parameters — see
   ``tests/scientific/test_missing_structure_mask_is_respected.py``. The same
   gating applies to an unseen ``study_id`` (e.g. the held-out study in a
   leave-study-out fold), so study context never leaks a bias for a study the
   model was never trained on.

2. **PU-correct training.** Unlabeled cysteines are never treated as hard
   negatives. Training follows the same Elkan-Noto two-step used by the
   traditional baselines (``plantpersulf.models.pu_risk``): a first pass fits a
   non-traditional classifier on (labeled-positive vs unlabeled), a held-out
   slice of positives estimates the label frequency ``c = P(s=1|y=1)``, and a
   second pass reweights every example before refitting.

Reproducibility: everything runs full-batch on CPU under a single
``torch.manual_seed(seed)``; two runs with the same inputs and seed produce
bit-identical scores and uncertainties.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from plantpersulf.models.pu_risk import estimate_label_frequency, pu_example_weights
from plantpersulf.models.traditional import TrainOnlyScaler

if TYPE_CHECKING:
    import torch


# ---------------------------------------------------------------------------
# public data types
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class BranchFeatures:
    """Per-row features for every branch, plus the structure missingness mask.

    ``structure`` columns are ``[accessibility_proxy, plddt]`` (matching
    ``features.structure``). ``structure_mask[i]`` is ``True`` iff row ``i`` has
    a real registered structure; when ``False`` the structure values are
    ignored entirely. ``study_ids`` supplies the study-context branch; ``None``
    means a single constant study.
    """

    sequence: list[list[float]]
    esm: list[list[float]]
    structure: list[list[float]]
    structure_mask: list[bool]
    study_ids: list[str] | None = None

    def n_rows(self) -> int:
        return len(self.sequence)


@dataclass(frozen=True)
class AblationConfig:
    """Toggles for the Codex "must ablate" matrix. Disabling a branch/feature
    zeroes that branch's contribution in *both* training and prediction."""

    use_esm: bool = True
    use_structure: bool = True
    use_plddt: bool = True
    use_accessibility: bool = True
    use_study_context: bool = True


@dataclass(frozen=True)
class RankerOutput:
    scores: list[float] = field(default_factory=list)
    uncertainty: list[float] = field(default_factory=list)


# ---------------------------------------------------------------------------
# helpers (validation + tensor assembly)
# ---------------------------------------------------------------------------


def _validate(features: BranchFeatures) -> None:
    n = len(features.sequence)
    if len(features.esm) != n or len(features.structure) != n:
        raise ValueError("all branches must have equal length (n_rows)")
    if len(features.structure_mask) != n:
        raise ValueError("structure_mask must have equal length (n_rows)")
    if features.study_ids is not None and len(features.study_ids) != n:
        raise ValueError("study_ids must have equal length (n_rows)")


def _study_vocab(train: BranchFeatures) -> dict[str, int]:
    if train.study_ids is None:
        return {"__const__": 0}
    return {study: i for i, study in enumerate(sorted(set(train.study_ids)))}


def _study_onehot(
    features: BranchFeatures,
    vocab: dict[str, int],
) -> tuple[list[list[float]], list[bool]]:
    """One-hot the study id against the *training* vocab. Rows whose study is
    unseen (or absent) get an all-zero vector and a ``False`` known-mask, so an
    unseen held-out study contributes nothing (no leaked bias)."""
    dim = len(vocab)
    rows: list[list[float]] = []
    known: list[bool] = []
    ids = features.study_ids or ["__const__"] * features.n_rows()
    for study in ids:
        vec = [0.0] * dim
        idx = vocab.get(study)
        if idx is None:
            known.append(False)
        else:
            vec[idx] = 1.0
            known.append(True)
        rows.append(vec)
    return rows, known


@dataclass(frozen=True)
class _BranchScalers:
    """Train-only standardizers, one per numeric branch. The structure scaler
    is fit on *present-structure rows only* so masked-out (junk) structure
    values never leak into its statistics."""

    seq: TrainOnlyScaler
    esm: TrainOnlyScaler
    struct: TrainOnlyScaler

    @classmethod
    def fit(cls, train: BranchFeatures, ablation: AblationConfig) -> _BranchScalers:
        present = [
            train.structure[i]
            for i, m in enumerate(train.structure_mask)
            if m and ablation.use_structure
        ]
        struct_rows = present or [[0.0] * len(train.structure[0])]
        return cls(
            seq=TrainOnlyScaler.fit(train.sequence),
            esm=TrainOnlyScaler.fit(train.esm),
            struct=TrainOnlyScaler.fit(struct_rows),
        )

    def apply(self, features: BranchFeatures) -> BranchFeatures:
        return BranchFeatures(
            sequence=self.seq.transform(features.sequence),
            esm=self.esm.transform(features.esm),
            structure=self.struct.transform(features.structure),
            structure_mask=features.structure_mask,
            study_ids=features.study_ids,
        )


def _apply_structure_ablation(
    structure: list[list[float]],
    ablation: AblationConfig,
) -> list[list[float]]:
    """Zero the accessibility (col 0) or pLDDT (col 1) column per the ablation.
    Column semantics match ``features.structure`` -> [contact_proxy, plddt]."""
    if ablation.use_accessibility and ablation.use_plddt:
        return structure
    out: list[list[float]] = []
    for row in structure:
        new = list(row)
        if not ablation.use_accessibility and len(new) > 0:
            new[0] = 0.0
        if not ablation.use_plddt and len(new) > 1:
            new[1] = 0.0
        out.append(new)
    return out


# ---------------------------------------------------------------------------
# the network
# ---------------------------------------------------------------------------


def _build_network(
    d_seq: int, d_esm: int, d_str: int, d_study: int, hidden: int, dropout: float
) -> torch.nn.Module:
    import torch
    from torch import nn

    class _Ranker(nn.Module):
        def __init__(self) -> None:
            super().__init__()
            self.enc_seq = nn.Linear(d_seq, hidden)
            self.enc_esm = nn.Linear(d_esm, hidden)
            self.enc_str = nn.Linear(d_str, hidden)
            self.enc_study = nn.Linear(d_study, hidden)
            self.gate = nn.Linear(4 * hidden, 4)
            self.head = nn.Linear(hidden, 1)
            self.dropout = nn.Dropout(dropout)

        def forward(
            self,
            seq: torch.Tensor,
            esm: torch.Tensor,
            struct: torch.Tensor,
            struct_active: torch.Tensor,
            study: torch.Tensor,
            study_active: torch.Tensor,
        ) -> torch.Tensor:
            e_seq = torch.relu(self.enc_seq(seq))
            e_esm = torch.relu(self.enc_esm(esm))
            e_str = torch.relu(self.enc_str(struct)) * struct_active
            e_study = torch.relu(self.enc_study(study)) * study_active
            cat = torch.cat([e_seq, e_esm, e_str, e_study], dim=1)
            gates = torch.sigmoid(self.gate(cat))
            fused = (
                gates[:, 0:1] * e_seq
                + gates[:, 1:2] * e_esm
                + gates[:, 2:3] * e_str
                + gates[:, 3:4] * e_study
            )
            fused = self.dropout(fused)
            out: torch.Tensor = self.head(fused).squeeze(-1)
            return out

    return _Ranker()


def _to_tensors(
    features: BranchFeatures,
    vocab: dict[str, int],
    ablation: AblationConfig,
) -> dict[str, torch.Tensor]:
    import torch

    struct = _apply_structure_ablation(features.structure, ablation)
    study_onehot, study_known = _study_onehot(features, vocab)

    esm_active = 1.0 if ablation.use_esm else 0.0
    struct_flag = 1.0 if ablation.use_structure else 0.0
    study_flag = 1.0 if ablation.use_study_context else 0.0

    struct_active = [
        [struct_flag if m else 0.0] for m in features.structure_mask
    ]
    study_active = [[study_flag if k else 0.0] for k in study_known]

    esm_rows = features.esm
    if not ablation.use_esm:
        esm_rows = [[v * esm_active for v in row] for row in esm_rows]

    return {
        "seq": torch.tensor(features.sequence, dtype=torch.float32),
        "esm": torch.tensor(esm_rows, dtype=torch.float32),
        "struct": torch.tensor(struct, dtype=torch.float32),
        "struct_active": torch.tensor(struct_active, dtype=torch.float32),
        "study": torch.tensor(study_onehot, dtype=torch.float32),
        "study_active": torch.tensor(study_active, dtype=torch.float32),
    }


def _train_network(
    net: torch.nn.Module,
    tensors: dict[str, torch.Tensor],
    s_labels: torch.Tensor,
    weights: torch.Tensor,
    seed: int,
    epochs: int,
    lr: float,
) -> None:
    import torch

    torch.manual_seed(seed)
    net.train()
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    loss_fn = torch.nn.BCEWithLogitsLoss(reduction="none")
    for _ in range(epochs):
        opt.zero_grad()
        logits = net(
            tensors["seq"], tensors["esm"], tensors["struct"],
            tensors["struct_active"], tensors["study"], tensors["study_active"],
        )
        per = loss_fn(logits, s_labels)
        loss = (per * weights).sum() / weights.sum()
        loss.backward()
        opt.step()


def _forward_scores(
    net: torch.nn.Module, tensors: dict[str, torch.Tensor]
) -> torch.Tensor:
    import torch

    with torch.no_grad():
        logits = net(
            tensors["seq"], tensors["esm"], tensors["struct"],
            tensors["struct_active"], tensors["study"], tensors["study_active"],
        )
        return torch.sigmoid(logits).to("cpu")


def _select_ranker_device(torch_mod: object) -> torch.device:
    """Device for ranker training. ``PLANTPERSULF_DEVICE`` (e.g. ``cuda``)
    opts in to the GPU; otherwise CPU is used by default so a model release
    stays bit-for-bit reproducible across machines (GPU reductions are not
    guaranteed bit-identical to CPU)."""
    import os

    forced = os.environ.get("PLANTPERSULF_DEVICE")
    if forced:
        return torch_mod.device(forced)  # type: ignore[attr-defined,no-any-return]
    return torch_mod.device("cpu")  # type: ignore[attr-defined,no-any-return]


# ---------------------------------------------------------------------------
# public entry point
# ---------------------------------------------------------------------------


def structure_ranker_scores(
    train: BranchFeatures,
    train_y: list[str],
    predict: BranchFeatures,
    seed: int,
    ablation: AblationConfig | None = None,
    *,
    hidden: int = 16,
    dropout: float = 0.2,
    epochs: int = 200,
    lr: float = 0.05,
    holdout_fraction: float = 0.2,
    n_mc_dropout: int = 16,
) -> RankerOutput:
    """Fit the gated-fusion PU ranker and score ``predict``.

    Returns per-row point scores (dropout off) and MC-dropout uncertainties
    (std over ``n_mc_dropout`` stochastic passes). Deterministic for a fixed
    seed. Requires >=2 labeled positives (Elkan-Noto needs a held-out positive
    slice to estimate the label frequency).
    """
    import random

    import torch

    if ablation is None:
        ablation = AblationConfig()
    _validate(train)
    _validate(predict)
    if len(train_y) != train.n_rows():
        raise ValueError("train_y must have equal length (n_rows)")

    positive_idx = [i for i, y in enumerate(train_y) if y == "positive"]
    if len(positive_idx) < 2:
        raise ValueError("ranker requires at least 2 labeled positives")

    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)
    device = _select_ranker_device(torch)

    vocab = _study_vocab(train)
    scalers = _BranchScalers.fit(train, ablation)
    train_s = scalers.apply(train)
    predict_s = scalers.apply(predict)
    train_t = {
        k: v.to(device) for k, v in _to_tensors(train_s, vocab, ablation).items()
    }
    predict_t = {
        k: v.to(device) for k, v in _to_tensors(predict_s, vocab, ablation).items()
    }
    s_labels = torch.tensor(
        [1.0 if y == "positive" else 0.0 for y in train_y], dtype=torch.float32
    ).to(device)

    d_seq = len(train.sequence[0])
    d_esm = len(train.esm[0])
    d_str = len(train.structure[0])
    d_study = len(vocab)

    # --- Elkan-Noto held-out slice of positives (only used to estimate c) ---
    rng = random.Random(seed)
    shuffled_pos = positive_idx[:]
    rng.shuffle(shuffled_pos)
    n_holdout = max(1, int(len(shuffled_pos) * holdout_fraction))
    holdout = set(shuffled_pos[:n_holdout])
    fit_idx = [i for i in range(train.n_rows()) if i not in holdout]

    def _subset(t: dict[str, torch.Tensor], idx: list[int]) -> dict[str, torch.Tensor]:
        sel = torch.tensor(idx, dtype=torch.long, device=device)
        return {k: v.index_select(0, sel) for k, v in t.items()}

    fit_t = _subset(train_t, fit_idx)
    fit_s = s_labels.index_select(
        0, torch.tensor(fit_idx, dtype=torch.long, device=device)
    )
    unit_w = torch.ones_like(fit_s)

    # --- Pass A: non-traditional classifier (positive vs unlabeled) ---
    net_a = _build_network(d_seq, d_esm, d_str, d_study, hidden, dropout).to(device)
    _train_network(net_a, fit_t, fit_s, unit_w, seed, epochs, lr)

    holdout_idx = sorted(holdout)
    holdout_t = _subset(train_t, holdout_idx)
    net_a.eval()
    holdout_scores = [float(x) for x in _forward_scores(net_a, holdout_t).tolist()]
    c = estimate_label_frequency(holdout_scores)

    fit_scores = [float(x) for x in _forward_scores(net_a, fit_t).tolist()]
    is_labeled_positive = [train_y[i] == "positive" for i in fit_idx]
    weights = pu_example_weights(fit_scores, is_labeled_positive, c)
    weight_t = torch.tensor(weights, dtype=torch.float32).to(device)

    # --- Pass B: reweighted refit -> final model ---
    net_b = _build_network(d_seq, d_esm, d_str, d_study, hidden, dropout).to(device)
    _train_network(net_b, fit_t, fit_s, weight_t, seed, epochs, lr)

    net_b.eval()
    point_scores = [float(x) for x in _forward_scores(net_b, predict_t).tolist()]

    # --- MC-dropout uncertainty (dropout ON), re-seeded for reproducibility ---
    torch.manual_seed(seed + 1)
    net_b.train()
    samples: list[list[float]] = []
    with torch.no_grad():
        for _ in range(n_mc_dropout):
            logits = net_b(
                predict_t["seq"], predict_t["esm"], predict_t["struct"],
                predict_t["struct_active"], predict_t["study"],
                predict_t["study_active"],
            )
            samples.append([float(x) for x in torch.sigmoid(logits).tolist()])

    uncertainty = _column_std(samples, len(point_scores))
    return RankerOutput(scores=point_scores, uncertainty=uncertainty)


def _column_std(samples: list[list[float]], n_cols: int) -> list[float]:
    if not samples:
        return [0.0] * n_cols
    out: list[float] = []
    n = len(samples)
    for j in range(n_cols):
        col = [row[j] for row in samples]
        mean = sum(col) / n
        var = sum((x - mean) ** 2 for x in col) / n
        out.append(var ** 0.5)
    return out
