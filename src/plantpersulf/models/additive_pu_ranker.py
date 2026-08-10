"""Deterministic additive PU ranker over frozen numeric feature rows."""

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
    device: str = "cpu"


@dataclass(frozen=True)
class AdditivePuModel:
    feature_names: tuple[str, ...]
    scaler: TrainOnlyScaler
    weights: tuple[float, ...]
    intercept: float

    def score(
        self,
        features: list[list[float]],
        device: str = "cpu",
        batch_size: int = 16384,
    ) -> tuple[float, ...]:
        """Score rows with an optional deterministic accelerator batch path."""
        scaled = self.scaler.transform(features)
        resolved_device = resolve_device(device)
        if resolved_device.type != "cuda":
            return tuple(
                self.intercept
                + sum(
                    weight * value
                    for weight, value in zip(self.weights, row, strict=True)
                )
                for row in scaled
            )
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        weights = torch.tensor(
            self.weights, dtype=torch.float32, device=resolved_device
        )
        intercept = torch.tensor(
            self.intercept, dtype=torch.float32, device=resolved_device
        )
        values: list[float] = []
        for start in range(0, len(scaled), batch_size):
            batch = torch.tensor(
                scaled[start : start + batch_size],
                dtype=torch.float32,
                device=resolved_device,
            )
            values.extend((batch @ weights + intercept).detach().cpu().tolist())
        return tuple(float(value) for value in values)

    def to_dict(self) -> dict[str, object]:
        return {
            "feature_names": list(self.feature_names),
            "scaler_mean": list(self.scaler.mean),
            "scaler_std": list(self.scaler.std),
            "weights": list(self.weights),
            "intercept": self.intercept,
        }

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> AdditivePuModel:
        return cls(
            feature_names=tuple(str(item) for item in value["feature_names"]),
            scaler=TrainOnlyScaler(
                mean=tuple(float(item) for item in value["scaler_mean"]),
                std=tuple(float(item) for item in value["scaler_std"]),
            ),
            weights=tuple(float(item) for item in value["weights"]),
            intercept=float(value["intercept"]),
        )


def resolve_device(device: str) -> torch.device:
    """Resolve the frozen execution policy without requiring a GPU."""
    if device == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if device in {"cpu", "cuda"}:
        if device == "cuda" and not torch.cuda.is_available():
            raise RuntimeError("CUDA was requested but is not available")
        return torch.device(device)
    raise ValueError("device must be one of: auto, cpu, cuda")


def fit_additive_pu_ranker(
    train_features: list[list[float]],
    train_labels: list[str],
    train_protein_ids: list[str],
    feature_names: tuple[str, ...],
    config: AdditivePuConfig,
) -> AdditivePuModel:
    """Fit a deterministic additive ranker using only the supplied train rows."""
    if not train_features:
        raise ValueError("train_features must not be empty")
    if not (len(train_features) == len(train_labels) == len(train_protein_ids)):
        raise ValueError("training row count mismatch")
    if len(feature_names) != len(train_features[0]) or any(
        len(row) != len(feature_names) for row in train_features
    ):
        raise ValueError("feature_names width mismatch")
    if config.epochs < 1 or config.learning_rate <= 0.0:
        raise ValueError("epochs and learning_rate must be positive")

    device = resolve_device(config.device)
    torch.manual_seed(config.seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(config.seed)
    torch.use_deterministic_algorithms(True)
    scaler = TrainOnlyScaler.fit(train_features)
    features = torch.tensor(
        scaler.transform(train_features), dtype=torch.float32, device=device
    )
    weights = torch.zeros(features.shape[1], requires_grad=True, device=device)
    intercept = torch.zeros(1, requires_grad=True, device=device)
    optimizer = torch.optim.Adam([weights, intercept], lr=config.learning_rate)
    for _ in range(config.epochs):
        optimizer.zero_grad()
        logits = features @ weights + intercept
        loss = combined_pu_ranking_loss(
            logits,
            train_labels,
            train_protein_ids,
            config.class_prior,
            config.pairwise_weight,
        )
        loss = loss + config.l1 * weights.abs().sum()
        loss = loss + config.l2 * weights.square().sum()
        loss.backward()  # type: ignore[no-untyped-call]
        optimizer.step()
    return AdditivePuModel(
        feature_names=feature_names,
        scaler=scaler,
        weights=tuple(float(value) for value in weights.detach()),
        intercept=float(intercept.detach().item()),
    )
