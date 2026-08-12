"""Auditable PyTorch port of the published Sul-BertGRU GRU/CNN head."""

from __future__ import annotations

import random
from typing import cast

import numpy as np
import torch
from numpy.typing import NDArray
from torch import nn


class _CnnFeature(nn.Module):
    def __init__(self, input_channels: int) -> None:
        super().__init__()
        self.layers = nn.Sequential(
            nn.Conv1d(input_channels, 32, kernel_size=2, padding=1),
            nn.ReLU(),
            nn.AvgPool1d(kernel_size=4),
            nn.Conv1d(32, 64, kernel_size=2, padding=1),
            nn.ReLU(),
            nn.AvgPool1d(kernel_size=4),
            nn.Conv1d(64, 128, kernel_size=2, padding=1),
            nn.ReLU(),
            nn.AvgPool1d(kernel_size=4),
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        output = self.layers(values)
        return cast(torch.Tensor, output.reshape(output.shape[0], -1))


class SulBertGruNetwork(nn.Module):
    """Published shared-GRU, three-attention, two-CNN classification head."""

    def __init__(
        self,
        *,
        input_size: int,
        hidden_size: int = 128,
        num_layers: int = 3,
    ) -> None:
        super().__init__()
        if hidden_size % 4:
            raise ValueError("Sul-BertGRU hidden size must be divisible by four")
        self.gru = nn.GRU(
            input_size,
            hidden_size,
            num_layers,
            batch_first=True,
        )
        self.attention = nn.ModuleList(
            [
                nn.MultiheadAttention(
                    hidden_size,
                    num_heads=4,
                    batch_first=True,
                )
                for _ in range(3)
            ]
        )
        self.attention_cnn = _CnnFeature(61)
        self.gru_cnn = _CnnFeature(61)
        with torch.no_grad():
            dummy = torch.zeros(1, 61, hidden_size)
            flattened = self.attention_cnn(dummy).shape[1]
        self.classifier = nn.Sequential(
            nn.Linear(2 * flattened, 512),
            nn.ReLU(),
            nn.Dropout(0.5),
            nn.Linear(512, 2),
        )

    def forward(self, windows: torch.Tensor) -> torch.Tensor:
        if windows.ndim != 3 or windows.shape[1] != 31:
            raise ValueError("Sul-BertGRU requires 31-residue BERT embeddings")
        segments = (windows[:, :15], windows[:, 16:], windows)
        gru_outputs = [self.gru(segment)[0] for segment in segments]
        attended = [
            attention(values, values, values, need_weights=False)[0]
            for attention, values in zip(
                self.attention, gru_outputs, strict=True
            )
        ]
        attention_features = self.attention_cnn(torch.cat(attended, dim=1))
        gru_features = self.gru_cnn(torch.cat(gru_outputs, dim=1))
        # The reference implementation applies sigmoid before softmax.
        # Retain that unusual published head rather than substituting BCE logits.
        logits = self.classifier(torch.cat((attention_features, gru_features), dim=1))
        return torch.sigmoid(logits)


def bert_window_embeddings(
    windows: list[str],
    *,
    tokenizer: object,
    model: object,
    device_name: str,
    batch_size: int = 128,
) -> NDArray[np.float32]:
    """Embed 31-aa character windows, dropping BERT's CLS and SEP tokens."""
    import torch

    if not windows or any(len(window) != 31 for window in windows):
        raise ValueError("Sul-BertGRU BERT adapter requires 31-aa windows")
    device = torch.device(device_name)
    outputs: list[NDArray[np.float32]] = []
    for start in range(0, len(windows), batch_size):
        batch = [list(window) for window in windows[start : start + batch_size]]
        encoded = tokenizer(  # type: ignore[operator]
            batch,
            is_split_into_words=True,
            padding=True,
            return_tensors="pt",
        )
        inputs = {
            name: value.to(device) if hasattr(value, "to") else value
            for name, value in encoded.items()
        }
        with torch.no_grad():
            hidden = model(**inputs).last_hidden_state  # type: ignore[operator]
        values = hidden[:, 1:32, :]
        if values.shape[1] != 31:
            raise RuntimeError("BERT adapter lost a Sul-BertGRU residue token")
        if hasattr(values, "to"):
            values = values.to("cpu", dtype=torch.float32).numpy()
        array = np.asarray(values, dtype=np.float32)
        outputs.append(array)
    return np.concatenate(outputs, axis=0)


def fit_sul_bertgru(
    embeddings: NDArray[np.floating],
    *,
    train_indices: NDArray[np.int64],
    train_labels: NDArray[np.int64],
    score_indices: NDArray[np.int64],
    seed: int,
    hidden_size: int,
    num_layers: int,
    epochs: int,
    batch_size: int,
    learning_rate: float,
    device_name: str,
) -> NDArray[np.float32]:
    """Fit only on adapter train rows and return positive-class scores."""
    if embeddings.ndim != 3 or embeddings.shape[1] != 31:
        raise ValueError("Sul-BertGRU embeddings must have shape (n, 31, d)")
    if len(train_indices) != len(train_labels):
        raise ValueError("Sul-BertGRU train indices and labels differ")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)
    device = torch.device(device_name)
    if device.type == "cuda":
        # ``MultiheadAttention`` otherwise selects a non-deterministic
        # memory-efficient CUDA kernel on this platform.
        torch.backends.cuda.enable_flash_sdp(False)
        torch.backends.cuda.enable_mem_efficient_sdp(False)
        torch.backends.cuda.enable_math_sdp(True)
    network = SulBertGruNetwork(
        input_size=embeddings.shape[2],
        hidden_size=hidden_size,
        num_layers=num_layers,
    ).to(device)
    optimizer = torch.optim.Adam(network.parameters(), lr=learning_rate)
    loss_fn = nn.CrossEntropyLoss()
    generator = np.random.default_rng(seed)
    for _ in range(epochs):
        network.train()
        order = generator.permutation(len(train_indices))
        for start in range(0, len(order), batch_size):
            selected = order[start : start + batch_size]
            row_indices = train_indices[selected]
            batch = torch.as_tensor(
                np.asarray(embeddings[row_indices], dtype=np.float32),
                device=device,
            )
            labels = torch.as_tensor(
                train_labels[selected], dtype=torch.long, device=device
            )
            optimizer.zero_grad()
            loss = loss_fn(network(batch), labels)
            loss.backward()
            optimizer.step()
    network.eval()
    output = np.empty(len(score_indices), dtype=np.float32)
    with torch.no_grad():
        for start in range(0, len(score_indices), batch_size):
            stop = min(start + batch_size, len(score_indices))
            selected = score_indices[start:stop]
            batch = torch.as_tensor(
                np.asarray(embeddings[selected], dtype=np.float32),
                device=device,
            )
            probabilities = torch.softmax(network(batch), dim=1)[:, 1]
            output[start:stop] = probabilities.to("cpu").numpy()
    return output
