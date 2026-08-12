"""Software-only tests for the isolated Sul-BertGRU architecture port."""

from __future__ import annotations

import numpy as np

from plantpersulf.competitors.sul_bertgru import (
    bert_window_embeddings,
    fit_sul_bertgru,
)


def test_port_trains_only_on_train_rows_and_scores_requested_rows() -> None:
    rng = np.random.default_rng(20260811)
    embeddings = rng.normal(size=(12, 31, 16)).astype(np.float32)
    labels = np.asarray([1, 0, 1, 0, 1, 0, 1, 0], dtype=np.int64)

    scores = fit_sul_bertgru(
        embeddings,
        train_indices=np.arange(8),
        train_labels=labels,
        score_indices=np.arange(8, 12),
        seed=0,
        hidden_size=64,
        num_layers=1,
        epochs=1,
        batch_size=4,
        learning_rate=0.002,
        device_name="cpu",
    )

    assert scores.shape == (4,)
    assert np.isfinite(scores).all()
    assert ((0.0 <= scores) & (scores <= 1.0)).all()


def test_bert_adapter_returns_exact_31_residue_windows() -> None:
    class Tokenizer:
        def __call__(self, batch, **kwargs):  # type: ignore[no-untyped-def]
            assert kwargs["is_split_into_words"] is True
            assert all(len(row) == 31 for row in batch)
            return {"input_ids": np.zeros((len(batch), 33), dtype=np.int64)}

    class Model:
        def __call__(self, **kwargs):  # type: ignore[no-untyped-def]
            size = kwargs["input_ids"].shape[0]
            return type(
                "Output",
                (),
                {"last_hidden_state": np.ones((size, 33, 4), dtype=np.float32)},
            )()

    values = bert_window_embeddings(
        ["X" * 15 + "C" + "M" * 15],
        tokenizer=Tokenizer(),
        model=Model(),
        device_name="cpu",
    )

    assert values.shape == (1, 31, 4)
    assert values.dtype == np.float32
