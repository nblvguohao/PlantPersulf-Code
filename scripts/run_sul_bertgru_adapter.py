#!/usr/bin/env python
"""Isolated, hash-bound Sul-BertGRU architecture adapter for Task 9.4.

The official repository lacks its ``Bert_config`` checkpoint.  This adapter
therefore retains the published GRU/attention/CNN head while using a frozen,
local ``google-bert/bert-base-cased`` snapshot.  Its output must be reported
as an architecture-adapted rerun, never as the original paper's checkpoint.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

import numpy as np

from plantpersulf.competitors.sul_bertgru import (
    bert_window_embeddings,
    fit_sul_bertgru,
)

_INPUT_FIELDS = (
    "partition",
    "global_protein_id",
    "cys_position",
    "sequence_window",
    "adapter_label",
)
_OUTPUT_FIELDS = ("partition", "global_protein_id", "cys_position", "score")


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--bert-snapshot", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--bert-batch-size", type=int, default=128)
    parser.add_argument("--train-batch-size", type=int, default=128)
    parser.add_argument("--hidden-size", type=int, default=128)
    parser.add_argument("--num-layers", type=int, default=3)
    parser.add_argument("--epochs", type=int, default=50)
    parser.add_argument("--learning-rate", type=float, default=0.002)
    return parser.parse_args()


def _read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != _INPUT_FIELDS:
            raise RuntimeError("Sul-BertGRU adapter input has invalid columns")
        rows = [dict(row) for row in reader]
    if not rows or any(
        row["partition"] not in {"train", "validation", "test"}
        or row["adapter_label"] not in {"positive", "putative_negative"}
        or len(row["sequence_window"]) != 31
        or row["sequence_window"][15] != "C"
        for row in rows
    ):
        raise RuntimeError("Sul-BertGRU adapter input is invalid")
    return rows


def _embed_all(
    rows: list[dict[str, str]],
    *,
    snapshot: Path,
    device: str,
    batch_size: int,
    destination: Path,
) -> np.memmap:
    from transformers import AutoModel, AutoTokenizer  # type: ignore[import-untyped]

    os.environ["TRANSFORMERS_OFFLINE"] = "1"
    tokenizer = AutoTokenizer.from_pretrained(snapshot, local_files_only=True)
    model = AutoModel.from_pretrained(snapshot, local_files_only=True).to(device)
    model.eval()
    matrix = np.lib.format.open_memmap(  # type: ignore[no-untyped-call]
        destination,
        mode="w+",
        dtype=np.float32,
        shape=(len(rows), 31, 768),
    )
    for start in range(0, len(rows), batch_size):
        stop = min(start + batch_size, len(rows))
        matrix[start:stop] = bert_window_embeddings(
            [row["sequence_window"] for row in rows[start:stop]],
            tokenizer=tokenizer,
            model=model,
            device_name=device,
            batch_size=batch_size,
        )
    matrix.flush()
    return matrix


def main() -> None:
    args = _arguments()
    if not args.bert_snapshot.is_dir():
        raise FileNotFoundError("Sul-BertGRU frozen BERT snapshot is missing")
    rows = _read_rows(args.input)
    embeddings_path = args.output.with_suffix(".bert_windows.npy")
    matrix = _embed_all(
        rows,
        snapshot=args.bert_snapshot,
        device=args.device,
        batch_size=args.bert_batch_size,
        destination=embeddings_path,
    )
    train_indices = np.asarray(
        [index for index, row in enumerate(rows) if row["partition"] == "train"],
        dtype=np.int64,
    )
    score_indices = np.asarray(
        [index for index, row in enumerate(rows) if row["partition"] != "train"],
        dtype=np.int64,
    )
    train_labels = np.asarray(
        [
            1 if rows[index]["adapter_label"] == "positive" else 0
            for index in train_indices
        ],
        dtype=np.int64,
    )
    if len(np.unique(train_labels)) != 2:
        raise RuntimeError("Sul-BertGRU adapter train partition needs both labels")
    scores = fit_sul_bertgru(
        matrix,
        train_indices=train_indices,
        train_labels=train_labels,
        score_indices=score_indices,
        seed=args.seed,
        hidden_size=args.hidden_size,
        num_layers=args.num_layers,
        epochs=args.epochs,
        batch_size=args.train_batch_size,
        learning_rate=args.learning_rate,
        device_name=args.device,
    )
    by_index = dict(zip(score_indices.tolist(), scores.tolist(), strict=True))
    with args.output.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=_OUTPUT_FIELDS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for index, row in enumerate(rows):
            if row["partition"] == "train":
                continue
            writer.writerow(
                {
                    "partition": row["partition"],
                    "global_protein_id": row["global_protein_id"],
                    "cys_position": row["cys_position"],
                    "score": by_index[index],
                }
            )


if __name__ == "__main__":
    main()
