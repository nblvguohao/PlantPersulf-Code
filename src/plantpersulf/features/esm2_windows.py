"""Frozen ESM-2 embeddings for the exact Task 9.4 comparison sites."""

from __future__ import annotations

import csv
import hashlib
import json
import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import cast

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class CysteineWindow:
    global_protein_id: str
    cys_position: int
    sequence_window: str

    @property
    def site_key(self) -> tuple[str, int]:
        return self.global_protein_id, self.cys_position


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_cysteine_windows(
    site_keys: set[tuple[str, int]],
    sequences: dict[str, str],
    *,
    flank: int = 15,
) -> tuple[CysteineWindow, ...]:
    """Build deterministic centered windows from registered reference sequences."""
    result: list[CysteineWindow] = []
    for global_protein_id, position in sorted(site_keys):
        sequence = sequences.get(global_protein_id)
        if sequence is None:
            raise RuntimeError(f"missing registered sequence: {global_protein_id}")
        if position < 1 or position > len(sequence) or sequence[position - 1] != "C":
            raise RuntimeError(
                f"comparison coordinate is not Cys: {(global_protein_id, position)}"
            )
        start = position - 1 - flank
        end = position + flank
        window = "X" * max(0, -start)
        window += sequence[max(0, start) : min(len(sequence), end)]
        window += "X" * max(0, end - len(sequence))
        if len(window) != 2 * flank + 1 or window[flank] != "C":
            raise RuntimeError(
                f"invalid comparison window: {(global_protein_id, position)}"
            )
        result.append(
            CysteineWindow(
                global_protein_id=global_protein_id,
                cys_position=position,
                sequence_window=window,
            )
        )
    return tuple(result)


def write_window_embedding_artifact(
    windows: tuple[CysteineWindow, ...],
    output_directory: Path,
    *,
    model_name: str,
    model_checkpoint: Path,
    input_sha256: dict[str, str],
    embedding_dim: int,
    batch_size: int,
    embed_batch: Callable[[list[str]], NDArray[np.float32]],
    device: str = "unspecified",
) -> None:
    """Atomically write a hash-bound, memory-mappable site feature artifact."""
    if output_directory.exists():
        raise FileExistsError(
            f"ESM feature artifact already exists: {output_directory}"
        )
    if not model_checkpoint.is_file():
        raise FileNotFoundError(f"ESM checkpoint is missing: {model_checkpoint}")
    if not windows or embedding_dim < 1 or batch_size < 1:
        raise ValueError("ESM artifact requires windows, dimension, and batch size")
    output_directory.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        dir=output_directory.parent, prefix=".esm2_windows_"
    ) as temporary_name:
        temporary = Path(temporary_name)
        keys_path = temporary / "site_keys.tsv"
        matrix_path = temporary / "embeddings.npy"
        with keys_path.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("row_index", "global_protein_id", "cys_position"))
            for index, row in enumerate(windows):
                writer.writerow((index, row.global_protein_id, row.cys_position))
        matrix = np.lib.format.open_memmap(  # type: ignore[no-untyped-call]
            matrix_path,
            mode="w+",
            dtype=np.float32,
            shape=(len(windows), embedding_dim),
        )
        for start in range(0, len(windows), batch_size):
            stop = min(start + batch_size, len(windows))
            embedded = np.asarray(
                embed_batch([row.sequence_window for row in windows[start:stop]]),
                dtype=np.float32,
            )
            if embedded.shape != (stop - start, embedding_dim):
                raise RuntimeError("ESM embedder returned an invalid batch shape")
            if not np.isfinite(embedded).all():
                raise RuntimeError("ESM embedder returned non-finite values")
            matrix[start:stop] = embedded
        matrix.flush()
        del matrix
        manifest = {
            "artifact_version": 1,
            "biological_values_modified": False,
            "device": device,
            "embedding_dim": embedding_dim,
            "input_sha256": dict(sorted(input_sha256.items())),
            "model_checkpoint": str(model_checkpoint),
            "model_checkpoint_sha256": _sha256(model_checkpoint),
            "model_name": model_name,
            "site_count": len(windows),
            "site_keys_sha256": _sha256(keys_path),
            "embeddings_sha256": _sha256(matrix_path),
            "window_length": len(windows[0].sequence_window),
        }
        (temporary / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(output_directory)


def load_window_embedding_artifact(
    directory: Path,
    *,
    required_site_keys: set[tuple[str, int]],
) -> dict[tuple[str, int], NDArray[np.float32]]:
    """Audit an ESM artifact and expose only exact, hash-verified rows."""
    try:
        manifest = json.loads((directory / "manifest.json").read_text("utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError("ESM feature manifest is invalid") from exc
    keys_path = directory / "site_keys.tsv"
    matrix_path = directory / "embeddings.npy"
    if _sha256(keys_path) != manifest.get("site_keys_sha256"):
        raise RuntimeError("ESM site-key SHA256 mismatch")
    if _sha256(matrix_path) != manifest.get("embeddings_sha256"):
        raise RuntimeError("ESM embedding SHA256 mismatch")
    matrix = np.load(matrix_path, mmap_mode="r")
    expected_shape = (manifest.get("site_count"), manifest.get("embedding_dim"))
    if matrix.shape != expected_shape or matrix.dtype != np.float32:
        raise RuntimeError("ESM embedding matrix metadata mismatch")
    indexed: dict[tuple[str, int], int] = {}
    with keys_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_fields = ("row_index", "global_protein_id", "cys_position")
        if tuple(reader.fieldnames or ()) != expected_fields:
            raise RuntimeError("ESM site-key table has invalid columns")
        for expected_index, row in enumerate(reader):
            if int(row["row_index"]) != expected_index:
                raise RuntimeError("ESM site-key row index is not contiguous")
            key = (row["global_protein_id"], int(row["cys_position"]))
            if key in indexed:
                raise RuntimeError(f"duplicate ESM site key: {key}")
            indexed[key] = expected_index
    if set(indexed) != required_site_keys:
        raise RuntimeError("ESM artifact does not cover the exact shared panels")
    return {key: matrix[index] for key, index in indexed.items()}


def extract_esm2_window_artifact(
    windows: tuple[CysteineWindow, ...],
    output_directory: Path,
    *,
    model_checkpoint: Path,
    input_sha256: dict[str, str],
    batch_size: int,
) -> None:
    """Run the frozen ESM-2 650M checkpoint on centered 31-aa windows."""
    import esm  # type: ignore[import-untyped]
    import torch

    from plantpersulf.features.esm2 import (
        EMBEDDING_DIM,
        ESM_MODEL_NAME,
        _get_esm_model,
        _select_device,
    )

    device = _select_device(torch)
    model, alphabet = _get_esm_model(esm, device)
    converter = alphabet.get_batch_converter()

    def embed_batch(sequences: list[str]) -> NDArray[np.float32]:
        batch = [(str(index), sequence) for index, sequence in enumerate(sequences)]
        _, _, tokens = converter(batch)
        with torch.no_grad():
            output = model(
                tokens.to(device), repr_layers=[33], return_contacts=False
            )
        center = output["representations"][33][:, 16, :]
        return cast(
            NDArray[np.float32],
            center.to("cpu", dtype=torch.float32).numpy(),
        )

    write_window_embedding_artifact(
        windows,
        output_directory,
        model_name=ESM_MODEL_NAME,
        model_checkpoint=model_checkpoint,
        input_sha256=input_sha256,
        embedding_dim=EMBEDDING_DIM,
        batch_size=batch_size,
        embed_batch=embed_batch,
        device=str(device),
    )
