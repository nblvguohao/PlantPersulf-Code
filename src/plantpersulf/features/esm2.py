"""Task 7 — frozen ESM-2 residue embeddings for benchmark cysteine positions.

Loads ESM-2 (esm2_t33_650M_UR50D, 1280-dim) once, embeds each benchmark
protein sequence, and records the per-position embedding at every labeled
cysteine coordinate. Embeddings are frozen and deterministic per model
checkpoint; no fine-tuning is performed.
"""

from __future__ import annotations

import csv
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ESM_MODEL_NAME = "esm2_t33_650M_UR50D"
EMBEDDING_DIM = 1280

BENCHMARK_FIELDS = (
    "protein_accession",
    "cys_position_in_protein",
    "label",
    "study_accession",
    "evidence_level",
    "source_sha256",
)


# Cache the 650M model+alphabet per device: the runner extracts embeddings in
# many small chunks (hundreds per experiment), and reloading the checkpoint from
# disk and re-moving it to the GPU each time would dominate the runtime.
_ESM_MODEL_CACHE: dict[str, tuple[Any, Any]] = {}


def _get_esm_model(esm_mod: Any, device: Any) -> tuple[Any, Any]:
    key = str(device)
    cached = _ESM_MODEL_CACHE.get(key)
    if cached is not None:
        return cached
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model, alphabet = esm_mod.pretrained.esm2_t33_650M_UR50D()
    model.eval()
    model = model.to(device)
    _ESM_MODEL_CACHE[key] = (model, alphabet)
    return model, alphabet


def _select_device(torch_mod: object) -> object:
    """Pick the compute device for ESM inference.

    ``PLANTPERSULF_DEVICE`` (e.g. ``cuda``, ``cuda:0``, ``cpu``) overrides
    auto-detection; otherwise CUDA is used when available, else CPU. The frozen
    embeddings are model-checkpoint deterministic; note that CPU and GPU results
    can differ in the last floating-point digits, so a benchmark should be
    embedded on a single device class.
    """
    import os

    forced = os.environ.get("PLANTPERSULF_DEVICE")
    if forced:
        return torch_mod.device(forced)  # type: ignore[attr-defined]
    if torch_mod.cuda.is_available():  # type: ignore[attr-defined]
        return torch_mod.device("cuda")  # type: ignore[attr-defined]
    return torch_mod.device("cpu")  # type: ignore[attr-defined]


@dataclass(frozen=True)
class ESM2FeatureRow:
    protein_accession: str
    cys_position: int
    label: str
    embedding: tuple[float, ...]


def _load_labels(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != BENCHMARK_FIELDS:
            raise RuntimeError(f"benchmark labels have invalid columns: {path}")
        return [dict(row) for row in reader]


def _load_proteome(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    cur_header = ""
    cur_lines: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if cur_header:
                acc = cur_header.strip().split("|")[1]
                sequences[acc] = "".join(cur_lines)
            cur_header = line
            cur_lines = []
        elif line:
            cur_lines.append(line)
    if cur_header:
        acc = cur_header.strip().split("|")[1]
        sequences[acc] = "".join(cur_lines)
    return sequences


def extract_esm2_embeddings(
    labels_path: Path,
    proteome_path: Path,
) -> tuple[ESM2FeatureRow, ...]:
    """Extract per-cysteine ESM-2 embeddings for all labeled proteins."""
    import esm  # type: ignore[import-untyped]
    import torch

    labels = _load_labels(labels_path)
    proteome = _load_proteome(proteome_path)

    # Deduplicate proteins needed
    needed = {row["protein_accession"] for row in labels}
    sequences = [
        (acc, proteome[acc]) for acc in sorted(needed) if acc in proteome
    ]
    if not sequences:
        return ()

    device = _select_device(torch)
    model, alphabet = _get_esm_model(esm, device)
    batch_converter = alphabet.get_batch_converter()

    # Batch all sequences at once
    batch_data = [(acc, seq) for acc, seq in sequences]
    _, _, batch_tokens = batch_converter(batch_data)

    def _forward(net: Any, tokens: Any, on_device: Any) -> Any:
        with torch.no_grad():
            out = net(tokens.to(on_device), repr_layers=[33], return_contacts=False)
        return out["representations"][33].to("cpu")

    token_representations: Any
    try:
        token_representations = _forward(model, batch_tokens, device)
    except torch.cuda.OutOfMemoryError:
        # ESM-2's O(L^2) attention can exceed even a 16GB GPU when several
        # very long proteins land in the same padded batch (the Arabidopsis
        # proteome has ~160 sequences >2000 aa, up to 5400 aa) — batching
        # multiplies the padded length's memory cost by the batch size, even
        # though any *one* of these sequences fits on the GPU alone. Rather
        # than losing the whole chunk (and silently degrading downstream rows
        # to a missing-embedding default), retry each sequence individually,
        # re-tokenized alone so it carries no other sequence's padding. Only
        # a sequence that is itself too long for the GPU falls further back
        # to CPU. Every embedding here is still real ESM-2 output, just
        # computed without the other sequences' padding overhead.
        torch.cuda.empty_cache()
        max_len = max(len(seq) for _, seq in sequences)
        print(
            f"    CUDA OOM on batch of {len(sequences)} (max protein "
            f"length={max_len} aa); retrying sequences individually ...",
            flush=True,
        )
        per_seq_reps = []
        for acc, seq in sequences:
            _, _, single_tokens = batch_converter([(acc, seq)])
            try:
                # [0] drops the batch-of-1 dim so each entry is (len+2, dim),
                # matching the per-item shape the batch path yields via [i].
                per_seq_reps.append(_forward(model, single_tokens, device)[0])
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                cpu_model, _ = _get_esm_model(esm, torch.device("cpu"))
                per_seq_reps.append(
                    _forward(cpu_model, single_tokens, torch.device("cpu"))[0]
                )
        token_representations = per_seq_reps

    # Map (protein, position) → embedding. token_representations is either a
    # single (batch, max_len+2, dim) tensor or a list of one-sequence
    # (len_i+2, dim) tensors; token_representations[i] is (len+2, dim) in
    # both cases, so the same slicing works for both.
    embed_map: dict[tuple[str, int], tuple[float, ...]] = {}
    for i, (acc, seq) in enumerate(sequences):
        rep = token_representations[i][1 : len(seq) + 1]  # strip BOS/EOS
        for pos in range(1, len(seq) + 1):
            embed_map[(acc, pos)] = tuple(float(v) for v in rep[pos - 1])

    rows: list[ESM2FeatureRow] = []
    for row in labels:
        acc = row["protein_accession"]
        pos = int(row["cys_position_in_protein"])
        emb = embed_map.get((acc, pos))
        if emb is None:
            continue
        rows.append(
            ESM2FeatureRow(
                protein_accession=acc,
                cys_position=pos,
                label=row["label"],
                embedding=emb,
            )
        )
    return tuple(rows)
