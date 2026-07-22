"""Task 7 — frozen ESM-2 residue embeddings for benchmark cysteine positions.

Loads ESM-2 (esm2_t33_650M_UR50D, 1280-dim) once, embeds each benchmark
protein sequence, and records the per-position embedding at every labeled
cysteine coordinate. Embeddings are frozen and deterministic per model
checkpoint; no fine-tuning is performed.
"""

from __future__ import annotations

import csv
import os
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


def _esm_cache_dir() -> Path:
    return Path(
        os.environ.get("ESM2_CACHE_DIR", "data/interim/esm2_cache")
    )


def _load_cached_embedding(accession: str) -> Any | None:
    """Return (seq_len, 1280) numpy array or None if not cached."""
    cache_file = _esm_cache_dir() / f"{accession}.npy"
    if not cache_file.is_file():
        return None
    import numpy as np

    arr = np.load(cache_file)
    if arr.ndim != 2 or arr.shape[1] != EMBEDDING_DIM:
        return None
    return arr


def _save_cached_embedding(accession: str, embedding: Any) -> None:
    import numpy as np

    cache_file = _esm_cache_dir() / f"{accession}.npy"
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    np.save(cache_file, embedding.astype(np.float32))


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

    # Deduplicate proteins needed. Filter out the rare very-long protein
    # (the Arabidopsis proteome has ~70 sequences >2500 aa, up to 5400 aa)
    # whose ESM-2 O(L^2) attention matrices can exceed 16 GB GPU VRAM even
    # for a single protein, and whose CPU fallback would be impractically
    # slow (>10 min/protein).  At most 3 positive cysteine sites (0.8% of
    # the benchmark's 390 positives) are affected; their embeddings default
    # to [0.0]*1280 — a small, documented, easily audited limitation.
    MAX_PROTEIN_LENGTH = 2500
    needed = {row["protein_accession"] for row in labels}
    sequences = [
        (acc, proteome[acc])
        for acc in sorted(needed)
        if acc in proteome and len(proteome[acc]) <= MAX_PROTEIN_LENGTH
    ]
    skipped_long = sorted(
        acc for acc in needed if acc in proteome
        and len(proteome[acc]) > MAX_PROTEIN_LENGTH
    )
    if skipped_long:
        lengths = [len(proteome[acc]) for acc in skipped_long]
        print(
            f"    Skipping {len(skipped_long)} protein(s) > "
            f"{MAX_PROTEIN_LENGTH} aa (max {max(lengths)} aa): "
            f"{', '.join(skipped_long[:5])}"
            + (f" ... ({len(skipped_long) - 5} more)" if len(skipped_long) > 5 else ""),
            flush=True,
        )
    if not sequences:
        return ()

    # ── cache lookup ──────────────────────────────────────────────────
    embed_map: dict[tuple[str, int], tuple[float, ...]] = {}
    cached_count = 0
    uncached_seqs: list[tuple[str, str]] = []
    for acc, seq in sequences:
        cached = _load_cached_embedding(acc)
        if cached is not None:
            for pos in range(len(seq)):
                embed_map[(acc, pos + 1)] = tuple(
                    float(v) for v in cached[pos]
                )
            cached_count += 1
        else:
            uncached_seqs.append((acc, seq))
    if cached_count:
        print(
            f"    ESM cache hit: {cached_count} / {len(sequences)} proteins "
            f"({len(embed_map)} residues loaded from disk)",
            flush=True,
        )

    if not uncached_seqs:
        return _rows_from_map(labels, embed_map)

    # ── compute uncached proteins ─────────────────────────────────────

    device = _select_device(torch)
    model, alphabet = _get_esm_model(esm, device)
    batch_converter = alphabet.get_batch_converter()

    batch_data = [(acc, seq) for acc, seq in uncached_seqs]
    _, _, batch_tokens = batch_converter(batch_data)

    def _forward(net: Any, tokens: Any, on_device: Any) -> Any:
        with torch.no_grad():
            out = net(tokens.to(on_device), repr_layers=[33], return_contacts=False)
        return out["representations"][33].to("cpu")

    token_representations: Any
    try:
        token_representations = _forward(model, batch_tokens, device)
    except torch.cuda.OutOfMemoryError:
        torch.cuda.empty_cache()
        max_len = max(len(seq) for _, seq in uncached_seqs)
        print(
            f"    CUDA OOM on batch of {len(uncached_seqs)} (max protein "
            f"length={max_len} aa); retrying sequences individually ...",
            flush=True,
        )
        per_seq_reps = []
        for acc, seq in uncached_seqs:
            _, _, single_tokens = batch_converter([(acc, seq)])
            try:
                per_seq_reps.append(_forward(model, single_tokens, device)[0])
            except torch.cuda.OutOfMemoryError:
                torch.cuda.empty_cache()
                cpu_model, _ = _get_esm_model(esm, torch.device("cpu"))
                per_seq_reps.append(
                    _forward(cpu_model, single_tokens, torch.device("cpu"))[0]
                )
        token_representations = per_seq_reps

    # ── save new embeddings to cache ──────────────────────────────────
    for i, (acc, seq) in enumerate(uncached_seqs):
        if isinstance(token_representations, list):
            rep = token_representations[i][1 : len(seq) + 1]
        else:
            rep = token_representations[i, 1 : len(seq) + 1]
        _save_cached_embedding(acc, rep.numpy())
        for pos in range(1, len(seq) + 1):
            embed_map[(acc, pos)] = tuple(float(v) for v in rep[pos - 1])

    return _rows_from_map(labels, embed_map)


def _rows_from_map(
    labels: list[dict[str, str]],
    embed_map: dict[tuple[str, int], tuple[float, ...]],
) -> tuple[ESM2FeatureRow, ...]:
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
