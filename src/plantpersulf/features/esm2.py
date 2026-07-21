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

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model.eval()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    batch_converter = alphabet.get_batch_converter()

    # Batch all sequences at once
    batch_data = [(acc, seq) for acc, seq in sequences]
    _, _, batch_tokens = batch_converter(batch_data)
    batch_tokens = batch_tokens.to(device)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=False)
    token_representations = results["representations"][33]

    # Map (protein, position) → embedding
    embed_map: dict[tuple[str, int], tuple[float, ...]] = {}
    for i, (acc, seq) in enumerate(sequences):
        rep = token_representations[i, 1 : len(seq) + 1]  # strip BOS/EOS
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
