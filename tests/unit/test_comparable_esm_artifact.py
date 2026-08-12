"""Software-policy tests for the shared-panel ESM feature artifact."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from plantpersulf.features.esm2_windows import (
    build_cysteine_windows,
    load_window_embedding_artifact,
    write_window_embedding_artifact,
)


def test_windows_keep_species_key_and_center_the_registered_cysteine() -> None:
    windows = build_cysteine_windows(
        {("species_a|P1", 3), ("species_b|P1", 7)},
        {
            "species_a|P1": "MMCMMMMMM",
            "species_b|P1": "MMMMMMCMM",
        },
    )

    assert [row.site_key for row in windows] == [
        ("species_a|P1", 3),
        ("species_b|P1", 7),
    ]
    assert all(len(row.sequence_window) == 31 for row in windows)
    assert all(row.sequence_window[15] == "C" for row in windows)
    assert windows[0].sequence_window != windows[1].sequence_window


def test_artifact_round_trip_binds_checkpoint_inputs_and_exact_sites(
    tmp_path: Path,
) -> None:
    windows = build_cysteine_windows(
        {("species_a|P1", 3), ("species_b|P2", 7)},
        {
            "species_a|P1": "MMCMMMMMM",
            "species_b|P2": "MMMMMMCMM",
        },
    )
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"policy-checkpoint")

    def embed_batch(sequences: list[str]) -> np.ndarray:
        return np.asarray(
            [
                [float(len(sequence)), float(sequence.count("X"))]
                for sequence in sequences
            ],
            dtype=np.float32,
        )

    output = tmp_path / "artifact"
    write_window_embedding_artifact(
        windows,
        output,
        model_name="POLICY_ESM",
        model_checkpoint=checkpoint,
        input_sha256={"species_a": "a" * 64, "species_b": "b" * 64},
        embedding_dim=2,
        batch_size=1,
        embed_batch=embed_batch,
    )

    loaded = load_window_embedding_artifact(
        output,
        required_site_keys={row.site_key for row in windows},
    )

    assert set(loaded) == {row.site_key for row in windows}
    assert loaded[("species_a|P1", 3)].tolist() == [31.0, 22.0]
    assert loaded[("species_b|P2", 7)].tolist() == [31.0, 22.0]
