"""Release audits for the multispecies v2 runtime manifest."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from plantpersulf.workflows.multispecies_v2 import (
    audit_v2_run_manifest,
    write_run_manifest,
)


def test_release_manifest_rejects_tampered_checkpoint_hash(tmp_path: Path) -> None:
    """A changed checkpoint makes a scientific execution manifest invalid."""
    checkpoint = tmp_path / "fold0.checkpoint.json"
    checkpoint.write_text('{"policy":"marker"}\n', encoding="utf-8")
    manifest_path = tmp_path / "manifest.json"
    write_run_manifest(
        manifest_path,
        config_sha256="a" * 64,
        split_sha256="b" * 64,
        code_revision="commit-marker",
        input_sha256={"development_sites": "c" * 64},
        command=["python", "script.py"],
        device="cpu",
        environment={"python": "3.10"},
        gpu_map=(),
        artifacts={"checkpoint": checkpoint},
        checkpoint_paths={"fold0": checkpoint},
    )
    checkpoint.write_text('{"policy":"changed"}\n', encoding="utf-8")

    with pytest.raises(RuntimeError, match="artifact SHA256 mismatch"):
        audit_v2_run_manifest(manifest_path)


def test_release_manifest_rejects_missing_required_field(tmp_path: Path) -> None:
    """A release manifest cannot silently omit its environment provenance."""
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema_version": 2,
                "config_sha256": "a" * 64,
                "split_sha256": "b" * 64,
                "code_revision": "commit-marker",
                "input_sha256": {"development_sites": "c" * 64},
                "command": ["python", "script.py"],
                "device": "cpu",
                "gpu_map": [],
                "dirty": False,
                "dirty_paths": [],
                "artifacts": {},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="required fields"):
        audit_v2_run_manifest(manifest)


def test_manifest_refuses_to_replace_nonidentical_existing_output(
    tmp_path: Path,
) -> None:
    """A second invocation cannot silently replace a different release record."""
    path = tmp_path / "manifest.json"
    path.write_text('{"policy":"old"}\n', encoding="utf-8")

    with pytest.raises(RuntimeError, match="non-identical manifest"):
        write_run_manifest(
            path,
            config_sha256="a" * 64,
            split_sha256="b" * 64,
            code_revision="commit-marker",
            input_sha256={"development_sites": "c" * 64},
            command=["python", "script.py"],
            device="cpu",
        )
