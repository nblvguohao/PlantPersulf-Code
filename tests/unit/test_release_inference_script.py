"""Gate 0 release package: inference_script input validation contract.

The blind-time scorer must reject any input that deviates from
feature_schema.json before a single row is scored — column set, unique site
keys, finite sequence features, and an explicit structure mask.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

_RELEASE_DIR = (
    Path(__file__).resolve().parents[2]
    / "results"
    / "candidates"
    / "multispecies_v2_candidate_release_v1"
)
_SCHEMA_PATH = _RELEASE_DIR / "feature_schema.json"
_BUNDLE_PATH = _RELEASE_DIR / "model_weights" / "structure_ranker_bundle.pt"

pytestmark = pytest.mark.skipif(
    not _SCHEMA_PATH.is_file() or not _BUNDLE_PATH.is_file(),
    reason="release package artifacts not present",
)


def _load_inference_script():
    spec = importlib.util.spec_from_file_location(
        "release_inference_script", _RELEASE_DIR / "inference_script.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def inference_script():
    return _load_inference_script()


def _valid_rows() -> list[dict[str, str]]:
    return [
        {
            "site_key": "tomato|SolycTest000001",
            "cys_position": "42",
            "hydrophobicity": "0.351",
            "cys_density": "0.0132",
            "local_positive_charge_density": "0.214",
            "contact_number_proxy": "11.7",
            "plddt": "84.2",
            "has_structure": "true",
        },
        {
            "site_key": "tomato|SolycTest000002",
            "cys_position": "21",
            "hydrophobicity": "0.076",
            "cys_density": "0.0041",
            "local_positive_charge_density": "0.091",
            "contact_number_proxy": "0.0",
            "plddt": "0.0",
            "has_structure": "false",
        },
    ]


def _write_input(
    tmp_path, rows: list[dict[str, str]], fieldnames: tuple[str, ...] | None = None
) -> Path:
    path = tmp_path / "candidates.tsv"
    import csv

    names = fieldnames or _load_inference_script().REQUIRED_INPUT_COLUMNS
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(names), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            {key: row[key] for key in names if key in row} for row in rows
        )
    return path


def test_schema_and_bundle_exist() -> None:
    schema = _SCHEMA_PATH.read_text(encoding="utf-8")
    assert '"schema_version": 1' in schema
    assert '"has_structure"' in schema


def test_accepts_valid_rows(inference_script, tmp_path) -> None:
    rows = inference_script._load_rows(_write_input(tmp_path, _valid_rows()))
    assert len(rows) == 2


def test_rejects_missing_column(inference_script, tmp_path) -> None:
    rows = [dict(row) for row in _valid_rows()]
    for row in rows:
        del row["plddt"]
    fieldnames = tuple(
        name
        for name in inference_script.REQUIRED_INPUT_COLUMNS
        if name != "plddt"
    )
    with pytest.raises(RuntimeError, match="columns"):
        inference_script._load_rows(
            _write_input(tmp_path, rows, fieldnames=fieldnames)
        )


def test_rejects_duplicate_site(inference_script, tmp_path) -> None:
    rows = _valid_rows()
    rows.append(dict(rows[0]))
    with pytest.raises(RuntimeError, match="duplicate site"):
        inference_script._load_rows(_write_input(tmp_path, rows))


def test_rejects_invalid_structure_mask(inference_script, tmp_path) -> None:
    rows = _valid_rows()
    rows[0]["has_structure"] = "maybe"
    with pytest.raises(RuntimeError, match="has_structure"):
        inference_script._load_rows(_write_input(tmp_path, rows))


def test_rejects_nonfinite_sequence_feature(inference_script, tmp_path) -> None:
    rows = _valid_rows()
    rows[0]["hydrophobicity"] = "nan"
    with pytest.raises(RuntimeError, match="non-finite"):
        inference_script._load_rows(_write_input(tmp_path, rows))


def test_rejects_nonpositive_position(inference_script, tmp_path) -> None:
    rows = _valid_rows()
    rows[0]["cys_position"] = "0"
    with pytest.raises(RuntimeError, match="cys_position"):
        inference_script._load_rows(_write_input(tmp_path, rows))


def test_rejects_empty_input(inference_script, tmp_path) -> None:
    with pytest.raises(RuntimeError, match="no data rows"):
        inference_script._load_rows(_write_input(tmp_path, []))


def test_tampered_bundle_rejected(inference_script, tmp_path) -> None:
    import hashlib
    import json

    bundle = _BUNDLE_PATH.read_bytes()
    assert len(bundle) > 100
    tampered = bytearray(bundle)
    tampered[100] ^= 0xFF

    tampered_path = tmp_path / "tampered.pt"
    tampered_path.write_bytes(bytes(tampered))
    manifest_path = _RELEASE_DIR / "model_weights" / "fit_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

    # Point the fit manifest at the tampered file, keeping the pinned hash.
    tampered_manifest = dict(manifest)
    tampered_manifest_path = tmp_path / "fit_manifest.json"
    tampered_manifest_path.write_text(
        json.dumps(tampered_manifest, indent=2), encoding="utf-8"
    )
    assert hashlib.sha256(bytes(tampered)).hexdigest() != manifest["bundle_sha256"]
    with pytest.raises(RuntimeError, match="SHA256 mismatch"):
        inference_script._verify_bundle(tampered_path, tampered_manifest_path)
