"""Provenance-preserving preparation of real global MMseqs2 inputs."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import ModuleType

import pytest

from plantpersulf.proteomics.global_mmseqs import (
    ReferenceProteomeInput,
    easy_cluster_command,
    global_cluster_rows_from_mmseqs_pairs,
    prepare_namespaced_fasta,
)


def _load_cluster_builder() -> ModuleType:
    path = Path("scripts/build_global_multispecies_clusters.py")
    spec = importlib.util.spec_from_file_location("global_mmseqs_builder", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_namespaced_fasta_retains_each_source_hash_and_avoids_accession_collision(
    tmp_path: Path,
) -> None:
    """Removing the species prefix would make these two source records collide."""
    arab = tmp_path / "arab.fasta"
    tomato = tmp_path / "tomato.fasta"
    arab.write_text(">sp|P_SHARED|marker\nACDE\n", encoding="utf-8")
    tomato.write_text(">sp|P_SHARED|marker\nCDEA\n", encoding="utf-8")
    output = tmp_path / "global.fasta"

    rows = prepare_namespaced_fasta(
        [
            ReferenceProteomeInput(
                "arabidopsis", arab, hashlib.sha256(arab.read_bytes()).hexdigest()
            ),
            ReferenceProteomeInput(
                "tomato", tomato, hashlib.sha256(tomato.read_bytes()).hexdigest()
            ),
        ],
        output,
    )

    assert [row.global_protein_id for row in rows] == [
        "arabidopsis|P_SHARED",
        "tomato|P_SHARED",
    ]
    assert output.read_text(encoding="utf-8").splitlines() == [
        ">arabidopsis|P_SHARED",
        "ACDE",
        ">tomato|P_SHARED",
        "CDEA",
    ]


def test_mmseqs_pair_parser_rejects_missing_global_protein(tmp_path: Path) -> None:
    """A partial MMseqs2 export must not become a silently incomplete split."""
    fasta = tmp_path / "arab.fasta"
    fasta.write_text(">sp|P1|marker\nACDE\n>sp|P2|marker\nCDEA\n", encoding="utf-8")
    namespaced = prepare_namespaced_fasta(
        [
            ReferenceProteomeInput(
                "arabidopsis", fasta, hashlib.sha256(fasta.read_bytes()).hexdigest()
            )
        ],
        tmp_path / "global.fasta",
    )
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("arabidopsis|P1\tarabidopsis|P1\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="missing proteins"):
        global_cluster_rows_from_mmseqs_pairs(namespaced, pairs, 0.3)


def test_mmseqs_command_locks_the_required_homology_thresholds(tmp_path: Path) -> None:
    """Changing coverage or cov-mode would change the meaning of a cluster."""
    command = easy_cluster_command(
        "mmseqs",
        tmp_path / "global.fasta",
        tmp_path / "result",
        tmp_path / "tmp",
        identity_threshold=0.3,
    )

    assert command[-8:] == [
        "--min-seq-id",
        "0.3",
        "-c",
        "0.5",
        "--cov-mode",
        "0",
        "--threads",
        "8",
    ]


def test_versioned_cluster_output_derives_matching_provenance_paths() -> None:
    """A v2 cluster table must never write its audit record under a v1 name."""
    builder = _load_cluster_builder()

    assert builder._provenance_path(
        Path("data/processed/clusters/global_multispecies_mmseqs2_30_v2.tsv"),
        "manifest.json",
    ) == Path("data/processed/clusters/global_multispecies_mmseqs2_v2.manifest.json")


def test_build_registers_mmseqs_version_commands_and_accession_mapping(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An audit record without version, commands, or mapping is not reproducible."""
    builder = _load_cluster_builder()
    proteomes: list[dict[str, str]] = []
    for species in ("arabidopsis", "tomato", "rice", "magnaporthe"):
        fasta = tmp_path / f"{species}.fasta"
        fasta.write_text(f">sp|P_{species}|marker\nACDE\n", encoding="utf-8")
        proteomes.append(
            {
                "species": species,
                "path": str(fasta),
                "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
            }
        )
    config = tmp_path / "multispecies_v2.yaml"
    config.write_text(
        "version: 2\n"
        "reference_proteomes:\n"
        + "".join(
            f"  - species: {row['species']}\n"
            f"    path: '{Path(row['path']).as_posix()}'\n"
            f"    sha256: {row['sha256']}\n"
            for row in proteomes
        )
        + "global_mmseqs2:\n"
        "  primary_identity_threshold: 0.3\n"
        "  sensitivity_identity_thresholds: [0.2, 0.4]\n"
        "  coverage: 0.5\n"
        "  cov_mode: 0\n"
        "  threads: 8\n"
        f"  cluster_table: '{(tmp_path / 'global_30.tsv').as_posix()}'\n",
        encoding="utf-8",
    )

    def fake_run(command: list[str], *, check: bool, **_: object) -> object:
        assert check
        if command[1:] == []:
            return type("Result", (), {"stdout": "MMseqs2 Version 15-6f452\n"})()
        prefix = Path(command[3])
        members = [
            line[1:]
            for line in Path(command[2]).read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]
        Path(str(prefix) + "_cluster.tsv").write_text(
            "".join(f"{member}\t{member}\n" for member in members),
            encoding="utf-8",
        )
        return object()

    monkeypatch.setattr(builder.subprocess, "run", fake_run)

    builder.build(config, "mmseqs")

    manifest = json.loads(
        (tmp_path / "global_30.manifest.json").read_text(encoding="utf-8")
    )
    assert manifest["mmseqs_version"] == "MMseqs2 Version 15-6f452"
    assert len(manifest["commands"]) == 3
    assert all("--min-seq-id" in command for command in manifest["commands"])
    mapping = tmp_path / "global_30.accession_mapping.tsv"
    assert mapping.is_file()
    assert manifest["accession_mapping"]["sha256"] == hashlib.sha256(
        mapping.read_bytes()
    ).hexdigest()
