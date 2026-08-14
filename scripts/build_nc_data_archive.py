#!/usr/bin/env python
"""Gate 6 (Task 7) — assemble the NC entry-gate data-transparency archive.

Builds ``results/data_archive/nc_entry_gates/`` so that any reviewer can
locate the raw source of any headline number within 24 hours:

- ``inventories/registered_inputs.tsv`` — every registered input file with
  SHA256 (straight from ``data/registry/files.tsv``);
- ``inventories/experiment_artifacts.tsv`` — key manifests/summaries of each
  experiment with computed SHA256;
- ``inventories/release_package.tsv`` — candidate-release package artifacts;
- ``exclusions/excluded_samples.tsv`` — every exclusion decision with reason
  and source reference (curated facts embedded below; the script is the
  provenance);
- ``negative_results/negative_results.md`` — every negative result, with its
  number, source file and status (negative results are never deleted);
- ``reproducibility/environment_lock.txt`` — ``pip freeze`` at build time;
- ``reproducibility/commands.md`` — exact commands to regenerate each
  artifact;
- ``README.md`` — navigation index mapping headline numbers to sources.

Deterministic: re-running rebuilds the archive from current sources.
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

ARCHIVE_DIR = _REPO_ROOT / "results" / "data_archive" / "nc_entry_gates"
RELEASE_DIR = _REPO_ROOT / "results" / "candidates" / "multispecies_v2_candidate_release_v1"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tsv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


# ---------------------------------------------------------------------------
# 1. registered inputs inventory
# ---------------------------------------------------------------------------


def _registered_inputs() -> list[dict[str, str]]:
    registry_files = _REPO_ROOT / "data" / "registry" / "files.tsv"
    rows: list[dict[str, str]] = []
    with registry_files.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for record in reader:
            rows.append(
                {
                    "dataset_accession": record.get("dataset_accession", ""),
                    "repository": record.get("repository", ""),
                    "record_type": record.get("record_type", ""),
                    "file_name": record.get("file_name", ""),
                    "file_category": record.get("file_category", ""),
                    "sha256": record.get("sha256", ""),
                    "status": record.get("status", ""),
                    "retrieved_at": record.get("retrieved_at", ""),
                }
            )
    return rows


# ---------------------------------------------------------------------------
# 2. experiment artifacts inventory
# ---------------------------------------------------------------------------


def _experiment_artifacts() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    experiments = _REPO_ROOT / "results" / "experiments"
    for manifest in sorted(experiments.glob("*/manifest.json")):
        experiment = manifest.parent.name
        rows.append(
            {
                "experiment": experiment,
                "artifact": "manifest.json",
                "sha256": _sha256(manifest),
                "note": "root run manifest (audit_v2_run_manifest audited)",
            }
        )
        literature = manifest.parent / "literature_random_protein"
        for name in ("summary.json", "manifest.json"):
            path = literature / name
            if path.is_file():
                rows.append(
                    {
                        "experiment": experiment,
                        "artifact": f"literature_random_protein/{name}",
                        "sha256": _sha256(path),
                        "note": "literature-track summary/manifest",
                    }
                )
        for report in sorted(manifest.parent.glob("*claim_gate_report.md")):
            rows.append(
                {
                    "experiment": experiment,
                    "artifact": report.name,
                    "sha256": _sha256(report),
                    "note": "claim gate report (wording discipline)",
                }
            )
    return rows


# ---------------------------------------------------------------------------
# 3. release package inventory
# ---------------------------------------------------------------------------


def _release_package() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not (RELEASE_DIR / "manifest.json").is_file():
        return rows
    manifest = json.loads((RELEASE_DIR / "manifest.json").read_text(encoding="utf-8"))
    for entry in manifest.get("produced_artifacts", []):
        artifact = entry.get("artifact")
        path = RELEASE_DIR / artifact
        if path and path.is_file():
            rows.append(
                {
                    "artifact": artifact,
                    "sha256": _sha256(path),
                    "status": entry.get("status", "complete"),
                }
            )
    for name in ("top_k_generation_manifest.json", "power_table.json"):
        path = RELEASE_DIR / name
        if path.is_file():
            rows.append(
                {"artifact": name, "sha256": _sha256(path), "status": "complete"}
            )
    return rows


# ---------------------------------------------------------------------------
# 4. exclusion ledger (curated facts; the script is the provenance)
# ---------------------------------------------------------------------------

_EXCLUSIONS = [
    {
        "id": "EX-001",
        "scope": "tomato candidate scan",
        "count": 72051,
        "reason": "site belongs to the release training panel (model was trained on it); excluded from the candidate table to avoid training contamination",
        "source": "results/candidates/multispecies_v2_candidate_release_v1/top_k_generation_manifest.json",
    },
    {
        "id": "EX-002",
        "scope": "structure feature branch",
        "count": 386009,
        "reason": "no registered AlphaFold structure; structure branch explicitly masked (never mean-imputed) - project missingness rule",
        "source": "results/candidates/multispecies_v2_candidate_release_v1/model_weights/fit_manifest.json",
    },
    {
        "id": "EX-003",
        "scope": "tomato candidate structure coverage",
        "count": 179736,
        "reason": "structure registry release v2 contains zero tomato accessions; every tomato candidate is structure-masked in this release",
        "source": "results/candidates/multispecies_v2_candidate_release_v1/top_k_generation_manifest.json",
    },
    {
        "id": "EX-004",
        "scope": "kiae271 119-peptide table",
        "count": 20,
        "reason": "rows dropped during coordinate+quality verification (8 no position alignment, 3 missing accession, 1 coordinate mismatch, 8 low localization)",
        "source": "results/known_controls/kiae271_differential_recovery_v1.json (parsed block)",
    },
    {
        "id": "EX-005",
        "scope": "blind primary analysis",
        "count": 99,
        "reason": "published kiae271 persulfidation sites excluded from the primary endpoint (published-positive contamination); kept as a pre-registered sensitivity stratum",
        "source": "results/candidates/multispecies_v2_candidate_release_v1/sap_protocol.json",
    },
    {
        "id": "EX-006",
        "scope": "Sul-BertGRU external dataset",
        "count": 1003,
        "reason": "1,003/2,705 external positive records have a non-C centre residue; dataset construction not reproducible as centred-Cys sites",
        "source": "docs/competitor_data_audit_2026-08-10.md",
    },
    {
        "id": "EX-007",
        "scope": "known-control registration",
        "count": 3,
        "reason": "APX1 (gene-name search yields two candidates, neither with Cys at the reported position), POD5 (zero Solanum hits), bZIP68 (no rice scoring path) - not registered as controls",
        "source": "src/plantpersulf/evaluation/known_controls.py (registration comments)",
    },
    {
        "id": "EX-008",
        "scope": "known-control scoring",
        "count": 2,
        "reason": "WRKY71 and ERF.D3 registered as position_shift; not scored against the model (coordinate offset vs published report)",
        "source": "src/plantpersulf/evaluation/known_controls.py",
    },
    {
        "id": "EX-009",
        "scope": "tomato AlphaFold structures (working registry)",
        "count": 76,
        "reason": (
            "76 tomato AlphaFold structures exist in the working registry "
            "(data/registry/alphafold_structures.tsv) but were deliberately "
            "excluded from the frozen candidate-release structure registry "
            "(release v2), and therefore never reach the tomato scoring "
            "path. These 76 structures were fetched to run the "
            "cross-species structural-context analysis and are drawn "
            "entirely from the 88 kiae271 published-positive proteins (76/88 "
            "resolved on AlphaFold DB) — not a proteome-wide or random "
            "sample. Quantified leakage risk if included: of the 179,736 "
            "tomato candidate sites, 580 (0.32%) sit on a protein with one "
            "of these 76 structures, and every one of those 73 proteins is "
            "a kiae271-positive protein; 'has AlphaFold structure' would "
            "therefore have been a near-deterministic proxy for 'protein "
            "carries a published site' among tomato candidates, leaking "
            "protein-level label information through a feature intended to "
            "be biologically orthogonal. Excluding the tomato branch from "
            "the frozen registry avoids this; recomputed 2026-08-14 as the "
            "quantitative justification for the choice already recorded in "
            "EX-003."
        ),
        "source": (
            "docs/cross_species_v2_findings_2026-08-14.md; "
            "data/registry/alphafold_structures.tsv; "
            "results/candidates/multispecies_v2_candidate_release_v1/top_k_candidates.tsv"
        ),
    },
]


# ---------------------------------------------------------------------------
# 5. negative results (curated; negative results are never deleted)
# ---------------------------------------------------------------------------

_NEGATIVE_RESULTS = [
    {
        "result": "kiae271 99-site recovery is statistically indistinguishable from chance",
        "numbers": "48/99 recovered (>50th percentile) = 48.5%, binomial p=0.84 vs 50% (95% CI 38.6-58.3%); mean percentile 49.8, median 47.7",
        "source": "results/known_controls/kiae271_differential_recovery_v1.json; docs/gate2_condition3_kiae271_expansion_2026-08-10.md",
        "status": "accepted negative; expands Gate 2 condition 3 power only, does not flip any gate",
    },
    {
        "result": "cross-species transfer of the Arabidopsis-trained model to tomato is weak",
        "numbers": "48/99 kiae271 sites recovered (48.5%, p=0.84); mean percentile 49.8",
        "source": "results/tomato_local_v1/summary.json (comparison_to_cross_species_baseline)",
        "status": "accepted negative; locked wording (weak 1.33x, no transfer claim)",
    },
    {
        "result": "multispecies v1 joint training shows no leave-one-species-out signal",
        "numbers": "pooled AP 0.0587 vs Arabidopsis-only 0.0095 (6.2x same-species gain); leave-one-species-out enrichment 0.90-1.00x for all four species",
        "source": "results/multispecies_v1/multispecies_summary.json; docs/multispecies_training_v1_2026-08-10.md",
        "status": "accepted negative; the 6x gain is same-species distribution memory, not transferable chemistry",
    },
    {
        "result": "v11 literature track: tomato remains the hardest species",
        "numbers": "tomato per-species test AP 0.0185 +/- 0.0060 (vs arabidopsis 0.4934, rice 0.6470)",
        "source": "results/experiments/multispecies_v2_global_clusters_v11/task9_6_claim_gate_report.md",
        "status": "accepted limitation; lockbox validation required for any external claim",
    },
    {
        "result": "tomato structure coverage in the frozen registry is zero",
        "numbers": "0/179,736 tomato candidate sites have a registered AlphaFold structure",
        "source": "results/candidates/multispecies_v2_candidate_release_v1/top_k_generation_manifest.json",
        "status": "accepted limitation; the structure branch cannot contribute to tomato ranking in this release",
    },
    {
        "result": "blind validation is underpowered at modest effect sizes",
        "numbers": "p0=0.02, K=200: power 0.386 at OR=2, 0.785 at OR=3",
        "source": "results/candidates/multispecies_v2_candidate_release_v1/power_table.json",
        "status": "pre-registered risk; full grid reported, K fixed before unblinding",
    },
    {
        "result": "Sul-BertGRU on the v11 same-benchmark comparison is the weakest model",
        "numbers": "v11 literature track mean macro AP: sul_bertgru 0.0213 vs structure_ranker 0.3863",
        "source": "results/experiments/multispecies_v2_global_clusters_v11/task9_6_claim_gate_report.md",
        "status": "internal same-benchmark comparison; never enters Gate 2 evidence",
    },
    {
        "result": (
            "the previously reported 'persulfidated cysteines are more "
            "buried' structural finding is a model-confidence artefact, not "
            "an independent structural signal"
        ),
        "numbers": (
            "pLDDT at persulfidated vs other cysteines in the same protein: "
            "Arabidopsis +1.8 (p=0.10), rice +3.5 (p=0.001), tomato -7.5 "
            "(p=0.005), Magnaporthe +3.5 (p=0.001); restricting the SASA "
            "comparison to pLDDT>=70 residues drops every species' "
            "significant SASA difference to non-significant (e.g. rice "
            "residue SASA diff -5.22 A2, p=0.001, unfiltered -> -0.57 A2, "
            "p=0.49, confident-only)"
        ),
        "source": "docs/cross_species_v2_findings_2026-08-14.md #3; results/cross_species_conservation/structural_context_v2.json (model_confidence_control)",
        "status": (
            "supersedes the structural-context finding in "
            "docs/phase_z_evidence_audit.md #5.5.2; that section carries an "
            "inline correction pointer and is retained for the record, not "
            "cited as a conclusion"
        ),
    },
    {
        "result": (
            "cross-kingdom conservation at the all-4-species level has zero "
            "statistical power given tomato's shallow ortholog coverage"
        ),
        "numbers": (
            ">=3-species cell: observed 8 vs expected 0.93 under "
            "independence (8.6x), permutation p=1.0e-04 (10,000 "
            "permutations); >=4-species cell: observed 0 vs expected 0.003, "
            "p=1.00 -- no attainable observation could have been "
            "significant, since tomato persulfidates only 8 of the 2,012 "
            "shared subfamilies"
        ),
        "source": "results/cross_species_conservation/conservation_v2.json (conservation_spectrum)",
        "status": (
            "accepted power limitation, not a biological null result; the "
            "3-species cell is reported as the informative one"
        ),
    },
]


def main() -> None:
    created = datetime.now(timezone.utc).isoformat()

    registered = _registered_inputs()
    _tsv(
        ARCHIVE_DIR / "inventories" / "registered_inputs.tsv",
        ["dataset_accession", "repository", "record_type", "file_name", "file_category", "sha256", "status", "retrieved_at"],
        registered,
    )

    experiments = _experiment_artifacts()
    _tsv(
        ARCHIVE_DIR / "inventories" / "experiment_artifacts.tsv",
        ["experiment", "artifact", "sha256", "note"],
        experiments,
    )

    release = _release_package()
    _tsv(
        ARCHIVE_DIR / "inventories" / "release_package.tsv",
        ["artifact", "sha256", "status"],
        release,
    )

    _tsv(
        ARCHIVE_DIR / "exclusions" / "excluded_samples.tsv",
        ["id", "scope", "count", "reason", "source"],
        _EXCLUSIONS,
    )

    negative_path = ARCHIVE_DIR / "negative_results" / "negative_results.md"
    negative_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# 阴性结果登记（negative results）— 永不删除",
        "",
        f"生成时间：{created}（脚本 `scripts/build_nc_data_archive.py`，脚本即 provenance）",
        "",
        "| 结果 | 数字 | 来源 | 状态 |",
        "|---|---|---|---|",
    ]
    for entry in _NEGATIVE_RESULTS:
        lines.append(
            f"| {entry['result']} | {entry['numbers']} | {entry['source']} | {entry['status']} |"
        )
    negative_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    env_lock = ARCHIVE_DIR / "reproducibility" / "environment_lock.txt"
    env_lock.parent.mkdir(parents=True, exist_ok=True)
    lock = subprocess.run(
        [sys.executable, "-m", "pip", "freeze"], capture_output=True, text=True, check=True
    ).stdout
    env_lock.write_text(f"# pip freeze at {created}\n{lock}", encoding="utf-8")

    commands_path = ARCHIVE_DIR / "reproducibility" / "commands.md"
    commands = f"""# 可复现命令清单（Gate 6）

生成时间：{created}

| 产物 | 命令 |
|---|---|
| 冻结模型 bundle + fit_manifest | `python scripts/freeze_structure_ranker_release.py --device cuda:0` |
| 番茄 Top-K 候选 + 对照 + 生成清单 | `python scripts/generate_tomato_topk_candidates.py --device cuda:0` |
| 盲法打分（Gate 3 用） | `python results/candidates/multispecies_v2_candidate_release_v1/inference_script.py --input <tsv> --output <tsv>` |
| 功效表 | `python scripts/compute_release_power_table.py` |
| 预冻结审计 | `python scripts/audit_candidate_release.py --expect-signed` |
| v11 双轨实验 | `python scripts/train_multispecies_v2.py --config configs/experiments/multispecies_v2_global_clusters_v11.yaml --prepare-development --prepare-literature-track --build-comparison-features --run-literature-baselines` |
| 本档案 | `python scripts/build_nc_data_archive.py` |
| 测试 | `python -m pytest tests/unit tests/scientific tests/release -q` |

所有 SHA256 以仓库内 manifest（`data/registry/files.tsv`、实验 manifest、
候选发布包 `manifest.json`）为权威来源；本档案仅是导航索引。
"""
    commands_path.write_text(commands, encoding="utf-8")

    readme = f"""# NC 入场门数据档案（Gate 6 骨架）

生成时间：{created}。任何数字的原始来源应可在 24 小时内从下表定位。

## 头条数字 → 来源

| 数字 | 来源文件 |
|---|---|
| v11 文献轨 structure_ranker mean macro AP 0.3863（内部基准，非外部证据） | `results/experiments/multispecies_v2_global_clusters_v11/literature_random_protein/summary.json` |
| v11 番茄单物种 test AP 0.0185 ± 0.0060 | `.../task9_6_claim_gate_report.md` |
| 冻结模型 bundle SHA256 `ab8a0353…` | `results/candidates/multispecies_v2_candidate_release_v1/model_weights/fit_manifest.json` |
| 训练面板 389,609 位点 SHA256 `30e432fb…` | 同上 |
| 番茄候选 179,736 位点、结构覆盖 0 | `.../top_k_generation_manifest.json` |
| 匹配对照 966 对 | `.../matched_controls.tsv` |
| 功效表（OR=2@K=200 → 0.386） | `.../power_table.json` |
| kiae271 恢复 48/99（48.5%，p=0.84） | `results/known_controls/kiae271_differential_recovery_v1.json` |
| multispecies v1 pooled AP 0.0587 | `results/multispecies_v1/multispecies_summary.json` |
| 注册输入 243 文件 | `inventories/registered_inputs.tsv`（源自 `data/registry/files.tsv`） |

## 目录

- `inventories/registered_inputs.tsv` — 全部注册输入（SHA256 逐文件）
- `inventories/experiment_artifacts.tsv` — 各实验 manifest/summary 哈希
- `inventories/release_package.tsv` — 候选发布包产物哈希
- `exclusions/excluded_samples.tsv` — 全部排除决定（含原因与来源）
- `negative_results/negative_results.md` — 全部阴性结果（永不删除）
- `reproducibility/environment_lock.txt` — pip freeze
- `reproducibility/commands.md` — 每项产物的复现命令

重建：`python scripts/build_nc_data_archive.py`
"""
    (ARCHIVE_DIR / "README.md").write_text(readme, encoding="utf-8")

    print(
        f"archive written: {ARCHIVE_DIR} "
        f"(registered={len(registered)} experiment_artifacts={len(experiments)} "
        f"release_artifacts={len(release)} exclusions={len(_EXCLUSIONS)} "
        f"negative_results={len(_NEGATIVE_RESULTS)})"
    )


if __name__ == "__main__":
    main()
