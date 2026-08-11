"""Content-addressed blind candidate-release writer."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from plantpersulf.provenance.hashing import hash_file

RELEASE_COLUMNS = (
    "site_key",
    "protein_accession",
    "arm",
    "selection_order",
    "matched_stratum",
)


@dataclass(frozen=True)
class ReleaseInputs:
    tomato_model_manifest: Path
    cross_crop_model_manifest: Path | None
    candidate_table: Path
    analysis_plan: Path | None
    total_noncontrol_k: int
    cross_crop_k: int
    matched_unlabeled_k: int
    process_control_k: int

    def validate(self) -> None:
        if self.analysis_plan is None:
            raise RuntimeError("frozen analysis plan is required")
        if (
            self.total_noncontrol_k <= 0
            or min(self.cross_crop_k, self.matched_unlabeled_k, self.process_control_k)
            < 0
        ):
            raise RuntimeError("invalid frozen arm sizes")
        if self.cross_crop_k / self.total_noncontrol_k > 0.25:
            raise RuntimeError("cross-crop arm exceeds the 25% cap")
        for path in (
            self.tomato_model_manifest,
            self.candidate_table,
            self.analysis_plan,
        ):
            if path is None or not path.is_file():
                raise RuntimeError(f"release input missing: {path}")
        plan = json.loads(self.analysis_plan.read_text(encoding="utf-8"))
        if (
            plan.get("status") != "frozen"
            or plan.get("power_passed") is not True
            or plan.get("collaborator_protocol_frozen") is not True
        ):
            raise RuntimeError("release arguments differ from the frozen analysis plan")
        if any(
            plan.get(k) != v
            for k, v in {
                "total_noncontrol_k": self.total_noncontrol_k,
                "cross_crop_k": self.cross_crop_k,
                "matched_unlabeled_k": self.matched_unlabeled_k,
                "process_control_k": self.process_control_k,
            }.items()
        ):
            raise RuntimeError("release arguments differ from the frozen analysis plan")
        model = json.loads(self.tomato_model_manifest.read_text(encoding="utf-8"))
        if (
            model.get("complete") is not True
            or model.get("selected_model") not in {"additive_pu", "pu_logistic"}
            or model.get("model_release_created") is not True
        ):
            raise RuntimeError("tomato model manifest has no selected release model")
        if self.cross_crop_k and self.cross_crop_model_manifest is None:
            raise RuntimeError("cross-crop candidates require a transfer manifest")
        if self.cross_crop_model_manifest is not None:
            if not self.cross_crop_model_manifest.is_file():
                raise RuntimeError("cross-crop transfer manifest is missing")
            transfer = json.loads(
                self.cross_crop_model_manifest.read_text(encoding="utf-8")
            )
            if (
                transfer.get("complete") is not True
                or transfer.get("claim_class") != "target_label_free_transfer_not_gate2"
                or transfer.get("source_gate_admitted") is not True
                or transfer.get("applicability_passed") is not True
                or transfer.get("wetlab_eligible") is not True
                or transfer.get("target_label_hash") is not None
            ):
                raise RuntimeError("cross-crop transfer manifest is not admitted")
            if (
                transfer.get("cross_source_domain_transfer_supported") is not True
                or transfer.get("pure_species_effect_supported") is not False
                or transfer.get("universal_cross_crop_generalization_supported")
                is not False
                or transfer.get("final_generalization_requires_tomato_blind_validation")
                is not True
                or not isinstance(
                    transfer.get("within_species_replication_supported"), bool
                )
                or not transfer.get("within_source_evidence")
                or not transfer.get("source_domain_evidence")
            ):
                raise RuntimeError(
                    "cross-crop transfer manifest lacks graded evidence safeguards"
                )


@dataclass(frozen=True)
class CandidateRelease:
    release_id: str
    output_dir: Path


def release_id_from_inputs(inputs: ReleaseInputs) -> str:
    inputs.validate()
    assert inputs.analysis_plan is not None
    material = "|".join(
        (
            hash_file(inputs.tomato_model_manifest, "sha256"),
            hash_file(inputs.cross_crop_model_manifest, "sha256")
            if inputs.cross_crop_model_manifest is not None
            else "",
            hash_file(inputs.candidate_table, "sha256"),
            hash_file(inputs.analysis_plan, "sha256"),
            str(inputs.total_noncontrol_k),
            str(inputs.cross_crop_k),
            str(inputs.matched_unlabeled_k),
            str(inputs.process_control_k),
        )
    )
    return "candidate-release-" + hashlib.sha256(material.encode()).hexdigest()[:16]


def blind_id(release_id: str, site_key: str) -> str:
    return (
        "BLIND-"
        + hashlib.sha256(f"{release_id}|{site_key}".encode()).hexdigest()[:12].upper()
    )


def build_candidate_release(
    inputs: ReleaseInputs, output_dir: Path
) -> CandidateRelease:
    if output_dir.exists():
        raise FileExistsError(output_dir)
    release_id = release_id_from_inputs(inputs)
    assert inputs.analysis_plan is not None
    with inputs.candidate_table.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != RELEASE_COLUMNS:
            raise RuntimeError("candidate table schema mismatch")
        rows = [dict(row) for row in reader]
    if {row["arm"] for row in rows} - {
        "T",
        "X",
        "matched_unlabeled",
        "process_control",
    }:
        raise RuntimeError("candidate table contains an unknown arm")
    by_arm = {
        arm: [row for row in rows if row["arm"] == arm]
        for arm in ("T", "X", "matched_unlabeled", "process_control")
    }
    if (
        len(by_arm["T"]) + len(by_arm["X"]) != inputs.total_noncontrol_k
        or len(by_arm["X"]) != inputs.cross_crop_k
        or len(by_arm["matched_unlabeled"]) != inputs.matched_unlabeled_k
        or len(by_arm["process_control"]) != inputs.process_control_k
    ):
        raise RuntimeError("candidate table arm counts differ from frozen K")
    noncontrol_proteins = [
        row["protein_accession"] for row in (*by_arm["T"], *by_arm["X"])
    ]
    if len(noncontrol_proteins) != len(set(noncontrol_proteins)):
        raise RuntimeError("more than one non-control site selected per protein")
    if len({row["site_key"] for row in rows}) != len(rows):
        raise RuntimeError("candidate table contains a duplicate site key")
    output_dir.mkdir(parents=True)
    private_rows = [
        {**row, "blind_id": blind_id(release_id, row["site_key"])}
        for row in sorted(
            rows, key=lambda row: (int(row["selection_order"]), row["site_key"])
        )
    ]
    with (output_dir / "private_key.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("blind_id", *RELEASE_COLUMNS), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(private_rows)
    with (output_dir / "blind_table.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("blind_id", "site_key", "protein_accession", "matched_stratum"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(
            {key: row[key] for key in writer.fieldnames} for row in private_rows
        )
    (output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "complete": True,
                "release_id": release_id,
                "candidate_table_sha256": hash_file(inputs.candidate_table, "sha256"),
                "analysis_plan_sha256": hash_file(inputs.analysis_plan, "sha256"),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (output_dir / "SHA256SUMS").write_text(
        "\n".join(
            f"{hash_file(output_dir / name, 'sha256')}  {name}"
            for name in ("private_key.tsv", "blind_table.tsv", "manifest.json")
        )
        + "\n",
        encoding="utf-8",
    )
    return CandidateRelease(release_id, output_dir)
