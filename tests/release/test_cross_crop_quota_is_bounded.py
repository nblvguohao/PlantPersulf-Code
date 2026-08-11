from pathlib import Path

import pytest

from plantpersulf.reporting.candidate_release import ReleaseInputs


def test_cross_crop_quota_cannot_exceed_one_quarter() -> None:
    with pytest.raises(RuntimeError, match="25%"):
        ReleaseInputs(
            Path("model"), Path("cross"), Path("scores"), Path("plan"), 20, 6, 20, 2
        ).validate()


def test_cross_crop_release_rejects_manifest_without_graded_claim_firewall(
    tmp_path: Path,
) -> None:
    model = tmp_path / "model.json"
    transfer = tmp_path / "transfer.json"
    candidates = tmp_path / "candidates.tsv"
    plan = tmp_path / "plan.json"
    model.write_text(
        '{"complete":true,"selected_model":"pu_logistic",'
        '"model_release_created":true}\n',
        encoding="utf-8",
    )
    transfer.write_text(
        '{"complete":true,'
        '"claim_class":"target_label_free_transfer_not_gate2",'
        '"source_gate_admitted":true,"applicability_passed":true,'
        '"wetlab_eligible":true,"target_label_hash":null}\n',
        encoding="utf-8",
    )
    candidates.write_text("policy-marker\n", encoding="utf-8")
    plan.write_text(
        '{"status":"frozen","power_passed":true,'
        '"collaborator_protocol_frozen":true,"total_noncontrol_k":4,'
        '"cross_crop_k":1,"matched_unlabeled_k":2,"process_control_k":1}\n',
        encoding="utf-8",
    )
    inputs = ReleaseInputs(model, transfer, candidates, plan, 4, 1, 2, 1)
    with pytest.raises(RuntimeError, match="graded evidence"):
        inputs.validate()
    transfer.write_text(
        '{"complete":true,'
        '"claim_class":"target_label_free_transfer_not_gate2",'
        '"source_gate_admitted":true,'
        '"cross_source_domain_transfer_supported":true,'
        '"within_species_replication_supported":false,'
        '"within_source_evidence":[{"study_accession":"source-1"}],'
        '"source_domain_evidence":[{"source_domain":"species|source-1"}],'
        '"pure_species_effect_supported":false,'
        '"universal_cross_crop_generalization_supported":false,'
        '"final_generalization_requires_tomato_blind_validation":true,'
        '"applicability_passed":true,"wetlab_eligible":true,'
        '"target_label_hash":null}\n',
        encoding="utf-8",
    )
    inputs.validate()
