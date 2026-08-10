"""RED (Phase S0, Task 4): explicit projection boundary for structure ablation
arms. Coverage-only must use only the mask column; it cannot consume contact
or pLDDT values regardless of what the raw structure array contains."""

from __future__ import annotations

import numpy as np
import pytest


class TestProjectStructureInputs:
    def test_coverage_only_cannot_consume_contact_or_plddt(self) -> None:
        from plantpersulf.models.structure_ranker import project_structure_inputs

        first = project_structure_inputs(
            structure=np.array([[11.0, 22.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_coverage_only",
        )
        second = project_structure_inputs(
            structure=np.array([[111.0, 222.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_coverage_only",
        )
        np.testing.assert_array_equal(first.values, second.values)
        np.testing.assert_array_equal(first.values, np.array([[1.0, 0.0]]))

    def test_sequence_only_exposes_no_structure_tensor(self) -> None:
        from plantpersulf.models.structure_ranker import project_structure_inputs

        result = project_structure_inputs(
            structure=np.array([[11.0, 22.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_only",
        )
        np.testing.assert_array_equal(result.values, np.array([[0.0, 0.0]]))
        np.testing.assert_array_equal(result.mask, np.array([[0.0]]))

    def test_sequence_contact_changes_with_contact_but_not_plddt(self) -> None:
        from plantpersulf.models.structure_ranker import project_structure_inputs

        base = project_structure_inputs(
            structure=np.array([[10.0, 50.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_contact",
        )
        different_contact = project_structure_inputs(
            structure=np.array([[20.0, 50.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_contact",
        )
        same_contact = project_structure_inputs(
            structure=np.array([[10.0, 99.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_contact",
        )
        assert base.values[0, 0] != different_contact.values[0, 0]
        np.testing.assert_array_equal(base.values, same_contact.values)

    def test_sequence_plddt_changes_with_plddt_but_not_contact(self) -> None:
        from plantpersulf.models.structure_ranker import project_structure_inputs

        base = project_structure_inputs(
            structure=np.array([[10.0, 50.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_plddt",
        )
        different_plddt = project_structure_inputs(
            structure=np.array([[10.0, 70.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_plddt",
        )
        same_plddt = project_structure_inputs(
            structure=np.array([[99.0, 50.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_plddt",
        )
        assert base.values[0, 1] != different_plddt.values[0, 1]
        np.testing.assert_array_equal(base.values, same_plddt.values)

    def test_sequence_contact_plddt_changes_with_both(self) -> None:
        from plantpersulf.models.structure_ranker import project_structure_inputs

        a = project_structure_inputs(
            structure=np.array([[10.0, 50.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_contact_plddt",
        )
        b = project_structure_inputs(
            structure=np.array([[20.0, 70.0]]),
            mask=np.array([[1.0]]),
            arm="sequence_contact_plddt",
        )
        assert a.values.shape[1] == 2
        assert b.values.shape[1] == 2
        assert not np.allclose(a.values, b.values)

    def test_missing_structure_produces_zero_projection_and_mask_0(self) -> None:
        from plantpersulf.models.structure_ranker import project_structure_inputs

        for arm in (
            "sequence_contact",
            "sequence_plddt",
            "sequence_contact_plddt",
            "sequence_coverage_only",
        ):
            result = project_structure_inputs(
                structure=np.array([[99.0, 99.0]]),
                mask=np.array([[0.0]]),
                arm=arm,
            )
            np.testing.assert_array_equal(result.values, np.zeros_like(result.values))
            np.testing.assert_array_equal(result.mask, np.array([[0.0]]))

    def test_unknown_arm_is_rejected(self) -> None:
        from plantpersulf.models.structure_ranker import project_structure_inputs

        with pytest.raises(ValueError, match="unknown arm"):
            project_structure_inputs(
                structure=np.array([[1.0, 2.0]]),
                mask=np.array([[1.0]]),
                arm="nonexistent",
            )

    def test_coverage_only_mask_is_preserved(self) -> None:
        from plantpersulf.models.structure_ranker import project_structure_inputs

        result = project_structure_inputs(
            structure=np.array([[11.0, 22.0], [33.0, 44.0]]),
            mask=np.array([[1.0], [0.0]]),
            arm="sequence_coverage_only",
        )
        np.testing.assert_array_equal(result.mask, np.array([[1.0], [0.0]]))
        assert result.values[0, 0] == 1.0
        assert result.values[1, 0] == 0.0
        assert result.values[0, 1] == 0.0
        assert result.values[1, 1] == 0.0
