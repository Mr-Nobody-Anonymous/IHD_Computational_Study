"""The frontier must stay honest about what simulation cannot decide.

The failure mode these tests guard against is a frontier map that quietly softens its own
conclusions: a gate that always passes, a ranking that reports a winner it cannot defend,
or a "deliverable" that restates the claim instead of bounding it.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from ihd_solutions import simulation_frontier as F

# ------------------------------------------------ LDL and NF-kB are one axis


def test_ldl_lowering_and_nfkb_inhibition_are_the_same_perturbation():
    result = F.ldl_and_nfkb_are_the_same_perturbation()
    assert result["identical_to_machine_precision"] is True
    assert result["max_absolute_difference_over_all_species"] == pytest.approx(0.0, abs=1e-12)


def test_the_identity_is_reported_as_a_defect_not_a_feature():
    result = F.ldl_and_nfkb_are_the_same_perturbation()
    assert result["status"] == "TWO_PARAMETERS_ONE_AXIS"
    assert "parameter against itself" in result["consequence"]


def test_the_effect_table_has_one_row_for_the_shared_axis():
    """Reporting LDL and NF-kB as two rows would imply a distinction the model lacks."""
    table = F.mechanism_effect_table()
    shared = table.loc[table["mechanism"] == "ldl_lowering_or_nfkb_inhibition"]
    assert len(shared) == 1
    assert "one row, not two" in shared["aliases_not_distinguished"].iloc[0]
    assert not any(
        "nfkb" in name for name in table["mechanism"] if name != "ldl_lowering_or_nfkb_inhibition"
    )


# ------------------------------------------------------- the ranking result


def test_four_downstream_mechanisms_cannot_be_ranked():
    result = F.ranking_is_unidentifiable()
    assert result["within_class_ranking_is_not"] is True
    assert result["structurally_blind"] is True


def test_the_separation_is_far_below_any_plausible_measurement():
    result = F.ranking_is_unidentifiable()
    assert result["max_terminal_crp_spread_mg_l"] < 1e-3
    assert result["spread_as_fraction_of_value"] < F.STRUCTURAL_BLINDNESS


def test_separating_the_mechanisms_needs_an_absurd_number_of_patients():
    result = F.ranking_is_unidentifiable()
    assert result["patients_per_arm_to_separate"] > 1e6


def test_the_ranking_result_does_not_claim_nothing_is_visible():
    """Upstream versus downstream IS distinguishable; only the within-class order is not."""
    result = F.ranking_is_unidentifiable()
    assert result["upstream_vs_downstream_is_visible"] is True


def test_the_ranking_result_explains_why_compute_cannot_fix_this_coordinate():
    text = F.ranking_is_unidentifiable()["why_more_compute_does_not_help_on_this_coordinate"]
    assert "not a numerical artefact" in text
    assert "multiplicative" in text


def test_the_ranking_result_is_labelled_as_coordinate_specific():
    """It is a statement about the day-180 terminal value, not about the model."""
    result = F.ranking_is_unidentifiable()
    assert result["status"] == "ranking_unidentifiable_on_terminal_value_only"
    assert "coordinate" in result["why_more_compute_does_not_help_on_this_coordinate"]


def test_the_report_carries_the_correction_so_both_claims_cannot_coexist():
    text = F.frontier_report()
    assert "CORRECTION TO THE PREVIOUS SECTION" in text
    assert "39,507x" in text
    assert "CAN rank the mechanisms through a time course" in text


def test_the_ranking_result_names_the_circularity():
    text = F.ranking_is_unidentifiable()["why_this_blocks_the_plan"]
    assert "circular" in text


def test_the_reduction_from_every_halving_is_about_half():
    table = F.mechanism_effect_table()
    treated = table.loc[table["mechanism"] != "untreated"]
    assert (treated["reduction_fraction"] > 0.49).all()
    assert (treated["reduction_fraction"] < 0.51).all()


def test_the_effect_table_rejects_a_nonsense_step():
    with pytest.raises(ValueError, match="strictly inside"):
        F.mechanism_effect_table(step=1.0)


def test_simulation_deliverables_are_all_caveated():
    frame = F.what_simulation_can_still_do()
    assert len(frame) == 5
    assert frame["caveat"].notna().all()
    assert (frame["caveat"].str.len() > 10).all()


# ------------------------------------------------------- validation targets


def test_no_validation_target_is_present():
    frame = F.validation_targets()
    assert len(frame) == 7
    assert not frame["currently_available"].any()


def test_every_validation_target_states_resolution_and_reason():
    frame = F.validation_targets()
    assert (frame["resolution_required"].str.len() > 10).all()
    assert (frame["why_it_is_the_target"].str.len() > 40).all()


def test_validation_target_dataclass_rejects_a_bare_claim():
    with pytest.raises(ValueError, match="why it is the right one"):
        F.ValidationTarget("m", "q", "r", "d", "", False)
    with pytest.raises(ValueError, match="resolution required"):
        F.ValidationTarget("m", "q", "", "d", "because", False)


def test_validation_arithmetic_states_the_observable_deficit():
    result = F.validation_arithmetic()
    assert result["observables_currently_available"] == 1
    assert result["minimum_observables_to_identify_all_rates"] == 14
    assert result["validation_targets_present"] == 0


# --------------------------------------------------------------- the gate


def test_the_transfer_gate_is_not_satisfied():
    summary = F.gate_summary()
    assert summary["transfer_authorised"] is False
    assert summary["checks_passed"] < summary["checks_total"]


def test_the_gate_names_defining_the_solution_as_the_first_fix():
    assert "define what 'the solution' means" in F.gate_summary()["first_item_to_fix"]


def test_every_gate_check_states_requirement_evidence_and_consequence():
    frame = F.transfer_gate()
    assert len(frame) == 8
    for row in frame.itertuples(index=False):
        assert len(row.requirement) > 10, row.key
        assert len(row.evidence) > 10, row.key
        assert len(row.if_failed) > 10, row.key


def test_a_passed_gate_check_would_flip_the_authorisation(monkeypatch):
    """The gate must be capable of passing; a gate that cannot pass is not a gate."""
    frame = F.transfer_gate().assign(satisfied=True)
    monkeypatch.setattr(F, "transfer_gate", lambda: frame)
    assert F.gate_summary()["transfer_authorised"] is True


def test_gate_check_dataclass_requires_a_consequence_for_failure():
    with pytest.raises(ValueError, match="what happens if it fails"):
        F.GateCheck("k", "requires something", False, "evidence", "")
    with pytest.raises(ValueError, match="must state its requirement"):
        F.GateCheck("k", "", True, "evidence", "n/a")


def test_gate_check_accepts_a_pass_without_a_consequence():
    check = F.GateCheck("k", "requires something", True, "evidence", "")
    assert check.satisfied is True


# ------------------------------------------------- the thing it CAN deliver


def test_the_specification_is_ordered_and_costs_are_stated():
    frame = F.experiment_specification()
    assert list(frame["priority"]) == [1, 2, 3, 4]
    assert frame["cost_class"].str.len().gt(5).all()
    assert "cheapest" in frame["why_this_first"].iloc[0]


def test_the_most_expensive_item_is_the_one_that_settles_efficacy():
    frame = F.experiment_specification()
    last = frame.iloc[-1]
    assert "randomised trial" in last["cost_class"]
    assert "efficacy" in last["resolves"]


# ---------------------------------------------------------------- report


def test_report_renders_and_states_no_transfer():
    text = F.frontier_report()
    assert "THE SIMULATION-FIRST FRONTIER" in text
    assert "TRANSFER AUTHORISED: NO" in text
    assert "CAN THE MODEL NAME A WINNER?" in text
    assert "[FAIL]" in text
    assert "[PASS]" not in text


def test_report_contains_no_nan():
    assert "nan" not in F.frontier_report().lower()


def test_writer_emits_seven_artefacts(tmp_path):
    summary = F.write_frontier_report(tmp_path)
    assert summary["artefacts"] == 7
    assert summary["transfer_authorised"] is False
    for name in (
        "tables/frontier_mechanism_effects.csv",
        "tables/frontier_validation_targets.csv",
        "tables/frontier_transfer_gate.csv",
        "tables/frontier_experiment_spec.csv",
        "tables/frontier_simulation_deliverables.csv",
        "DISCLAIMER.txt",
        "frontier_report.txt",
    ):
        assert (tmp_path / name).exists(), name


def test_written_gate_round_trips(tmp_path):
    F.write_frontier_report(tmp_path)
    frame = pd.read_csv(tmp_path / "tables" / "frontier_transfer_gate.csv")
    assert len(frame) == 8
    assert not frame["satisfied"].any()


def test_disclaimer_says_no_result_is_improved():
    assert "No result here is improved" in F.DISCLAIMER


def test_power_arithmetic_is_finite_and_positive():
    result = F.ranking_is_unidentifiable()
    assert math.isfinite(result["patients_per_arm_to_separate"])
    assert result["patients_per_arm_to_separate"] > 0
    assert np.isfinite(result["between_patient_sd_mg_l"])
