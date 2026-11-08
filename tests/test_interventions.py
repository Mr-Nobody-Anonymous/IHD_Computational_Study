"""Each gate item must map to a real instrument, and the verdict must stay honest.

The failure mode: a register that cites a paper for every item and implies that citing it
closes the item. These tests pin the outcome counts, force the register to agree with the
live gate, and keep the correction to the repository's own earlier claim in place.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ihd_solutions import interventions as I
from ihd_solutions import simulation_frontier as F

# ---------------------------------------------------------------- the register


def test_every_gate_item_has_a_response():
    assert len(I.GATE_RESPONSES) == 8
    keys = {response.gate_key for response in I.GATE_RESPONSES}
    assert keys == set(F.transfer_gate()["key"])


def test_every_response_cites_a_primary_source():
    for response in I.GATE_RESPONSES:
        assert response.citations, response.gate_key
        assert len(response.citations[0]) > 20, response.gate_key


def test_every_response_states_remaining_work():
    for response in I.GATE_RESPONSES:
        assert len(response.remaining_work) > 40, response.gate_key


def test_no_item_is_reported_resolved():
    """If an item ever passes, this test should be the thing that changes."""
    assert all(response.outcome != "RESOLVED" for response in I.GATE_RESPONSES)


def test_register_outcomes_agree_with_the_live_gate():
    frame = I.gate_matrix()
    assert frame["outcome_agrees_with_gate"].all()
    assert not frame["gate_passes_now"].any()


def test_register_rejects_an_unknown_outcome():
    with pytest.raises(ValueError, match="unknown outcome"):
        I.GateResponse("k", "instrument", ("a citation",), "verdict", "SOLVED", "does", "work")


def test_register_rejects_an_uncited_response():
    with pytest.raises(ValueError, match="must cite"):
        I.GateResponse("k", "instrument", (), "verdict", "STILL_FAILED", "does", "work")


def test_register_rejects_a_response_with_no_remaining_work():
    with pytest.raises(ValueError, match="must state what work remains"):
        I.GateResponse("k", "instrument", ("a citation",), "verdict", "STILL_FAILED", "does", "")


def test_register_refuses_to_answer_a_gate_item_that_does_not_exist(monkeypatch):
    bogus = I.GateResponse(
        "not_a_gate_item", "estimand", ("x" * 30,), "v", "STILL_FAILED", "d", "w"
    )
    monkeypatch.setattr(I, "GATE_RESPONSES", (bogus,))
    with pytest.raises(ValueError, match="no matching gate check"):
        I.gate_matrix()


def test_summary_counts_are_self_consistent():
    summary = I.instrument_summary()
    assert summary["gate_items"] == 8
    assert summary["resolved"] + summary["applicable_but_not_done"] + summary["still_failed"] == 8
    assert summary["verdicts_consistent_with_live_gate"] is True


def test_instruments_cover_every_response():
    for response in I.GATE_RESPONSES:
        assert response.instrument in I.INSTRUMENTS


def test_each_instrument_names_its_definition_and_source():
    for key, (definition, citations, verdict) in I.INSTRUMENTS.items():
        assert len(definition) > 60, key
        assert citations, key
        assert len(verdict) > 40, key


# ------------------------------------------------------------- the correction


def test_the_correction_reverses_the_repository_s_own_claim():
    correction = I.terminal_value_was_the_wrong_coordinate()
    assert "cannot rank" in correction["superseded_claim"]
    assert "CAN" in correction["correct_claim"] or "ranks them" in correction["correct_claim"]


def test_the_time_course_separation_dwarfs_the_terminal_reading():
    power = I.discrimination_power()
    assert power["separation_mg_l"] > 1.0
    assert power["contrast_with_terminal_value"]["fold_improvement"] > 1000


def test_the_time_course_design_is_resolvable_in_a_small_trial():
    power = I.discrimination_power()
    assert power["resolvable"] is True
    assert power["patients_per_arm"] < I.RESOLVABLE_TRIAL_SIZE


def test_the_correction_keeps_what_is_still_true():
    text = I.terminal_value_was_the_wrong_coordinate()["what_is_still_true"]
    assert "no time course" in text
    assert "terminal CRP" in text


def test_the_correction_names_the_three_errors_of_the_original_reading():
    items = I.terminal_value_was_the_wrong_coordinate()["what_was_wrong_with_the_original_reading"]
    assert len(items) == 3
    assert any("not converged" in item for item in items)
    assert any("design choice" in item for item in items)
    assert any("unobserved" in item for item in items)


def test_the_chosen_observable_is_a_readout_not_a_mechanism_node():
    """Ranking mechanisms by a node they are themselves perturbing is circular."""
    design = I.t_optimal_design()
    readouts = design.loc[design["is_readout"]]
    assert set(readouts["observable"]) == set(I.READOUT_SPECIES)
    chosen = I.discrimination_power()["observable"]
    assert chosen in I.READOUT_SPECIES


def test_mechanism_nodes_are_flagged_as_circular_evidence():
    design = I.t_optimal_design()
    nodes = design.loc[~design["is_readout"]]
    assert len(nodes) == 4
    # IL-6 is the one node that carries only part of the information, because two of the
    # four perturbations act upstream of it; it decays anyway. The three strict nodes are
    # exact restatements of the perturbation and are flagged.
    assert (nodes["circular_warning"].str.contains("circular")).all()


def test_the_returned_observable_is_downstream_of_every_knob():
    """CRP therefore separates the mechanisms even though it is not itself perturbed."""
    design = I.t_optimal_design()
    crp = design.loc[design["observable"] == "crp"].iloc[0]
    assert crp["is_readout"] is True or bool(crp["is_readout"])
    assert crp["terminal_ratio_to_best"] < 0.001
    assert crp["best_day"] < 30


def test_strict_mechanism_nodes_separate_at_day_180_only_because_they_are_perturbed():
    """Not evidence: the node IS the turned-down quantity. Reported, and flagged."""
    design = I.t_optimal_design()
    strict = design.loc[design["observable"].isin(["il1b", "casp1", "nlrp3"])]
    assert (strict["terminal_ratio_to_best"] == 1.0).all()
    assert (strict["circular_warning"] != "").all()


# ------------------------------------------------------------------- estimand


def test_the_estimand_is_not_fully_specified():
    result = I.estimand_specification()
    assert result["fully_specified"] is False
    assert result["attributes_stated"] == 2
    assert result["attributes_total"] == 5


def test_the_blocking_attribute_is_intercurrent_events():
    result = I.estimand_specification()
    assert result["blocking_attribute"] == "intercurrent_events"
    assert result["draft"]["intercurrent_events"].startswith("NOT SPECIFIED")


def test_the_estimand_names_the_comparator_as_guideline_therapy():
    assert "guideline therapy" in I.estimand_specification()["draft"]["treatment"]


def test_gate_item_one_still_fails():
    assert "STILL_FAILED" in I.estimand_specification()["gate_item_1_status"]


# ------------------------------------------------------------------ reporting


def test_report_renders_every_section_and_item():
    text = I.intervention_report()
    assert "WHAT THE LITERATURE SAYS ABOUT EACH GATE ITEM" in text
    assert "THE CORRECTION" in text
    assert "THE T-OPTIMAL DESIGN" in text
    assert "STILL FAILED                  8" in text
    for response in I.GATE_RESPONSES:
        assert response.gate_key in text


def test_report_contains_no_nan():
    assert "nan" not in I.intervention_report().lower()


def test_report_does_not_claim_the_gate_passes():
    text = I.intervention_report()
    assert "TRANSFER AUTHORISED: YES" not in text


def test_writer_emits_six_artefacts(tmp_path):
    summary = I.write_intervention_report(str(tmp_path))
    assert summary["artefacts"] == 6
    assert summary["still_failed"] == 8
    for name in (
        "tables/interventions_gate_matrix.csv",
        "tables/interventions_t_optimal_design.csv",
        "tables/interventions_estimand_draft.csv",
        "tables/interventions_instruments.csv",
        "DISCLAIMER.txt",
        "intervention_report.txt",
    ):
        assert (tmp_path / name).exists(), name


def test_written_matrix_round_trips(tmp_path):
    I.write_intervention_report(str(tmp_path))
    frame = pd.read_csv(tmp_path / "tables" / "interventions_gate_matrix.csv")
    assert len(frame) == 8
    assert (frame["outcome"] == "STILL_FAILED").all()
    assert not frame["gate_passes_now"].any()


def test_disclaimer_separates_methods_from_data():
    assert "The methods." in I.DISCLAIMER
    assert "This study's data" in I.DISCLAIMER
