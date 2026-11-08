"""Tests for the simulation-only cure-discovery engine.

The most important assertions here are negative: the engine must never be able to
call its output a cure, and the excluded criterion must stay excluded.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ihd_solutions.cure import CRITERION_THRESHOLD, CureCriterion
from ihd_solutions.discovery import (
    CAPABILITY_NAMES,
    DISCLAIMER,
    DISCOVERY_STAGES,
    DISCOVERY_TAG,
    EXCLUDED_CRITERION,
    NORMALISATION_LVEF_GAIN_PP,
    SIMULATION_ADMISSIBLE,
    BoundSet,
    Specification,
    best_achievable,
    bounds_for,
    capability_ablation,
    capability_gap_table,
    decision_value_of_information,
    discovery_report,
    excess_risk_remaining,
    feasibility_search,
    feasible_volume_fraction,
    required_capabilities,
    run_discovery_study,
    sample_capabilities,
    stages_frame,
)


@pytest.fixture(scope="module")
def search():
    return feasibility_search(draws=400, seed=5)


# --------------------------------------------------------------------------
# The definitional guard: simulation cannot produce a cure.
# --------------------------------------------------------------------------


def test_randomised_proof_is_excluded_and_named():
    assert CureCriterion.RANDOMISED_PROOF.value == EXCLUDED_CRITERION
    assert EXCLUDED_CRITERION not in SIMULATION_ADMISSIBLE
    assert len(SIMULATION_ADMISSIBLE) == 4
    assert set(SIMULATION_ADMISSIBLE) | {EXCLUDED_CRITERION} == {c.value for c in CureCriterion}


def test_specification_exposes_no_is_cure_attribute():
    spec = best_achievable(BoundSet.CONCEIVABLE)
    assert spec.is_candidate is True
    assert not hasattr(spec, "is_cure")
    assert not hasattr(spec, "cure")


def test_a_passing_specification_still_blocks_on_proof():
    """Even a perfect capability vector leaves the proof criterion unaddressed."""
    spec = best_achievable(BoundSet.CONCEIVABLE)
    assert spec.is_candidate
    assert EXCLUDED_CRITERION not in spec.scores()


def test_disclaimer_says_it_is_not_a_cure():
    assert "never a cure" in DISCLAIMER
    assert "randomised_proof" in DISCLAIMER


def test_report_states_the_exclusion_explicitly(search):
    samples, summary = search
    text = discovery_report(
        summary,
        capability_gap_table(),
        decision_value_of_information(samples),
        capability_ablation(),
        required_capabilities(),
    )
    assert "cannot return a cure" in text.lower()
    assert "definitional trick" in text
    assert DISCLAIMER in text


# --------------------------------------------------------------------------
# Required vector and risk arithmetic.
# --------------------------------------------------------------------------


def test_required_vector_meets_every_admissible_criterion():
    required = required_capabilities()
    spec = Specification(**required)
    assert spec.is_candidate
    for name, value in spec.scores().items():
        assert value == pytest.approx(CRITERION_THRESHOLD, abs=0.02), name


def test_required_lvef_scales_with_the_threshold():
    assert required_capabilities(threshold=0.8)["lvef_gain_pp"] == pytest.approx(20.0)
    assert required_capabilities(threshold=0.5)["lvef_gain_pp"] == pytest.approx(12.5)


def test_tighter_threshold_demands_a_lower_hazard_ratio():
    loose = required_capabilities(threshold=0.5)["hazard_ratio"]
    tight = required_capabilities(threshold=0.9)["hazard_ratio"]
    assert tight < loose


def test_excess_risk_remaining_is_monotone_in_hazard_ratio():
    values = [excess_risk_remaining(hr, 0.1393, 0.0198) for hr in (0.2, 0.4, 0.6, 1.0)]
    assert values == sorted(values)
    assert values[-1] == pytest.approx(1.0)


def test_no_treatment_leaves_all_excess_risk():
    assert excess_risk_remaining(1.0, 0.20, 0.02) == pytest.approx(1.0)


def test_full_normalisation_leaves_no_excess_risk():
    required = required_capabilities(0.20, 0.02, threshold=1.0)
    assert excess_risk_remaining(required["hazard_ratio"], 0.20, 0.02) == pytest.approx(
        0.0, abs=1e-9
    )


def test_risk_arithmetic_rejects_an_inverted_reference():
    with pytest.raises(ValueError, match="reference risk"):
        excess_risk_remaining(0.5, 0.02, 0.20)


def test_required_capabilities_rejects_an_impossible_target():
    with pytest.raises(ValueError):
        required_capabilities(baseline_risk=0.5, reference_risk=0.5)


# --------------------------------------------------------------------------
# Feasibility: the headline result.
# --------------------------------------------------------------------------


def test_demonstrated_bounds_cannot_reach_the_target():
    corner = best_achievable(BoundSet.DEMONSTRATED)
    assert not corner.is_candidate
    assert "myocardial_restoration" in corner.blocking()
    assert "burden_reversal" in corner.blocking()


def test_conceivable_bounds_can_reach_the_target_at_the_corner():
    corner = best_achievable(BoundSet.CONCEIVABLE)
    assert corner.is_candidate
    assert corner.blocking() == ()


def test_demonstrated_feasible_volume_is_exactly_zero():
    volume, windows = feasible_volume_fraction(BoundSet.DEMONSTRATED)
    assert volume == 0.0
    blocked = windows.loc[windows["passing_fraction"].eq(0.0), "capability"]
    assert set(blocked) == {"substrate_removal", "lvef_gain_pp", "hazard_ratio"}


def test_conceivable_feasible_volume_is_positive_but_tiny():
    volume, _ = feasible_volume_fraction(BoundSet.CONCEIVABLE)
    assert 0.0 < volume < 0.01


def test_analytic_volume_equals_the_product_of_windows():
    volume, windows = feasible_volume_fraction(BoundSet.CONCEIVABLE)
    assert volume == pytest.approx(float(windows["passing_fraction"].prod()))


def test_substrate_clearance_is_the_tightest_window():
    """The modelled architecture's ceiling only just clears the bar."""
    _, windows = feasible_volume_fraction(BoundSet.CONCEIVABLE)
    tightest = windows.loc[windows["passing_fraction"].idxmin(), "capability"]
    assert tightest == "substrate_removal"


def test_search_summary_reports_both_boxes(search):
    _, summary = search
    assert set(summary["bound_set"]) == {b.value for b in BoundSet}
    demonstrated = summary.loc[summary["bound_set"].eq("demonstrated")].iloc[0]
    assert demonstrated["p_candidate"] == 0.0
    assert not demonstrated["corner_is_candidate"]


def test_search_samples_stay_inside_their_bounds(search):
    samples, _ = search
    for which in BoundSet:
        subset = samples.loc[samples["bound_set"].eq(which.value)]
        for name, (low, high) in bounds_for(which).as_dict().items():
            assert subset[name].between(low, high).all(), name


def test_search_is_reproducible():
    first, _ = feasibility_search(draws=60, seed=99)
    second, _ = feasibility_search(draws=60, seed=99)
    pd.testing.assert_frame_equal(first, second)


def test_every_sample_row_is_tagged(search):
    samples, summary = search
    assert samples["status"].eq(DISCOVERY_TAG).all()
    assert summary["status"].eq(DISCOVERY_TAG).all()


def test_criteria_met_count_matches_the_scores(search):
    samples, _ = search
    scores = samples[list(SIMULATION_ADMISSIBLE)]
    expected = scores.ge(CRITERION_THRESHOLD).sum(axis=1)
    assert samples["criteria_met"].eq(expected).all()


def test_sample_capabilities_rejects_inverted_bounds():
    bounds = bounds_for(BoundSet.CONCEIVABLE)
    broken = type(bounds)(
        driver_durability=(1.0, 0.0),
        substrate_removal=bounds.substrate_removal,
        lvef_gain_pp=bounds.lvef_gain_pp,
        hazard_ratio=bounds.hazard_ratio,
    )
    with pytest.raises(ValueError, match="upper bound"):
        sample_capabilities(10, broken, 1)


# --------------------------------------------------------------------------
# Gaps, value of information, ablation.
# --------------------------------------------------------------------------


def test_gap_table_ranks_myocardial_restoration_first():
    gaps = capability_gap_table()
    assert gaps.iloc[0]["capability"] == "lvef_gain_pp"
    assert gaps.iloc[0]["gap_factor"] > 5.0


def test_driver_durability_is_already_reachable():
    gaps = capability_gap_table()
    row = gaps.loc[gaps["capability"].eq("driver_durability")].iloc[0]
    assert row["reachable_within_demonstrated_bounds"]
    assert row["gap_factor"] < 1.0


def test_every_gap_carries_a_provenance():
    gaps = capability_gap_table()
    assert gaps["provenance"].str.len().gt(10).all()
    assert gaps["status"].eq(DISCOVERY_TAG).all()


def test_hazard_ratio_gap_is_direction_correct():
    gaps = capability_gap_table()
    row = gaps.loc[gaps["capability"].eq("hazard_ratio")].iloc[0]
    assert row["lower_is_better"]
    assert row["best_demonstrated"] > row["required"]
    assert row["gap_factor"] > 1.0


def test_value_of_information_shares_are_bounded(search):
    samples, _ = search
    voi = decision_value_of_information(samples)
    assert voi["variance_share"].between(0.0, 1.0).all()
    assert set(voi["capability"]) == set(CAPABILITY_NAMES)
    assert voi["variance_share"].is_monotonic_decreasing


def test_value_of_information_handles_a_constant_output():
    frame = pd.DataFrame(
        {
            "driver_durability": [0.1, 0.5, 0.9, 0.4],
            "worst_score": [0.3, 0.3, 0.3, 0.3],
        }
    )
    voi = decision_value_of_information(frame)
    assert voi["variance_share"].eq(0.0).all()


def test_value_of_information_handles_a_constant_input():
    frame = pd.DataFrame(
        {
            "driver_durability": [0.5, 0.5, 0.5, 0.5],
            "worst_score": [0.1, 0.4, 0.7, 0.9],
        }
    )
    voi = decision_value_of_information(frame)
    assert voi["variance_share"].eq(0.0).all()


def test_value_of_information_rejects_a_missing_outcome():
    with pytest.raises(ValueError, match="worst_score"):
        decision_value_of_information(pd.DataFrame({"driver_durability": [0.1, 0.2]}))


def test_binary_outcome_is_supported_when_positives_exist():
    frame = pd.DataFrame(
        {
            "driver_durability": [0.1, 0.2, 0.8, 0.9],
            "is_candidate": [False, False, True, True],
        }
    )
    voi = decision_value_of_information(frame, outcome="is_candidate")
    assert voi.iloc[0]["variance_share"] > 0.5


def test_ablation_identifies_exactly_three_indispensable_capabilities():
    ablation = capability_ablation()
    indispensable = set(ablation.loc[ablation["indispensable"], "reverted_capability"])
    assert indispensable == {"substrate_removal", "lvef_gain_pp", "hazard_ratio"}
    solved = ablation.loc[~ablation["indispensable"], "reverted_capability"]
    assert list(solved) == ["driver_durability"]


def test_ablation_records_what_breaks():
    ablation = capability_ablation()
    row = ablation.loc[ablation["reverted_capability"].eq("lvef_gain_pp")].iloc[0]
    assert row["blocking_after_revert"] == "myocardial_restoration"
    assert row["worst_score"] == pytest.approx(3.9 / NORMALISATION_LVEF_GAIN_PP)


# --------------------------------------------------------------------------
# Stages, outputs, CLI.
# --------------------------------------------------------------------------


def test_stages_are_ordered_and_have_exit_gates():
    frame = stages_frame()
    assert list(frame["stage"]) == list(range(len(DISCOVERY_STAGES)))
    assert frame["exit_gate"].str.len().gt(20).all()
    assert frame["failure_means"].str.len().gt(20).all()
    assert frame["executable_now"].any()
    assert not frame["executable_now"].all()


def test_outputs_are_quarantined_and_tagged(tmp_path):
    target = tmp_path / "disc"
    summary = run_discovery_study(target, draws=60, seed=2)
    assert summary["status"] == "specification_only"
    assert summary["excluded_criterion"] == EXCLUDED_CRITERION
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 10
    assert not any("solutions" in name for name in written)
    for path in (target / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].eq(DISCOVERY_TAG).all(), path.name


def test_summary_volumes_match_the_analytic_computation(tmp_path):
    summary = run_discovery_study(tmp_path / "disc", draws=40, seed=2)
    assert summary["feasible_volume_demonstrated"] == 0.0
    assert summary["feasible_volume_conceivable"] == pytest.approx(
        feasible_volume_fraction(BoundSet.CONCEIVABLE)[0]
    )


def test_scores_are_bounded_for_extreme_inputs():
    spec = Specification(
        driver_durability=1.0, substrate_removal=5.0, lvef_gain_pp=100.0, hazard_ratio=0.01
    )
    for value in spec.scores().values():
        assert 0.0 <= value <= 1.0


def test_cli_discover_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    target = tmp_path / "cli-disc"
    assert main(["discover", "--output-dir", str(target), "--draws", "40"]) == 0
    out = capsys.readouterr().out
    assert DISCLAIMER in out
    assert "Wrote 10 artefacts" in out


def test_cli_discover_rejects_an_inverted_reference(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(
        [
            "discover",
            "--output-dir",
            str(tmp_path / "bad"),
            "--draws",
            "20",
            "--baseline-risk",
            "0.01",
            "--reference-risk",
            "0.20",
        ]
    )
    assert code == 1
    assert "Error:" in capsys.readouterr().err


def test_worst_score_is_the_minimum_not_an_average():
    spec = Specification(
        driver_durability=1.0, substrate_removal=1.0, lvef_gain_pp=0.0, hazard_ratio=0.1
    )
    assert spec.worst_score() == 0.0
    assert not spec.is_candidate
    assert np.isclose(spec.worst_score(), min(spec.scores().values()))
