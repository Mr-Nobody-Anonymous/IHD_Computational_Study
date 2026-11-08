"""Plaque regression: calibration recovery, monotonicity, and no eradication claim."""

from __future__ import annotations

import numpy as np
import pytest

from ihd_solutions.catalog import validate_combination
from ihd_solutions.evidence import EndpointClass
from ihd_solutions.regression import (
    CALIBRATION_HORIZON_YEARS,
    DEFAULT_CALIBRATION,
    NONCALCIFIED_POINTS,
    PAV_POINTS,
    calibration_frame,
    calibration_summary,
    fit_calibration,
    imaging_anchor_frame,
    regimen_regression,
    regression_grid,
    simulate_regression,
    years_to_noncalcified_depletion,
)


def test_calibration_fit_is_reported_not_hidden():
    calibration = fit_calibration()
    assert calibration.pav_slope > 0  # higher LDL-C means faster growth
    assert 0.5 < calibration.pav_r_squared <= 1.0
    summary = calibration_summary(calibration)
    assert summary["pav_r_squared"] == pytest.approx(calibration.pav_r_squared)
    assert summary["caveats"]
    frame = calibration_frame(calibration)
    assert len(frame) == len(PAV_POINTS)
    assert "residual_pp_per_year" in frame.columns


def test_fitted_regression_threshold_matches_the_imaging_literature():
    # Serial-imaging trials place the switch from progression to regression around
    # 70-80 mg/dL. The fit must land there without being told to.
    assert 60.0 < DEFAULT_CALIBRATION.threshold_ldl_mg_dl < 90.0


@pytest.mark.parametrize("ldl,expected_loss,years", NONCALCIFIED_POINTS)
def test_model_reproduces_its_noncalcified_calibration(ldl, expected_loss, years):
    result = simulate_regression(ldl, years=years, step_years=years / 4.0)
    assert result.noncalcified_depletion_fraction == pytest.approx(expected_loss, rel=1e-6)


def test_lower_achieved_ldl_gives_more_regression():
    high = simulate_regression(90.0, years=2.0)
    low = simulate_regression(35.0, years=2.0)
    assert low.noncalcified_depletion_fraction > high.noncalcified_depletion_fraction
    assert low.percent_atheroma_volume[-1] < high.percent_atheroma_volume[-1]


def test_above_threshold_plaque_still_grows():
    result = simulate_regression(130.0, years=2.0)
    assert result.percent_atheroma_volume[-1] > result.percent_atheroma_volume[0]
    assert result.noncalcified_depletion_fraction == pytest.approx(0.0)
    assert any("above the fitted regression threshold" in note for note in result.notes)


def test_calcified_burden_rises_and_total_burden_is_never_eradicated():
    result = simulate_regression(20.0, years=30.0)
    assert result.calcified_increased
    assert not result.eradicated
    assert float(result.total_mm3[-1]) > 0.0
    # Even with the lipid pool essentially gone, dense calcium remains.
    assert result.noncalcified_depletion_fraction > 0.9
    assert result.total_reduction_fraction < 0.7


def test_zero_calcification_is_the_most_optimistic_allowed_case():
    optimistic = simulate_regression(20.0, years=30.0, calcification_fraction=0.0)
    default = simulate_regression(20.0, years=30.0)
    assert optimistic.total_reduction_fraction > default.total_reduction_fraction
    assert not optimistic.eradicated


def test_extrapolation_beyond_the_calibration_horizon_is_flagged():
    within = simulate_regression(40.0, years=CALIBRATION_HORIZON_YEARS)
    beyond = simulate_regression(40.0, years=10.0)
    assert not within.extrapolated_beyond_evidence
    assert beyond.extrapolated_beyond_evidence
    assert any("extrapolation" in note for note in beyond.notes)


def test_time_to_depletion_is_none_when_there_is_no_net_loss():
    assert years_to_noncalcified_depletion(130.0) is None
    years = years_to_noncalcified_depletion(40.0)
    assert years is not None and years > CALIBRATION_HORIZON_YEARS
    assert years_to_noncalcified_depletion(20.0) < years


def test_time_to_depletion_rejects_impossible_targets():
    for target in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="target_fraction"):
            years_to_noncalcified_depletion(40.0, target_fraction=target)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"achieved_ldl_mg_dl": 0.0},
        {"achieved_ldl_mg_dl": float("nan")},
        {"achieved_ldl_mg_dl": 40.0, "years": 0.0},
        {"achieved_ldl_mg_dl": 40.0, "step_years": 0.0},
        {"achieved_ldl_mg_dl": 40.0, "initial_noncalcified_mm3": 0.0},
        {"achieved_ldl_mg_dl": 40.0, "calcification_fraction": 1.0},
        {"achieved_ldl_mg_dl": 40.0, "calcification_fraction": -0.1},
    ],
)
def test_invalid_scenarios_are_rejected(kwargs):
    with pytest.raises(ValueError):
        simulate_regression(**kwargs)


def test_trajectory_is_on_a_regular_grid_and_starts_at_baseline():
    result = simulate_regression(45.0, years=2.0, step_years=0.5)
    assert result.years[0] == 0.0
    assert result.years[-1] == pytest.approx(2.0)
    assert np.all(np.diff(result.years) > 0)
    assert result.noncalcified_mm3[0] == pytest.approx(120.0)
    assert result.calcified_mm3[0] == pytest.approx(60.0)


def test_mass_balance_of_the_calcification_assumption():
    fraction = 0.4
    result = simulate_regression(30.0, years=3.0, calcification_fraction=fraction)
    lost = 120.0 - float(result.noncalcified_mm3[-1])
    gained = float(result.calcified_mm3[-1]) - 60.0
    assert gained == pytest.approx(fraction * lost)


def test_regimen_regression_uses_the_achieved_ldl(catalog, populations, with_background):
    population = populations["stable_ccs"]
    plain = regimen_regression(with_background(), catalog, population)
    intensified = regimen_regression(with_background("evolocumab"), catalog, population)
    assert intensified.achieved_ldl_mg_dl < plain.achieved_ldl_mg_dl
    assert intensified.noncalcified_depletion_fraction > plain.noncalcified_depletion_fraction


def test_imaging_anchors_are_reported_with_zero_mace_credit(anchors):
    frame = imaging_anchor_frame(anchors)
    assert {"ASTEROID", "GLAGOV", "FLOWPROMOTE"} <= set(frame["anchor_id"])
    assert (frame["mace_credit"] == 0).all()
    for anchor_id in frame["anchor_id"]:
        assert anchors[anchor_id].endpoint_class is EndpointClass.PLAQUE_IMAGING


def test_imaging_evidence_never_moves_the_event_curve(catalog, populations, anchors, days):
    from ihd_solutions.risk import simulate_events

    population = populations["stable_ccs"]
    regimen = validate_combination(catalog, ("statin_high_intensity",))
    simulation = simulate_events(regimen, catalog, population, anchors, days)
    # The statin still earns its CTT per-unit credit; ASTEROID must not shadow it.
    assert simulation.events_avoided_per_1000 > 0
    assert all(
        anchors[decision.anchor_id].endpoint_class is not EndpointClass.PLAQUE_IMAGING
        for decision in simulation.decisions
        if decision.credited
    )


def test_regression_grid_ranks_by_structural_effect_and_claims_no_events(catalog, populations):
    frame = regression_grid(
        catalog,
        populations["stable_ccs"],
        ("statin_high_intensity", "aspirin_secondary_prevention"),
        ("evolocumab", "ezetimibe", "obicetrapib", "verve102_base_editing"),
    )
    assert not frame.empty
    assert (frame["mace_credit_from_imaging"] == 0).all()
    assert not frame["eradicated"].any()
    assert frame["noncalcified_depletion_pct"].is_monotonic_decreasing
    assert np.isfinite(frame["achieved_ldl_mg_dl"]).all()
