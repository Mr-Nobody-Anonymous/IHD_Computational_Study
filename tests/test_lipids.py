"""Biomarker trajectories and multiplicative combination logic."""

from __future__ import annotations

import numpy as np
import pytest

from ihd_solutions.lipids import (
    APOB_FLOOR_MG_DL,
    HSCRP_FLOOR_MG_L,
    LDL_FLOOR_MG_DL,
    LPA_FLOOR_NMOL_L,
    simulate_biomarkers,
    trajectory_frame,
)
from ihd_solutions.pharmacology import Adherence
from ihd_solutions.units import MG_DL_PER_MMOL_L


def test_single_agent_reproduces_its_reported_reduction(catalog, population, make_regimen):
    regimen = make_regimen("evolocumab")
    # Day 364 is a dosing day (26 x 14), so the sawtooth is at its peak.
    trajectory = simulate_biomarkers(regimen, catalog, population, np.linspace(0, 364, 365))
    expected = population.ldl_c_mg_dl * (1.0 - catalog["evolocumab"].ldl_pct)
    assert trajectory.ldl_mg_dl[-1] == pytest.approx(expected, rel=1e-6)
    assert trajectory.ldl_reduction_fraction[-1] == pytest.approx(0.59, rel=1e-6)


def test_combination_is_multiplicative_not_additive(catalog, population, make_regimen):
    regimen = make_regimen("statin_high_intensity", "ezetimibe")
    trajectory = simulate_biomarkers(regimen, catalog, population, np.linspace(0, 365, 366))
    expected = population.ldl_c_mg_dl * (1 - 0.50) * (1 - 0.18)
    assert trajectory.ldl_mg_dl[-1] == pytest.approx(expected, rel=1e-6)
    additive = population.ldl_c_mg_dl * (1 - 0.50 - 0.18)
    assert trajectory.ldl_mg_dl[-1] > additive


def test_stacking_never_goes_below_the_physiological_floor(catalog, population, make_regimen):
    regimen = make_regimen(
        "statin_high_intensity", "ezetimibe", "bempedoic_acid", "evolocumab", "diet_ldl_lowering"
    )
    trajectory = simulate_biomarkers(regimen, catalog, population, np.linspace(0, 730, 200))
    assert trajectory.ldl_mg_dl.min() >= LDL_FLOOR_MG_DL
    assert trajectory.apob_mg_dl.min() >= APOB_FLOOR_MG_DL


def test_lpa_and_hscrp_floors_are_respected(catalog, populations, make_regimen):
    population = populations["ascvd_high_lpa"]
    regimen = make_regimen("olpasiran", "ziltivekimab")
    trajectory = simulate_biomarkers(regimen, catalog, population, np.linspace(0, 730, 200))
    assert trajectory.lpa_nmol_l.min() >= LPA_FLOOR_NMOL_L
    assert trajectory.hscrp_mg_l.min() >= HSCRP_FLOOR_MG_L


def test_no_ldl_therapy_leaves_ldl_untouched(catalog, population, make_regimen):
    regimen = make_regimen("colchicine")
    trajectory = simulate_biomarkers(regimen, catalog, population, np.linspace(0, 365, 100))
    assert np.allclose(trajectory.ldl_mg_dl, population.ldl_c_mg_dl)
    years = 365.0 / 365.25
    assert trajectory.total_ldl_exposure_mmol_years == pytest.approx(
        population.ldl_c_mg_dl / MG_DL_PER_MMOL_L * years
    )


def test_cumulative_exposure_is_monotone_and_lower_with_therapy(
    catalog, population, make_regimen, days
):
    untreated = simulate_biomarkers(make_regimen("colchicine"), catalog, population, days)
    treated = simulate_biomarkers(
        make_regimen("statin_high_intensity", "evolocumab"), catalog, population, days
    )
    assert np.all(np.diff(treated.cumulative_ldl_mmol_years) >= 0)
    assert treated.total_ldl_exposure_mmol_years < untreated.total_ldl_exposure_mmol_years
    assert treated.time_averaged_ldl_mg_dl < untreated.time_averaged_ldl_mg_dl


def test_time_below_goal_increases_with_potency(catalog, population, make_regimen, days):
    moderate = simulate_biomarkers(
        make_regimen("statin_moderate_intensity"), catalog, population, days
    )
    intensive = simulate_biomarkers(
        make_regimen("statin_high_intensity"), catalog, population, days
    )
    assert intensive.fraction_time_below(55.0) >= moderate.fraction_time_below(55.0)
    assert 0.0 <= intensive.fraction_time_below(55.0) <= 1.0


def test_adherence_reduces_the_achieved_reduction(catalog, population, make_regimen, days):
    regimen = make_regimen("statin_high_intensity")
    full = simulate_biomarkers(regimen, catalog, population, days)
    partial = simulate_biomarkers(
        regimen,
        catalog,
        population,
        days,
        {"statin_high_intensity": Adherence(taking_fraction=0.5)},
    )
    assert partial.time_averaged_ldl_mg_dl > full.time_averaged_ldl_mg_dl


def test_blood_pressure_reduction_accumulates_across_agents(catalog, population, make_regimen):
    regimen = make_regimen("bp_intensification", "exercise_cardiac_rehab")
    trajectory = simulate_biomarkers(regimen, catalog, population, np.linspace(0, 730, 300))
    expected_drop = 10.0 + 4.0
    assert trajectory.systolic_bp_mmhg[-1] == pytest.approx(
        population.systolic_bp_mmhg - expected_drop
    )


def test_two_sglt2_agents_cannot_be_combined(catalog, populations, make_regimen):
    with pytest.raises(ValueError, match="exclusive group"):
        make_regimen("empagliflozin", "sglt2_inhibitor_hf")


def test_pressure_cannot_fall_below_the_safety_floor(catalog, populations, make_regimen):
    from dataclasses import replace

    population = replace(populations["ischemic_cm_heart_failure"], systolic_bp_mmhg=95.0)
    regimen = make_regimen("bp_intensification", "semaglutide_2_4mg", "smoking_cessation")
    trajectory = simulate_biomarkers(regimen, catalog, population, np.linspace(0, 730, 300))
    assert trajectory.systolic_bp_mmhg.min() >= 80.0


def test_components_expose_each_agents_contribution(catalog, population, make_regimen, days):
    regimen = make_regimen("statin_high_intensity", "ezetimibe")
    trajectory = simulate_biomarkers(regimen, catalog, population, days)
    assert set(trajectory.components) == {"statin_high_intensity", "ezetimibe"}
    assert trajectory.components["statin_high_intensity"][-1] == pytest.approx(0.50)


def test_trajectory_frame_columns_and_length(catalog, population, make_regimen, days):
    trajectory = simulate_biomarkers(make_regimen("ezetimibe"), catalog, population, days)
    frame = trajectory_frame(trajectory, "ezetimibe", population.population_id)
    assert len(frame) == len(days)
    assert {"ldl_c_mg_dl", "ldl_c_mmol_l", "cumulative_ldl_mmol_years"} <= set(frame.columns)


@pytest.mark.parametrize(
    "grid", [np.array([1.0]), np.array([0.0, 0.0, 1.0]), np.array([0.0, np.nan])]
)
def test_invalid_grids_are_rejected(catalog, population, make_regimen, grid):
    with pytest.raises(ValueError):
        simulate_biomarkers(make_regimen("ezetimibe"), catalog, population, grid)


def test_adherence_for_a_non_member_is_rejected(catalog, population, make_regimen, days):
    with pytest.raises(ValueError, match="outside the regimen"):
        simulate_biomarkers(
            make_regimen("ezetimibe"), catalog, population, days, {"colchicine": Adherence()}
        )


def test_all_therapies_produce_finite_biomarkers(catalog, population, make_regimen, days):
    for therapy_id in catalog:
        trajectory = simulate_biomarkers(make_regimen(therapy_id), catalog, population, days)
        for array in (
            trajectory.ldl_mg_dl,
            trajectory.apob_mg_dl,
            trajectory.lpa_nmol_l,
            trajectory.hscrp_mg_l,
            trajectory.systolic_bp_mmhg,
        ):
            assert np.isfinite(array).all(), therapy_id
