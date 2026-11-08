"""Plaque remodelling scenarios: geometry, thresholds and horizon bookkeeping."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ihd_solutions.plaque import (
    CAP_THRESHOLD_UM,
    STRESS_THRESHOLD_KPA,
    PlaqueParameters,
    simulate_plaque,
)


def test_initial_stenosis_is_derived_from_geometry():
    params = PlaqueParameters()
    expected = 100.0 * 2.2 * 0.40 / (np.pi * 1.5**2)
    assert params.initial_stenosis_pct == pytest.approx(expected)
    assert params.initial_stenosis_pct == pytest.approx(12.449453, abs=1e-5)


def test_lumen_radius_shrinks_as_the_core_grows():
    params = PlaqueParameters()
    assert params.lumen_radius_mm(0.4) > params.lumen_radius_mm(0.8)
    with pytest.raises(ValueError, match="occludes"):
        params.lumen_radius_mm(4.0)


def test_thinner_cap_raises_hoop_stress():
    params = PlaqueParameters()
    assert params.hoop_stress_kpa(0.4, 60.0) > params.hoop_stress_kpa(0.4, 120.0)


def test_untreated_scenario_crosses_before_treated_scenarios(catalog, population, make_regimen):
    untreated = simulate_plaque(make_regimen("colchicine"), catalog, population, months=60)
    statin = simulate_plaque(make_regimen("statin_high_intensity"), catalog, population, months=60)
    intensive = simulate_plaque(
        make_regimen("statin_high_intensity", "evolocumab"), catalog, population, months=60
    )
    assert untreated.crossed
    assert untreated.first_crossing_month < 60
    assert intensive.min_cap_um > statin.min_cap_um > untreated.min_cap_um
    assert intensive.max_stress_kpa < statin.max_stress_kpa < untreated.max_stress_kpa


def test_intact_horizon_stops_at_the_first_crossing(catalog, population, make_regimen):
    outcome = simulate_plaque(make_regimen("colchicine"), catalog, population, months=60)
    assert outcome.intact_months == outcome.first_crossing_month
    subset = outcome.timeline.loc[~outcome.timeline["post_threshold_counterfactual"]]
    assert len(subset) == outcome.intact_months + 1
    assert subset["cap_thickness_um"].min() == pytest.approx(outcome.min_cap_um)
    assert subset["hoop_stress_kpa"].max() == pytest.approx(outcome.max_stress_kpa)


def test_counterfactual_rows_are_labelled(catalog, population, make_regimen):
    outcome = simulate_plaque(make_regimen("colchicine"), catalog, population, months=60)
    if outcome.crossed:
        assert outcome.timeline["post_threshold_counterfactual"].any()
        after = outcome.timeline.loc[outcome.timeline["post_threshold_counterfactual"]]
        assert after["month"].min() == outcome.first_crossing_month + 1


def test_stable_scenario_reports_no_crossing(catalog, population, make_regimen):
    outcome = simulate_plaque(
        make_regimen("statin_high_intensity", "evolocumab", "colchicine"),
        catalog,
        population,
        months=60,
    )
    assert not outcome.crossed
    assert outcome.first_crossing_month is None
    assert outcome.intact_months == 60
    assert not outcome.timeline["post_threshold_counterfactual"].any()


def test_threshold_constants_match_the_literature_heuristics():
    assert CAP_THRESHOLD_UM == 65.0
    assert STRESS_THRESHOLD_KPA == 300.0


def test_core_area_never_falls_below_its_floor(catalog, population, make_regimen):
    params = PlaqueParameters(lipid_influx_coeff=0.0)
    outcome = simulate_plaque(
        make_regimen("statin_high_intensity"), catalog, population, months=60, params=params
    )
    assert outcome.timeline["core_area_mm2"].min() >= params.core_floor_mm2


def test_cap_thickness_is_bounded(catalog, population, make_regimen):
    outcome = simulate_plaque(
        make_regimen("statin_high_intensity", "evolocumab"), catalog, population, months=120
    )
    assert outcome.timeline["cap_thickness_um"].max() <= PlaqueParameters().cap_ceiling_um
    assert outcome.timeline["cap_thickness_um"].min() >= 5.0


def test_mmp_suppression_slows_cap_loss(catalog, population, make_regimen):
    with_colchicine = simulate_plaque(
        make_regimen("statin_high_intensity", "colchicine"), catalog, population, months=60
    )
    without = simulate_plaque(make_regimen("statin_high_intensity"), catalog, population, months=60)
    assert with_colchicine.min_cap_um >= without.min_cap_um


@pytest.mark.parametrize(
    "changes",
    [
        {"reference_radius_mm": 0.0},
        {"core_floor_mm2": 1.0},
        {"initial_core_area_mm2": 10.0},
        {"lipid_influx_coeff": -1.0},
        {"cap_repair_coeff": -1.0},
    ],
)
def test_parameter_validation(changes):
    with pytest.raises(ValueError):
        PlaqueParameters(**changes)


def test_zero_months_is_rejected(catalog, population, make_regimen):
    with pytest.raises(ValueError, match="months must be positive"):
        simulate_plaque(make_regimen("colchicine"), catalog, population, months=0)


def test_timeline_shape_and_finiteness(catalog, population, make_regimen):
    outcome = simulate_plaque(make_regimen("ezetimibe"), catalog, population, months=36)
    assert len(outcome.timeline) == 37
    numeric = outcome.timeline.select_dtypes(include=[np.number])
    assert np.isfinite(numeric.to_numpy()).all()


def test_extreme_influx_fails_loudly(catalog, population, make_regimen):
    params = replace(PlaqueParameters(), lipid_influx_coeff=5.0)
    with pytest.raises(ValueError):
        simulate_plaque(make_regimen("colchicine"), catalog, population, months=60, params=params)
