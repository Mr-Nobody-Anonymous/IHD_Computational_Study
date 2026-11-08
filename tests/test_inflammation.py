"""Inflammatory cascade: integrator accuracy and the CIRT negative control."""

from __future__ import annotations

import numpy as np
import pytest

from ihd_solutions.inflammation import (
    STATE_NAMES,
    UNTREATED,
    InflammationParameters,
    InflammatoryIntervention,
    initial_state,
    intervention_for,
    simulate_inflammation,
    steady_state,
    with_parameters,
)


def test_long_run_converges_to_the_analytic_equilibrium():
    frame = simulate_inflammation(days=4000.0, step_days=1.0)
    final = frame.iloc[-1][list(STATE_NAMES)].to_numpy(dtype=float)
    assert np.allclose(final, steady_state(), rtol=1e-6)


def test_inhibited_equilibrium_matches_the_integrator():
    intervention = InflammatoryIntervention("il1", il1_inhibition=0.7)
    frame = simulate_inflammation(intervention, days=4000.0, step_days=1.0)
    final = frame.iloc[-1][list(STATE_NAMES)].to_numpy(dtype=float)
    assert np.allclose(final, steady_state(intervention=intervention), rtol=1e-6)


def test_states_stay_non_negative_and_finite():
    frame = simulate_inflammation(days=365.0, step_days=0.5)
    values = frame[list(STATE_NAMES)].to_numpy(dtype=float)
    assert np.isfinite(values).all()
    assert (values >= 0).all()


def test_methotrexate_is_an_exact_negative_control(catalog):
    intervention = intervention_for(catalog["methotrexate_low_dose"])
    treated = simulate_inflammation(intervention, days=180.0, step_days=0.5)
    untreated = simulate_inflammation(UNTREATED, days=180.0, step_days=0.5)
    assert np.allclose(
        treated[list(STATE_NAMES)].to_numpy(), untreated[list(STATE_NAMES)].to_numpy()
    )


def test_il6_blockade_lowers_crp_more_than_il1_blockade(catalog):
    il6 = simulate_inflammation(intervention_for(catalog["ziltivekimab"]), days=365.0)
    il1 = simulate_inflammation(intervention_for(catalog["canakinumab"]), days=365.0)
    assert il6["hsCRP"].iloc[-1] < il1["hsCRP"].iloc[-1]


def test_stronger_inhibition_monotonically_lowers_crp():
    values = [
        simulate_inflammation(InflammatoryIntervention("x", il1_inhibition=strength), days=365.0)[
            "hsCRP"
        ].iloc[-1]
        for strength in (0.0, 0.2, 0.5, 0.8)
    ]
    assert all(later < earlier for earlier, later in zip(values, values[1:], strict=False))


def test_statin_reduces_the_oxldl_driver(catalog):
    intervention = intervention_for(catalog["statin_high_intensity"])
    assert intervention.oxldl_reduction == pytest.approx(0.50)
    frame = simulate_inflammation(intervention, days=365.0)
    untreated = simulate_inflammation(UNTREATED, days=365.0)
    assert frame["oxLDL"].iloc[-1] < untreated["oxLDL"].iloc[-1]


def test_colchicine_maps_to_nlrp3_and_mmp_targets(catalog):
    intervention = intervention_for(catalog["colchicine"])
    assert intervention.nlrp3_inhibition > 0
    assert intervention.mmp_suppression > 0
    assert intervention.il6_inhibition == 0.0


def test_unknown_mechanism_group_engages_nothing(catalog):
    intervention = intervention_for(catalog["pelacarsen"])
    assert intervention.nlrp3_inhibition == 0.0
    assert intervention.il1_inhibition == 0.0
    assert intervention.il6_inhibition == 0.0


def test_integrator_is_step_size_convergent():
    coarse = simulate_inflammation(days=180.0, step_days=1.0)["hsCRP"].iloc[-1]
    fine = simulate_inflammation(days=180.0, step_days=0.1)["hsCRP"].iloc[-1]
    assert coarse == pytest.approx(fine, rel=1e-6)


def test_overflow_is_reported_instead_of_written():
    params = with_parameters(oxldl_influx=1e308)
    with pytest.raises(ValueError, match="non-finite"):
        simulate_inflammation(UNTREATED, params, days=50.0, step_days=1.0)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"days": 0.0},
        {"step_days": 0.0},
        {"days": 1.0, "step_days": 2.0},
        {"state0": np.zeros(3)},
        {"state0": np.array([1.0, -1.0, 0.0, 0.0, 0.0, 0.0])},
        {"state0": np.array([np.nan, 0.0, 0.0, 0.0, 0.0, 0.0])},
    ],
)
def test_invalid_arguments_are_rejected(kwargs):
    with pytest.raises(ValueError):
        simulate_inflammation(**kwargs)


def test_parameter_and_intervention_validation():
    with pytest.raises(ValueError):
        InflammationParameters(crp_decay=0.0)
    with pytest.raises(ValueError):
        InflammatoryIntervention("x", il1_inhibition=1.0)
    assert initial_state().shape == (len(STATE_NAMES),)
