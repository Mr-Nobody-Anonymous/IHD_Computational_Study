"""Neutral trials must stay neutral.

ZEUS, CIRT, CLEAR SYNERGY and ISCHEMIA published intervals that span 1.0. Their
point estimates (0.99, 1.01, 0.99, 0.93) are noise around the null, so the default
policy applies HR = 1 exactly. A neutral drug must not create benefit, must not
create harm, and must not change a regimen's event curve in any credit mode.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ihd_solutions.evidence import CreditMode
from ihd_solutions.risk import (
    NEUTRAL_CI_POLICIES,
    RiskParameters,
    effective_hazard_ratio,
    simulate_events,
)
from ihd_solutions.uncertainty import monte_carlo_regimen

NEUTRAL_ANCHORS = ("ZEUS", "CIRT", "CLEAR-SYNERGY", "ISCHEMIA")
NEUTRAL_THERAPIES = {
    "ZEUS": "ziltivekimab",
    "CIRT": "methotrexate_low_dose",
    "CLEAR-SYNERGY": "colchicine",
    "ISCHEMIA": "pci_stable_cad",
}


def test_registry_records_reported_neutrality(anchors):
    for anchor_id in NEUTRAL_ANCHORS:
        assert anchors[anchor_id].neutral_at_report, anchor_id
        assert anchors[anchor_id].is_neutral, anchor_id
    for anchor_id in ("LODOCO2", "CANTOS-150", "SELECT", "FOURIER", "CTT-STATIN"):
        assert not anchors[anchor_id].neutral_at_report, anchor_id


def test_reported_neutrality_survives_resampling(anchors):
    anchor = anchors["ZEUS"]
    resampled = replace(anchor, hr=0.62, ci_low=0.55, ci_high=0.70)
    assert not resampled.is_neutral  # the sampled interval no longer spans 1
    assert resampled.neutral_at_report  # but the published verdict is unchanged


def test_effective_hazard_ratio_holds_neutral_anchors_at_the_null(anchors):
    params = RiskParameters()
    for anchor_id in NEUTRAL_ANCHORS:
        assert effective_hazard_ratio(anchors[anchor_id], params) == 1.0
    assert effective_hazard_ratio(anchors["LODOCO2"], params) == pytest.approx(0.69)


def test_point_estimate_policy_is_available_as_a_sensitivity_analysis(anchors):
    params = RiskParameters(neutral_ci_policy="point_estimate")
    assert effective_hazard_ratio(anchors["ZEUS"], params) == pytest.approx(0.99)
    assert effective_hazard_ratio(anchors["CIRT"], params) == pytest.approx(1.01)


def test_missing_effect_estimate_cannot_be_applied(anchors):
    pending = next(a for a in anchors.values() if not a.has_effect_estimate)
    with pytest.raises(ValueError, match="no effect estimate"):
        effective_hazard_ratio(pending, RiskParameters())


def test_unknown_policy_is_rejected():
    assert {"null_effect", "point_estimate"} == NEUTRAL_CI_POLICIES
    with pytest.raises(ValueError, match="neutral_ci_policy"):
        RiskParameters(neutral_ci_policy="optimistic")


@pytest.mark.parametrize("mode", list(CreditMode))
@pytest.mark.parametrize("anchor_id", NEUTRAL_ANCHORS)
def test_adding_a_neutral_therapy_changes_nothing(
    anchor_id, mode, catalog, populations, anchors, with_background, days
):
    anchor = anchors[anchor_id]
    therapy_id = NEUTRAL_THERAPIES[anchor_id]
    population = populations[anchor.population_id]
    if therapy_id not in catalog:
        pytest.skip(f"{therapy_id} not catalogued")
    try:
        regimen = with_background(therapy_id)
    except ValueError:
        pytest.skip(f"{therapy_id} not admissible in {population.population_id}")
    reference = simulate_events(with_background(), catalog, population, anchors, days, mode=mode)
    simulation = simulate_events(regimen, catalog, population, anchors, days, mode=mode)
    assert simulation.events_avoided_per_1000 == pytest.approx(
        reference.events_avoided_per_1000, abs=1e-12
    )
    assert np.allclose(simulation.relative_risk, reference.relative_risk, atol=1e-12)
    assert any(note.startswith("neutral_trials_held_at_null:") for note in simulation.notes)


def test_neutral_therapy_can_look_beneficial_only_under_the_point_estimate_policy(
    catalog, populations, anchors, with_background, days
):
    population = populations["ascvd_ckd_inflammation"]
    params = RiskParameters(neutral_ci_policy="point_estimate")
    reference = simulate_events(
        with_background(), catalog, population, anchors, days, params=params
    )
    simulation = simulate_events(
        with_background("ziltivekimab"), catalog, population, anchors, days, params=params
    )
    # ZEUS point estimate 0.99 leaks a small, spurious benefit under this policy.
    assert simulation.events_avoided_per_1000 > reference.events_avoided_per_1000
    assert simulation.events_avoided_per_1000 - reference.events_avoided_per_1000 < 5.0


def test_monte_carlo_keeps_uncertainty_but_no_direction_for_neutral_trials(
    catalog, populations, anchors, with_background, short_days
):
    population = populations["ascvd_ckd_inflammation"]
    kwargs = {"draws": 80, "seed": 23}
    treated = monte_carlo_regimen(
        with_background("ziltivekimab"), catalog, population, anchors, short_days, **kwargs
    )
    reference = monte_carlo_regimen(
        with_background(), catalog, population, anchors, short_days, **kwargs
    )
    difference = treated.events_avoided_per_1000 - reference.events_avoided_per_1000
    assert abs(float(np.median(difference))) < 5.0
    # Both directions occur: the interval is informative, the direction is not.
    assert np.any(difference > 0)
    assert np.any(difference < 0)
