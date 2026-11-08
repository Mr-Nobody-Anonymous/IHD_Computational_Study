"""Cross-cutting scientific guards applied to every therapy and population."""

from __future__ import annotations

import numpy as np
import pytest

from ihd_solutions.evidence import (
    CreditMode,
    EndpointClass,
    EvidenceTier,
    assert_no_biomarker_surrogate_credit,
)
from ihd_solutions.risk import simulate_events


@pytest.fixture(scope="module")
def all_simulations(catalog, populations, anchors, make_regimen):
    from ihd_solutions.units import time_grid

    days = time_grid(3.0, 90.0)
    simulations = []
    for population in populations.values():
        for therapy_id in catalog:
            for mode in CreditMode:
                simulations.append(
                    simulate_events(
                        make_regimen(therapy_id), catalog, population, anchors, days, mode
                    )
                )
    return simulations


def test_no_simulation_credits_a_biomarker_pathway(all_simulations):
    decisions = [d for simulation in all_simulations for d in simulation.decisions]
    assert_no_biomarker_surrogate_credit(decisions)
    assert all(d.pathway != "biomarker" for d in decisions)


def test_every_credited_decision_points_at_a_real_anchor(all_simulations, anchors):
    for simulation in all_simulations:
        for decision in simulation.decisions:
            if not decision.credited:
                continue
            assert decision.anchor_id in anchors
            anchor = anchors[decision.anchor_id]
            assert anchor.has_effect_estimate
            assert anchor.evidence_tier is not EvidenceTier.MECHANISTIC
            assert anchor.endpoint_class is not EndpointClass.SURROGATE
            assert anchor.endpoint_class is not EndpointClass.FUNCTIONAL
            assert anchor.endpoint_class is not EndpointClass.PENDING


def test_strict_mode_only_uses_randomised_hard_endpoints(all_simulations, anchors):
    for simulation in all_simulations:
        if simulation.mode is not CreditMode.STRICT:
            continue
        for decision in simulation.decisions:
            if not decision.credited:
                continue
            anchor = anchors[decision.anchor_id]
            assert anchor.evidence_tier in {
                EvidenceTier.META_ANALYSIS,
                EvidenceTier.RANDOMIZED_OUTCOMES,
                EvidenceTier.RANDOMIZED_TOPLINE,
            }
            assert anchor.endpoint_class is EndpointClass.MACE


def test_no_scenario_claims_more_than_the_configured_cap(all_simulations):
    for simulation in all_simulations:
        assert simulation.relative_risk.min() >= 0.15 - 1e-9
        assert simulation.relative_risk_reduction < 0.85 + 1e-9
        assert 0.0 <= simulation.event_risk < 1.0


def test_cumulative_incidence_is_monotone_in_time(all_simulations):
    for simulation in all_simulations:
        assert np.all(np.diff(simulation.cumulative_incidence) >= -1e-12)
        assert np.all(np.diff(simulation.baseline_incidence) >= -1e-12)


def test_inflammatory_therapies_alone_never_reduce_events_in_strict_mode(
    catalog, populations, anchors, make_regimen, short_days
):
    """Any event change must come from an anchor, never from the hsCRP trajectory."""
    for therapy_id in ("ziltivekimab", "nlrp3_inhibitor", "methotrexate_low_dose"):
        for population in populations.values():
            simulation = simulate_events(
                make_regimen(therapy_id),
                catalog,
                population,
                anchors,
                short_days,
                CreditMode.STRICT,
            )
            assert abs(simulation.events_avoided_per_1000) < 25.0
            if therapy_id == "nlrp3_inhibitor":
                assert simulation.events_avoided_per_1000 == pytest.approx(0.0, abs=1e-9)


def test_regenerative_and_lpa_therapies_are_never_credited(
    catalog, populations, anchors, make_regimen, short_days
):
    for therapy_id in ("ipsc_cm_sheet", "pelacarsen", "olpasiran", "muvalaplin_oral"):
        for population in populations.values():
            for mode in CreditMode:
                simulation = simulate_events(
                    make_regimen(therapy_id), catalog, population, anchors, short_days, mode
                )
                assert simulation.events_avoided_per_1000 == pytest.approx(0.0, abs=1e-9), (
                    therapy_id,
                    population.population_id,
                    mode,
                )


def test_deterministic_replay_of_the_same_scenario(
    catalog, population, anchors, with_background, short_days
):
    regimen = with_background("evolocumab", "colchicine")
    first = simulate_events(regimen, catalog, population, anchors, short_days)
    second = simulate_events(regimen, catalog, population, anchors, short_days)
    assert np.array_equal(first.cumulative_incidence, second.cumulative_incidence)
    assert first.summary_row() == second.summary_row()
