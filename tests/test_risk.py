"""Event-risk scenarios: the scientific guards that prevent overclaiming."""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from ihd_solutions.evidence import CreditMode, EndpointClass
from ihd_solutions.pharmacology import Adherence
from ihd_solutions.risk import (
    RiskParameters,
    benefit_ramp,
    classify_therapies,
    decisions_frame,
    per_unit_relative_risk,
    simulate_events,
)


def _unmatched_population(population):
    """A population id with no declared anchor, so nothing can be credited."""
    return replace(population, population_id="no_evidence_population")


def test_background_only_matches_its_own_baseline_exactly(
    catalog, population, anchors, make_regimen, days
):
    """A regimen with no creditable effect must reproduce the untreated curve."""
    simulation = simulate_events(
        make_regimen("colchicine"), catalog, _unmatched_population(population), anchors, days
    )
    assert simulation.events_avoided_per_1000 == pytest.approx(0.0, abs=1e-9)
    assert np.allclose(simulation.relative_risk, 1.0)


def test_baseline_curve_is_reproduced_when_nothing_is_credited(
    catalog, population, anchors, make_regimen, days
):
    simulation = simulate_events(
        make_regimen("methotrexate_low_dose"),
        catalog,
        replace(population, population_id="unknown_pop"),
        anchors,
        days,
    )
    assert np.allclose(simulation.cumulative_incidence, simulation.baseline_incidence)


def test_strict_mode_denies_unproven_ldl_modalities(
    catalog, population, anchors, with_background, days
):
    for therapy_id in ("enlicitide_oral", "azd0780_oral", "verve102_base_editing", "inclisiran"):
        simulation = simulate_events(
            with_background(therapy_id), catalog, population, anchors, days, CreditMode.STRICT
        )
        reference = simulate_events(
            with_background(), catalog, population, anchors, days, CreditMode.STRICT
        )
        assert simulation.events_avoided_per_1000 == pytest.approx(
            reference.events_avoided_per_1000
        ), therapy_id
        decision = next(d for d in simulation.decisions if d.therapy_id == therapy_id)
        assert not decision.credited


def test_surrogate_mode_credits_unproven_ldl_modalities_with_a_flag(
    catalog, population, anchors, with_background, days
):
    reference = simulate_events(
        with_background(), catalog, population, anchors, days, CreditMode.SURROGATE
    )
    for therapy_id in ("enlicitide_oral", "azd0780_oral", "verve102_base_editing", "inclisiran"):
        simulation = simulate_events(
            with_background(therapy_id), catalog, population, anchors, days, CreditMode.SURROGATE
        )
        assert simulation.events_avoided_per_1000 > reference.events_avoided_per_1000
        assert any("transported_to" in note for note in simulation.notes), therapy_id


def test_neutral_il6_blockade_produces_no_benefit_despite_large_crp_drop(
    catalog, populations, anchors, with_background, days
):
    population = populations["ascvd_ckd_inflammation"]
    reference = simulate_events(
        with_background(), catalog, population, anchors, days, CreditMode.STRICT
    )
    simulation = simulate_events(
        with_background("ziltivekimab"), catalog, population, anchors, days, CreditMode.STRICT
    )
    assert simulation.biomarkers.hscrp_mg_l[-1] < reference.biomarkers.hscrp_mg_l[-1] * 0.5
    delta = simulation.events_avoided_per_1000 - reference.events_avoided_per_1000
    assert abs(delta) < 5.0
    decision = next(d for d in simulation.decisions if d.therapy_id == "ziltivekimab")
    assert "neutral_confidence_interval" in decision.warnings


def test_methotrexate_negative_control_does_not_improve_events(
    catalog, population, anchors, with_background, days
):
    reference = simulate_events(with_background(), catalog, population, anchors, days)
    simulation = simulate_events(
        with_background("methotrexate_low_dose"), catalog, population, anchors, days
    )
    assert simulation.events_avoided_per_1000 <= reference.events_avoided_per_1000 + 1e-9
    assert abs(simulation.events_avoided_per_1000 - reference.events_avoided_per_1000) < 3.0


def test_colchicine_benefit_is_population_specific(
    catalog, populations, anchors, with_background, days
):
    regimen = with_background("colchicine")
    stable = simulate_events(regimen, catalog, populations["stable_ccs"], anchors, days)
    early_post_mi = simulate_events(regimen, catalog, populations["post_mi_pci"], anchors, days)
    stable_reference = simulate_events(
        with_background(), catalog, populations["stable_ccs"], anchors, days
    )
    post_mi_reference = simulate_events(
        with_background(), catalog, populations["post_mi_pci"], anchors, days
    )
    stable_gain = stable.events_avoided_per_1000 - stable_reference.events_avoided_per_1000
    post_mi_gain = early_post_mi.events_avoided_per_1000 - post_mi_reference.events_avoided_per_1000
    assert stable_gain > 10.0
    assert abs(post_mi_gain) < 5.0
    assert stable_gain > post_mi_gain


def test_heart_failure_endpoint_is_never_credited_as_mace(
    catalog, populations, anchors, with_background, days
):
    population = populations["ischemic_cm_heart_failure"]
    reference = simulate_events(with_background(), catalog, population, anchors, days)
    simulation = simulate_events(
        with_background("sglt2_inhibitor_hf"), catalog, population, anchors, days
    )
    assert simulation.events_avoided_per_1000 == pytest.approx(reference.events_avoided_per_1000)
    decision = next(d for d in simulation.decisions if d.therapy_id == "sglt2_inhibitor_hf")
    assert not decision.credited
    assert anchors[decision.anchor_id].endpoint_class is EndpointClass.HF_COMPOSITE


def test_regenerative_sheet_never_receives_event_credit(
    catalog, populations, anchors, make_regimen, days
):
    population = populations["ischemic_cm_heart_failure"]
    reference = simulate_events(
        make_regimen("bp_intensification"), catalog, population, anchors, days
    )
    simulation = simulate_events(
        make_regimen("bp_intensification", "ipsc_cm_sheet"), catalog, population, anchors, days
    )
    assert simulation.events_avoided_per_1000 == pytest.approx(reference.events_avoided_per_1000)


def test_lpa_lowering_is_not_credited_without_outcomes(
    catalog, populations, anchors, with_background, days
):
    population = populations["ascvd_high_lpa"]
    reference = simulate_events(with_background(), catalog, population, anchors, days)
    for therapy_id in ("pelacarsen", "olpasiran", "muvalaplin_oral"):
        simulation = simulate_events(
            with_background(therapy_id), catalog, population, anchors, days, CreditMode.SURROGATE
        )
        assert simulation.events_avoided_per_1000 == pytest.approx(
            reference.events_avoided_per_1000
        )
        assert simulation.biomarkers.lpa_nmol_l[-1] < population.lp_a_nmol_l
        # Pelacarsen now has a *negative* outcomes trial (Lp(a)HORIZON 2026); the
        # others are still unproven. Either way the event curve must not move.
        assert (
            "lpa_lowering_not_outcomes_validated" in simulation.notes
            or f"neutral_trials_held_at_null:{therapy_id}" in simulation.notes
            or f"negative_outcomes_trial_no_benefit:{therapy_id}" in simulation.notes
        )


def test_pci_for_stable_disease_does_not_reduce_events(
    catalog, population, anchors, with_background, days
):
    reference = simulate_events(with_background(), catalog, population, anchors, days)
    simulation = simulate_events(
        with_background("pci_stable_cad"), catalog, population, anchors, days
    )
    assert abs(simulation.events_avoided_per_1000 - reference.events_avoided_per_1000) < 8.0
    assert catalog["pci_stable_cad"].symptom_benefit_only


def test_outcomes_proven_pcsk9_antibody_is_credited_in_strict_mode(
    catalog, population, anchors, with_background, days
):
    reference = simulate_events(with_background(), catalog, population, anchors, days)
    simulation = simulate_events(with_background("evolocumab"), catalog, population, anchors, days)
    assert simulation.events_avoided_per_1000 > reference.events_avoided_per_1000
    decision = next(d for d in simulation.decisions if d.therapy_id == "evolocumab")
    assert decision.credited
    assert decision.anchor_id == "FOURIER"


def test_fixed_anchor_and_per_unit_pathway_are_not_double_counted(
    catalog, population, anchors, with_background, days
):
    partition = classify_therapies(
        with_background("evolocumab"), population, anchors, CreditMode.STRICT
    )
    assert "evolocumab" in partition.fixed
    assert "statin_high_intensity" in partition.per_unit
    assert "evolocumab" not in partition.ldl_participants(catalog, CreditMode.STRICT)


def test_number_needed_to_treat_is_the_inverse_of_arr(
    catalog, population, anchors, with_background, days
):
    simulation = simulate_events(with_background("evolocumab"), catalog, population, anchors, days)
    assert simulation.number_needed_to_treat == pytest.approx(
        1.0 / simulation.absolute_risk_reduction
    )
    assert simulation.events_avoided_per_1000 == pytest.approx(
        1000.0 * simulation.absolute_risk_reduction
    )


def test_nnt_is_infinite_when_nothing_is_credited(catalog, population, anchors, make_regimen, days):
    simulation = simulate_events(
        make_regimen("colchicine"),
        catalog,
        replace(population, population_id="no_evidence_population"),
        anchors,
        days,
    )
    assert simulation.number_needed_to_treat == float("inf")
    assert simulation.relative_risk_reduction == pytest.approx(0.0, abs=1e-9)


def test_combined_benefit_is_capped(catalog, population, anchors, make_regimen, days):
    regimen = make_regimen(
        "statin_high_intensity",
        "ezetimibe",
        "bempedoic_acid",
        "evolocumab",
        "colchicine",
        "bp_intensification",
        "aspirin_secondary_prevention",
        "icosapent_ethyl",
    )
    params = RiskParameters(max_combined_rrr=0.50)
    simulation = simulate_events(
        regimen, catalog, population, anchors, days, CreditMode.EXPLORATORY, params
    )
    assert simulation.relative_risk.min() >= 0.50 - 1e-12
    assert simulation.capped
    assert "combined_relative_risk_reduction_capped" in simulation.notes


def test_better_adherence_never_reduces_credited_benefit(
    catalog, population, anchors, with_background, days
):
    regimen = with_background("evolocumab")
    full = simulate_events(regimen, catalog, population, anchors, days)
    half = simulate_events(
        regimen,
        catalog,
        population,
        anchors,
        days,
        adherence={tid: Adherence(taking_fraction=0.5) for tid in regimen.therapy_ids},
    )
    assert full.events_avoided_per_1000 > half.events_avoided_per_1000


def test_longer_horizon_accumulates_more_absolute_benefit(
    catalog, population, anchors, with_background, short_days, days
):
    regimen = with_background("evolocumab")
    short = simulate_events(regimen, catalog, population, anchors, short_days)
    long = simulate_events(regimen, catalog, population, anchors, days)
    assert long.events_avoided_per_1000 > short.events_avoided_per_1000


def test_higher_baseline_risk_gives_larger_absolute_benefit(
    catalog, populations, anchors, with_background, days
):
    regimen = with_background("evolocumab")
    low = simulate_events(regimen, catalog, populations["stable_ccs"], anchors, days)
    high = simulate_events(regimen, catalog, populations["post_mi_hscrp_high"], anchors, days)
    assert high.events_avoided_per_1000 > low.events_avoided_per_1000
    assert high.relative_risk_reduction == pytest.approx(low.relative_risk_reduction, rel=0.05)


def test_first_year_benefit_is_attenuated_relative_to_later_years():
    years = np.array([0.0, 0.5, 1.0, 5.0])
    ramp = benefit_ramp(years, RiskParameters())
    assert ramp[0] == pytest.approx(0.5)
    assert ramp[2] == pytest.approx(1.0)
    assert np.all(np.diff(ramp) >= 0)


def test_per_unit_relative_risk_matches_ctt_arithmetic():
    years = np.array([5.0])
    rr = per_unit_relative_risk(np.array([1.0]), 0.78, years)
    assert rr[0] == pytest.approx(0.78)
    rr_two = per_unit_relative_risk(np.array([2.0]), 0.78, years)
    assert rr_two[0] == pytest.approx(0.78**2)
    early = per_unit_relative_risk(np.array([1.0]), 0.78, np.array([0.0]))
    assert early[0] > rr[0]


def test_per_unit_relative_risk_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="non-negative"):
        per_unit_relative_risk(np.array([-1.0]), 0.78, np.array([1.0]))
    with pytest.raises(ValueError, match="must not exceed 1"):
        per_unit_relative_risk(np.array([1.0]), 1.2, np.array([1.0]))
    with pytest.raises(ValueError):
        per_unit_relative_risk(np.array([1.0]), 0.0, np.array([1.0]))


@pytest.mark.parametrize(
    "kwargs",
    [
        {"year1_benefit_fraction": 1.5},
        {"benefit_ramp_years": 0.0},
        {"max_combined_rrr": 1.0},
        {"sbp_unit_mmhg": 0.0},
    ],
)
def test_risk_parameter_validation(kwargs):
    with pytest.raises(ValueError):
        RiskParameters(**kwargs)


def test_invalid_grid_is_rejected(catalog, population, anchors, make_regimen):
    with pytest.raises(ValueError):
        simulate_events(make_regimen("ezetimibe"), catalog, population, anchors, np.array([1.0]))


def test_summary_row_and_decisions_frame_are_consistent(
    catalog, population, anchors, with_background, days
):
    simulation = simulate_events(with_background("colchicine"), catalog, population, anchors, days)
    row = simulation.summary_row()
    assert row["population_id"] == population.population_id
    assert "colchicine" in row["credited_therapies"]
    frame = decisions_frame([simulation])
    assert len(frame) == len(simulation.decisions)
    assert set(frame["therapy_id"]) == {d.therapy_id for d in simulation.decisions}


def test_every_therapy_and_population_combination_stays_finite(
    catalog, populations, anchors, make_regimen, short_days
):
    for population in populations.values():
        for therapy_id in catalog:
            simulation = simulate_events(
                make_regimen(therapy_id), catalog, population, anchors, short_days
            )
            assert np.isfinite(simulation.cumulative_incidence).all()
            assert 0.0 <= simulation.event_risk < 1.0
            assert simulation.relative_risk.min() > 0.0
