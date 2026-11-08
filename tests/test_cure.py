"""Cure criteria: the framework must refuse to call anything a cure.

These tests encode the scientific guard rails of the cure search: animal data
cannot satisfy human proof, partial regression cannot score as eradication, two
half-effective procedures cannot be summed into a restored heart, and the
lifetime model must stay labelled as a genetic scenario.
"""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from ihd_solutions import cure
from ihd_solutions.catalog import DosingProfile, EvidenceTier
from ihd_solutions.evidence import CreditMode
from ihd_solutions.risk import simulate_events


@pytest.fixture(scope="module")
def assessments(catalog, anchors, populations, days):
    return cure.assess_catalog(catalog, anchors, populations["stable_ccs"], days)


def test_weights_and_criteria_are_explicit():
    assert math.isclose(sum(cure.CRITERION_WEIGHTS.values()), 1.0)
    assert set(cure.CRITERION_WEIGHTS) == {c.value for c in cure.CureCriterion}
    for criterion in cure.CureCriterion:
        assert len(cure.CRITERION_DESCRIPTIONS[criterion]) > 40


def test_every_therapy_is_scored_on_every_criterion(assessments, catalog):
    assert len(assessments) == len(catalog)
    for criterion in cure.CureCriterion:
        assert assessments[criterion.value].between(0.0, 1.0).all()
    assert assessments["cure_readiness_index"].between(0.0, 1.0).all()
    assert assessments["rationale"].str.len().gt(50).all()


def test_no_therapy_satisfies_every_criterion(assessments):
    assert not assessments["criteria_blocking"].eq("").any()
    assert assessments["cure_readiness_index"].max() < cure.CRITERION_THRESHOLD


def test_animal_only_modalities_cannot_earn_human_proof(assessments):
    for therapy_id in ("pcsk9_vaccine", "car_treg_oxldl", "yap101_sav1_gene_therapy"):
        row = assessments.loc[assessments["therapy_id"].eq(therapy_id)].iloc[0]
        assert row["randomised_proof"] == 0.0
        assert row["strict_rrr_pct"] == pytest.approx(0.0)


def test_failed_and_neutral_trials_do_not_score_as_proven_benefit(assessments):
    for therapy_id in ("pelacarsen", "ziltivekimab", "methotrexate_low_dose", "pci_stable_cad"):
        row = assessments.loc[assessments["therapy_id"].eq(therapy_id)].iloc[0]
        assert row["randomised_proof"] <= 0.5
        assert row["strict_rrr_pct"] == pytest.approx(0.0, abs=1e-9)


def test_only_one_time_editing_reaches_durable_driver_control(assessments, catalog):
    satisfied = assessments.loc[
        assessments["durable_driver_control"] >= cure.CRITERION_THRESHOLD, "therapy_id"
    ].tolist()
    assert set(satisfied) == {"verve102_base_editing", "ctx310_angptl3_editing"}
    for therapy_id in satisfied:
        assert catalog[therapy_id].profile is DosingProfile.ONE_TIME_EDIT


def test_burden_reversal_is_capped_because_regression_is_partial(assessments):
    assert assessments["burden_reversal"].max() <= 0.6
    scored = assessments.loc[assessments["burden_reversal"] > 0, "therapy_id"].tolist()
    assert set(scored) == set(cure.BURDEN_REVERSAL_EVIDENCE)


def test_myocardial_restoration_is_capped_by_the_human_data(assessments):
    assert assessments["myocardial_restoration"].max() <= 0.5
    best = assessments.sort_values("myocardial_restoration", ascending=False).iloc[0]
    assert best["therapy_id"] == "biovat_hf_patch"


def test_non_ischaemic_gene_therapy_earns_no_restoration_credit(assessments):
    row = assessments.loc[assessments["therapy_id"].eq("ab1002_gene_therapy")].iloc[0]
    assert row["myocardial_restoration"] == 0.0
    assert row["strict_rrr_pct"] == pytest.approx(0.0)


def test_search_reports_no_covering_combination_and_names_the_ceilings(assessments):
    result = cure.search_cure_combination(assessments)
    assert not result.satisfied
    assert result.selected == ()
    assert set(result.binding_constraints) == {
        "burden_reversal",
        "myocardial_restoration",
        "risk_normalisation",
    }
    assert result.missing
    assert any("no combination" in note for note in result.notes)


def test_relaxing_the_threshold_still_cannot_normalise_risk(assessments):
    relaxed = cure.search_cure_combination(assessments, threshold=0.5)
    assert not relaxed.satisfied
    assert "risk_normalisation" in relaxed.missing


def test_coverage_combines_by_maximum_never_by_summation():
    frame = pd.DataFrame(
        [
            {
                "therapy_id": "a",
                "durable_driver_control": 1.0,
                "burden_reversal": 0.4,
                "myocardial_restoration": 0.4,
                "risk_normalisation": 0.0,
                "randomised_proof": 1.0,
            },
            {
                "therapy_id": "b",
                "durable_driver_control": 0.0,
                "burden_reversal": 0.4,
                "myocardial_restoration": 0.4,
                "risk_normalisation": 1.0,
                "randomised_proof": 1.0,
            },
        ]
    )
    result = cure.search_cure_combination(frame, threshold=0.8)
    assert not result.satisfied
    # 0.4 + 0.4 must not become 0.8.
    assert result.coverage["burden_reversal"] == pytest.approx(0.4)
    assert result.coverage["myocardial_restoration"] == pytest.approx(0.4)


def test_search_finds_a_cover_when_one_genuinely_exists():
    frame = pd.DataFrame(
        [
            {
                "therapy_id": "driver",
                "durable_driver_control": 1.0,
                "burden_reversal": 0.9,
                "myocardial_restoration": 0.0,
                "risk_normalisation": 0.9,
                "randomised_proof": 1.0,
            },
            {
                "therapy_id": "muscle",
                "durable_driver_control": 0.0,
                "burden_reversal": 0.0,
                "myocardial_restoration": 0.9,
                "risk_normalisation": 0.0,
                "randomised_proof": 0.0,
            },
        ]
    )
    result = cure.search_cure_combination(frame, threshold=0.8)
    assert result.satisfied
    assert result.selected == ("driver", "muscle")
    assert result.missing == ()


def test_search_rejects_invalid_arguments(assessments):
    with pytest.raises(ValueError):
        cure.search_cure_combination(assessments, threshold=0.0)
    with pytest.raises(ValueError):
        cure.search_cure_combination(assessments, max_size=0)


def test_assessment_rejects_incomplete_scores():
    with pytest.raises(ValueError, match="scores must cover"):
        cure.CureAssessment("x", {"burden_reversal": 0.5}, {})


# --------------------------------------------------------------------------- lifetime


def test_earlier_intervention_removes_more_lifetime_risk():
    scenarios = [cure.lifetime_exposure_scenario(age, 0.62) for age in (30, 40, 50, 60, 70)]
    risks = [scenario.scenario_lifetime_risk for scenario in scenarios]
    assert risks == sorted(risks)
    assert scenarios[0].relative_risk_reduction > scenarios[-1].relative_risk_reduction


def test_lifetime_scenario_is_labelled_as_genetic_not_randomised():
    scenario = cure.lifetime_exposure_scenario(40, 0.62)
    assert any("Mendelian-randomisation scenario" in note for note in scenario.notes)
    assert any("does not model" in note for note in scenario.notes)
    assert 0.0 < scenario.scenario_lifetime_risk < scenario.baseline_lifetime_risk
    assert scenario.events_avoided_per_1000 > 0


def test_no_reduction_means_no_benefit():
    scenario = cure.lifetime_exposure_scenario(40, 0.0)
    assert scenario.relative_risk == pytest.approx(1.0)
    assert scenario.scenario_lifetime_risk == pytest.approx(scenario.baseline_lifetime_risk)
    assert scenario.events_avoided_per_1000 == pytest.approx(0.0)


def test_full_window_coverage_reproduces_the_genetic_contrast():
    scenario = cure.lifetime_exposure_scenario(
        20, 38.67 / 160.0, exposure_start_age=20, baseline_ldl_mg_dl=160.0
    )
    assert scenario.exposure_fraction_covered == pytest.approx(1.0)
    assert scenario.delta_mmol_l == pytest.approx(1.0)
    assert scenario.relative_risk == pytest.approx(0.46)


def test_lifetime_relative_risk_is_floored():
    scenario = cure.lifetime_exposure_scenario(20, 0.95, baseline_ldl_mg_dl=400.0)
    assert scenario.relative_risk >= cure.MIN_SCENARIO_RR
    assert any("floor" in note for note in scenario.notes)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"start_age": 10, "ldl_reduction_fraction": 0.5},
        {"start_age": 80, "ldl_reduction_fraction": 0.5},
        {"start_age": 40, "ldl_reduction_fraction": 1.0},
        {"start_age": 40, "ldl_reduction_fraction": -0.1},
        {"start_age": 40, "ldl_reduction_fraction": 0.5, "genetic_hr_per_mmol_l": 1.0},
        {"start_age": 40, "ldl_reduction_fraction": 0.5, "baseline_lifetime_risk": 1.0},
    ],
)
def test_lifetime_scenario_rejects_invalid_inputs(kwargs):
    with pytest.raises(ValueError):
        cure.lifetime_exposure_scenario(**kwargs)


def test_lifetime_grid_is_complete_and_labelled():
    frame = cure.lifetime_grid(start_ages=(30, 50))
    assert len(frame) == 2 * len(cure.LIFETIME_MODALITIES)
    assert (frame["evidence_class"] == "mendelian_randomisation_scenario").all()
    assert frame["lifetime_rrr_pct"].between(0.0, 100.0).all()


# --------------------------------------------------------------------------- residual


def test_residual_risk_shows_that_excess_risk_survives(
    catalog, populations, anchors, with_background, days
):
    population = populations["stable_ccs"]
    simulation = simulate_events(
        with_background("colchicine", "bp_intensification"),
        catalog,
        population,
        anchors,
        days,
        mode=CreditMode.STRICT,
    )
    residual = cure.residual_risk(simulation, population)
    assert 0.0 < residual.excess_risk_remaining_fraction < 1.0
    assert not residual.normalised
    assert residual.treated_risk_pct < residual.baseline_risk_pct
    assert residual.reference_risk_pct < residual.treated_risk_pct


def test_residual_risk_rejects_a_nonpositive_reference(
    catalog, populations, anchors, with_background, days
):
    population = populations["stable_ccs"]
    simulation = simulate_events(
        with_background(), catalog, population, anchors, days, mode=CreditMode.STRICT
    )
    with pytest.raises(ValueError):
        cure.residual_risk(simulation, population, reference_annual_rate=0.0)


def test_report_states_the_verdict_and_the_non_claims(assessments):
    cover = cure.search_cure_combination(assessments)
    lifetime = cure.lifetime_grid(start_ages=(30, 60))
    report = cure.cure_search_report(assessments, cover, lifetime, [])
    assert "VERDICT: no catalogued intervention" in report
    assert "BINDING CONSTRAINTS" in report
    assert "No cure is demonstrated" in report
    assert "Mendelian-randomisation scenario" in report


def test_durable_score_is_reduced_for_animal_only_evidence(catalog):
    therapy = catalog["verve102_base_editing"]
    human, _ = cure._durable_driver_score(therapy)
    animal, reason = cure._durable_driver_score(
        replace(therapy, evidence_tier=EvidenceTier.MECHANISTIC)
    )
    assert animal < human
    assert "animal" in reason


def test_waning_edit_scores_below_a_permanent_one(catalog):
    therapy = catalog["verve102_base_editing"]
    permanent, _ = cure._durable_driver_score(therapy)
    waning, _ = cure._durable_driver_score(replace(therapy, waning_half_life_days=1825.0))
    assert permanent == 1.0
    assert waning < permanent


def test_therapies_that_do_not_touch_the_driver_score_zero(catalog):
    for therapy_id in ("pci_stable_cad", "coronary_sinus_reducer", "biovat_hf_patch"):
        score, _ = cure._durable_driver_score(catalog[therapy_id])
        assert score == 0.0


def test_numpy_grid_is_accepted_by_assess_catalog(catalog, anchors, populations):
    days = np.array([0.0, 365.25, 730.5])
    frame = cure.assess_catalog(catalog, anchors, populations["stable_ccs"], days)
    assert len(frame) == len(catalog)


def test_neutral_randomised_anchor_supersedes_biomarker_transport(
    catalog, populations, anchors, with_background, days
):
    """A looser mode can score *lower*, and must say why.

    Obicetrapib has only biomarker evidence in chronic coronary syndrome, so
    surrogate mode transports the CTT relationship to it. Exploratory mode admits
    its own neutral BROADWAY event analysis instead, which is held at the null and
    displaces the transport (no double counting). The inversion is intentional.
    """
    population = populations["stable_ccs"]
    regimen = with_background("obicetrapib")
    reference = simulate_events(
        with_background(), catalog, population, anchors, days, mode=CreditMode.SURROGATE
    )
    surrogate = simulate_events(
        regimen, catalog, population, anchors, days, mode=CreditMode.SURROGATE
    )
    exploratory = simulate_events(
        regimen, catalog, population, anchors, days, mode=CreditMode.EXPLORATORY
    )
    assert surrogate.events_avoided_per_1000 > reference.events_avoided_per_1000
    assert exploratory.events_avoided_per_1000 == pytest.approx(
        reference.events_avoided_per_1000, abs=1e-9
    )
    assert any(
        note.startswith("neutral_anchor_supersedes_biomarker_transport:")
        for note in exploratory.notes
    )


def test_animal_evidence_never_transports_in_any_mode(
    catalog, populations, anchors, with_background, days
):
    population = populations["stable_ccs"]
    reference = simulate_events(
        with_background(), catalog, population, anchors, days, mode=CreditMode.EXPLORATORY
    )
    for therapy_id in ("pcsk9_vaccine", "car_treg_oxldl"):
        for mode in CreditMode:
            simulation = simulate_events(
                with_background(therapy_id), catalog, population, anchors, days, mode=mode
            )
            assert simulation.events_avoided_per_1000 == pytest.approx(
                reference.events_avoided_per_1000, abs=1e-9
            ), (therapy_id, mode)
