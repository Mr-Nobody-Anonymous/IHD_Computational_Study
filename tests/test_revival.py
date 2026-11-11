"""Tests for the control-arm correction.

The module under test exists to LOWER a number. These tests therefore pin four
properties that an optimistic model would fail: a correction may never raise an
estimate, the zero-comparator rung must reproduce the published figure exactly, no
amount of re-searching may manufacture a robust design, and stratifying a population
may never be reported as a gain.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ihd_solutions import revival as R
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solver import SCAR_RECOVERABLE_LVEF_PP
from ihd_solutions.recombine import ATLAS5, Composite, cost_ledger, evaluate_composite

FAST = 600


@pytest.fixture(scope="module")
def ladder():
    return R.correction_ladder(ATLAS5, 1200)


@pytest.fixture(scope="module")
def sensitivity():
    return R.natural_history_sensitivity(ATLAS5, 800)


@pytest.fixture(scope="module")
def search():
    return R.revival_search(800)


# ------------------------------------------------------------------ documented sources


def test_every_deficiency_cites_a_source():
    table = R.correction_table()
    assert len(table) == len(R.CORRECTIONS) >= 5
    assert (table["source"].str.len() > 20).all()
    assert (table["finding"].str.len() > 40).all()


def test_the_control_arm_value_is_documented_not_invented():
    """3.4 pp is read off REVIVED-BCIS2's OMT arm, not chosen for convenience."""
    assert pytest.approx(3.4) == R.NATURAL_HISTORY_PP
    low, high = R.NATURAL_HISTORY_RANGE
    assert low == pytest.approx(1.1) and high == pytest.approx(3.4)
    assert low <= R.NATURAL_HISTORY_PP <= high


def test_lvef_change_sd_reuses_the_gap4_family():
    """The SD must come from the existing sweep, not a second invented number."""
    from ihd_solutions.gap_solutions import cato_power_analysis

    swept = set(cato_power_analysis()["assumed_sd_pp"].astype(float))
    assert R.LVEF_CHANGE_SD_PP in swept


def test_stich_and_revived_are_both_cited():
    joined = " ".join(item.source + item.finding for item in R.CORRECTIONS)
    assert "STICH" in joined
    assert "REVIVED" in joined


# ------------------------------------------------------- direction: only ever downward


def test_correction_only_ever_lowers_the_estimate():
    """The single most important property: no comparator can raise the claim."""
    uncorrected = evaluate_composite(ATLAS5, 900)
    values = [
        R.evaluate_revival(ATLAS5, 900, natural_history_pp=value).probability_of_target
        for value in (1.1, 2.25, 3.4)
    ]
    assert all(value <= uncorrected.probability_of_target + 1e-12 for value in values)


def test_corrected_risk_adjustment_uses_each_simulated_patients_dose_count():
    draws, seed = 900, 123
    gain, used = R.revival_draws(ATLAS5, draws, seed)
    target = gain >= CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP
    survival = (1.0 - ATLAS5.per_dose_risk) ** used
    result = R.evaluate_revival(ATLAS5, draws, seed)
    assert result.risk_adjusted_probability == pytest.approx(float(np.mean(target * survival)))
    assert result.cumulative_adverse_risk == pytest.approx(float(np.mean(1.0 - survival)))


def test_more_correction_never_yields_more_probability():
    values = [0.0, 1.1, 2.25, 3.4]
    probabilities = [
        R.evaluate_revival(ATLAS5, 900, natural_history_pp=v).probability_of_target for v in values
    ]
    assert probabilities == sorted(probabilities, reverse=True)


def test_zero_comparator_reproduces_the_published_figure_exactly():
    """Rung 0 must BE the published number, not a noisy draw near it."""
    published = evaluate_composite(ATLAS5, 900).probability_of_target
    corrected = R.evaluate_revival(ATLAS5, 900, natural_history_pp=0.0)
    assert corrected.probability_of_target == pytest.approx(published, abs=1e-12)


def test_comparator_is_signed_and_not_floored():
    """Flooring N(3.4, 5) would rent the mean up to 4.16 pp and over-subtract.

    The physical bound belongs on the outcome, so the comparator must be allowed to be
    negative for patients whose function worsened on guideline therapy.
    """
    draws = R.natural_history_draws(4000)
    assert (draws < 0.0).any()


def test_outcome_never_goes_negative():
    gain, _ = R.revival_draws(ATLAS5, 600)
    assert (gain >= 0.0).all()


def test_hypothetical_uplift_crosses_threshold_only_as_an_explicit_assumption():
    table = R.hypothetical_uniform_uplift_table(R.ATLAS5_FRUGAL, draws=5000)
    current, target_only, risk_adjusted = table.itertuples(index=False)
    assert current.uniform_additional_lvef_pp == 0.0
    assert not current.target_80_met
    assert target_only.target_80_met
    assert not target_only.risk_adjusted_80_met
    assert risk_adjusted.risk_adjusted_80_met
    assert (table["evidence_status"].str.contains("no molecule or mechanism supports it")).all()
    assert (table["status"] == "hypothetical_uniform_uplift_sensitivity_unvalidated").all()


def test_uplift_evidence_gate_fails_closed_when_no_trial_record_is_supplied():
    decision = R.adjudicate_lvef_uplift_evidence()
    assert decision["evidence_eligible"] is False
    assert decision["credited_mean_incremental_lvef_pp"] == 0.0
    assert "randomized, concurrent controlled comparison" in decision["missing_requirements"]
    assert "positive lower 95% confidence bound" in decision["missing_requirements"]
    assert (
        "independently validated patient-level response model" in decision["missing_requirements"]
    )


def test_evidence_gate_requires_validated_patient_level_effect_draws():
    evidence = {
        "candidate": "test candidate",
        "randomized_controlled": True,
        "same_population": True,
        "same_intervention": True,
        "lvef_endpoint": True,
        "validated_response_model": True,
        "source_id": "doi:10.test/example",
        "lower_95_ci_pp": 1.0,
        "patient_level_incremental_lvef_pp": [2.0, 4.0, -1.0],
    }
    decision = R.adjudicate_lvef_uplift_evidence(evidence, expected_draws=3)
    assert decision["evidence_eligible"] is True
    assert decision["credited_mean_incremental_lvef_pp"] == pytest.approx(5.0 / 3.0)
    assert decision["missing_requirements"] == "none"
    assert "independently verified" in decision["evidence_status"]
    assert not R.adjudicate_lvef_uplift_evidence(evidence, expected_draws=4)["evidence_eligible"]


def test_evidence_gate_rejects_class_analogy_and_nonpositive_confidence_bound():
    evidence = {
        "randomized_controlled": True,
        "same_population": True,
        "same_intervention": False,
        "lvef_endpoint": True,
        "validated_response_model": True,
        "source_id": "doi:10.test/example",
        "lower_95_ci_pp": 2.0,
        "patient_level_incremental_lvef_pp": [2.0, 4.0, -1.0],
    }
    decision = R.adjudicate_lvef_uplift_evidence(evidence, expected_draws=3)
    assert not decision["evidence_eligible"]
    assert decision["credited_mean_incremental_lvef_pp"] == 0.0
    evidence["same_intervention"] = True
    evidence["lower_95_ci_pp"] = 0.0
    assert (
        R.adjudicate_lvef_uplift_evidence(evidence, expected_draws=3)[
            "credited_mean_incremental_lvef_pp"
        ]
        == 0.0
    )


def test_evidence_gated_table_applies_only_aligned_validated_effect_draws():
    draws = 5000
    evidence = {
        "candidate": "synthetic test fixture only",
        "randomized_controlled": True,
        "same_population": True,
        "same_intervention": True,
        "lvef_endpoint": True,
        "validated_response_model": True,
        "source_id": "doi:10.test/synthetic-fixture",
        "lower_95_ci_pp": 1.0,
        "patient_level_incremental_lvef_pp": np.full(draws, 5.0),
    }
    row = R.evidence_gated_uplift_table(draws=draws, evidence=evidence).iloc[0]
    assert row["evidence_eligible"]
    assert row["credited_mean_incremental_lvef_pp"] == 5.0
    assert row["target_80_met"]
    assert "independently verified" in row["evidence_status"]


def test_evidence_gated_table_does_not_credit_the_counterfactual_uplift():
    row = R.evidence_gated_uplift_table(draws=5000).iloc[0]
    assert not row["evidence_eligible"]
    assert row["credited_mean_incremental_lvef_pp"] == 0.0
    assert not row["target_80_met"]
    assert not row["risk_adjusted_80_met"]
    assert row["additional_uplift_needed_for_80pct_target_pp"] > 0.0


def test_zero_comparator_is_exactly_zero_not_noise():
    assert (R.natural_history_draws(500, mean_pp=0.0) == 0.0).all()


def test_natural_history_mean_matches_the_documented_value():
    draws = R.natural_history_draws(40_000)
    assert float(np.mean(draws)) == pytest.approx(R.NATURAL_HISTORY_PP, abs=0.1)


# ------------------------------------------------------------------- the ladder


def test_ladder_has_a_rung_for_every_deficiency(ladder):
    assert list(ladder["rung"])[:3] == [0, 1, 2]
    assert ladder["fraction_of_original"].iloc[0] == pytest.approx(1.0)


def test_ladder_shows_the_loss(ladder):
    as_published = float(ladder.iloc[0]["probability_of_target"])
    corrected = float(ladder.iloc[1]["probability_of_target"])
    assert corrected < as_published
    assert as_published - corrected > 0.2


def test_ladder_survives_below_a_quarter_at_the_upper_comparator(ladder):
    """At +3.4 pp less than 75% of the claim survives; the exact value is pinned
    loosely because the rung sits on simulated draws."""
    assert float(ladder.iloc[1]["fraction_of_original"]) < 0.75


# ------------------------------------------------------------- ATLAS-5 does not survive


def test_atlas5_is_no_longer_robust(sensitivity):
    at_upper = sensitivity.loc[sensitivity["natural_history_pp"] == 3.4].iloc[0]
    assert not bool(at_upper["robust"])
    assert not bool(at_upper["survives_correction"])


def test_no_composition_is_robust_after_correction(search):
    frame, selection = search
    assert selection["any_robust"] is False
    assert not frame["corrected_robust"].any()


def test_the_search_no_longer_selects_atlas5(search):
    _, selection = search
    assert selection["atlas5_still_selected"] is False


def test_full_search_at_every_documented_comparator_still_fails_the_80pct_bar():
    frontier = R.comparator_frontier(draws=2000)
    assert list(frontier["natural_history_pp"]) == [1.1, 2.25, 3.4]
    assert (frontier["composition_dose_pairs_searched"] == 256).all()
    assert not frontier["target_80_met"].any()
    assert not frontier["risk_adjusted_80_met"].any()
    assert (
        frontier.iloc[0]["best_probability_of_target"]
        > frontier.iloc[-1]["best_probability_of_target"]
    )


def test_search_covers_every_precedented_dose_count(search):
    frame, selection = search
    maximum = R.max_precedented_doses()
    assert selection["dose_counts_searched"] == list(range(1, maximum + 1))
    assert len(frame) == 64 * maximum
    assert set(frame["doses"]) == set(range(1, maximum + 1))


def test_cross_dose_search_selects_the_highest_precedented_retry_budget(search):
    _, selection = search
    assert selection["best_corrected"]["doses"] == R.max_precedented_doses()
    assert selection["best_corrected"]["mechanisms"] == "retry_redundancy|maturation"
    assert selection["best_corrected"]["label"] == "ATLAS-5F (simulation-only frugal variant)"


def test_raw_search_still_finds_a_high_number(search):
    """The raw number must remain high: the point is that it is UNEARNED."""
    frame, selection = search
    assert selection["best_raw_probability"] > 0.95
    assert selection["best_corrected_probability"] < 0.65


def test_corrected_search_still_respects_precedent(search):
    frame, _ = search
    admissible = frame.loc[frame["admissible"]]
    assert not admissible.empty


# ------------------------------------------------ stratifying is not a population gain


def test_stratification_claims_no_population_gain():
    table = R.transmurality_stratification()
    assert (table["population_level_gain_claimed"] == 0.0).all()


def test_transmural_bands_match_the_published_recovery_rates():
    table = R.transmurality_stratification().set_index("transmurality_band")
    assert float(table.loc["no_lge", "conditional_recovery"]) == pytest.approx(0.91)
    assert float(table.loc["lge_over_75", "conditional_recovery"]) == pytest.approx(0.08)


def test_stratification_is_labelled_a_stratification():
    table = R.transmurality_stratification()
    assert table["reason"].str.contains("stratification").all()


# --------------------------------------------------------------- surrogate validity


def test_surrogate_table_carries_the_absent_mortality_link():
    table = R.surrogate_translation(800)
    mortality = table.loc[table["endpoint"] == "all cause mortality"].iloc[0]
    assert not bool(mortality["link_established_in_this_population"])


def test_saturated_mean_is_labelled_unusable_rather_than_a_pass():
    """A mean sitting at the scar ceiling is a property of the ceiling.

    The model clamps every increment at SCAR_RECOVERABLE_LVEF_PP, so once most draws
    saturate the mean carries no information about expected clinical benefit. The table
    must say so instead of reporting the STE comparison as met.
    """
    table = R.surrogate_translation(900)
    assert (~table["comparison_usable"]).all()
    assert table["model_exceeds_observed_range"].all()
    # The corrected output is still several times the largest human observation.
    assert (table["corrected_median_gain_pp"] > 2 * R.LARGEST_HUMAN_OBSERVATION_PP).all()


def test_the_observed_ceiling_is_the_documented_value():
    assert pytest.approx(6.9) == R.LARGEST_HUMAN_OBSERVATION_PP
    from ihd_solutions.gap_solutions import BIOVAT_LVEF_GAIN_17MO_PP

    assert pytest.approx(BIOVAT_LVEF_GAIN_17MO_PP) == R.LARGEST_HUMAN_OBSERVATION_PP


def test_below_the_ceiling_the_ste_comparison_is_live():
    table = R.surrogate_translation(900)
    row = table.loc[table["endpoint"] == "all cause mortality"].iloc[0]
    assert float(row["corrected_median_gain_pp"]) < SCAR_RECOVERABLE_LVEF_PP


# -------------------------------------------------------------- verdict and report


def test_verdict_prices_the_attribution_loss(ladder, sensitivity, search):
    frame, selection = search
    verdict = R.revival_verdict(selection, ladder, sensitivity)
    assert verdict["attribution_loss_pp"] > 20.0
    assert verdict["no_longer_robust"] is True
    assert verdict["atlas5_still_best_after_correction"] is False


def test_report_states_the_failure_plainly(ladder, sensitivity, search):
    frame, selection = search
    text = R.revival_report(
        ladder,
        sensitivity,
        R.surrogate_translation(800),
        R.transmurality_stratification(),
        frame,
        selection,
        cost_ledger(ATLAS5),
    )
    assert "NO LONGER ROBUST" in text
    assert "ATTRIBUTION LOSS" in text
    assert "CONTROL ARM" in text
    assert "0.00%" in text


def test_report_never_claims_a_cure(ladder, sensitivity, search):
    frame, selection = search
    text = R.revival_report(
        ladder,
        sensitivity,
        R.surrogate_translation(800),
        R.transmurality_stratification(),
        frame,
        selection,
        cost_ledger(ATLAS5),
    )
    lowered = text.lower()
    assert "not a therapy" in lowered
    assert "no cure is declared" in lowered
    assert "cure achieved" not in lowered


def test_study_writes_twelve_artefacts(tmp_path):
    summary = R.run_revival_study(tmp_path / "revival", draws=FAST)
    files = sorted((tmp_path / "revival").rglob("*"))
    written = [path for path in files if path.is_file()]
    assert len(written) == summary["artefacts"] == 12
    assert (tmp_path / "revival" / "DISCLAIMER.txt").exists()


def test_study_summary_reports_the_failure(tmp_path):
    summary = R.run_revival_study(tmp_path / "revival", draws=FAST)
    assert summary["no_longer_robust"] is True
    assert summary["any_robust_after_correction"] is False
    assert summary["attribution_loss_pp"] > 20.0
    assert summary["composition_dose_pairs_searched"] == 256
    assert summary["draws_per_pair"] == 2000
    assert summary["best_corrected_doses"] == 4
    assert summary["any_robust_at_documented_comparator"] is False
    assert summary["best_probability_at_any_documented_comparator"] < CRITERION_THRESHOLD
    assert summary["status"] == R.REVIVAL_TAG


def test_every_table_carries_the_tag(tmp_path):
    R.run_revival_study(tmp_path / "revival", draws=FAST)
    for path in (tmp_path / "revival" / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        allowed_statuses = {
            R.REVIVAL_TAG,
            "hypothetical_uniform_uplift_sensitivity_unvalidated",
            "evidence_gated_uplift_assessment_unvalidated",
        }
        assert frame["status"].isin(allowed_statuses).all(), path.name


def test_cli_exposes_the_revival_command():
    from ihd_solutions.cli import build_parser

    parser = build_parser()
    args = parser.parse_args(["revival", "--draws", "100"])
    assert args.command == "revival"
    assert args.draws == 100
    assert parser.parse_args(["revival"]).draws == 20_000


# ------------------------------------------------------------------ structural guards


def test_correction_does_not_change_the_design_only_its_attribution():
    """The mechanisms kept are the same; what changes is what they are credited with."""
    design = Composite(enabled=ATLAS5.enabled, doses=ATLAS5.doses)
    corrected = R.evaluate_revival(design, 800, natural_history_pp=3.4)
    assert corrected.mechanisms == len(ATLAS5.enabled)
    assert corrected.saturated_fraction <= 1.0


def test_scar_recoverable_ceiling_is_unmoved():
    """The correction subtracts attribution; it does not rescale the model."""
    assert pytest.approx(25.0) == SCAR_RECOVERABLE_LVEF_PP
    assert pytest.approx(0.8) == CRITERION_THRESHOLD


# ------------------------------------------------- alternative report branches


def _render(search, **overrides):
    """Render the report against a supplied verdict so every branch is reachable."""
    frame, selection = search
    ladder = R.correction_ladder(ATLAS5, 600)
    sensitivity = R.natural_history_sensitivity(ATLAS5, 600)
    verdict = R.revival_verdict(selection, ladder, sensitivity)
    verdict.update(overrides)
    # Re-vote the verdict with the overrides applied: the report re-derives it, so the
    # branch is chosen by the frames, not by the dict.
    return R.revival_report(
        ladder,
        sensitivity,
        R.surrogate_translation(600),
        R.transmurality_stratification(),
        frame,
        selection,
        cost_ledger(ATLAS5),
    ), verdict


def test_report_branch_still_robust(search):
    """The 'still robust' branch must exist and must not be reachable by accident."""
    frame, selection = search
    ladder = R.correction_ladder(ATLAS5, 300).copy()
    ladder.loc[1, "probability_of_target"] = 0.85
    sensitivity = R.natural_history_sensitivity(ATLAS5, 300)
    text = R.revival_report(
        ladder,
        sensitivity,
        R.surrogate_translation(600),
        R.transmurality_stratification(),
        frame,
        selection,
        cost_ledger(ATLAS5),
    )
    assert "still clears the bar" in text
    assert "NO LONGER ROBUST" not in text


def test_report_includes_the_expanded_search_result(search):
    frame, selection = search
    ladder = R.correction_ladder(ATLAS5, 300)
    sensitivity = R.natural_history_sensitivity(ATLAS5, 300)
    text = R.revival_report(
        ladder,
        sensitivity,
        R.surrogate_translation(600),
        R.transmurality_stratification(),
        frame,
        selection,
        cost_ledger(ATLAS5),
    )
    assert "Search: 256 composition-dose pairs across" in text
    assert "retry_redundancy|maturation at 4 doses" in text
    assert "NO composition is robust after the correction" in text


def test_report_branch_a_composition_survives(search):
    frame, selection = search
    robust = dict(selection)
    robust["any_robust"] = True
    robust["best_corrected"] = {
        "mechanisms": "retry_redundancy|maturation",
        "count": 2,
        "doses": 4,
        "corrected_probability_of_target": 0.83,
        "corrected_risk_adjusted": 0.79,
    }
    ladder = R.correction_ladder(ATLAS5, 300)
    sensitivity = R.natural_history_sensitivity(ATLAS5, 300)
    text = R.revival_report(
        ladder,
        sensitivity,
        R.surrogate_translation(600),
        R.transmurality_stratification(),
        frame,
        robust,
        cost_ledger(ATLAS5),
    )
    assert "At least one design clears the model's robustness threshold." in text


def test_report_branch_robust_at_every_comparator(search):
    frame, selection = search
    sensitivity = R.natural_history_sensitivity(ATLAS5, 300).copy()
    sensitivity["survives_correction"] = True
    text = R.revival_report(
        R.correction_ladder(ATLAS5, 300),
        sensitivity,
        R.surrogate_translation(600),
        R.transmurality_stratification(),
        frame,
        selection,
        cost_ledger(ATLAS5),
    )
    assert "Robust across the entire documented comparator range" in text


def test_report_branch_robust_at_no_comparator(search):
    frame, selection = search
    sensitivity = R.natural_history_sensitivity(ATLAS5, 300).copy()
    sensitivity["survives_correction"] = False
    text = R.revival_report(
        R.correction_ladder(ATLAS5, 300),
        sensitivity,
        R.surrogate_translation(600),
        R.transmurality_stratification(),
        frame,
        selection,
        cost_ledger(ATLAS5),
    )
    assert "Robust at NO documented comparator value" in text


def test_report_names_the_saturation_caveat_when_it_applies(search):
    frame, selection = search
    surrogate = R.surrogate_translation(600).copy()
    surrogate["corrected_saturated_fraction"] = 0.9
    text = R.revival_report(
        R.correction_ladder(ATLAS5, 300),
        R.natural_history_sensitivity(ATLAS5, 300),
        surrogate,
        R.transmurality_stratification(),
        frame,
        selection,
        cost_ledger(ATLAS5),
    )
    assert "of draws also sit at the" in text


def test_revival_study_ledger_is_tagged(tmp_path):
    R.run_revival_study(tmp_path / "revival", draws=FAST)
    ledger = pd.read_csv(tmp_path / "revival" / "tables" / "revival_cost_ledger.csv")
    assert (ledger["status"] == R.REVIVAL_TAG).all()
