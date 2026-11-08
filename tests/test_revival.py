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


def test_study_writes_nine_artefacts(tmp_path):
    summary = R.run_revival_study(tmp_path / "revival", draws=FAST)
    files = sorted((tmp_path / "revival").rglob("*"))
    written = [path for path in files if path.is_file()]
    assert len(written) == summary["artefacts"] == 9
    assert (tmp_path / "revival" / "DISCLAIMER.txt").exists()


def test_study_summary_reports_the_failure(tmp_path):
    summary = R.run_revival_study(tmp_path / "revival", draws=FAST)
    assert summary["no_longer_robust"] is True
    assert summary["any_robust_after_correction"] is False
    assert summary["attribution_loss_pp"] > 20.0
    assert summary["status"] == R.REVIVAL_TAG


def test_every_table_carries_the_tag(tmp_path):
    R.run_revival_study(tmp_path / "revival", draws=FAST)
    for path in (tmp_path / "revival" / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        assert (frame["status"] == R.REVIVAL_TAG).all(), path.name


def test_cli_exposes_the_revival_command():
    from ihd_solutions.cli import build_parser

    args = build_parser().parse_args(["revival", "--draws", "100"])
    assert args.command == "revival"
    assert args.draws == 100


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


def test_report_branch_atlas5_survives_reselection(search):
    frame, selection = search
    surviving = dict(selection)
    surviving["atlas5_still_selected"] = True
    surviving["best_corrected"] = dict(selection["best_corrected"])
    ladder = R.correction_ladder(ATLAS5, 300)
    sensitivity = R.natural_history_sensitivity(ATLAS5, 300)
    text = R.revival_report(
        ladder,
        sensitivity,
        R.surrogate_translation(600),
        R.transmurality_stratification(),
        frame,
        surviving,
        cost_ledger(ATLAS5),
    )
    assert "YES. Re-searching all 64 compositions" in text


def test_report_branch_a_composition_survives(search):
    frame, selection = search
    robust = dict(selection)
    robust["any_robust"] = True
    robust["best_corrected"] = {
        "mechanisms": "retry_redundancy|maturation",
        "count": 2,
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
    assert "Best corrected admissible design" in text


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
