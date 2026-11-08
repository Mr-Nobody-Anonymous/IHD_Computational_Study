"""Tests for the precedent registry and ATLAS-4.

Four things are pinned. First, the evidence discipline: animal work and ongoing
trials can support a dose count but can never earn human efficacy credit. Second,
the correction: cumulative dosing beats best-of-k on the same draws, and the
independence requirement was a property of the wrong operator. Third, the honesty of
the 100% question: the rodent-optimistic design does read exactly 1.0 and must be
flagged and excluded, and the criteria that cannot reach 1.0 are capped by rules
rather than by an unexplored design space. Fourth, the artefact contract.
"""

from __future__ import annotations

from dataclasses import replace
from itertools import product

import numpy as np
import pandas as pd
import pytest

from ihd_solutions import precedent as P
from ihd_solutions.atlas3 import MAX_PRECEDENTED_ATTEMPTS as ATLAS3_LIMIT
from ihd_solutions.atlas4 import (
    ANIMAL_OPTIMISTIC,
    ATLAS4,
    ATLAS4_TAG,
    DISCLAIMER,
    NON_VIABILITY_SIGNAL_QUALITY,
    SPLIT_RETENTION,
    DosingMode,
    GraftDesign,
    correlation_by_operator,
    criterion_ceiling_audit,
    design_draws,
    evaluate_design,
    exhaustive_design_search,
    measurement_modality_falsifier,
    operator_contrast,
    run_atlas4_study,
    search_summary,
    translation_sensitivity,
)
from ihd_solutions.cure import CRITERION_THRESHOLD, MAX_SINGLE_AGENT_RISK_NORMALISATION
from ihd_solutions.gap_solver import SCAR_RECOVERABLE_LVEF_PP

# ----------------------------------------------------------- the registry


def test_registry_is_populated_and_sourced():
    assert len(P.PRECEDENTS) >= 10
    for entry in P.PRECEDENTS:
        assert entry.url.startswith("https://")
        assert entry.citation
        assert entry.finding
        assert entry.supports


def test_keys_are_unique():
    keys = [entry.key for entry in P.PRECEDENTS]
    assert len(keys) == len(set(keys))


def test_invalid_design_rejected():
    with pytest.raises(ValueError, match="unknown design"):
        P.Precedent(
            key="x",
            citation="c",
            url="https://example.org",
            design="press_release",
            population="p",
            doses_tested=1,
            interval_months=0.0,
            finding="f",
            supports=(P.SUPPORTS_DOSES,),
            efficacy_direction="positive",
        )


def test_invalid_direction_rejected():
    with pytest.raises(ValueError, match="unknown efficacy_direction"):
        P.Precedent(
            key="x",
            citation="c",
            url="https://example.org",
            design=P.HUMAN_RCT,
            population="p",
            doses_tested=1,
            interval_months=0.0,
            finding="f",
            supports=(P.SUPPORTS_DOSES,),
            efficacy_direction="promising",
        )


def test_precedent_without_a_lever_rejected():
    with pytest.raises(ValueError, match="at least one lever"):
        P.Precedent(
            key="x",
            citation="c",
            url="https://example.org",
            design=P.RODENT,
            population="p",
            doses_tested=1,
            interval_months=0.0,
            finding="f",
            supports=(),
            efficacy_direction="positive",
        )


@pytest.mark.parametrize("key", ["bolli_rat_2016", "bolli_rat_2018", "bolli_pig_2024"])
def test_animal_work_never_earns_human_efficacy_credit(key):
    entry = next(e for e in P.PRECEDENTS if e.key == key)
    assert not entry.is_human
    assert not entry.earns_human_efficacy_credit


def test_ongoing_trial_supports_doses_but_not_efficacy():
    cato = next(e for e in P.PRECEDENTS if e.key == "cato")
    assert cato.is_human
    assert cato.efficacy_direction == "ongoing"
    assert not cato.earns_human_efficacy_credit
    assert P.SUPPORTS_DOSES in cato.supports


def test_null_human_trial_earns_no_efficacy_credit():
    trial = next(e for e in P.PRECEDENTS if e.key == "vrtovec_remedium")
    assert trial.efficacy_direction == "null"
    assert not trial.earns_human_efficacy_credit


def test_wrong_population_positive_trial_is_still_labelled():
    hope = next(e for e in P.PRECEDENTS if e.key == "hope2")
    assert "NOT ischaemic heart disease" in hope.population


def test_precedented_dose_limit_is_four():
    assert P.max_precedented_doses() == 4


def test_atlas3_limit_has_been_corrected_to_the_registry():
    """The correction is now applied in code, not just documented."""
    assert ATLAS3_LIMIT == P.max_precedented_doses() == 4


def test_measurement_anchor_is_the_published_correlation():
    assert P.measurement_quality_anchor() == (0.76, 0.80)


def test_responder_fraction_is_the_human_observation():
    assert P.retry_responder_fraction() == pytest.approx(10.0 / 30.0)


def test_increment_shape_is_normalised_not_absolute():
    ratios = P.cumulative_increment_ratios()
    assert ratios[0] == pytest.approx(1.0)
    assert all(r > 0 for r in ratios)
    # Normalised so no rodent magnitude is imported as a human claim.
    assert max(ratios) < 2.0


def test_verdict_states_the_null_trial():
    text = P.human_repeat_efficacy_verdict()
    assert "NULL" in text
    assert "not established" in text


def test_precedent_frame_round_trips():
    frame = P.precedent_frame()
    assert len(frame) == len(P.PRECEDENTS)
    assert "earns_human_efficacy_credit" in frame.columns


# ------------------------------------------------------- the correction


def test_split_retention_matches_the_rat_arithmetic():
    assert pytest.approx((13.0 / 3.0) / 4.7) == SPLIT_RETENTION
    assert 0.9 < SPLIT_RETENTION < 1.0


def test_default_design_is_within_precedent():
    assert ATLAS4.doses == 4
    assert not ATLAS4.beyond_precedent


def test_default_uses_the_cumulative_operator():
    assert ATLAS4.mode is DosingMode.CUMULATIVE


def test_default_temperes_repeat_response_to_the_human_fraction():
    assert ATLAS4.repeat_responder_probability == pytest.approx(P.retry_responder_fraction())


def test_zero_doses_rejected():
    with pytest.raises(ValueError, match="at least one dose"):
        GraftDesign(doses=0)


def test_bad_mode_rejected():
    with pytest.raises(ValueError, match="mode must be a DosingMode"):
        GraftDesign(mode="cumulative")


@pytest.mark.parametrize(
    "field", ["failure_correlation", "repeat_responder_probability", "measurement_quality"]
)
def test_fraction_fields_validated(field):
    with pytest.raises(ValueError):
        GraftDesign(**{field: 1.2})


def test_draws_must_be_positive():
    with pytest.raises(ValueError):
        evaluate_design(ATLAS4, draws=0)


def test_operators_are_identical_at_one_dose():
    frame = operator_contrast(draws=1200).set_index(["doses", "mode"])
    assert frame.loc[(1, "best_of"), "probability_of_target"] == pytest.approx(
        frame.loc[(1, "cumulative"), "probability_of_target"], abs=1e-12
    )


def test_cumulative_beats_best_of_at_every_higher_dose_count():
    frame = operator_contrast(draws=1500).set_index(["doses", "mode"])
    for doses in (2, 3, 4, 6):
        assert (
            frame.loc[(doses, "cumulative"), "probability_of_target"]
            > frame.loc[(doses, "best_of"), "probability_of_target"]
        )


def test_human_tempered_design_is_robust_but_not_certain():
    result = evaluate_design(ATLAS4, 2500)
    assert result.robust_with_margin
    assert result.probability_of_target < 1.0
    assert not result.saturation_artefact


def test_gain_is_capped_by_the_recoverable_ceiling():
    gain, _, _ = design_draws(ATLAS4, 1200)
    assert gain.max() <= SCAR_RECOVERABLE_LVEF_PP + 1e-9
    assert gain.min() >= 0.0


# ------------------------------------------- the 100% question, answered


def test_rodent_optimism_produces_certainty_and_is_flagged():
    """This is the answer to 'keep going until everything is 100%'.

    Removing the human tempering does produce exactly 1.0. It is an artefact of
    assuming patients respond like rats, and the model says so.
    """
    result = evaluate_design(ANIMAL_OPTIMISTIC, 2000)
    assert result.probability_of_target == pytest.approx(1.0)
    assert result.saturation_artefact
    assert result.saturated_fraction > 0.0


def test_certainty_survives_even_a_harsh_translation_penalty():
    """The certainty is driven by the responder assumption, not by the biology."""
    harsh = replace(ANIMAL_OPTIMISTIC, translation_floor=0.4)
    assert evaluate_design(harsh, 1500).probability_of_target >= 0.99


def test_human_tempering_is_what_removes_the_certainty():
    tempered = evaluate_design(ATLAS4, 2000)
    optimistic = evaluate_design(ANIMAL_OPTIMISTIC, 2000)
    assert tempered.probability_of_target < optimistic.probability_of_target


def test_no_criterion_set_is_entirely_at_one():
    ceilings = criterion_ceiling_audit()
    assert not ceilings["already_at_100pct"].all()


def test_two_criteria_already_reach_one():
    ceilings = criterion_ceiling_audit().set_index("criterion")
    assert ceilings.loc["durable_driver_control", "catalogue_ceiling"] == pytest.approx(1.0)
    assert ceilings.loc["randomised_proof", "catalogue_ceiling"] == pytest.approx(1.0)


def test_the_unreachable_criteria_are_capped_by_rules():
    ceilings = criterion_ceiling_audit().set_index("criterion")
    assert ceilings.loc["myocardial_restoration", "catalogue_ceiling"] <= 0.5
    assert (
        ceilings.loc["risk_normalisation", "catalogue_ceiling"]
        <= MAX_SINGLE_AGENT_RISK_NORMALISATION
    )
    assert ceilings.loc["burden_reversal", "catalogue_ceiling"] <= 0.6


def test_every_criterion_carries_a_blocker():
    ceilings = criterion_ceiling_audit()
    assert len(ceilings) == 5
    for row in ceilings.itertuples(index=False):
        assert isinstance(row.blocker, str) and len(row.blocker) > 30


def test_threshold_flags_agree_with_the_numbers():
    ceilings = criterion_ceiling_audit()
    for row in ceilings.itertuples(index=False):
        assert row.meets_threshold == (row.catalogue_ceiling >= CRITERION_THRESHOLD)


# ------------------------------------------------------ the falsifiers


def test_non_viability_signals_gain_exactly_nothing():
    frame = measurement_modality_falsifier(draws=1500).set_index("modality")
    for name in ("iron_oxide_mri", "acute_retention_only"):
        assert frame.loc[name, "gain_vs_no_feedback_pp"] == pytest.approx(0.0, abs=1e-12)
        assert not bool(frame.loc[name, "viability_specific"])


def test_a_viability_signal_does_gain_something():
    frame = measurement_modality_falsifier(draws=2000).set_index("modality")
    assert frame.loc["reporter_gene_pet", "gain_vs_no_feedback_pp"] > 0.0
    assert bool(frame.loc["reporter_gene_pet", "viability_specific"])


def test_non_viability_quality_is_zero():
    assert NON_VIABILITY_SIGNAL_QUALITY == 0.0


def test_independence_requirement_is_operator_specific():
    frame = correlation_by_operator(draws=2000)
    retained = {}
    for mode in DosingMode:
        rows = frame.loc[frame["mode"].eq(mode.value)]
        free = float(rows.loc[rows["failure_correlation"].eq(0.0)].iloc[0]["gain_vs_single_pp"])
        shared = float(rows.loc[rows["failure_correlation"].eq(1.0)].iloc[0]["gain_vs_single_pp"])
        retained[mode] = shared / free
    assert retained[DosingMode.CUMULATIVE] > retained[DosingMode.BEST_OF]
    assert retained[DosingMode.CUMULATIVE] > 0.8


def test_cumulative_stays_robust_even_when_doses_are_correlated():
    frame = correlation_by_operator(draws=2000).set_index(["mode", "failure_correlation"])
    assert bool(frame.loc[("cumulative", 1.0), "robust"])


def test_headroom_gating_reduces_the_claim():
    gated = evaluate_design(replace(ATLAS4, headroom_gating=True), 2000)
    ungated = evaluate_design(replace(ATLAS4, headroom_gating=False), 2000)
    assert gated.probability_of_target <= ungated.probability_of_target


def test_translation_penalty_is_swept_not_assumed():
    frame = translation_sensitivity(draws=1500)
    fixed = frame.loc[frame["repeat_responder_probability"].eq(P.retry_responder_fraction())]
    ordered = fixed.sort_values("translation_floor")["probability_of_target"].to_numpy()
    assert (np.diff(ordered) >= -1e-12).all()


def test_more_procedural_risk_lowers_risk_adjusted_success():
    cheap = evaluate_design(replace(ATLAS4, per_dose_risk=0.0), 1500)
    dear = evaluate_design(replace(ATLAS4, per_dose_risk=0.15), 1500)
    assert dear.risk_adjusted_probability < cheap.risk_adjusted_probability
    assert cheap.cumulative_adverse_risk == pytest.approx(0.0)


# --------------------------------------------------- the exhaustive search


@pytest.fixture(scope="module")
def grid():
    return exhaustive_design_search(draws=400)


def test_grid_covers_the_full_cartesian_product(grid):
    assert len(grid) == 2 * 5 * 3 * 3 * 2 * 2
    keys = [
        "mode",
        "doses",
        "failure_correlation",
        "feedback_learning",
        "headroom_gating",
        "repeat_responder_probability",
    ]
    assert not grid.duplicated(subset=keys).any()
    expected = set(
        product(
            [m.value for m in DosingMode],
            (1, 2, 3, 4, 6),
            (0.0, 0.5, 1.0),
            (0.0, 0.5, 1.0),
            (True, False),
            (1.0 / 3.0, 1.0),
        )
    )
    assert set(grid[keys].itertuples(index=False, name=None)) == expected


def test_grid_keeps_the_saturated_points_visible(grid):
    assert grid["saturation_artefact"].any()
    assert (grid.loc[grid["saturation_artefact"], "probability_of_target"] >= 1.0).all()


def test_summary_counts_are_consistent(grid):
    summary = search_summary(grid)
    assert summary["configurations_tested"] == len(grid)
    assert summary["robust"] == int(grid["robust"].sum())
    assert summary["saturation_artefacts"] == int(grid["saturation_artefact"].sum())
    assert summary["admissible"] <= summary["robust_with_margin"]


def test_admissible_points_exclude_saturation_and_optimism(grid):
    summary = search_summary(grid)
    assert summary["admissible"] > 0
    best = summary["best_admissible"]
    assert best is not None
    assert best["probability_of_target"] < 1.0
    assert best["doses"] <= P.max_precedented_doses()
    assert best["repeat_responder_probability"] <= P.retry_responder_fraction() + 1e-9


def test_the_best_point_is_not_the_highest_probability_point(grid):
    """The filter must actually bind, or it is decoration."""
    summary = search_summary(grid)
    assert summary["best_admissible"]["probability_of_target"] < grid["probability_of_target"].max()


# ------------------------------------------------------------- artefacts


@pytest.fixture(scope="module")
def study(tmp_path_factory):
    target = tmp_path_factory.mktemp("atlas4") / "atlas4"
    return target, run_atlas4_study(target, draws=400)


def test_outputs_are_counted(study):
    target, summary = study
    written = [p for p in target.rglob("*") if p.is_file()]
    assert len(written) == summary["outputs"] == 11


def test_disclaimer_is_written(study):
    target, summary = study
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER.strip()
    assert summary["cures_declared"] == 0
    assert "null" in DISCLAIMER


def test_every_table_is_tagged(study):
    target, _ = study
    tables = list((target / "tables").glob("*.csv"))
    assert len(tables) == 9
    for path in tables:
        frame = pd.read_csv(path)
        assert "status" in frame.columns, path.name
        assert (frame["status"] == ATLAS4_TAG).all(), path.name


def test_report_answers_the_hundred_percent_question(study):
    target, _ = study
    text = (target / "reports" / "atlas4_report.txt").read_text()
    assert "KEEP GOING UNTIL ALL ARE 100%" in text
    assert "saturation artefact" in text
    assert "No cure is declared." in text
    assert "PROPOSED in print in 2019, still NOT TESTED" in text


def test_summary_reports_the_correction(study):
    _, summary = study
    assert summary["precedented_dose_limit"] == 4
    assert summary["configurations_tested"] == 360
    assert summary["saturation_artefacts"] > 0
    assert summary["animal_optimistic_probability"] > summary["probability_of_target"]


def test_cli_atlas4_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(["atlas4", "--output-dir", str(tmp_path / "cli-atlas4"), "--draws", "200"])
    assert code == 0
    out = capsys.readouterr().out
    assert "ATLAS-4" in out
    assert "Wrote 11 artefacts" in out
    assert (tmp_path / "cli-atlas4" / "DISCLAIMER.txt").exists()
