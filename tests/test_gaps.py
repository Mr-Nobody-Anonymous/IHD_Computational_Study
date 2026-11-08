"""The gap register must classify honestly and quantify only what it can.

These tests exist to stop the register decaying into a list of excuses. A gap must be
classified, must name what would close it, and where a cost is quoted it must be the
value the arithmetic actually produces.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from ihd_solutions import gaps as G

# --------------------------------------------------------------- the register itself


def test_register_covers_every_class():
    frame = G.gap_register()
    assert set(frame["gap_class"]) == set(G.GAP_CLASSES)


def test_register_has_thirteen_gaps_and_no_duplicates():
    frame = G.gap_register()
    assert len(frame) == 13
    assert frame["key"].is_unique


def test_no_gap_is_left_unclassified_or_unclosable():
    frame = G.gap_register()
    assert frame["closes_with"].notna().all()
    assert (frame["closes_with"].str.len() >= 10).all()


def test_every_gap_states_a_question_not_a_topic():
    frame = G.gap_register()
    for question in frame["question"]:
        assert question.rstrip().endswith("?"), question


def test_only_epistemic_gaps_explain_why_a_number_cannot_settle_them():
    frame = G.gap_register()
    for row in frame.itertuples(index=False):
        if row.gap_class == "epistemic":
            assert isinstance(row.cannot_be_known_because, str), row.key
            assert len(row.cannot_be_known_because) > 30
        else:
            assert pd.isna(row.cannot_be_known_because), row.key


def test_no_gap_reports_an_infinite_cost():
    """Infinity is not a cost; it is a refusal to quantify, and must be None."""
    frame = G.gap_register()
    finite = frame["cost_value"].dropna()
    assert np.isfinite(finite).all()


def test_gap_dataclass_rejects_an_unknown_class():
    with pytest.raises(ValueError, match="unknown gap class"):
        G.Gap("k", "a question long enough to pass", "made_up", "why", "cost", 1.0, "u", "x")


def test_gap_dataclass_rejects_a_missing_closer():
    with pytest.raises(ValueError, match="must name what would close it"):
        G.Gap("k", "a question long enough to pass", "analytic", "why", "cost", 1.0, "u", "")


def test_gap_dataclass_rejects_a_topic_instead_of_a_question():
    with pytest.raises(ValueError, match="must state the question"):
        G.Gap("k", "molecular dynamics", "analytic", "why", "cost", 1.0, "u", "a gpu")


def test_gap_dataclass_rejects_infinite_cost():
    with pytest.raises(ValueError, match="not infinity"):
        G.Gap("k", "a question long enough to pass", "analytic", "why", "c", math.inf, "u", "x")


def test_summary_counts_match_the_register():
    frame = G.gap_register()
    summary = G.gap_summary(frame)
    assert summary["gaps_total"] == len(frame)
    assert sum(summary["by_class"].values()) == len(frame)
    assert summary["gaps_that_no_computation_can_close"] == len(frame) - 6


def test_summary_names_external_validation_as_the_first_thing_to_build():
    assert "external validation" in G.gap_summary()["first_thing_to_build"]


# ------------------------------------------------------------------ the open defect


def test_the_viability_table_defect_is_detected_not_papered_over():
    defect = G.transmurality_reviewer_check()
    assert defect["contradiction_detected"] is True
    assert defect["status"] == "OPEN_DEFECT_IN_THIS_STUDY"
    # 0.44 -> 0.43 across the mid-wall boundary is a 0.01 step against a 0.091 floor.
    assert defect["step_across_midwall"] == pytest.approx(0.01, abs=1e-12)
    assert defect["step_across_midwall"] < defect["step_required_exceeds"]


def test_the_defect_names_the_direction_of_repair_not_a_new_number():
    defect = G.transmurality_reviewer_check()
    assert "LOWER" in defect["required_direction_of_repair"]
    # It must not supply replacement values - that would be invention.
    assert "recovery_51_75" in defect
    assert defect["recovery_51_75"] == pytest.approx(0.43)


def test_the_defect_says_who_it_flatters():
    assert "intervention" in G.transmurality_reviewer_check()["who_this_flatters"]


# ------------------------------------------------------------- analytic quantities


def test_composite_adverse_event_risk_is_not_the_single_attempt_risk():
    result = G.composite_adverse_event_probability()
    assert result["composite_any_event_risk"] == pytest.approx(1.0 - (1.0 - 0.0550) ** 3, abs=1e-12)
    assert result["composite_any_event_risk"] > 3.0 * 0.0550 / 2.0
    assert result["inflation_vs_single"] == pytest.approx(2.83, abs=0.02)


def test_composite_adverse_event_risk_rejects_impossible_inputs():
    with pytest.raises(ValueError, match="at least one exposure"):
        G.composite_adverse_event_probability(exposures=0)
    with pytest.raises(ValueError, match="strictly inside"):
        G.composite_adverse_event_probability(per_attempt=1.0)


def test_composite_risk_is_monotone_in_exposures():
    risks = [
        G.composite_adverse_event_probability(exposures=n)["composite_any_event_risk"]
        for n in (1, 2, 3, 5)
    ]
    assert risks == sorted(risks)
    assert risks[0] == pytest.approx(0.0550)


def test_selection_produces_no_population_level_gain():
    result = G.response_guided_selection_bias(n=20_000)
    assert result["population_level_gain"] == 0.0
    assert result["status"] == "selection_is_not_a_population_gain"


def test_selection_readout_is_inflated_by_roughly_an_order_of_magnitude():
    result = G.response_guided_selection_bias(n=20_000)
    assert 5.0 < result["inflated_readout_factor"] < 15.0
    assert result["recoverable_room_selected_pp"] > result["recoverable_room_all_pp"]


def test_selection_effect_is_capped_by_the_oracle():
    result = G.response_guided_selection_bias(n=20_000)
    assert result["effect_in_selected_patients_pp"] <= result["oracle_selectable_ceiling_pp"]


def test_selection_is_seeded_and_reproducible():
    a = G.response_guided_selection_bias(n=5_000)
    b = G.response_guided_selection_bias(n=5_000)
    assert a["selection_gain_pp"] == b["selection_gain_pp"]


def test_selection_rejects_an_impossible_fraction():
    with pytest.raises(ValueError, match="strictly inside"):
        G.response_guided_selection_bias(n=100, selection_fraction=1.0)


def test_crispr_enumerable_counts_grow_with_mismatch_depth():
    frame = G.crispr_offtarget_enumerable_count()
    assert frame["expected_genome_wide_sites"].is_monotonic_increasing
    assert frame.loc[frame["max_mismatches"] == 0, "expected_genome_wide_sites"].iloc[0] < 1.0


def test_crispr_exact_match_is_essentially_unique_in_the_genome():
    """One 20-mer target site: the match itself, with nothing else expected to exist."""
    frame = G.crispr_offtarget_enumerable_count()
    exact = frame.loc[frame["max_mismatches"] == 0, "expected_genome_wide_sites"].iloc[0]
    assert exact < 0.01


def test_crispr_screen_covers_a_small_fraction_of_the_search_space():
    result = G.crispr_screen_completeness()
    assert result["coverage_fraction_at_3mm"] < 0.05
    assert result["can_it_report_SAFETY"] is False
    assert result["jaccard_of_two_lists"] == pytest.approx(0.4)


def test_crispr_guide_too_short_is_rejected():
    with pytest.raises(ValueError, match="guide too short"):
        G.crispr_offtarget_enumerable_count(guide_length=3)


def test_winner_curse_inflation_grows_with_the_grid():
    small = G.pick_best_of_k_measurement(10, measurement_sd=2.0, true_spread=8.0)
    large = G.pick_best_of_k_measurement(400, measurement_sd=2.0, true_spread=8.0)
    assert large["winner_inflation_pp"] > small["winner_inflation_pp"]
    assert large["winner_inflation_pp"] > 3.0


def test_winner_curse_is_zero_without_measurement_noise():
    result = G.pick_best_of_k_measurement(400, measurement_sd=0.0, true_spread=8.0)
    assert result["winner_inflation_pp"] == pytest.approx(0.0)


def test_winner_curse_rejects_bad_inputs():
    with pytest.raises(ValueError, match="at least one candidate"):
        G.pick_best_of_k_measurement(0, 1.0, 1.0)
    with pytest.raises(ValueError, match="cannot be negative"):
        G.pick_best_of_k_measurement(5, -1.0, 1.0)


def test_external_validation_needs_no_laboratory():
    result = G.external_benchmark_gap()
    assert result["current_rung"] == "convergence"
    assert result["cheapest_rung_above"] == "external"
    assert "GPU" in result["does_not_require"]
    assert "NOT DONE" in result["honest_status"]


# ------------------------------------------------------- the trial resolution floor


def test_each_patch_trial_observed_less_than_it_could_have_detected():
    frame = G.trial_resolution_floor()
    studies = frame.loc[frame["study"].notna()]
    patch = studies.loc[studies["study"].str.contains("BioVAT|AATS", na=False)]
    assert patch["effect_below_own_floor"].all()


def test_the_large_trial_is_not_flagged():
    frame = G.trial_resolution_floor()
    big = frame.loc[frame["study"].str.contains("REVIVED", na=False)].iloc[0]
    assert big["effect_below_own_floor"] is False or not big["effect_below_own_floor"]


def test_minimum_detectable_difference_follows_the_formula():
    frame = G.trial_resolution_floor(sd=5.0)
    row = frame.loc[frame["n_per_arm"] == 8].iloc[0]
    expected = G.POWER_Z_TOTAL * 5.0 * math.sqrt(2.0 / 8.0)
    assert row["min_detectable_difference_pp"] == pytest.approx(expected, rel=1e-9)


def test_required_trial_size_matches_the_floor():
    n = G.required_trial_size(3.9, sd=5.0)
    assert n == 26
    assert G.required_trial_size(3.4, sd=5.0) == 34


def test_required_trial_size_rejects_impossible_targets():
    with pytest.raises(ValueError, match="must be positive"):
        G.required_trial_size(0.0)
    with pytest.raises(ValueError, match="sd must be positive"):
        G.required_trial_size(3.0, sd=0.0)


def test_resolution_floor_rejects_bad_inputs():
    with pytest.raises(ValueError, match="sd must be positive"):
        G.trial_resolution_floor(sd=0.0)
    with pytest.raises(ValueError, match="power z"):
        G.trial_resolution_floor(power_z=0.0)


# ------------------------------------------------------------- the epistemic gaps


def test_adverse_event_blind_spot_reports_the_absent_denominator():
    result = G.adverse_event_blind_spot()
    assert result["reported_per_trial_counts"] is True
    assert result["reported_pooled_count"] is False
    assert result["total_patients_implanted"] == 28


def test_pooled_event_fraction_is_not_called_an_incidence():
    result = G.adverse_event_blind_spot()
    assert 0.0 <= result["pooled_any_event_fraction"] <= 1.0
    assert "not an incidence" in result["note"]


def test_causal_chain_identifies_the_missing_link_as_the_outcome():
    result = G.causal_gap()
    assert result["missing_link"] == "outcome"
    assert "does not supply the missing link" in result["why_simulation_cannot_repair_it"]
    assert "randomised" in result["what_would_repair_it"]


def test_causal_chain_only_licenses_a_hypothesis():
    assert "hypothesis" in G.causal_gap()["current_evidence_licenses"]


def test_prior_sensitivity_shows_the_recommendation_flipping():
    frame = G.prior_sensitivity()
    assert frame["conclusion_depends_on_prior"].all()
    assert frame["recommendation"].nunique() == 2


def test_prior_sensitivity_covers_the_full_grid():
    frame = G.prior_sensitivity()
    assert len(frame) == 4 * 4 * 4
    assert frame["prior_continue_weak"].min() == pytest.approx(0.20)
    assert frame["prior_continue_strong"].max() == pytest.approx(0.95)


def test_prior_weights_are_normalised_inside_each_cell():
    """Raw grid probabilities do not form a distribution; normalising stops the
    expected value from being inflated by whichever cell has the largest raw sum."""
    frame = G.prior_sensitivity()
    total = frame["normalised_weight_weak"] + frame["normalised_weight_strong"]
    assert (total < 1.0).all()
    assert (total > 0.0).all()


def test_continuing_carries_a_cost_so_the_sweep_can_flip():
    """Without a cost term 'more evidence is never bad' and the sweep is decorative."""
    assert G.PROGRAMME_COST_UNITS > 1.0
    assert G.PROGRAMME_COST_UNITS < 3.0


# ---------------------------------------------------------------------- reporting


def test_report_renders_every_section():
    text = G.gap_report()
    assert "WHAT THIS STUDY CANNOT ANSWER" in text
    assert "AN OPEN DEFECT IN THIS STUDY" in text
    assert "BELOW ITS OWN FLOOR" in text
    assert "patients per arm needed to resolve a given effect" in text
    # The discredited prompt: a count of four classes must all appear.
    for label in ("ANALYTIC", "IN PRACTICE", "EPISTEMIC", "IN PRINCIPLE"):
        assert label in text


def test_report_contains_no_nan_rows():
    text = G.gap_report()
    assert "nan" not in text.lower()


def test_disclaimer_says_nothing_is_improved_by_running_it():
    assert "No result in this repository is improved" in G.DISCLAIMER


def test_writer_emits_six_artefacts(tmp_path):
    summary = G.write_gap_report(tmp_path)
    assert summary["artefacts"] == 6
    assert summary["register_rows"] == 13
    for name in (
        "tables/gaps_register.csv",
        "tables/gaps_trial_resolution_floor.csv",
        "tables/gaps_crispr_enumerable_sites.csv",
        "tables/gaps_prior_sensitivity.csv",
        "tables/gaps_transmurality_defect.csv",
        "DISCLAIMER.txt",
        "gap_report.txt",
    ):
        assert (tmp_path / name).exists(), name


def test_written_register_round_trips(tmp_path):
    G.write_gap_report(tmp_path)
    frame = pd.read_csv(tmp_path / "tables" / "gaps_register.csv")
    assert len(frame) == 13
    assert int(frame["cost_value"].notna().sum()) == 7


def test_prior_sensitivity_refuses_to_report_a_single_recommendation(monkeypatch):
    """If the grid ever stops flipping, the sweep is decorative and must fail loudly."""
    monkeypatch.setattr(G, "PROGRAMME_COST_UNITS", 0.0)
    with pytest.raises(AssertionError, match="prior-dominated"):
        G.prior_sensitivity()


def test_report_skips_a_gap_class_with_no_members(monkeypatch):
    """The empty-class guard must exist even though all four classes are populated."""
    frame = G.gap_register()
    frame = pd.concat([frame.assign(gap_class="analytic")], ignore_index=True)
    monkeypatch.setattr(G, "gap_register", lambda: frame)
    monkeypatch.setattr(
        G,
        "gap_summary",
        lambda f=None: {
            "gaps_total": len(frame),
            "by_class": {"analytic": len(frame)},
            "gaps_with_a_quantified_cost": 7,
            "gaps_that_no_computation_can_close": 0,
            "first_thing_to_build": "x",
            "status": "y",
        },
    )
    text = G.gap_report()
    assert "in practice" not in text.lower()
    assert "epistemic" not in text.lower()
    assert "ANALYTIC" in text
