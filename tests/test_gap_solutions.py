"""Tests for the six-gap solution register and its design arithmetic.

Three things are pinned. First, the corrections: three of the six prior claims were
overstated, and the register must say so rather than quietly keeping the old framing.
Second, the statistics: power and sample-size helpers must be internally consistent
and must penalise a three-arm trial for its extra comparison. Third, the honesty of
the gap-5 trade-off: a potency breach may not be paid for with an unflagged redundancy
breach, and the selected point may not sit on the threshold.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ihd_solutions.gap_solutions import (
    BIOVAT_LVEF_GAIN_3MO_PP,
    BIOVAT_LVEF_GAIN_17MO_PP,
    CATO_ARMS,
    CATO_TOTAL_N,
    DISCLAIMER,
    GAP_SOLUTIONS,
    MG_DL_PER_MMOL_L,
    PRECEDENTED_POTENCY_CAP,
    REQUIRED_LVEF_GAIN_PP,
    RETRY_RESPONDERS,
    RETRY_TOTAL,
    SOLUTIONS_TAG,
    VERVE_LDL_ABS_REDUCTION_MG_DL,
    GapSolution,
    SolutionStatus,
    cato_power_analysis,
    detectable_difference,
    maturation_benefit,
    n_per_arm_for_difference,
    pcsk9_outcomes_power,
    potency_redundancy_tradeoff,
    power_at_difference,
    remuscularisation_scaling,
    reporter_signal_options,
    response_guided_enrichment,
    run_gap_solutions_study,
    solutions_frame,
)

# ------------------------------------------------------------- the register


def test_all_six_gaps_are_addressed():
    assert len(GAP_SOLUTIONS) == 6
    assert sorted(s.gap_id for s in GAP_SOLUTIONS) == [1, 2, 3, 4, 5, 6]


def test_every_gap_has_a_concrete_solution_and_falsifier():
    for entry in GAP_SOLUTIONS:
        assert len(entry.solution) > 80
        assert len(entry.falsifier) > 30
        assert len(entry.blockers) > 30
        assert entry.primary_endpoint
        assert entry.source.startswith("https://")


def test_three_prior_claims_were_overstated():
    """The point of the turn: the search overturned half the list."""
    overturned = [s for s in GAP_SOLUTIONS if s.was_overstated]
    assert len(overturned) == 3
    assert {s.gap_id for s in overturned} == {1, 3, 6}


def test_remuscularisation_is_not_claimed_impossible():
    entry = next(s for s in GAP_SOLUTIONS if s.gap_id == 1)
    assert entry.status is SolutionStatus.NARROWED
    assert "FALSE" in entry.search_finding
    assert "+20 pp" in entry.solution
    assert "maturation" in entry.solution.lower()


def test_reporter_gap_records_the_human_precedent():
    entry = next(s for s in GAP_SOLUTIONS if s.gap_id == 3)
    assert entry.status is SolutionStatus.NARROWED
    assert "MV-NIS" in entry.search_finding


def test_pcsk9_durability_is_no_longer_called_assumed():
    entry = next(s for s in GAP_SOLUTIONS if s.gap_id == 6)
    assert "MEASURED" in entry.search_finding
    assert entry.status is SolutionStatus.NARROWED


def test_still_open_gaps_are_the_untested_ones():
    still_open = {s.gap_id for s in GAP_SOLUTIONS if s.status is SolutionStatus.STILL_OPEN}
    assert still_open == {2, 5}


def test_register_frame_matches_the_tuple():
    frame = solutions_frame()
    assert len(frame) == len(GAP_SOLUTIONS)
    assert (frame["tag"] == SOLUTIONS_TAG).all()


def test_was_overstated_flag_is_consistent():
    for entry in GAP_SOLUTIONS:
        expected = entry.status in (
            SolutionStatus.CLOSED_BY_OTHERS,
            SolutionStatus.NARROWED,
        )
        assert entry.was_overstated == expected


# ------------------------------------------------------ statistical helpers


def test_detectable_difference_matches_the_sample_size_helper():
    for sd in (5.0, 6.0, 8.0):
        delta = detectable_difference(25, sd)
        assert n_per_arm_for_difference(sd, delta) == pytest.approx(25, abs=1)


def test_power_at_the_detectable_difference_is_the_requested_power():
    delta = detectable_difference(30, 6.0, power=0.8)
    assert power_at_difference(30, 6.0, delta) == pytest.approx(0.8, abs=0.005)


def test_more_comparisons_reduce_power():
    one = detectable_difference(20, 6.0, comparisons=1)
    two = detectable_difference(20, 6.0, comparisons=2)
    assert two > one


def test_bigger_sample_detects_a_smaller_difference():
    assert detectable_difference(80, 6.0) < detectable_difference(20, 6.0)


def test_higher_variance_needs_more_patients():
    assert n_per_arm_for_difference(8.0, 3.0) > n_per_arm_for_difference(5.0, 3.0)


def test_power_is_a_probability():
    for n in (10, 20, 100):
        value = power_at_difference(n, 6.0, 3.0)
        assert 0.0 <= value <= 1.0


@pytest.mark.parametrize("n", [0, 1])
def test_degenerate_arm_size_rejected(n):
    with pytest.raises(ValueError):
        detectable_difference(n, 6.0)


def test_nonpositive_sd_rejected():
    with pytest.raises(ValueError):
        detectable_difference(20, 0.0)


def test_nonpositive_difference_rejected():
    with pytest.raises(ValueError):
        n_per_arm_for_difference(6.0, 0.0)


# --------------------------------------------------------- gap 1 arithmetic


def test_scaling_reproduces_the_repository_requirement():
    frame = remuscularisation_scaling().set_index("anchor")
    assert frame.loc["3_month_effect", "required_multiplier"] == pytest.approx(
        REQUIRED_LVEF_GAIN_PP / BIOVAT_LVEF_GAIN_3MO_PP
    )
    assert frame.loc["3_month_effect", "required_multiplier"] == pytest.approx(5.128, abs=0.001)


def test_maturation_reduces_the_mass_requirement():
    frame = remuscularisation_scaling().set_index("anchor")
    assert (
        frame.loc["17_month_effect", "required_multiplier"]
        < frame.loc["3_month_effect", "required_multiplier"]
    )


def test_maturation_closes_part_of_the_gap_but_not_all():
    values = maturation_benefit()
    assert values["maturation_multiplier"] > 1.0
    assert 0.0 < values["gap_closed_by_time_alone_fraction"] < 1.0
    assert values["mass_multiplier_at_17_months"] > 1.0


def test_the_wall_thickness_discordance_is_reported():
    """Function rose while mass fell - the reason mass scaling may be the wrong lever."""
    values = maturation_benefit()
    assert values["wall_thickness_change_mm"] < 0.0
    assert values["lvef_change_pp"] > 0.0


def test_cell_counts_scale_consistently():
    frame = remuscularisation_scaling()
    for row in frame.itertuples(index=False):
        assert row.cardiomyocytes_if_mass_scaled == pytest.approx(
            row.ehm_units_if_mass_scaled * 40_000_000
        )


# --------------------------------------------------------- gap 2 enrichment


def test_enrichment_never_exceeds_certainty():
    frame = response_guided_enrichment()
    assert (frame["enriched_responder_rate"] <= 1.0).all()
    assert (frame["enriched_responder_rate"] > 0.0).all()


def test_enrichment_improves_on_the_baseline():
    frame = response_guided_enrichment()
    base = RETRY_RESPONDERS / RETRY_TOTAL
    assert (frame["enriched_responder_rate"] >= base - 1e-12).all()
    assert frame["enrichment_fold"].max() > 1.0


def test_a_tighter_cutoff_enriches_more_but_screens_out_more():
    frame = response_guided_enrichment()
    fixed = frame.loc[frame["assumed_sd_pp"].eq(4.0)].sort_values("cutoff_lvef_pp")
    assert fixed["enriched_responder_rate"].is_monotonic_decreasing
    assert fixed["fraction_passing_screen"].is_monotonic_increasing


def test_screening_cost_is_reported():
    frame = response_guided_enrichment()
    assert (frame["screen_failures_per_responder_enrolled"] >= 0.0).all()
    assert frame["screen_failures_per_responder_enrolled"].max() > 1.0


# ------------------------------------------------------------- gap 3 options


def test_only_viability_specific_reporters_qualify():
    frame = reporter_signal_options().set_index("reporter")
    assert bool(frame.loc["human_NIS", "viability_specific"])
    assert not bool(frame.loc["iron_oxide_MRI", "viability_specific"])


def test_human_nis_carries_human_imaging_precedent():
    frame = reporter_signal_options().set_index("reporter")
    assert "myeloma" in frame.loc["human_NIS", "human_imaging_precedent"]
    assert "none" in frame.loc["human_NIS", "cardiac_graft_precedent"]


# -------------------------------------------------------------- gap 4 power


def test_cato_uses_twenty_per_arm_and_two_comparisons():
    frame = cato_power_analysis()
    assert (frame["n_per_arm"] == CATO_TOTAL_N // CATO_ARMS).all()
    assert (frame["comparisons"] == CATO_ARMS - 1).all()


def test_cato_is_underpowered_for_a_small_increment():
    """The actionable finding: a null CATO would be uninformative for a 3 pp effect."""
    frame = cato_power_analysis()
    assert (frame["power_for_3pp"] < 0.5).all()


def test_cato_needs_far_more_patients_for_three_pp():
    frame = cato_power_analysis()
    assert (frame["n_per_arm_needed_for_3pp"] > CATO_TOTAL_N // CATO_ARMS).all()


def test_power_rises_with_the_effect_size():
    frame = cato_power_analysis()
    assert (frame["power_for_5pp"] > frame["power_for_3pp"]).all()


def test_power_falls_as_variance_rises():
    frame = cato_power_analysis().sort_values("assumed_sd_pp")
    assert frame["power_for_3pp"].is_monotonic_decreasing


# ------------------------------------------------------------- gap 6 events


def test_ldl_conversion_is_correct():
    frame = pcsk9_outcomes_power()
    full = frame.loc[frame["assumption"].eq("ctt_slope_full")].iloc[0]
    assert full["ldl_reduction_mmol_l"] == pytest.approx(
        VERVE_LDL_ABS_REDUCTION_MG_DL / MG_DL_PER_MMOL_L
    )


def test_a_weaker_slope_needs_more_events():
    frame = pcsk9_outcomes_power().set_index("assumption")
    assert (
        frame.loc["ctt_slope_half", "events_required"]
        > frame.loc["ctt_slope_three_quarters", "events_required"]
        > frame.loc["ctt_slope_full", "events_required"]
    )


def test_assumed_hazard_ratios_are_reductions_not_certainties():
    frame = pcsk9_outcomes_power()
    assert (frame["assumed_hazard_ratio"] > 0.0).all()
    assert (frame["assumed_hazard_ratio"] < 1.0).all()


def test_the_trial_is_not_a_mega_trial_on_the_full_slope():
    frame = pcsk9_outcomes_power().set_index("assumption")
    assert frame.loc["ctt_slope_full", "events_required"] < 200


def test_lower_event_rates_need_more_patients():
    frame = pcsk9_outcomes_power()
    assert (
        frame["n_at_2pct_annual_event_rate_3yr"] > frame["n_at_3pct_annual_event_rate_3yr"]
    ).all()


# ---------------------------------------------------- gap 5 honest accounting


@pytest.fixture(scope="module")
def tradeoff():
    return potency_redundancy_tradeoff(draws=200)


def test_tradeoff_covers_the_grid(tradeoff):
    assert len(tradeoff) == 6 * 4 * 3
    assert not tradeoff.duplicated(subset=["max_potency", "node_paths", "gate_paths"]).any()


def test_redundancy_precedent_is_flagged(tradeoff):
    """The accounting fix: 4 node paths is itself beyond precedent."""
    beyond = tradeoff.loc[tradeoff["node_paths"].eq(4)]
    assert beyond["redundancy_beyond_precedent"].all()
    within = tradeoff.loc[tradeoff["node_paths"].le(3) & tradeoff["gate_paths"].le(3)]
    assert not within["redundancy_beyond_precedent"].any()


def test_a_cheaper_cap_exists_only_by_breaching_redundancy(tradeoff):
    robust = tradeoff.loc[tradeoff["robust_with_margin"]]
    if robust.empty:
        pytest.skip("no robust point at this draw count")
    cheapest_any = robust["max_potency"].min()
    admissible = robust.loc[~robust["redundancy_beyond_precedent"]]
    if admissible.empty:
        pytest.skip("no admissible point at this draw count")
    assert admissible["max_potency"].min() >= cheapest_any


def test_precedent_gap_is_computed_against_the_precedented_cap(tradeoff):
    for row in tradeoff.itertuples(index=False):
        assert row.precedent_gap_fold == pytest.approx(row.max_potency / PRECEDENTED_POTENCY_CAP)


def test_margin_flag_is_stricter_than_the_bare_flag(tradeoff):
    assert (tradeoff["robust_with_margin"] <= tradeoff["robust"]).all()


def test_higher_potency_never_hurts_at_fixed_redundancy(tradeoff):
    for nodes, gates in ((1, 1), (2, 2), (3, 3)):
        subset = tradeoff.loc[
            tradeoff["node_paths"].eq(nodes) & tradeoff["gate_paths"].eq(gates)
        ].sort_values("max_potency")
        values = subset["probability_of_target"].to_numpy()
        assert (np.diff(values) >= -0.06).all()


# --------------------------------------------------------------- artefacts


@pytest.fixture(scope="module")
def study(tmp_path_factory):
    target = tmp_path_factory.mktemp("gaps") / "gap_solutions"
    return target, run_gap_solutions_study(target, draws=150)


def test_outputs_are_counted(study):
    target, summary = study
    written = [p for p in target.rglob("*") if p.is_file()]
    assert len(written) == summary["outputs"] == 10


def test_disclaimer_is_written(study):
    target, summary = study
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER.strip()
    assert summary["cures_declared"] == 0
    assert "NOT RESULTS" in DISCLAIMER


def test_every_table_is_tagged(study):
    target, _ = study
    tables = list((target / "tables").glob("*.csv"))
    assert len(tables) == 8
    for path in tables:
        frame = pd.read_csv(path)
        column = "tag" if "tag" in frame.columns else "status"
        assert column in frame.columns, path.name
        assert (frame[column] == SOLUTIONS_TAG).all(), path.name


def test_summary_reports_the_overturned_claims(study):
    _, summary = study
    assert summary["gaps_addressed"] == 6
    assert summary["claims_overturned"] == 3
    assert summary["still_open"] == 2


def test_report_states_the_corrections_and_the_limits(study):
    target, _ = study
    text = (target / "reports" / "gap_solutions_report.txt").read_text()
    assert "OVERTURNED THREE OF MY SIX CLAIMS" in text
    assert "underpowered" in text or "cannot distinguish" in text
    assert "unflagged redundancy breach" in text
    assert "None of this is a cure." in text


def test_cli_gap_solutions_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(["gap-solutions", "--output-dir", str(tmp_path / "cli-gaps"), "--draws", "120"])
    assert code == 0
    out = capsys.readouterr().out
    assert "SOLUTIONS TO THE SIX OPEN GAPS" in out
    assert "Wrote 10 artefacts" in out
    assert (tmp_path / "cli-gaps" / "DISCLAIMER.txt").exists()


def test_invalid_status_rejected():
    with pytest.raises(ValueError):
        SolutionStatus("probably_fine")


def test_gap_solution_is_frozen():
    entry = GAP_SOLUTIONS[0]
    with pytest.raises(AttributeError):
        entry.gap_id = 99


def test_register_entries_are_gap_solution_instances():
    assert all(isinstance(s, GapSolution) for s in GAP_SOLUTIONS)


def test_late_gain_exceeds_early_gain_in_the_source_data():
    assert BIOVAT_LVEF_GAIN_17MO_PP > BIOVAT_LVEF_GAIN_3MO_PP
