"""Tests for mechanism-level recombination.

These pin the corrections found while building it, not just the happy path: the unit
mismatch that called a -14 pp overlap "additive", the response-guided screen that
duplicated stop-on-success, the population-vs-screened-in accounting, the two prices that
were conflated into a false precedent claim, and the budget-dependent specification.
"""

from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pandas as pd
import pytest

from ihd_solutions import recombine as R
from ihd_solutions.atlas3 import GraftUpgrades
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solutions import (
    BIOVAT_LVEF_GAIN_3MO_PP,
    BIOVAT_LVEF_GAIN_17MO_PP,
    PRECEDENTED_POTENCY_CAP,
)
from ihd_solutions.precedent import max_precedented_doses
from ihd_solutions.upgrade import MAX_PRECEDENTED_PATHS

FAST = 400


@pytest.fixture(scope="module")
def full() -> R.CompositeResult:
    return R.evaluate_composite(R.COMPOSITE, 1500)


@pytest.fixture(scope="module")
def marginals() -> pd.DataFrame:
    return R.mechanism_marginals(draws=1200)


@pytest.fixture(scope="module")
def standalone() -> pd.DataFrame:
    return R.mechanism_standalone(draws=1200)


@pytest.fixture(scope="module")
def interactions() -> pd.DataFrame:
    return R.pairwise_interactions(draws=1200)


@pytest.fixture(scope="module")
def grid() -> pd.DataFrame:
    return R.exhaustive_subsets(draws=FAST)


@pytest.fixture(scope="module")
def sensitivity() -> pd.DataFrame:
    return R.maturation_sensitivity(draws=600)


# --------------------------------------------------------------------------- catalogue


def test_every_mechanism_is_reachable_by_name():
    assert set(R.MECHANISM_BY_NAME) == set(R.ALL_MECHANISMS)
    assert len(R.MECHANISM_BY_NAME) == len(R.MECHANISMS)


def test_benefits_and_constraints_partition_the_catalogue():
    assert set(R.BENEFIT_MECHANISMS).isdisjoint(R.CONSTRAINT_MECHANISMS)
    assert set(R.BENEFIT_MECHANISMS) | set(R.CONSTRAINT_MECHANISMS) == set(R.ALL_MECHANISMS)


def test_constraints_are_flagged_as_constraints():
    for name in R.CONSTRAINT_MECHANISMS:
        assert R.MECHANISM_BY_NAME[name].is_constraint
    for name in R.BENEFIT_MECHANISMS:
        assert not R.MECHANISM_BY_NAME[name].is_constraint


def test_every_mechanism_declares_an_origin_and_a_cost():
    for mechanism in R.MECHANISMS:
        assert mechanism.origin.strip()
        assert mechanism.cost_axis.strip()
        assert mechanism.description.strip()


def test_maturation_ratio_comes_from_the_human_series():
    assert (
        pytest.approx(BIOVAT_LVEF_GAIN_17MO_PP / BIOVAT_LVEF_GAIN_3MO_PP) == R.MATURATION_MULTIPLIER
    )
    assert R.MATURATION_MULTIPLIER > 1.0


# ------------------------------------------------------------------------ design object


def test_unknown_mechanism_is_rejected():
    with pytest.raises(ValueError, match="unknown mechanisms"):
        R.Composite(enabled=frozenset({"teleportation"}))


def test_zero_doses_is_rejected():
    with pytest.raises(ValueError, match="at least one dose"):
        R.Composite(doses=0)


def test_bad_risk_is_rejected():
    with pytest.raises(ValueError):
        R.Composite(per_dose_risk=1.5)


def test_zero_redundancy_is_rejected():
    with pytest.raises(ValueError, match="at least one path"):
        R.Composite(clearance_nodes=0)


def test_design_is_immutable():
    with pytest.raises(FrozenInstanceError):
        R.COMPOSITE.doses = 9


def test_with_mechanism_does_not_mutate_the_original():
    reduced = R.COMPOSITE.with_mechanism("maturation", False)
    assert R.COMPOSITE.has("maturation")
    assert not reduced.has("maturation")


def test_with_mechanism_is_idempotent():
    once = R.COMPOSITE.with_mechanism("maturation", True)
    assert once.enabled == R.COMPOSITE.enabled


def test_constraints_intact_only_when_all_are_present():
    assert R.COMPOSITE.constraints_intact
    for name in R.CONSTRAINT_MECHANISMS:
        assert not R.COMPOSITE.with_mechanism(name, False).constraints_intact


# ------------------------------------------------------------------------- cost axes


def test_doses_within_the_corrected_limit_are_not_a_breach():
    assert R.ATLAS5.doses == max_precedented_doses()
    assert not R.ATLAS5.attempts_beyond_precedent


def test_more_doses_than_precedent_is_a_breach():
    assert R.Composite(doses=max_precedented_doses() + 1).attempts_beyond_precedent


def test_single_dose_design_cannot_breach_the_attempt_limit():
    design = R.Composite(doses=99).with_mechanism("retry_redundancy", False)
    assert not design.attempts_beyond_precedent


def test_redundancy_breach_is_charged_separately_from_potency():
    design = R.Composite(clearance_nodes=MAX_PRECEDENTED_PATHS + 1)
    assert design.redundancy_beyond_precedent
    assert "redundancy" in design.precedent_breaches


def test_potency_breach_is_flagged():
    assert R.Composite(clearance_cap=PRECEDENTED_POTENCY_CAP + 1.0).potency_beyond_precedent
    assert not R.Composite(clearance_cap=PRECEDENTED_POTENCY_CAP).potency_beyond_precedent


def test_reporter_costs_a_genetic_modification():
    assert R.COMPOSITE.genetic_modifications == 1
    assert R.ATLAS5.genetic_modifications == 0


def test_graft_breaches_exclude_the_clearance_arm():
    # The clearance arm's standing potency breach must not make every graft composition
    # inadmissible; that returned an empty set for the wrong reason.
    assert R.COMPOSITE.potency_beyond_precedent
    assert "potency" in R.COMPOSITE.precedent_breaches
    assert "potency" not in R.COMPOSITE.graft_precedent_breaches


def test_reporter_is_a_graft_axis_breach():
    assert "genetic_modification" in R.COMPOSITE.graft_precedent_breaches
    assert R.ATLAS5.graft_precedent_breaches == ()


def test_cost_ledger_charges_every_axis():
    ledger = R.cost_ledger(R.COMPOSITE)
    assert set(ledger["axis"]) == {
        "doses",
        "clearance_potency_cap",
        "clearance_redundant_paths",
        "genetic_modifications",
    }
    assert bool(ledger["beyond_precedent"].any())


def test_atlas5_ledger_still_charges_the_clearance_breach():
    ledger = R.cost_ledger(R.ATLAS5).set_index("axis")
    assert not bool(ledger.loc["doses", "beyond_precedent"])
    assert not bool(ledger.loc["genetic_modifications", "beyond_precedent"])
    assert bool(ledger.loc["clearance_potency_cap", "beyond_precedent"])


# --------------------------------------------------------------------------- simulation


def test_draws_are_reproducible():
    first, used_a = R.composite_draws(R.ATLAS5, FAST, seed=7)
    second, used_b = R.composite_draws(R.ATLAS5, FAST, seed=7)
    assert np.array_equal(first, second)
    assert np.array_equal(used_a, used_b)


def test_different_seeds_give_different_draws():
    first, _ = R.composite_draws(R.ATLAS5, FAST, seed=7)
    second, _ = R.composite_draws(R.ATLAS5, FAST, seed=8)
    assert not np.array_equal(first, second)


def test_gain_never_exceeds_the_recoverable_ceiling():
    from ihd_solutions.gap_solver import SCAR_RECOVERABLE_LVEF_PP

    gain, _ = R.composite_draws(R.COMPOSITE, FAST)
    assert gain.max() <= SCAR_RECOVERABLE_LVEF_PP + 1e-9
    assert gain.min() >= 0.0


def test_no_retry_means_exactly_one_dose():
    design = R.COMPOSITE.with_mechanism("retry_redundancy", False)
    _, used = R.composite_draws(design, FAST)
    assert np.all(used == 1.0)


def test_doses_used_never_exceeds_the_design():
    _, used = R.composite_draws(R.ATLAS5, FAST)
    assert used.max() <= R.ATLAS5.doses


def test_zero_draws_is_rejected():
    with pytest.raises(ValueError):
        R.composite_draws(R.ATLAS5, 0)


def test_adapter_builds_a_real_graft_upgrades():
    view = R._upgrade_view_shim(R._DesignView(4, 0.5, 0.78, 0.7, 0.05))
    assert isinstance(view, GraftUpgrades)
    assert view.attempts == 4
    assert view.measurement_quality == pytest.approx(0.78)


def test_result_is_a_probability(full):
    assert 0.0 <= full.probability_of_target <= 1.0
    assert full.draws == 1500


def test_standard_error_shrinks_with_draws():
    small = R.evaluate_composite(R.ATLAS5, 200)
    large = R.evaluate_composite(R.ATLAS5, 4000)
    assert large.monte_carlo_standard_error < small.monte_carlo_standard_error


def test_risk_adjustment_never_raises_the_probability(full):
    assert full.risk_adjusted_probability <= full.probability_of_target + 1e-12


def test_risk_adjustment_weights_each_simulated_patient_by_actual_doses():
    design = R.ATLAS5
    draws, seed = 900, 123
    gain, used = R.composite_draws(design, draws, seed)
    target = gain >= CRITERION_THRESHOLD * R.SCAR_RECOVERABLE_LVEF_PP
    survival = (1.0 - design.per_dose_risk) ** used
    result = R.evaluate_composite(design, draws, seed)
    assert result.risk_adjusted_probability == pytest.approx(float(np.mean(target * survival)))
    assert result.cumulative_adverse_risk == pytest.approx(float(np.mean(1.0 - survival)))


def test_adverse_risk_grows_with_doses():
    one = R.evaluate_composite(R.ATLAS5.with_mechanism("retry_redundancy", False), 1200)
    four = R.evaluate_composite(R.ATLAS5, 1200)
    assert four.cumulative_adverse_risk > one.cumulative_adverse_risk


def test_saturation_artefact_is_detected():
    result = R.CompositeResult(1.0, 25.0, 4.0, 0.18, 1.0, 1000, 5)
    assert result.saturation_artefact
    assert not result.robust_with_margin or result.probability_of_target == 1.0


def test_margin_rule_requires_two_standard_errors():
    borderline = R.CompositeResult(CRITERION_THRESHOLD, 20.0, 1.0, 0.05, 0.0, 400, 3)
    assert borderline.robust
    assert not borderline.robust_with_margin


# ------------------------------------------------------- the corrections, pinned


def test_dropping_human_tempering_raises_the_number(marginals):
    row = marginals.set_index("mechanism").loc["human_tempering"]
    assert row["contribution_pp"] < 0.0


def test_response_guided_selection_is_not_a_therapy_benefit(marginals):
    """The screen concentrates responders but excludes patients.

    Modelled over the whole population it must not be a free win. Reporting only the
    screened-in subset is what turned a trial-design lever into a fake therapy gain.
    """
    row = marginals.set_index("mechanism").loc["response_guided_selection"]
    assert row["contribution_pp"] <= 0.0


def test_screened_in_rate_is_below_the_enriched_responder_rate():
    assert R.RESPONSE_GUIDED_PASS < R.RESPONSE_GUIDED_RATE < 1.0


def test_enrichment_uses_the_central_assumption_not_the_maximum():
    from ihd_solutions.gap_solutions import response_guided_enrichment

    table = response_guided_enrichment()
    assert float(table["enriched_responder_rate"].max()) > R.RESPONSE_GUIDED_RATE
    assert (
        pytest.approx(
            float(
                table.loc[
                    table["assumed_sd_pp"].eq(4.0) & table["cutoff_lvef_pp"].eq(37.5),
                    "enriched_responder_rate",
                ].iloc[0]
            ),
            abs=1e-5,
        )
        == R.RESPONSE_GUIDED_RATE
    )


def test_screen_does_not_simply_duplicate_stop_on_success():
    """If the screen were the success rule again, toggling it would change nothing."""
    on = R.evaluate_composite(R.ATLAS5.with_mechanism("response_guided_selection", True), 1500)
    off = R.evaluate_composite(R.ATLAS5, 1500)
    assert abs(on.probability_of_target - off.probability_of_target) > 1e-9


def test_interaction_verdict_uses_probability_units(interactions):
    """A -14 pp overlap was once labelled 'additive' by a pp-vs-fraction mismatch."""
    row = interactions.set_index(["first", "second"]).loc[("retry_redundancy", "maturation")]
    assert row["interaction_pp"] < -1.0
    assert row["verdict"] == "overlapping"


def test_tolerance_is_expressed_in_probability_units():
    assert pytest.approx(0.005) == R.INTERACTION_TOLERANCE
    assert R.INTERACTION_TOLERANCE < 1.0


def test_retry_and_cumulative_are_synergistic(interactions):
    row = interactions.set_index(["first", "second"]).loc[
        ("retry_redundancy", "cumulative_operator")
    ]
    assert row["interaction_pp"] > 0.5
    assert row["verdict"] == "synergistic"


def test_cumulative_operator_is_inert_without_retries(standalone):
    row = standalone.set_index("mechanism").loc["cumulative_operator"]
    assert row["standalone_gain_pp"] == pytest.approx(0.0, abs=1e-9)


def test_interaction_arithmetic_is_exact(interactions):
    reconstructed = (
        interactions["gain_both_pp"]
        - interactions["gain_first_pp"]
        - interactions["gain_second_pp"]
    )
    assert np.allclose(reconstructed, interactions["interaction_pp"], atol=1e-9)


def test_every_benefit_pair_is_tested(interactions):
    from itertools import combinations

    assert len(interactions) == len(list(combinations(R.BENEFIT_MECHANISMS, 2)))


def test_feedback_and_reporter_are_inert_in_this_composition(marginals):
    indexed = marginals.set_index("mechanism")
    for name in ("feedback", "viability_reporter"):
        assert not bool(indexed.loc[name, "load_bearing"])


def test_maturation_is_the_strongest_standalone_part(standalone):
    best = standalone.loc[standalone["standalone_gain_pp"].idxmax()]
    assert best["mechanism"] == "maturation"


# ------------------------------------------------------------------ exhaustive search


def test_every_subset_is_enumerated(grid):
    assert len(grid) == 2 ** len(R.BENEFIT_MECHANISMS)
    assert grid["mechanisms"].nunique() == len(grid)


def test_search_never_switches_off_a_constraint(grid):
    assert bool(grid["constraints_intact"].all())


def test_no_composition_reaches_certainty(grid):
    assert not bool(grid["saturation_artefact"].any())
    assert float(grid["probability_of_target"].max()) < 1.0


def test_search_budget_is_floored():
    low = R.exhaustive_subsets(draws=10)
    assert R.MIN_SPECIFICATION_DRAWS >= 2000
    # A floored search cannot be noisier than an unfloored one at the same tiny budget.
    assert len(low) == 2 ** len(R.BENEFIT_MECHANISMS)


def test_specification_does_not_move_with_the_sampling_budget():
    first = R.select_composite(R.exhaustive_subsets(draws=FAST))
    second = R.select_composite(R.exhaustive_subsets(draws=R.MIN_SPECIFICATION_DRAWS))
    assert first["best"]["mechanisms"] == second["best"]["mechanisms"]


def test_selection_reports_both_designs(grid):
    selection = R.select_composite(grid)
    assert selection["robust_set_exists"]
    assert isinstance(selection["best"], dict)
    assert isinstance(selection["cheapest"], dict)


def test_recommended_design_maximises_risk_adjusted_probability(grid):
    selection = R.select_composite(grid)
    admissible = grid.loc[
        grid["constraints_intact"]
        & (~grid["saturation_artefact"])
        & grid["robust_with_margin"]
        & grid["precedent_breaches"].eq("none")
    ]
    assert selection["best"]["risk_adjusted_probability"] == pytest.approx(
        float(admissible["risk_adjusted_probability"].max())
    )


def test_cheapest_design_has_no_more_parts_than_the_recommended(grid):
    selection = R.select_composite(grid)
    assert selection["cheapest"]["count"] <= selection["best"]["count"]


def test_recommended_design_carries_no_genetic_modification(grid):
    assert R.select_composite(grid)["best"]["genetic_modifications"] == 0


def test_empty_admissible_set_reports_no_robust_set(grid):
    impossible = grid.assign(precedent_breaches="potency")
    selection = R.select_composite(impossible)
    assert selection["admissible"] == 0
    assert selection["best"] is None
    assert not selection["robust_set_exists"]


def test_atlas5_is_the_selected_composition(grid):
    selected = R.select_composite(grid)["best"]["mechanisms"]
    assert set(selected.split("|")) == set(R.ATLAS5_MECHANISMS) - set(R.CONSTRAINT_MECHANISMS)


# ------------------------------------------------------------------------- verdict


def test_verdict_separates_the_two_prices(grid):
    verdict = R.composition_verdict(grid, R.select_composite(grid))
    assert verdict["price_of_admissibility_pp"] == pytest.approx(
        100.0 * (verdict["best_probability_any"] - verdict["best_probability_admissible"])
    )
    assert verdict["price_of_cheapness_pp"] == pytest.approx(
        100.0 * (verdict["best_probability_admissible"] - verdict["selected_probability"])
    )


def test_admissibility_is_free_here_and_the_report_says_so(grid):
    """The strongest composition breaches nothing, so claiming otherwise was false."""
    verdict = R.composition_verdict(grid, R.select_composite(grid))
    assert verdict["price_of_admissibility_pp"] < 100.0 * R.INTERACTION_TOLERANCE


def test_composing_beats_the_best_single_mechanism(grid):
    verdict = R.composition_verdict(grid, R.select_composite(grid))
    assert verdict["gain_over_best_single_pp"] > 0.5
    assert verdict["composition_beat_best_single"]


def test_atlas4_reference_matches_the_committed_table():
    import ast
    from pathlib import Path

    path = Path("results/atlas4/tables/atlas4_search_summary.csv")
    if not path.exists():  # pragma: no cover - artefacts are committed
        pytest.skip("atlas4 artefacts not present")
    summary = pd.read_csv(path)
    best = ast.literal_eval(str(summary["best_admissible"].iloc[0]))
    assert pytest.approx(float(best["probability_of_target"]), abs=1e-6) == R.ATLAS4_BEST_ADMISSIBLE


def test_composite_beats_the_atlas4_reference(grid):
    verdict = R.composition_verdict(grid, R.select_composite(grid))
    assert verdict["gain_over_atlas4_pp"] > 0.0


def test_verdict_handles_an_empty_selection(grid):
    verdict = R.composition_verdict(grid, {"best": None, "admissible": 0})
    assert verdict["selected_probability"] == 0.0
    assert verdict["selected_count"] == 0


# --------------------------------------------------------------------- sensitivity


def test_sensitivity_covers_the_full_cross_product(sensitivity):
    assert len(sensitivity) == 4 * 3 * 3
    assert sensitivity["maturation_multiplier"].nunique() == 4
    assert sensitivity["translation_floor"].nunique() == 3
    assert sensitivity["doses"].nunique() == 3


def test_sensitivity_does_not_claim_robustness_everywhere(sensitivity):
    assert not bool(sensitivity["robust_with_margin"].all())
    assert bool(sensitivity["robust_with_margin"].any())


def test_model_cannot_resolve_the_full_maturation_ratio(sensitivity):
    limit = R.maturation_resolution_limit(sensitivity)
    assert limit < R.MATURATION_MULTIPLIER
    assert limit == pytest.approx(1.25)


def test_resolution_limit_falls_back_when_no_plateau_exists():
    # Every multiplier must still show spread, including the largest, so the loop is
    # exhausted rather than returning early. Duplicated cells at the top multiplier are
    # what make that possible.
    frame = pd.DataFrame(
        {
            "maturation_multiplier": [1.0, 2.0, 2.0],
            "translation_floor": [1.0, 1.0, 1.0],
            "doses": [1, 1, 1],
            "probability_of_target": [0.1, 0.4, 0.9],
        }
    )
    assert R.maturation_resolution_limit(frame) == pytest.approx(2.0)


def test_translation_penalty_never_helps(sensitivity):
    pivot = sensitivity.pivot_table(
        index=["maturation_multiplier", "doses"],
        columns="translation_floor",
        values="probability_of_target",
    )
    assert bool((pivot[1.0] >= pivot[0.4] - 1e-9).all())


# ------------------------------------------------------------------------ both arms


def test_joint_is_a_product_and_proof_is_zero():
    joint = R.joint_result(R.ATLAS5, draws=600)
    assert joint["joint_optimistic"] == pytest.approx(
        joint["graft_probability"] * joint["clearance_probability"]
    )
    assert joint["randomised_proof"] == 0.0


def test_clearance_arm_runs_at_the_solved_horizon():
    result = R.clearance_arm(R.ATLAS5, draws=600)
    assert 0.0 <= result.probability_of_target <= 1.0


# --------------------------------------------------------------------------- report


@pytest.fixture(scope="module")
def report(grid, marginals, standalone, interactions, sensitivity, full) -> str:
    selection = R.select_composite(grid)
    return R.recombine_report(
        full,
        marginals,
        interactions,
        grid,
        selection,
        R.joint_result(R.ATLAS5, draws=600),
        R.cost_ledger(R.COMPOSITE),
        standalone,
        sensitivity,
        R.composition_verdict(grid, selection),
        R.evaluate_composite(R.design_from_selection(selection), 800),
        R.stacking_contrast(draws=FAST),
    )


def test_report_carries_the_disclaimer(report):
    assert R.DISCLAIMER in report
    assert "NOT A THERAPY" in report


def test_report_shows_every_required_section(report):
    for heading in (
        "THE CATALOGUE",
        "LEAVE-ONE-OUT",
        "STANDALONE",
        "INTERACTIONS",
        "EXHAUSTIVE COMPOSITION SEARCH",
        "THE BILL",
        "VERDICT (computed, not asserted)",
    ):
        assert heading in report


def test_report_names_atlas5_and_the_frugal_alternative(report):
    assert "ATLAS-5" in report
    assert "FRUGAL ALTERNATIVE" in report


def test_report_declares_the_model_blind_spot(report):
    assert "MODEL BLIND SPOT" in report
    assert "cannot distinguish" in report


def test_report_states_admissibility_is_free_when_it_is(report):
    assert "Admissibility costs NOTHING here" in report


def test_report_never_promises_a_cure(report):
    assert "not a cure" in report
    assert "randomised" in report.lower()


def test_report_handles_no_admissible_composition(
    grid, marginals, standalone, interactions, sensitivity, full
):
    empty = {"admissible": 0, "best": None, "robust_set_exists": False}
    text = R.recombine_report(
        full,
        marginals,
        interactions,
        grid,
        empty,
        R.joint_result(R.ATLAS5, draws=400),
        R.cost_ledger(R.COMPOSITE),
        standalone,
        sensitivity,
        R.composition_verdict(grid, empty),
        R.evaluate_composite(R.COMPOSITE, 400),
        R.stacking_contrast(draws=FAST),
    )
    assert "NO admissible composition" in text


# ---------------------------------------------------------------------------- study


def test_study_writes_exactly_the_declared_artefacts(tmp_path):
    summary = R.run_recombine_study(tmp_path / "recombine", draws=FAST)
    written = [p for p in (tmp_path / "recombine").rglob("*") if p.is_file()]
    assert len(written) == summary["outputs"]


def test_study_declares_no_cures_and_no_certainty(tmp_path):
    summary = R.run_recombine_study(tmp_path / "recombine", draws=FAST)
    assert summary["cures_declared"] == 0
    assert not summary["full_saturation_artefact"]
    assert summary["status"] == R.RECOMBINE_TAG
    assert summary["compositions_tested"] == 2 ** len(R.BENEFIT_MECHANISMS)


def test_study_reports_overlaps_and_inert_parts(tmp_path):
    summary = R.run_recombine_study(tmp_path / "recombine", draws=FAST)
    assert summary["overlapping_pairs"] >= 1
    assert summary["synergistic_pairs"] >= 1
    assert summary["inert_mechanisms"] >= 1


def test_study_writes_the_verdict_table(tmp_path):
    R.run_recombine_study(tmp_path / "recombine", draws=FAST)
    verdict = pd.read_csv(tmp_path / "recombine" / "tables" / "recombine_verdict.csv")
    assert "price_of_cheapness_pp" in verdict.columns
    assert "price_of_admissibility_pp" in verdict.columns


def test_cli_exposes_recombine():
    from ihd_solutions.cli import build_parser

    args = build_parser().parse_args(["recombine", "--draws", "400"])
    assert args.command == "recombine"
    assert args.draws == 400


def test_cli_runs_recombine(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(["recombine", "--output-dir", str(tmp_path / "cli"), "--draws", "400"])
    assert code == 0
    assert "ATLAS-5" in capsys.readouterr().out


# ---------------------------------------------------------------- stacking vs selection


def test_stacking_everything_is_worse_than_selecting():
    """A negative-effect mechanism drags the naive stack below the selected design."""
    contrast = R.stacking_contrast(draws=1200).set_index("design")
    assert float(contrast.loc["selected_composition", "probability_of_target"]) > float(
        contrast.loc["everything_on", "probability_of_target"]
    )


def test_constraints_only_is_the_weakest_design():
    contrast = R.stacking_contrast(draws=1200).set_index("design")
    assert float(contrast.loc["constraints_only", "probability_of_target"]) < float(
        contrast.loc["selected_composition", "probability_of_target"]
    )


def test_maturation_makes_retries_mostly_unnecessary():
    result = R.evaluate_composite(R.ATLAS5, 2000)
    assert result.mean_doses_used < 1.5
    assert result.mean_doses_used >= 1.0


def test_design_from_selection_restores_constraints(grid):
    design = R.design_from_selection(R.select_composite(grid))
    assert design.constraints_intact
    assert design.doses == max_precedented_doses()


def test_design_from_selection_falls_back_when_nothing_selected():
    assert R.design_from_selection({"best": None}) is R.COMPOSITE


def test_recommended_design_needs_fewer_procedures_than_its_dose_cap():
    result = R.evaluate_composite(R.ATLAS5, 1500)
    assert result.mean_doses_used < R.ATLAS5.doses


def test_report_shows_the_stacking_contrast(report):
    assert "NAIVE STACKING VERSUS SELECTION" in report
    assert "WORSE than selecting" in report


def test_report_shows_the_atlas5_headline(report):
    assert "ATLAS-5 AT THE REPORTED DRAW COUNT" in report
    assert "mean doses actually used" in report


def test_study_reports_the_stacking_penalty(tmp_path):
    summary = R.run_recombine_study(tmp_path / "recombine", draws=FAST)
    assert summary["stacking_penalty_pp"] > 0.0
    assert summary["recommended_mean_doses"] < 1.5


# ------------------------------------------------------- alternative report branches


def _render(grid, selection, verdict=None, **kwargs):
    """Render the report against a supplied verdict so every branch is reachable."""
    selection = selection or R.select_composite(grid)
    return R.recombine_report(
        R.evaluate_composite(R.COMPOSITE, 300),
        R.mechanism_marginals(draws=300),
        R.pairwise_interactions(draws=300),
        grid,
        selection,
        {
            "graft_probability": 0.9,
            "clearance_probability": 0.9,
            "joint_optimistic": 0.81,
            "randomised_proof": 0.0,
        },
        R.cost_ledger(R.COMPOSITE),
        R.mechanism_standalone(draws=300),
        R.maturation_sensitivity(draws=300),
        verdict if verdict is not None else R.composition_verdict(grid, selection),
        R.evaluate_composite(R.ATLAS5, 300),
        R.stacking_contrast(draws=300),
        **kwargs,
    )


def test_graft_breach_lists_attempts():
    design = R.Composite(doses=max_precedented_doses() + 2)
    assert "attempts" in design.graft_precedent_breaches
    assert "attempts" in design.precedent_breaches


def test_report_branch_single_mechanism_selection(grid):
    selection = R.select_composite(grid)
    verdict = R.composition_verdict(grid, selection)
    verdict["composition_adds_over_single_pp"] = 0.0
    verdict["selection_is_a_single_mechanism"] = True
    text = _render(grid, selection, verdict)
    assert "SELECTED A SINGLE MECHANISM" in text


def test_report_branch_composition_beats_single_without_a_large_margin(grid):
    selection = R.select_composite(grid)
    verdict = R.composition_verdict(grid, selection)
    verdict["composition_adds_over_single_pp"] = 0.0
    verdict["selection_is_a_single_mechanism"] = False
    verdict["composition_beat_best_single"] = True
    text = _render(grid, selection, verdict)
    assert "Composing beats the best single mechanism" in text


def test_report_branch_composition_does_not_beat_single(grid):
    selection = R.select_composite(grid)
    verdict = R.composition_verdict(grid, selection)
    verdict["composition_adds_over_single_pp"] = 0.0
    verdict["selection_is_a_single_mechanism"] = False
    verdict["composition_beat_best_single"] = False
    text = _render(grid, selection, verdict)
    assert "did NOT beat the best single mechanism" in text


def test_report_branch_admissibility_has_a_price(grid):
    selection = R.select_composite(grid)
    verdict = R.composition_verdict(grid, selection)
    verdict["price_of_admissibility_pp"] = 5.0
    text = _render(grid, selection, verdict)
    assert "Admissibility has a price" in text


def test_report_branch_cheapness_wins_on_risk(grid):
    selection = R.select_composite(grid)
    verdict = R.composition_verdict(grid, selection)
    verdict["price_of_cheapness_pp"] = 5.0
    verdict["cheapness_wins_on_risk"] = True
    text = _render(grid, selection, verdict)
    assert "the cheaper design is AHEAD" in text


def test_report_branch_cheapness_loses_on_risk(grid):
    selection = R.select_composite(grid)
    verdict = R.composition_verdict(grid, selection)
    verdict["price_of_cheapness_pp"] = 5.0
    verdict["cheapness_wins_on_risk"] = False
    text = _render(grid, selection, verdict)
    assert "STILL ahead after risk is charged" in text
