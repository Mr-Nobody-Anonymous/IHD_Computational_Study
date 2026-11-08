"""Tests for ATLAS-3, the graft arm given a variance architecture.

Three things are pinned here. First, the diagnosis: the published factors already
exceed the requirement, so the graft shortfall is variance, not capability. Second,
the falsifiers: retry redundancy must contribute exactly nothing when attempts
share a failure mode, and feedback exactly nothing when it cannot learn. Third,
the honesty of the bill: every extra attempt is charged a risk, three attempts is
beyond precedent, and nothing here is called a cure.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from ihd_solutions.atlas3 import (
    ATLAS3,
    ATLAS3_TAG,
    DISCLAIMER,
    FEEDBACK_MAX_UPLIFT,
    LEVER_NAMES,
    MAX_PRECEDENTED_ATTEMPTS,
    REQUIRED_MULTIPLIER,
    GraftUpgrades,
    attempt_sweep,
    current_state_voi,
    evaluate_graft,
    feedback_sweep,
    graft_draws,
    independence_sweep,
    joint_capability_probability,
    lever_ablation,
    maximum_tolerable_correlation,
    minimum_robust_attempts,
    precedented_option,
    run_atlas3_study,
    safety_sweep,
)
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solver import (
    BASELINE_SCAR_FRACTION,
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
)


@pytest.fixture(scope="module")
def single():
    return evaluate_graft(GraftUpgrades(attempts=1, translation_floor=1.0), draws=1500)


@pytest.fixture(scope="module")
def staged():
    return evaluate_graft(ATLAS3, draws=1500)


# ---------------------------------------------------------------- validation


def test_zero_attempts_rejected():
    with pytest.raises(ValueError, match="at least one attempt"):
        GraftUpgrades(attempts=0)


@pytest.mark.parametrize(
    "field",
    [
        "failure_correlation",
        "feedback_learning",
        "measurement_quality",
        "per_attempt_risk",
        "translation_floor",
    ],
)
def test_fraction_fields_reject_out_of_range(field):
    with pytest.raises(ValueError):
        GraftUpgrades(**{field: 1.5})
    with pytest.raises(ValueError):
        GraftUpgrades(**{field: -0.1})


def test_fraction_fields_accept_one():
    config = GraftUpgrades(failure_correlation=1.0, feedback_learning=1.0, measurement_quality=1.0)
    assert config.failure_correlation == 1.0


def test_unknown_lever_rejected():
    with pytest.raises(ValueError, match="unknown lever"):
        ATLAS3.without("potency")


def test_draws_must_be_positive():
    with pytest.raises(ValueError):
        evaluate_graft(ATLAS3, draws=0)


@pytest.mark.parametrize("lever", LEVER_NAMES)
def test_every_named_lever_is_ablatable(lever):
    assert ATLAS3.without(lever) != ATLAS3


# ------------------------------------------------------- the core diagnosis


def test_requirement_matches_the_criterion():
    assert pytest.approx(CRITERION_THRESHOLD / BASELINE_SCAR_FRACTION) == REQUIRED_MULTIPLIER


def test_published_capability_already_exceeds_the_requirement():
    product = float(np.prod([f.demonstrated_multiplier for f in REMUSCULARISATION_FACTORS]))
    assert product > REQUIRED_MULTIPLIER * 3.0


def test_capability_is_sufficient_but_one_graft_is_still_a_coin_flip(single):
    """The nominal-product trap, pinned: capability exceeds the bar, one shot fails."""
    assert 0.40 < single.probability_of_target < 0.62
    assert not single.robust


def test_staged_retry_beats_a_single_graft(single, staged):
    assert staged.probability_of_target > single.probability_of_target + 0.30


def test_staged_graft_is_robust_with_margin(staged):
    assert staged.robust
    assert staged.robust_with_margin


def test_no_configuration_reaches_certainty():
    """Any P of 1.0 in this model is an artefact, not a result."""
    assert evaluate_graft(ATLAS3, draws=1200).probability_of_target < 1.0
    ceiling = evaluate_graft(
        GraftUpgrades(attempts=6, feedback_learning=1.0, measurement_quality=1.0), draws=1200
    )
    assert ceiling.probability_of_target < 1.0


def test_feedback_uplift_is_bounded():
    assert 0.0 < FEEDBACK_MAX_UPLIFT < 1.0


# ------------------------------------------------------------- falsifier one


def test_correlated_attempts_gain_exactly_nothing():
    frame = independence_sweep(draws=900)
    worst = frame.loc[frame["failure_correlation"].eq(1.0)].iloc[0]
    assert worst["redundancy_gain_pp"] == pytest.approx(0.0, abs=1e-12)


def test_correlated_attempts_equal_the_first_attempt():
    best, _, first = graft_draws(GraftUpgrades(attempts=4, failure_correlation=1.0), 600)
    assert np.array_equal(best, first)


def test_independence_gain_decreases_as_failure_becomes_shared():
    frame = independence_sweep(draws=1200).set_index("failure_correlation")
    gains = frame["redundancy_gain_pp"]
    assert gains.loc[0.0] > gains.loc[0.5] > gains.loc[1.0]


def test_correlation_sweep_is_not_robust_when_failure_is_shared():
    frame = independence_sweep(draws=1200).set_index("failure_correlation")
    assert not bool(frame.loc[1.0, "robust"])
    assert not bool(frame.loc[0.75, "robust"])


# ------------------------------------------------------------- falsifier two


def test_zero_learning_gains_exactly_nothing():
    frame = feedback_sweep(draws=900)
    zero = frame.loc[frame["feedback_learning"].eq(0.0)].iloc[0]
    assert zero["feedback_gain_pp"] == pytest.approx(0.0, abs=1e-12)


def test_feedback_gain_increases_with_learning():
    frame = feedback_sweep(draws=1200).set_index("feedback_learning")
    gains = frame["feedback_gain_pp"]
    assert gains.loc[0.0] < gains.loc[0.5] < gains.loc[1.0]


def test_unmeasurable_engraftment_disables_feedback():
    """A clinical signal gates learning: quality 0 must equal learning 0."""
    blind = evaluate_graft(ATLAS3.without("measurement"), 1200)
    deaf = evaluate_graft(ATLAS3.without("feedback"), 1200)
    assert blind.probability_of_target == pytest.approx(deaf.probability_of_target, abs=1e-12)


# --------------------------------------------------------------- ablation


def test_every_lever_is_load_bearing():
    frame = lever_ablation(ATLAS3, draws=1500)
    assert set(frame["removed_lever"]) == set(LEVER_NAMES)
    assert frame["load_bearing"].all()
    assert (frame["probability_loss_pp"] > 0.0).all()


def test_retry_redundancy_is_the_dominant_lever():
    frame = lever_ablation(ATLAS3, draws=1500).set_index("removed_lever")
    assert frame["probability_loss_pp"].idxmax() == "retry_redundancy"


def test_removing_redundancy_returns_to_the_single_graft(single):
    frame = lever_ablation(ATLAS3, draws=1500).set_index("removed_lever")
    assert frame.loc["retry_redundancy", "probability_of_target"] == pytest.approx(
        single.probability_of_target, abs=0.02
    )


# ------------------------------------------------------------- the safety bill


def test_three_attempts_is_within_precedent(staged):
    """CORRECTED: the limit is 4 (CATO, HOPE-2), so three was never beyond precedent."""
    assert MAX_PRECEDENTED_ATTEMPTS == 4
    assert ATLAS3.attempts <= MAX_PRECEDENTED_ATTEMPTS
    assert not staged.beyond_precedent


def test_the_limit_comes_from_the_registry():
    from ihd_solutions.precedent import max_precedented_doses

    assert max_precedented_doses() == MAX_PRECEDENTED_ATTEMPTS


def test_six_attempts_is_still_beyond_precedent():
    assert GraftUpgrades(attempts=6).beyond_precedent


def test_two_attempts_stays_within_precedent():
    assert not GraftUpgrades(attempts=2).beyond_precedent


def test_every_extra_attempt_is_charged_a_risk():
    frame = attempt_sweep(draws=1200).set_index("attempts")
    risks = frame["cumulative_adverse_risk"]
    assert (np.diff(risks.values) > 0.0).all()


def test_risk_adjusted_success_never_exceeds_raw_success(staged):
    assert staged.risk_adjusted_probability < staged.probability_of_target
    assert staged.cumulative_adverse_risk > 0.0


def test_responders_are_not_exposed_to_further_attempts(staged):
    """Stopping on success keeps the median number of procedures at one."""
    assert 1.0 <= staged.median_attempts_used <= ATLAS3.attempts
    assert staged.median_attempts_used == pytest.approx(1.0)


def test_success_is_monotone_in_attempts():
    frame = attempt_sweep(draws=1500).sort_values("attempts")
    probabilities = frame["probability_of_target"].to_numpy()
    assert (np.diff(probabilities) >= -1e-12).all()


def test_more_procedural_risk_lowers_risk_adjusted_success():
    frame = safety_sweep(per_attempt_risks=(0.0, 0.05, 0.20), draws=1200)
    three = frame.loc[frame["attempts"].eq(3)].set_index("per_attempt_risk")
    adjusted = three["risk_adjusted_probability"]
    assert adjusted.loc[0.0] > adjusted.loc[0.05] > adjusted.loc[0.20]


def test_zero_risk_leaves_success_undiscounted():
    free = evaluate_graft(GraftUpgrades(attempts=3, per_attempt_risk=0.0), 900)
    assert free.cumulative_adverse_risk == pytest.approx(0.0)
    assert free.risk_adjusted_probability == pytest.approx(free.probability_of_target)


def test_translation_penalty_is_swept_not_asserted():
    perfect = evaluate_graft(GraftUpgrades(attempts=3, translation_floor=1.0), 1200)
    penalised = evaluate_graft(GraftUpgrades(attempts=3, translation_floor=0.4), 1200)
    assert penalised.probability_of_target < perfect.probability_of_target


def test_gain_is_capped_by_the_recoverable_ceiling():
    gains, _, _ = graft_draws(GraftUpgrades(attempts=4), 900)
    assert gains.max() <= SCAR_RECOVERABLE_LVEF_PP + 1e-9
    assert gains.min() >= 0.0


# ------------------------------------------------------- the specification


def test_two_attempts_already_clear_the_bar():
    assert minimum_robust_attempts(ATLAS3, draws=1500) == 2


def test_feasibility_does_not_require_breaking_precedent():
    """The precedent breach buys headroom; it is not what makes the design work."""
    minimum = minimum_robust_attempts(ATLAS3, draws=1500)
    assert minimum is not None and minimum <= MAX_PRECEDENTED_ATTEMPTS


def test_the_specification_does_not_move_with_the_sampling_budget():
    """At small N the 2-SE margin wrongly rejected two attempts; the floor fixes it."""
    assert minimum_robust_attempts(ATLAS3, draws=200) == minimum_robust_attempts(ATLAS3, draws=4000)


def test_precedented_option_is_robust_and_cheaper(staged):
    option = precedented_option(draws=4000)
    # The cheapest ROBUST design, not the attempt limit itself.
    assert option.attempts < ATLAS3.attempts
    assert option.attempts <= MAX_PRECEDENTED_ATTEMPTS
    assert not option.beyond_precedent
    assert option.robust_with_margin
    assert option.cumulative_adverse_risk < staged.cumulative_adverse_risk
    assert option.probability_of_target < staged.probability_of_target


def test_the_third_attempt_is_headroom_not_feasibility(staged):
    option = precedented_option(draws=4000)
    assert staged.probability_of_target - option.probability_of_target > 0.05


def test_no_robust_design_when_attempts_are_fully_correlated():
    assert minimum_robust_attempts(replace(ATLAS3, failure_correlation=1.0), draws=900) is None


def test_tolerable_correlation_is_interior():
    rho = maximum_tolerable_correlation(ATLAS3, draws=900)
    assert rho is not None and 0.0 < rho < 1.0


def test_tolerable_correlation_brackets_the_bar():
    rho = maximum_tolerable_correlation(ATLAS3, draws=2000)
    assert rho is not None
    below = evaluate_graft(replace(ATLAS3, failure_correlation=rho - 0.10), 2000)
    above = evaluate_graft(replace(ATLAS3, failure_correlation=min(rho + 0.10, 1.0)), 2000)
    assert below.robust
    assert not above.robust


# -------------------------------------------------------- the recomputed VOI


def test_joint_is_the_product_and_proof_is_absent():
    joint = joint_capability_probability(draws=900)
    assert joint["joint_probability_optimistic"] == pytest.approx(
        joint["clearance_probability"] * joint["graft_probability"], abs=1e-12
    )
    assert joint["randomised_proof"] == 0.0


def test_voi_is_ranked_and_arithmetically_consistent():
    frame = current_state_voi(draws=900)
    assert len(frame) == 2
    assert frame["voi_gain"].is_monotonic_decreasing
    for row in frame.itertuples():
        assert row.voi_gain == pytest.approx(row.joint_if_resolved - row.joint_now, abs=1e-12)
        assert 0.0 < row.joint_if_resolved <= 1.0


def test_clearance_remains_the_higher_value_uncertainty():
    frame = current_state_voi(draws=1200)
    assert frame.iloc[0]["capability"] == "substrate_removal"


def test_voi_resolving_both_is_not_claimed():
    frame = current_state_voi(draws=900)
    assert (frame["joint_if_resolved"] < 1.0).all()


# --------------------------------------------------------------- artefacts


@pytest.fixture(scope="module")
def study(tmp_path_factory):
    target = tmp_path_factory.mktemp("atlas3") / "atlas3"
    return target, run_atlas3_study(target, draws=400)


def test_outputs_are_counted(study):
    target, summary = study
    written = [p for p in target.rglob("*") if p.is_file()]
    assert len(written) == summary["outputs"] == 10


def test_disclaimer_is_written_and_unambiguous(study):
    target, summary = study
    text = (target / "DISCLAIMER.txt").read_text()
    assert "INVENTED" in text
    assert text.strip() == DISCLAIMER.strip()
    assert summary["disclaimer"] == DISCLAIMER
    assert summary["cures_declared"] == 0


def test_every_table_is_tagged(study):
    target, _ = study
    tables = list((target / "tables").glob("*.csv"))
    assert len(tables) == 8
    for path in tables:
        frame = pd.read_csv(path)
        assert "status" in frame.columns, path.name
        assert (frame["status"] == ATLAS3_TAG).all(), path.name


def test_report_states_the_limits(study):
    target, _ = study
    text = (target / "reports" / "atlas3_report.txt").read_text()
    assert "not a cure" in text
    assert "MAXIMUM TOLERABLE CORRELATION" in text
    assert "BEYOND PRECEDENT" in text
    assert "THE CHEAPEST PRECEDENTED OPTION" in text
    assert "CORRECTED from 2" in text
    assert "Randomised proof remains zero." in text


def test_summary_exposes_the_price_and_the_spec(study):
    _, summary = study
    assert summary["status"] == ATLAS3_TAG
    assert summary["minimum_robust_attempts"] == 2
    assert summary["precedented_option_robust"] is True
    assert 0.8 < summary["precedented_option_probability"] < 0.9
    assert summary["cumulative_adverse_risk"] > 0.0
    assert 0.0 < summary["maximum_tolerable_correlation"] < 1.0
    assert 0.0 < summary["joint_probability_optimistic"] < 1.0


def test_cli_atlas3_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(["atlas3", "--output-dir", str(tmp_path / "cli-atlas3"), "--draws", "200"])
    assert code == 0
    out = capsys.readouterr().out
    assert "ATLAS-3" in out
    assert "Wrote 10 artefacts" in out
    assert (tmp_path / "cli-atlas3" / "DISCLAIMER.txt").exists()
