"""Tests for ATLAS-2, the upgraded architecture.

Two things are pinned here. First, the solve: the clearance arm must be robust
across the whole population at the solved operating point, with every active
upgrade load-bearing. Second, the honesty of the bill: the precedent gap, the
redundancy limit, and the fact that a conditioning gate cannot create success.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ihd_solutions.atlas import ATLAS_POTENCY
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.upgrade import (
    ATLAS2,
    DISCLAIMER,
    MAX_PRECEDENTED_PATHS,
    SOLVED_MONTHS,
    UPGRADE_NAMES,
    Upgrades,
    achieved_clearance_draws,
    coupled_graft_with_upgrades,
    evaluate_upgrades,
    horizon_sweep,
    minimum_robust_cap,
    potency_ceiling_sweep,
    redundancy_sweep,
    run_upgrade_study,
    source_gate_tradeoff,
    source_required_for_target,
    upgrade_ablation,
    upgrade_report,
)


@pytest.fixture(scope="module")
def solved():
    return evaluate_upgrades(ATLAS2, months=SOLVED_MONTHS, draws=600)


@pytest.fixture(scope="module")
def baseline():
    return evaluate_upgrades(Upgrades(), months=24.0, draws=600)


@pytest.fixture(scope="module")
def ablation() -> pd.DataFrame:
    return upgrade_ablation(ATLAS2, draws=500, months=SOLVED_MONTHS)


# --------------------------------------------------------------------------
# The solve.
# --------------------------------------------------------------------------


def test_atlas2_is_robust_across_the_whole_population(solved):
    assert solved.robust
    assert solved.probability_of_target >= 0.8


def test_the_solve_has_margin_not_a_boundary_estimate(solved):
    """A point estimate sitting on the bar would be the borderline case again."""
    assert solved.robust_with_margin
    assert solved.probability_of_target - 2 * solved.monte_carlo_standard_error >= 0.8
    assert solved.eligible_fraction == pytest.approx(1.0)
    assert solved.unconditional_probability == pytest.approx(solved.probability_of_target)


def test_the_unupgraded_baseline_is_not_robust(baseline):
    assert not baseline.robust
    assert baseline.probability_of_target < 0.05


def test_the_upgrade_is_a_large_real_improvement(solved, baseline):
    assert solved.probability_of_target - baseline.probability_of_target > 0.7


def test_median_clearance_clears_the_target(solved):
    assert solved.median_clearance >= CRITERION_THRESHOLD


def test_no_patient_is_excluded_to_buy_the_result(solved):
    """The earlier gate-based version excluded 44%. This one must not."""
    assert not ATLAS2.verified_source
    assert solved.eligible_fraction == 1.0


# --------------------------------------------------------------------------
# Every active upgrade is load-bearing; the unused one is reported as unused.
# --------------------------------------------------------------------------


def test_all_active_upgrades_are_load_bearing(ablation):
    active = ablation.loc[ablation["active_in_construct"]]
    assert len(active) == 3
    assert active["load_bearing"].all()


def test_titration_is_the_largest_contributor(ablation):
    ranked = ablation.loc[ablation["active_in_construct"]].sort_values(
        "probability_loss_pp", ascending=False
    )
    assert ranked.iloc[0]["removed_upgrade"] == "titration"
    assert ranked.iloc[0]["probability_loss_pp"] > 50.0


def test_verified_source_is_reported_as_unused_not_as_inert(ablation):
    row = ablation.set_index("removed_upgrade").loc["verified_source"]
    assert not row["active_in_construct"]
    assert not row["load_bearing"]


def test_ablation_covers_every_upgrade(ablation):
    assert set(ablation["removed_upgrade"]) == set(UPGRADE_NAMES)


def test_unknown_upgrade_is_rejected():
    with pytest.raises(ValueError, match="unknown upgrade"):
        ATLAS2.without("not_an_upgrade")


# --------------------------------------------------------------------------
# The horizon correction to the earlier audit.
# --------------------------------------------------------------------------


def test_source_is_binding_at_24_months_but_not_beyond_36():
    """The 55.4% ceiling was horizon-specific, not a permanent bound."""
    assert source_required_for_target(months=24.0) == pytest.approx(0.2230, abs=1e-3)
    assert source_required_for_target(months=36.0) == pytest.approx(0.0, abs=1e-6)
    assert source_required_for_target(months=60.0) == pytest.approx(0.0, abs=1e-6)


def test_precedented_potency_alone_does_not_fix_it():
    config = Upgrades(
        gate_paths=ATLAS2.gate_paths,
        node_paths=ATLAS2.node_paths,
        titrate=True,
        max_potency=ATLAS_POTENCY,
    )
    result = evaluate_upgrades(config, months=SOLVED_MONTHS, draws=500)
    assert not result.robust
    assert result.probability_of_target < 0.25


def test_time_and_potency_are_substitutes():
    grid = horizon_sweep(months_grid=(24.0, 120.0), caps=(2.0,), draws=350)
    short = grid.loc[grid["horizon_months"].eq(24.0)].iloc[0]
    long = grid.loc[grid["horizon_months"].eq(120.0)].iloc[0]
    assert long["probability_of_target"] > short["probability_of_target"]


# --------------------------------------------------------------------------
# The bill: precedent gap and redundancy limit.
# --------------------------------------------------------------------------


def test_the_solved_cap_is_beyond_precedent_and_reported_as_so():
    assert ATLAS2.max_potency > ATLAS_POTENCY
    gap = ATLAS2.max_potency / ATLAS_POTENCY
    assert 1.5 < gap < 3.0


def test_redundancy_stays_within_the_precedented_count():
    assert ATLAS2.gate_paths <= MAX_PRECEDENTED_PATHS
    assert not ATLAS2.redundancy_beyond_precedent


def test_redundancy_beyond_precedent_is_flagged():
    assert Upgrades(gate_paths=MAX_PRECEDENTED_PATHS + 1).redundancy_beyond_precedent


def test_redundancy_helps_and_then_saturates():
    sweep = redundancy_sweep(paths=(1, 2, 3), draws=350, months=SOLVED_MONTHS)
    values = list(sweep["probability_of_target"])
    assert values[1] > values[0]
    assert values[2] >= values[1]


def test_potency_ceiling_sweep_marks_the_robust_boundary():
    sweep = potency_ceiling_sweep(ceilings=(1.284, 3.0), draws=350, months=SOLVED_MONTHS)
    low = sweep.loc[sweep["max_potency"].eq(1.284)].iloc[0]
    high = sweep.loc[sweep["max_potency"].eq(3.0)].iloc[0]
    assert not low["robust"]
    assert high["robust"]


def test_minimum_robust_cap_is_found_and_is_minimal():
    cap = minimum_robust_cap(months=SOLVED_MONTHS, draws=300)
    assert cap is not None
    fields = {f: getattr(ATLAS2, f) for f in ATLAS2.__dataclass_fields__}
    at_cap = evaluate_upgrades(
        Upgrades(**{**fields, "max_potency": cap}), months=SOLVED_MONTHS, draws=400
    )
    below = evaluate_upgrades(
        Upgrades(**{**fields, "max_potency": cap * 0.8}), months=SOLVED_MONTHS, draws=400
    )
    assert at_cap.probability_of_target >= 0.75
    assert below.probability_of_target < at_cap.probability_of_target


def test_minimum_robust_cap_returns_none_when_unreachable():
    config = Upgrades(gate_paths=1, node_paths=1, titrate=False)
    assert minimum_robust_cap(months=24.0, upgrades=config, draws=200, upper=2.0) is None


# --------------------------------------------------------------------------
# A conditioning gate cannot create success.
# --------------------------------------------------------------------------


def test_conditioning_does_not_change_the_unconditional_probability():
    """A run-in selects who is reported on; it does not improve outcomes."""
    sweep = source_gate_tradeoff(thresholds=(0.0, 0.10), draws=500)
    values = list(sweep["unconditional_probability"])
    assert values[0] == pytest.approx(values[1], abs=0.02)


def test_a_tighter_gate_excludes_more_and_raises_only_the_conditional_number():
    sweep = source_gate_tradeoff(thresholds=(0.0, 0.35), draws=500).set_index("source_required")
    loose, tight = sweep.loc[0.0], sweep.loc[0.35]
    assert tight["eligible_fraction"] < loose["eligible_fraction"]
    assert tight["conditional_probability"] >= loose["conditional_probability"]


def test_unconditional_is_never_above_conditional():
    sweep = source_gate_tradeoff(thresholds=(0.10, 0.35), draws=400)
    assert (sweep["unconditional_probability"] <= sweep["conditional_probability"] + 1e-9).all()


def test_an_impossible_gate_excludes_everyone():
    result = evaluate_upgrades(
        Upgrades(titrate=True, verified_source=True, source_required=0.99),
        months=SOLVED_MONTHS,
        draws=300,
    )
    assert result.eligible_fraction == 0.0
    assert result.probability_of_target == 0.0


# --------------------------------------------------------------------------
# The sequencing idea is recovered, not invented anew.
# --------------------------------------------------------------------------


def test_zero_coupling_still_contributes_exactly_nothing():
    row = coupled_graft_with_upgrades(
        ATLAS2, months=SOLVED_MONTHS, bed_coupling=0.0, draws=500
    ).iloc[0]
    assert row["sequencing_gain_pp"] == pytest.approx(0.0, abs=1e-9)


def test_the_bed_is_now_actually_prepared():
    row = coupled_graft_with_upgrades(
        ATLAS2, months=SOLVED_MONTHS, bed_coupling=1.0, draws=500
    ).iloc[0]
    assert row["median_bed_quality"] > 0.9


def test_sequencing_delivers_its_gain_again():
    row = coupled_graft_with_upgrades(
        ATLAS2, months=SOLVED_MONTHS, bed_coupling=1.0, draws=500
    ).iloc[0]
    assert row["sequencing_gain_pp"] > 15.0
    assert row["graft_probability"] > row["graft_probability_unsequenced"]


def test_the_graft_is_still_not_robust():
    row = coupled_graft_with_upgrades(
        ATLAS2, months=SOLVED_MONTHS, bed_coupling=1.0, draws=500
    ).iloc[0]
    assert not row["graft_robust"]
    assert row["graft_probability"] < 0.8


def test_coupled_graft_rejects_bad_coupling():
    with pytest.raises(ValueError):
        coupled_graft_with_upgrades(ATLAS2, bed_coupling=1.5, draws=200)


# --------------------------------------------------------------------------
# Engine, validation, report, quarantine.
# --------------------------------------------------------------------------


def test_per_draw_engine_shapes_and_eligibility():
    achieved, potency, eligible = achieved_clearance_draws(ATLAS2, months=SOLVED_MONTHS, draws=200)
    assert achieved.shape == potency.shape == (200,)
    assert eligible == 1.0
    assert (potency >= 1.0).all()
    assert (potency <= ATLAS2.max_potency + 1e-9).all()


def test_titration_never_exceeds_its_cap():
    achieved, potency, _ = achieved_clearance_draws(ATLAS2, months=SOLVED_MONTHS, draws=150)
    assert (potency <= ATLAS2.max_potency + 1e-9).all()


def test_results_are_reproducible():
    first = evaluate_upgrades(ATLAS2, months=SOLVED_MONTHS, draws=200, seed=3)
    second = evaluate_upgrades(ATLAS2, months=SOLVED_MONTHS, draws=200, seed=3)
    assert first.probability_of_target == second.probability_of_target


def test_invalid_upgrade_configs_are_rejected():
    with pytest.raises(ValueError):
        Upgrades(gate_paths=0)
    with pytest.raises(ValueError):
        Upgrades(max_potency=0.5)
    with pytest.raises(ValueError):
        Upgrades(titration_steps=1)
    with pytest.raises(ValueError):
        Upgrades(source_required=1.5)


def test_report_states_the_solve_the_price_and_what_remains(solved, baseline, ablation):
    text = upgrade_report(
        baseline,
        solved,
        ablation,
        redundancy_sweep(paths=(1, 3), draws=250, months=SOLVED_MONTHS),
        potency_ceiling_sweep(ceilings=(1.284, 3.0), draws=250, months=SOLVED_MONTHS),
        horizon_sweep(months_grid=(60.0,), caps=(2.59,), draws=250),
        coupled_graft_with_upgrades(ATLAS2, months=SOLVED_MONTHS, draws=250),
        ATLAS2.max_potency,
    )
    assert DISCLAIMER in text
    assert "precedent gap" in text
    assert "A CORRECTION TO THE EARLIER AUDIT" in text
    assert "not used at this horizon" in text
    assert "is not a cure." in text


def test_outputs_are_counted_and_quarantined(tmp_path):
    target = tmp_path / "upgrade"
    summary = run_upgrade_study(target, draws=250)
    assert summary["status"] == "invented_architecture_unvalidated"
    assert summary["cures_declared"] == 0
    assert summary["robust"] is True
    assert summary["eligible_fraction"] == 1.0
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 10
    assert not any(name.startswith("solutions") for name in written)


def test_every_table_is_tagged(tmp_path):
    target = tmp_path / "upgrade"
    run_upgrade_study(target, draws=250)
    for path in (target / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].eq("atlas2_invented_architecture_unvalidated").all(), path.name


def test_cli_upgrade_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(["upgrade", "--output-dir", str(tmp_path / "cli-upgrade"), "--draws", "200"])
    assert code == 0
    out = capsys.readouterr().out
    assert DISCLAIMER in out
    assert "Wrote 10 artefacts" in out
