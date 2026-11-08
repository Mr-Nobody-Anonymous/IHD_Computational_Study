"""Tests for the gap-closure engine.

Key guards: the analytic ceiling must match the simulated one (that agreement is what
identified the structural cause), the source-reduction default must not perturb any
existing artefact, factors must saturate rather than multiply without bound, and no
verdict may be reported as evidence.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.gap_solver import (
    BASELINE_SCAR_FRACTION,
    DISCLAIMER,
    GAP_TAG,
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
    SOURCE_REDUCTION_PRECEDENT_MAX,
    Factor,
    GapVerdict,
    analytic_clearance_ceiling,
    assess_gaps,
    factor_frame,
    gap_report,
    hazard_gap_table,
    hazard_ratio_from_risks,
    minimum_remuscularisation_set,
    minimum_source_reduction,
    remuscularisation_scenario,
    remuscularisation_search,
    run_gap_study,
    solve_substrate_gap,
)
from ihd_solutions.hypothesis import (
    _ZERO_STRENGTHS,
    ClearanceParameters,
    simulate_clearance,
)

PERFECT = ClearanceParameters(gate_specificity=1.0, **{n: 1.0 for n in _ZERO_STRENGTHS})


@pytest.fixture(scope="module")
def matrix() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "population_id": "stable_ccs",
                "mode": "strict",
                "baseline_risk_pct": 13.9292,
                "composite_risk_pct": 4.1386,
            },
            {
                "population_id": "secondary_post_mi",
                "mode": "strict",
                "baseline_risk_pct": 20.148,
                "composite_risk_pct": 6.714,
            },
            {
                "population_id": "lifetime_exposure",
                "mode": "strict",
                "baseline_risk_pct": 1.98013,
                "composite_risk_pct": 1.98013,
            },
        ]
    )


# --------------------------------------------------------------------------
# Gap 3: hazard ratio, recomputed per population.
# --------------------------------------------------------------------------


def test_hazard_ratio_round_trips():
    baseline, hr = 0.20, 0.4
    treated = 1.0 - (1.0 - baseline) ** hr
    assert hazard_ratio_from_risks(baseline, treated) == pytest.approx(hr)


def test_hazard_ratio_rejects_a_zero_baseline():
    with pytest.raises(ValueError, match="positive"):
        hazard_ratio_from_risks(0.0, 0.0)


def test_required_hr_is_stricter_for_higher_baselines(matrix):
    table = hazard_gap_table(matrix).set_index("population_id")
    assert table.loc["secondary_post_mi", "required_hr"] < table.loc["stable_ccs", "required_hr"]


def test_the_composite_closes_the_gap_in_stable_ccs(matrix):
    row = hazard_gap_table(matrix).set_index("population_id").loc["stable_ccs"]
    assert row["closed"]
    assert row["achieved_hr"] <= 0.298
    assert row["extra_rrr_needed_pct"] == 0.0


def test_higher_risk_populations_stay_open(matrix):
    row = hazard_gap_table(matrix).set_index("population_id").loc["secondary_post_mi"]
    assert not row["closed"]
    assert row["extra_rrr_needed_pct"] > 0.0


def test_a_population_at_reference_is_not_a_gap(matrix):
    row = hazard_gap_table(matrix).set_index("population_id").loc["lifetime_exposure"]
    assert row["at_reference_risk"]
    assert row["closed"]


def test_extra_rrr_needed_is_the_reduction_on_current_risk(matrix):
    row = hazard_gap_table(matrix).set_index("population_id").loc["secondary_post_mi"]
    baseline = row["baseline_risk_pct"] / 100.0
    target = 1.0 - (1.0 - baseline) ** row["required_hr"]
    expected = 100.0 * (1.0 - target / (row["composite_risk_pct"] / 100.0))
    assert row["extra_rrr_needed_pct"] == pytest.approx(expected)


# --------------------------------------------------------------------------
# Gap 2: the plateau is structural.
# --------------------------------------------------------------------------


def test_source_reduction_defaults_to_zero_and_changes_nothing():
    """Guards every pre-existing hypothesis and discovery artefact."""
    implicit = simulate_clearance(PERFECT, "x", months=24.0, step_months=0.5)
    explicit = simulate_clearance(PERFECT, "x", months=24.0, step_months=0.5, source_reduction=0.0)
    np.testing.assert_allclose(implicit.substrate, explicit.substrate)


def test_analytic_ceiling_matches_the_simulated_plateau():
    """This agreement is what identified the constant source as the cause."""
    simulated = simulate_clearance(
        PERFECT, "x", months=600.0, step_months=0.5
    ).substrate_removal_fraction
    assert analytic_clearance_ceiling(0.0) == pytest.approx(simulated, abs=1e-3)


@pytest.mark.parametrize("strength", [0.0, 0.2, 0.35, 0.5])
def test_analytic_ceiling_matches_simulation_at_every_strength(strength):
    simulated = simulate_clearance(
        PERFECT, "x", months=600.0, step_months=0.5, source_reduction=strength
    ).substrate_removal_fraction
    assert analytic_clearance_ceiling(strength) == pytest.approx(simulated, abs=2e-3)


def test_source_reduction_monotonically_raises_clearance():
    values = [analytic_clearance_ceiling(s) for s in (0.0, 0.2, 0.35, 0.5, 0.8)]
    assert values == sorted(values)
    assert values[0] < values[-1]


def test_eradication_would_require_abolishing_the_apoptotic_source():
    """Full suppression eradicates in the model - and has no precedent whatsoever."""
    assert analytic_clearance_ceiling(SOURCE_REDUCTION_PRECEDENT_MAX) < 1.0
    assert analytic_clearance_ceiling(1.0) == pytest.approx(1.0)
    assert SOURCE_REDUCTION_PRECEDENT_MAX < 1.0


def test_source_reduction_rejects_out_of_range():
    with pytest.raises(ValueError):
        analytic_clearance_ceiling(1.5)
    with pytest.raises(ValueError):
        simulate_clearance(PERFECT, months=12.0, source_reduction=-0.1)


def test_the_precedented_strength_closes_the_gap_with_perfect_modules():
    minimum = minimum_source_reduction(horizon_months=24.0)
    assert minimum is not None
    assert minimum <= SOURCE_REDUCTION_PRECEDENT_MAX


def test_central_modules_cannot_close_the_gap_in_two_years():
    """The honest caveat: delivery and time are jointly binding."""
    minimum = minimum_source_reduction(horizon_months=24.0, perfect_modules=False)
    assert minimum is None or minimum > SOURCE_REDUCTION_PRECEDENT_MAX


def test_central_modules_can_close_it_given_a_long_horizon():
    minimum = minimum_source_reduction(horizon_months=240.0, perfect_modules=False)
    assert minimum is not None
    assert minimum <= SOURCE_REDUCTION_PRECEDENT_MAX


def test_clearance_table_flags_unprecedented_strengths():
    table = solve_substrate_gap()
    strong = table.loc[table["source_reduction"].gt(SOURCE_REDUCTION_PRECEDENT_MAX)]
    assert not strong.empty
    assert not strong["within_precedent"].any()


def test_impossible_clearance_target_returns_none():
    assert minimum_source_reduction(required=0.999999, horizon_months=6.0) is None


def test_clearance_table_is_tagged():
    assert solve_substrate_gap()["status"].eq(GAP_TAG).all()


# --------------------------------------------------------------------------
# Gap 1: remuscularisation, saturating.
# --------------------------------------------------------------------------


def test_baseline_reproduces_the_biovat_result():
    scenario = remuscularisation_scenario(())
    assert scenario.lvef_gain_pp == pytest.approx(3.9)
    assert scenario.scar_fraction == pytest.approx(BASELINE_SCAR_FRACTION)
    assert not scenario.meets_target


def test_baseline_scar_fraction_matches_primate_graft_measurements():
    """Independent check: primate grafts occupy 10-17% of scar."""
    assert 0.10 <= BASELINE_SCAR_FRACTION <= 0.17


def test_scar_fraction_saturates_at_full_replacement():
    scenario = remuscularisation_scenario(REMUSCULARISATION_FACTORS)
    assert scenario.saturated
    assert scenario.scar_fraction == 1.0
    assert scenario.lvef_gain_pp == pytest.approx(SCAR_RECOVERABLE_LVEF_PP)


def test_factors_cannot_multiply_past_the_scar():
    huge = Factor("absurd", 1000.0, "synthetic", "https://example.com", "rodent")
    scenario = remuscularisation_scenario((huge,))
    assert scenario.lvef_gain_pp <= SCAR_RECOVERABLE_LVEF_PP
    assert scenario.saturated


def test_a_minimum_two_factor_set_reaches_the_target():
    minimum = minimum_remuscularisation_set()
    assert len(minimum) == 2
    chosen = tuple(f for f in REMUSCULARISATION_FACTORS if f.name in minimum)
    scenario = remuscularisation_scenario(chosen)
    assert scenario.meets_target
    assert scenario.lvef_gain_pp >= CRITERION_THRESHOLD * SCAR_RECOVERABLE_LVEF_PP


def test_no_single_factor_reaches_the_target():
    for factor in REMUSCULARISATION_FACTORS:
        assert not remuscularisation_scenario((factor,)).meets_target


def test_no_factor_is_a_human_efficacy_result():
    scenario = remuscularisation_scenario(REMUSCULARISATION_FACTORS)
    assert not scenario.any_human_efficacy
    assert all(f.species != "human_outcomes" for f in REMUSCULARISATION_FACTORS)


def test_search_is_ordered_by_fewest_factors_first():
    search = remuscularisation_search()
    passing = search.loc[search["meets_target"]]
    assert passing.iloc[0]["n_factors"] == passing["n_factors"].min()
    assert len(search) == 2 ** len(REMUSCULARISATION_FACTORS)


def test_every_factor_carries_a_url_and_a_species():
    frame = factor_frame()
    assert frame["source_url"].str.startswith("https://").all()
    assert frame["species"].str.len().gt(3).all()
    assert frame["demonstrated_multiplier"].gt(1.0).all()


def test_factor_requires_an_https_source():
    with pytest.raises(ValueError, match="https"):
        Factor("x", 2.0, "p", "http://example.com", "rodent")


def test_scenario_rejects_an_impossible_baseline_fraction():
    with pytest.raises(ValueError):
        remuscularisation_scenario((), baseline_scar_fraction=1.5)


# --------------------------------------------------------------------------
# Verdicts, report, outputs.
# --------------------------------------------------------------------------


def test_all_three_gaps_get_a_verdict(matrix):
    results, detail = assess_gaps(matrix)
    assert len(results) == 3
    assert {r.verdict for r in results} <= set(GapVerdict)
    assert detail["hazard_populations_tested"] == 3


def test_verdicts_record_the_joint_delivery_time_constraint(matrix):
    _, detail = assess_gaps(matrix)
    assert detail["clearance_closes_in_24m_within_precedent"]
    assert detail["clearance_needs_near_perfect_delivery_for_24m"]


def test_report_states_it_is_not_evidence(matrix):
    text = gap_report(matrix)
    assert DISCLAIMER in text
    assert "None of this is evidence" in text
    assert "structural" in text
    assert "chronic-injury efficacy is" in text


def test_outputs_are_tagged_and_quarantined(tmp_path, matrix):
    target = tmp_path / "gaps"
    summary = run_gap_study(target, matrix=matrix)
    assert summary["status"] == "scenario_only"
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 9
    for path in (target / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].eq(GAP_TAG).all(), path.name


def test_missing_composite_matrix_is_a_clear_error(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError, match="composite"):
        run_gap_study(tmp_path / "gaps")


def test_cli_gaps_subcommand(tmp_path, capsys, matrix):
    from ihd_solutions.cli import main

    source = tmp_path / "matrix.csv"
    matrix.to_csv(source, index=False)
    target = tmp_path / "cli-gaps"
    assert main(["gaps", "--output-dir", str(target), "--composite-matrix", str(source)]) == 0
    out = capsys.readouterr().out
    assert DISCLAIMER in out
    assert "Wrote 9 artefacts" in out
