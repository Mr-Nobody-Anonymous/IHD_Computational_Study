"""Tests for the quarantined REVERT-1 hypothesis model.

Two jobs here. First, prove the invention can never contaminate the
evidence-graded outputs. Second, prove the coupling claim is actually being
tested - including the possibility that it fails.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest

from ihd_solutions.catalog import load_catalog
from ihd_solutions.cure import assess_catalog, search_cure_combination
from ihd_solutions.evidence import load_anchors
from ihd_solutions.hypothesis import (
    ARMS,
    DISCLAIMER,
    HYPOTHESIS_TAG,
    KILL_CRITERIA,
    REVERT1_MODULES,
    ClearanceParameters,
    LesionState,
    ModuleSpec,
    PrecedentLevel,
    ablation_study,
    architecture_ceiling,
    arm_parameters,
    best_case_horizon_sweep,
    design_report,
    kill_criteria_frame,
    module_frame,
    monte_carlo_hypothesis,
    required_hazard_ratio,
    run_hypothesis_study,
    sensitivity_summary,
    simulate_clearance,
    target_product_profile,
)


@pytest.fixture(scope="module")
def ablation() -> pd.DataFrame:
    return ablation_study()


@pytest.fixture(scope="module")
def draws() -> pd.DataFrame:
    return monte_carlo_hypothesis(draws=60, seed=7)


def arm_row(ablation: pd.DataFrame, arm: str):
    return ablation.loc[ablation["arm"].eq(arm)].iloc[0]


# --------------------------------------------------------------------------
# Quarantine: the invention must never reach the evidence tables.
# --------------------------------------------------------------------------


def test_no_module_id_appears_in_the_therapy_catalog():
    catalog_ids = set(load_catalog())
    for spec in REVERT1_MODULES:
        assert spec.module_id not in catalog_ids
    assert not any("revert" in therapy_id.lower() for therapy_id in catalog_ids)


def test_no_hypothesis_id_appears_in_the_anchor_registry():
    therapy_ids = {anchor.therapy_id for anchor in load_anchors().values()}
    for spec in REVERT1_MODULES:
        assert spec.module_id not in therapy_ids
    assert not any("revert" in therapy_id.lower() for therapy_id in therapy_ids)


def test_cure_search_never_selects_the_invention(catalog, anchors, population, days):
    assessments = assess_catalog(catalog, anchors, population, days)
    cover = search_cure_combination(assessments)
    payload = json.dumps(
        {"assessments": assessments.to_dict("records"), "selected": list(cover.selected)},
        default=str,
    ).lower()
    assert "revert" not in payload
    for spec in REVERT1_MODULES:
        assert spec.module_id.lower() not in payload


def test_hypothesis_outputs_are_tagged_on_every_row(tmp_path):
    run_hypothesis_study(tmp_path / "hyp", draws=20, seed=3)
    tables = sorted((tmp_path / "hyp" / "tables").glob("*.csv"))
    assert len(tables) == 8
    for path in tables:
        frame = pd.read_csv(path)
        assert "status" in frame.columns, path.name
        assert frame["status"].eq(HYPOTHESIS_TAG).all(), path.name


def test_hypothesis_writes_only_inside_its_own_directory(tmp_path):
    target = tmp_path / "hyp"
    summary = run_hypothesis_study(target, draws=20, seed=3)
    assert summary["status"] == "hypothesis_only"
    assert summary["disclaimer"] == DISCLAIMER
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 10
    assert not any("solutions" in name for name in written)


def test_report_leads_with_the_disclaimer(ablation, draws):
    text = design_report(
        ablation,
        target_product_profile(),
        draws,
        sensitivity_summary(draws),
        architecture_ceiling(),
    )
    assert DISCLAIMER in text.split("\n")[2]
    assert "not a therapy" in text.lower()
    assert "does not regenerate myocardium" in text


# --------------------------------------------------------------------------
# The coupling claim, tested so that it could fail.
# --------------------------------------------------------------------------


def test_efferocytosis_alone_is_worse_than_untreated(ablation):
    """The core design claim: uptake without export grows the lesion."""
    alone = arm_row(ablation, "efferocytosis_only")
    untreated = arm_row(ablation, "untreated")
    assert alone.substrate_removal_pct < untreated.substrate_removal_pct
    assert alone.phagocyte_depletion
    assert not untreated.phagocyte_depletion


def test_efferocytosis_alone_looks_good_at_twelve_months_then_fails(ablation):
    """A short study would be misled - this is why the horizon matters."""
    alone = arm_row(ablation, "efferocytosis_only")
    assert alone.necrotic_core_change_12m_pct < 0.0
    assert alone.necrotic_core_change_pct > 0.0


def test_full_architecture_beats_every_partial_arm(ablation):
    full = arm_row(ablation, "revert1_full")
    others = ablation.loc[~ablation["arm"].eq("revert1_full")]
    assert (full.substrate_removal_pct > others.substrate_removal_pct).all()
    assert full.necrotic_core_change_pct < 0.0
    assert not full.phagocyte_depletion


def test_co_delivery_beats_either_module_alone(ablation):
    combined = arm_row(ablation, "efflux_plus_efferocytosis")
    assert combined.substrate_removal_pct > arm_row(ablation, "efflux_only").substrate_removal_pct
    assert (
        combined.substrate_removal_pct
        > arm_row(ablation, "efferocytosis_only").substrate_removal_pct
    )


def test_the_gate_is_what_makes_the_dose_reachable(ablation):
    gated = arm_row(ablation, "revert1_full")
    ungated = arm_row(ablation, "revert1_ungated")
    assert gated.substrate_removal_pct > ungated.substrate_removal_pct


def test_driver_lowering_alone_does_not_clear_the_existing_core(ablation):
    """Mirrors the real clinical picture: lipid lowering stabilises, not erases."""
    driver = arm_row(ablation, "driver_only")
    assert driver.necrotic_core_change_pct > 0.0
    assert driver.substrate_removal_pct < arm_row(ablation, "revert1_full").substrate_removal_pct


def test_untreated_lesion_progresses(ablation):
    untreated = arm_row(ablation, "untreated")
    assert untreated.substrate_removal_pct < 0.0
    assert untreated.necrotic_core_change_pct > 0.0


@pytest.mark.parametrize("module", ["gate_specificity", "efflux_gain", "plasma_driver_reduction"])
def test_stronger_modules_never_reduce_clearance(module):
    weak = ClearanceParameters(**{module: 0.2})
    strong = ClearanceParameters(**{module: 0.9})
    low = simulate_clearance(weak, months=24.0, step_months=0.5).substrate_removal_fraction
    high = simulate_clearance(strong, months=24.0, step_months=0.5).substrate_removal_fraction
    assert high >= low


def test_zero_strength_modules_reproduce_the_untreated_arm():
    zero = ClearanceParameters(
        gate_specificity=1.0,
        crystal_solubilisation=0.0,
        efflux_gain=0.0,
        efferocytosis_gain=0.0,
        influx_block=0.0,
        resolution=0.0,
        nucleation_block=0.0,
        plasma_driver_reduction=0.0,
    )
    direct = simulate_clearance(zero, months=12.0, step_months=0.5)
    untreated = simulate_clearance(arm_parameters("untreated"), months=12.0, step_months=0.5)
    np.testing.assert_allclose(direct.substrate, untreated.substrate)


def test_substrate_excludes_stabilising_macrocalcification():
    result = simulate_clearance(arm_parameters("revert1_full"), months=12.0, step_months=0.5)
    expected = result.free_cholesterol + result.necrotic_debris + result.microcalcification
    np.testing.assert_allclose(result.substrate, expected)


def test_trajectory_is_nonnegative_and_viability_bounded():
    result = simulate_clearance(arm_parameters("efferocytosis_only"), months=60.0, step_months=0.5)
    for series in (
        result.free_cholesterol,
        result.necrotic_debris,
        result.microcalcification,
        result.macrophage_viability,
    ):
        assert np.all(series >= 0.0)
        assert np.all(np.isfinite(series))
    assert np.all(result.macrophage_viability <= 1.0)


def test_integrator_is_step_size_stable():
    fine = simulate_clearance(arm_parameters("revert1_full"), months=12.0, step_months=0.125)
    coarse = simulate_clearance(arm_parameters("revert1_full"), months=12.0, step_months=0.5)
    assert abs(fine.substrate_removal_fraction - coarse.substrate_removal_fraction) < 0.05


# --------------------------------------------------------------------------
# Back-calculated specification.
# --------------------------------------------------------------------------


def test_required_hazard_ratio_lands_between_reference_and_baseline():
    baseline, reference = 0.20, 0.02
    hr = required_hazard_ratio(baseline, reference, allowed_excess_remaining=0.2)
    assert 0.0 < hr < 1.0
    implied = 1.0 - (1.0 - baseline) ** hr
    assert implied == pytest.approx(reference + 0.2 * (baseline - reference))


def test_zero_allowed_excess_demands_full_normalisation():
    hr = required_hazard_ratio(0.20, 0.02, allowed_excess_remaining=0.0)
    implied = 1.0 - (1.0 - 0.20) ** hr
    assert implied == pytest.approx(0.02)


def test_tighter_criteria_demand_a_lower_hazard_ratio():
    loose = required_hazard_ratio(0.20, 0.02, allowed_excess_remaining=0.5)
    tight = required_hazard_ratio(0.20, 0.02, allowed_excess_remaining=0.1)
    assert tight < loose


def test_hazard_ratio_rejects_a_reference_above_baseline():
    with pytest.raises(ValueError, match="reference risk"):
        required_hazard_ratio(0.02, 0.20)


def test_target_product_profile_reports_real_gaps():
    profile = target_product_profile()
    assert profile.required_substrate_removal == pytest.approx(0.8)
    assert profile.substrate_gap_factor > 1.0
    assert profile.myocardial_gap_factor > 6.0
    assert profile.required_lvef_gain_pp == pytest.approx(25.0)


def test_profile_handles_a_zero_current_benchmark():
    profile = target_product_profile(best_current_lvef_gain_pp=0.0001)
    assert profile.myocardial_gap_factor > 1000


def test_profile_rejects_implausible_ejection_fractions():
    with pytest.raises(ValueError, match="ejection fraction"):
        target_product_profile(normal_lvef_pp=20.0, advanced_hf_lvef_pp=25.0)


# --------------------------------------------------------------------------
# Ceiling, Monte Carlo, sensitivity.
# --------------------------------------------------------------------------


def test_architecture_never_eradicates_the_lesion():
    ceiling = architecture_ceiling()
    assert ceiling["eradication"] is False
    assert ceiling["asymptotic_removal_fraction"] < 1.0
    assert ceiling["target_reachable"] is True
    assert ceiling["months_to_target"] > 24.0


def test_an_impossible_target_is_reported_as_unreachable():
    ceiling = architecture_ceiling(required_removal=0.99)
    assert ceiling["target_reachable"] is False
    assert ceiling["months_to_target"] is None


def test_perfect_modules_bound_the_central_assumptions():
    sweep = best_case_horizon_sweep(horizons=(12.0, 24.0))
    for horizon in (12.0, 24.0):
        rows = sweep.loc[sweep["horizon_months"].eq(horizon)]
        perfect = rows.loc[rows["scenario"].eq("perfect_modules")].iloc[0]
        central = rows.loc[rows["scenario"].eq("central_assumptions")].iloc[0]
        assert perfect.substrate_removal_pct >= central.substrate_removal_pct


def test_monte_carlo_is_reproducible_and_bounded():
    first = monte_carlo_hypothesis(draws=25, seed=11)
    second = monte_carlo_hypothesis(draws=25, seed=11)
    pd.testing.assert_frame_equal(first, second)
    assert first["substrate_removal"].le(1.0).all()
    assert first["status"].eq(HYPOTHESIS_TAG).all()


def test_monte_carlo_does_not_claim_the_cure_bar_is_met(draws):
    """Even generous uniform priors do not buy a cure - the honest headline."""
    assert float(draws["meets_substrate_target"].mean()) < 0.5


def test_sensitivity_ranks_the_gate_as_a_top_driver(draws):
    ranked = sensitivity_summary(draws)
    assert ranked["status"].eq(HYPOTHESIS_TAG).all()
    assert "gate_specificity" in set(ranked.head(3)["module"])
    assert ranked["spearman_rho"].abs().is_monotonic_decreasing


def test_sensitivity_handles_a_constant_output():
    frame = pd.DataFrame(
        {
            "gate_specificity": [0.1, 0.5, 0.9, 0.3],
            "substrate_removal": [0.2, 0.2, 0.2, 0.2],
        }
    )
    ranked = sensitivity_summary(frame)
    assert ranked["spearman_rho"].eq(0.0).all()


# --------------------------------------------------------------------------
# Precedents, kill criteria, validation.
# --------------------------------------------------------------------------


def test_every_module_declares_a_precedent_and_a_failure_mode():
    frame = module_frame()
    assert len(frame) == len(REVERT1_MODULES)
    for spec in REVERT1_MODULES:
        assert spec.source_url.startswith("https://")
        assert spec.known_failure_mode
        assert spec.precedent_reference
    assert frame["source_url"].str.startswith("https://").all()


def test_no_module_claims_human_outcomes_evidence():
    allowed = {
        PrecedentLevel.HUMAN_MECHANISM,
        PrecedentLevel.MOUSE_GENETIC,
        PrecedentLevel.MOUSE_PHARMACOLOGIC,
        PrecedentLevel.IN_VITRO,
        PrecedentLevel.ENGINEERING_ONLY,
    }
    assert {spec.precedent for spec in REVERT1_MODULES} <= allowed


def test_kill_criteria_are_ordered_and_actionable():
    frame = kill_criteria_frame()
    assert list(frame["order"]) == list(range(1, len(KILL_CRITERIA) + 1))
    assert frame["cheapest_experiment"].str.len().gt(20).all()
    assert frame["falsifying_observation"].str.len().gt(20).all()
    assert "ex vivo" in frame.iloc[1]["cheapest_experiment"].lower()


def test_calcium_module_does_not_attempt_total_decalcification():
    mineral = next(spec for spec in REVERT1_MODULES if spec.module_id == "M6_mineral")
    assert "NOT" in mineral.known_failure_mode
    assert "nucleation" in mineral.target
    assert "micro-deposits" in mineral.role


def test_module_spec_requires_an_https_source():
    with pytest.raises(ValueError, match="https"):
        ModuleSpec("x", "t", "r", PrecedentLevel.IN_VITRO, "p", "http://example.com", "f")


def test_unknown_arm_is_rejected():
    with pytest.raises(ValueError, match="unknown arm"):
        arm_parameters("revert2_magic")


def test_every_arm_zeroes_the_modules_it_excludes():
    params = arm_parameters("efferocytosis_only")
    assert params.efferocytosis_gain > 0.0
    assert params.efflux_gain == 0.0
    assert params.crystal_solubilisation == 0.0
    assert params.plasma_driver_reduction == 0.0
    assert set(ARMS) == set(ablation_study()["arm"])


@pytest.mark.parametrize("field", ["gate_specificity", "efflux_gain", "nucleation_block"])
def test_parameters_reject_out_of_range_strengths(field):
    with pytest.raises(ValueError):
        ClearanceParameters(**{field: 1.5})
    with pytest.raises(ValueError):
        ClearanceParameters(**{field: -0.1})


def test_integrator_rejects_an_unstable_step():
    with pytest.raises(ValueError, match="step_months"):
        simulate_clearance(ClearanceParameters(), months=12.0, step_months=1.0)


def test_integrator_rejects_a_nonpositive_horizon():
    with pytest.raises(ValueError):
        simulate_clearance(ClearanceParameters(), months=0.0)


def test_lesion_state_rejects_an_empty_or_dead_lesion():
    with pytest.raises(ValueError):
        LesionState(free_cholesterol=0.0)
    with pytest.raises(ValueError):
        LesionState(macrophage_viability=0.0)


def test_delivery_scales_payloads_but_not_the_systemic_arm():
    params = ClearanceParameters(gate_specificity=0.5, plasma_driver_reduction=0.7)
    delivered = params.delivered()
    assert delivered.efflux_gain == pytest.approx(0.4)
    assert delivered.plasma_driver_reduction == pytest.approx(0.7)
    assert delivered.gate_specificity == 1.0


def test_core_change_at_month_zero_is_zero():
    result = simulate_clearance(arm_parameters("untreated"), months=12.0, step_months=0.5)
    assert result.necrotic_core_change_at(0.0) == pytest.approx(0.0)
    with pytest.raises(ValueError):
        result.necrotic_core_change_at(-1.0)


def test_cli_hypothesis_subcommand_is_quarantined_and_loud(tmp_path, capsys):
    from ihd_solutions.cli import main

    target = tmp_path / "cli-hyp"
    assert main(["hypothesis", "--output-dir", str(target), "--draws", "15"]) == 0
    out = capsys.readouterr().out
    assert DISCLAIMER in out
    assert "INVENTED HYPOTHESIS" in out
    assert "Wrote 10 artefacts" in out
    assert (target / "reports" / "revert1_design_report.txt").exists()


def test_cli_hypothesis_rejects_a_bad_horizon(tmp_path, capsys):
    from ihd_solutions.cli import main

    assert main(["hypothesis", "--output-dir", str(tmp_path / "bad"), "--months", "0"]) == 1
    assert "Error:" in capsys.readouterr().err
