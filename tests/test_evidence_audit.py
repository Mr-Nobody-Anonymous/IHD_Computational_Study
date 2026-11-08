"""Tests for the evidence-layer uncertainty audit.

These pin a correction to the layer every construct inherits, so a future edit
cannot quietly restore point-estimate scoring of risk normalisation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ihd_solutions.catalog import load_catalog, validate_combination
from ihd_solutions.composite import assemble_composite
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.evidence import CreditMode, load_anchors
from ihd_solutions.evidence_audit import (
    DISCLAIMER,
    EFFECT_ESTIMATE_CRITERIA,
    ROBUSTNESS_BAR,
    criterion_ledger,
    evidence_audit_matrix,
    evidence_audit_report,
    mode_summary,
    risk_normalisation_robustness,
    run_evidence_audit_study,
)
from ihd_solutions.population import load_populations
from ihd_solutions.uncertainty import MonteCarloResult, monte_carlo_regimen
from ihd_solutions.units import time_grid

IDS = ("stable_ccs", "lifetime_exposure", "secondary_post_mi")


@pytest.fixture(scope="module")
def matrix() -> pd.DataFrame:
    return evidence_audit_matrix(population_ids=IDS, draws=200)


@pytest.fixture(scope="module")
def days() -> np.ndarray:
    return time_grid(5.0, 30.0)


# --------------------------------------------------------------------------
# The Monte Carlo now records what the criterion is computed from.
# --------------------------------------------------------------------------


def test_monte_carlo_records_per_draw_risks(days):
    catalog, anchors, populations = load_catalog(), load_anchors(), load_populations()
    population = populations["stable_ccs"]
    composite = assemble_composite(population, catalog, anchors, days, mode=CreditMode.STRICT)
    regimen = validate_combination(catalog, composite.therapy_ids)
    result = monte_carlo_regimen(
        regimen, catalog, population, anchors, days, mode=CreditMode.STRICT, draws=40
    )
    assert result.baseline_risk is not None
    assert result.treated_risk is not None
    assert result.baseline_risk.shape == (40,)
    assert result.treated_risk.shape == (40,)
    assert (result.baseline_risk > 0.0).all()
    assert (result.treated_risk <= result.baseline_risk + 1e-9).all()


def test_draws_helper_returns_none_when_risks_were_not_recorded():
    legacy = MonteCarloResult(
        regimen_label="x",
        population_id="p",
        mode=CreditMode.STRICT,
        events_avoided_per_1000=np.zeros(3),
        time_averaged_ldl_mg_dl=np.zeros(3),
    )
    assert legacy.risk_normalisation_draws() is None


def test_draws_formula_matches_the_point_estimate_formula():
    """Same quantity, so the nominal score and the distribution are comparable."""
    result = MonteCarloResult(
        regimen_label="x",
        population_id="p",
        mode=CreditMode.STRICT,
        events_avoided_per_1000=np.zeros(3),
        time_averaged_ldl_mg_dl=np.zeros(3),
        baseline_risk=np.array([0.10, 0.0198, 0.01]),
        treated_risk=np.array([0.0198, 0.01, 0.005]),
    )
    values = result.risk_normalisation_draws(0.0198)
    assert values is not None
    # baseline 0.10 treated at reference -> fully normalised
    assert values[0] == pytest.approx(1.0)
    # baseline at or below reference -> score is zero, not negative or undefined
    assert values[1] == 0.0
    assert values[2] == 0.0


def test_draws_are_clipped_to_the_unit_interval():
    result = MonteCarloResult(
        regimen_label="x",
        population_id="p",
        mode=CreditMode.STRICT,
        events_avoided_per_1000=np.zeros(2),
        time_averaged_ldl_mg_dl=np.zeros(2),
        baseline_risk=np.array([0.10, 0.10]),
        treated_risk=np.array([0.20, 0.001]),
    )
    values = result.risk_normalisation_draws()
    assert values is not None
    assert (values >= 0.0).all()
    assert (values <= 1.0).all()


# --------------------------------------------------------------------------
# The correction itself.
# --------------------------------------------------------------------------


def test_strict_credit_has_no_robust_criterion_left(matrix):
    strict = matrix.loc[matrix["mode"].eq("strict")]
    assert strict["nominal_criteria_met"].max() == 1
    assert strict["robust_criteria_met"].max() == 0
    assert strict["robust_pass"].sum() == 0


def test_the_strict_nominal_pass_is_retracted(matrix):
    row = matrix.loc[matrix["population_id"].eq("stable_ccs") & matrix["mode"].eq("strict")].iloc[0]
    assert row["nominal_pass"]
    assert not row["robust_pass"]
    assert row["nominal_is_misleading"]
    assert row["probability_of_target"] < 0.35


def test_exploratory_credit_mostly_survives(matrix):
    exploratory = matrix.loc[matrix["mode"].eq("exploratory")]
    assert exploratory["robust_pass"].sum() >= 2
    assert exploratory["robust_criteria_met"].max() == 2


def test_corrected_count_never_exceeds_the_nominal_count(matrix):
    assert (matrix["robust_criteria_met"] <= matrix["nominal_criteria_met"]).all()


def test_no_population_is_declared_a_cure(matrix):
    assert not matrix["is_cure"].any()
    assert matrix["robust_criteria_met"].lt(5).all()


def test_matrix_covers_every_requested_mode(matrix):
    assert set(matrix["mode"]) == {m.value for m in CreditMode}
    assert len(matrix) == len(IDS) * len(CreditMode)


def test_unknown_population_is_rejected():
    with pytest.raises(ValueError, match="unknown population"):
        evidence_audit_matrix(population_ids=("not_real",), draws=50)


# --------------------------------------------------------------------------
# Median instability and borderline honesty.
# --------------------------------------------------------------------------


def test_lifetime_exposure_is_bimodal_not_a_clean_pass(matrix):
    row = matrix.loc[
        matrix["population_id"].eq("lifetime_exposure") & matrix["mode"].eq("exploratory")
    ].iloc[0]
    assert row["nominal_risk_normalisation"] >= CRITERION_THRESHOLD
    assert row["mass_at_floor"] > 0.2
    assert row["mass_at_ceiling"] > 0.2
    assert row["median_uninformative"]
    assert not row["robust_pass"]


def test_borderline_cases_are_not_counted_as_retractions(matrix):
    for row in matrix.itertuples(index=False):
        if row.borderline:
            assert not row.nominal_is_misleading


def test_borderline_band_widens_with_fewer_draws(days):
    catalog, anchors, populations = load_catalog(), load_anchors(), load_populations()
    population = populations["secondary_post_mi"]
    few = risk_normalisation_robustness(
        population, CreditMode.SURROGATE, catalog, anchors, days, draws=60
    )
    many = risk_normalisation_robustness(
        population, CreditMode.SURROGATE, catalog, anchors, days, draws=400
    )
    assert few.monte_carlo_standard_error > 0.0
    assert many.monte_carlo_standard_error < few.monte_carlo_standard_error * 2.0


def test_standard_error_is_zero_at_the_poles(days):
    catalog, anchors, populations = load_catalog(), load_anchors(), load_populations()
    result = risk_normalisation_robustness(
        populations["lifetime_exposure"], CreditMode.STRICT, catalog, anchors, days, draws=60
    )
    assert result.probability_of_target == 0.0
    assert result.monte_carlo_standard_error == 0.0
    assert not result.borderline


def test_robustness_is_reproducible(days):
    catalog, anchors, populations = load_catalog(), load_anchors(), load_populations()
    args = (populations["stable_ccs"], CreditMode.STRICT, catalog, anchors, days)
    first = risk_normalisation_robustness(*args, draws=80, seed=5)
    second = risk_normalisation_robustness(*args, draws=80, seed=5)
    assert first.probability_of_target == second.probability_of_target


# --------------------------------------------------------------------------
# Scope: only the effect-estimate criterion is audited.
# --------------------------------------------------------------------------


def test_only_risk_normalisation_is_an_effect_estimate():
    assert {"risk_normalisation"} == EFFECT_ESTIMATE_CRITERIA
    ledger = criterion_ledger().set_index("criterion")
    assert ledger.loc["risk_normalisation", "audited"]
    assert ledger["audited"].sum() == 1
    assert len(ledger) == 5


def test_ledger_explains_why_the_others_are_not_audited():
    ledger = criterion_ledger()
    assert ledger["reason"].str.len().gt(10).all()
    assert ledger["caveat"].str.len().gt(10).all()


def test_driver_control_is_flagged_as_structural_not_durable_proof():
    row = criterion_ledger().set_index("criterion").loc["durable_driver_control"]
    assert row["kind"] == "structural_flag"
    assert "lifetime durability" in row["caveat"]


# --------------------------------------------------------------------------
# Report, summary and quarantine.
# --------------------------------------------------------------------------


def test_mode_summary_counts_add_up(matrix):
    summary = mode_summary(matrix).set_index("mode")
    for mode in summary.index:
        group = matrix.loc[matrix["mode"].eq(mode)]
        assert summary.loc[mode, "nominal_risk_passes"] == int(group["nominal_pass"].sum())
        assert summary.loc[mode, "robust_risk_passes"] == int(group["robust_pass"].sum())


def test_report_leads_with_the_disclaimer_and_states_the_correction(matrix):
    text = evidence_audit_report(matrix, mode_summary(matrix))
    assert DISCLAIMER in text
    assert "The evidence layer was partly nominal too" in text
    assert "A MEDIAN IS NOT A PROBABILITY EITHER" in text
    assert "zero of five" in text


def test_outputs_are_counted_and_quarantined(tmp_path):
    target = tmp_path / "evidence_audit"
    summary = run_evidence_audit_study(target, population_ids=("stable_ccs",), draws=120)
    assert summary["status"] == "evidence_layer_partly_nominal_corrected"
    assert summary["cures_declared"] == 0
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 6
    assert not any(name.startswith("solutions") for name in written)


def test_every_table_is_tagged(tmp_path):
    target = tmp_path / "evidence_audit"
    run_evidence_audit_study(target, population_ids=("stable_ccs",), draws=120)
    for path in (target / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        assert "status" in frame.columns, path.name
        assert frame["status"].eq("evidence_layer_uncertainty_audit").all(), path.name


def test_cli_audit_evidence_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    target = tmp_path / "cli-evidence-audit"
    code = main(
        [
            "audit-evidence",
            "--output-dir",
            str(target),
            "--populations",
            "stable_ccs",
            "--draws",
            "100",
        ]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert DISCLAIMER in out
    assert "Wrote 6 artefacts" in out


def test_robustness_bar_matches_the_cure_threshold():
    assert ROBUSTNESS_BAR == CRITERION_THRESHOLD
