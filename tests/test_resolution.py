"""Tests for the step-by-step residual resolution.

Two findings are load-bearing and pinned here: no single clearance module suffices,
and the nominal two-factor remuscularisation answer is not robust. The second is a
deliberate retraction of an earlier conclusion, so it must not silently regress.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ihd_solutions.catalog import load_catalog
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.evidence import EvidenceTier
from ihd_solutions.gap_solver import (
    REMUSCULARISATION_FACTORS,
    SCAR_RECOVERABLE_LVEF_PP,
    SOURCE_REDUCTION_PRECEDENT_MAX,
)
from ihd_solutions.resolution import (
    DISCLAIMER,
    MODE_ORDER,
    OUTCOMES_TIERS,
    RESOLUTION_TAG,
    cheapest_clearance_point,
    chronic_penalty_sensitivity,
    clearance_at,
    clearance_operating_points,
    evidence_upgrade_levers,
    minimum_potency,
    minimum_robust_set,
    module_leverage_table,
    population_closure,
    potency_scaled,
    remuscularisation_robustness,
    resolution_report,
    robust_remuscularisation_search,
    run_resolution_study,
)


@pytest.fixture(scope="module")
def matrix() -> pd.DataFrame:
    rows = []
    for pop, baseline, risks in (
        ("closed_pop", 13.93, (4.14, 2.84, 2.45)),
        ("surrogate_pop", 20.15, (6.71, 5.10, 4.00)),
        ("exploratory_pop", 25.92, (17.11, 12.00, 4.78)),
        ("open_pop", 32.97, (27.97, 24.00, 22.00)),
    ):
        for mode, treated in zip(MODE_ORDER, risks, strict=True):
            rows.append(
                {
                    "population_id": pop,
                    "mode": mode,
                    "baseline_risk_pct": baseline,
                    "composite_risk_pct": treated,
                    "therapy_ids": "statin_high_intensity;colchicine",
                }
            )
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Step 1: evidence maturity, not biology.
# --------------------------------------------------------------------------


def test_closure_reports_the_first_mode_that_works(matrix):
    closure = population_closure(matrix).set_index("population_id")
    assert closure.loc["closed_pop", "first_closing_mode"] == "strict"
    assert closure.loc["surrogate_pop", "first_closing_mode"] == "surrogate"
    assert closure.loc["exploratory_pop", "first_closing_mode"] == "exploratory"
    assert closure.loc["open_pop", "first_closing_mode"] == "none"


def test_closure_is_monotone_in_permissiveness(matrix):
    closure = population_closure(matrix).set_index("population_id")
    rank = {mode: index for index, mode in enumerate(MODE_ORDER)}
    for _, row in closure.iterrows():
        statuses = [row[f"{mode}_status"] for mode in MODE_ORDER]
        if "closed" in statuses:
            first = statuses.index("closed")
            assert all(s == "closed" for s in statuses[first:])
    assert rank["strict"] < rank["exploratory"]


def test_a_closed_population_has_zero_residual(matrix):
    closure = population_closure(matrix).set_index("population_id")
    assert closure.loc["closed_pop", "strict_residual_pct"] == 0.0
    assert closure.loc["open_pop", "strict_residual_pct"] > 0.0


def test_levers_exclude_agents_that_already_have_outcomes_proof():
    """Aspirin and bempedoic acid must never be listed as needing a trial."""
    matrix = pd.read_csv("results/composite/tables/composite_test_matrix.csv")
    levers = evidence_upgrade_levers(matrix)
    genuine = levers.loc[~levers["already_outcomes_tier"]]
    catalog = load_catalog()
    for therapy_id in genuine["therapy_id"]:
        assert catalog[therapy_id].evidence_tier not in OUTCOMES_TIERS
    assert "aspirin_secondary_prevention" not in set(genuine["therapy_id"])
    assert "bempedoic_acid" not in set(genuine["therapy_id"])


def test_reordering_is_labelled_distinctly_from_an_evidence_gap():
    matrix = pd.read_csv("results/composite/tables/composite_test_matrix.csv")
    levers = evidence_upgrade_levers(matrix)
    reordering = levers.loc[levers["already_outcomes_tier"]]
    assert not reordering.empty
    assert reordering["reason"].eq("composition_reordering").all()
    genuine = levers.loc[~levers["already_outcomes_tier"]]
    assert genuine["reason"].eq("evidence_maturation_required").all()


def test_the_top_lever_is_an_editing_agent_without_outcomes_data():
    matrix = pd.read_csv("results/composite/tables/composite_test_matrix.csv")
    levers = evidence_upgrade_levers(matrix)
    top = levers.iloc[0]
    assert not top["already_outcomes_tier"]
    assert top["populations_unlocked"] == levers["populations_unlocked"].max()
    assert EvidenceTier(top["current_tier"]) not in OUTCOMES_TIERS


def test_levers_handle_a_matrix_that_closes_everywhere(matrix):
    simple = matrix.loc[matrix["population_id"].eq("closed_pop")]
    assert evidence_upgrade_levers(simple).empty


# --------------------------------------------------------------------------
# Step 2: no single module suffices; the operating point is joint.
# --------------------------------------------------------------------------


def test_potency_one_reproduces_the_central_parameters():
    from ihd_solutions.resolution import CENTRAL

    scaled = potency_scaled(1.0)
    assert scaled == CENTRAL


def test_potency_is_capped_at_one_per_module():
    scaled = potency_scaled(10.0)
    for name in scaled.__dataclass_fields__:
        assert getattr(scaled, name) <= 1.0


def test_potency_rejects_a_nonpositive_multiplier():
    with pytest.raises(ValueError):
        potency_scaled(0.0)


def test_no_single_module_reaches_the_target():
    """The correction to 'near-perfect delivery': it is joint performance."""
    leverage = module_leverage_table()
    assert not leverage["reaches_target_alone"].any()
    assert leverage["removal_alone_pct"].max() < 100.0 * CRITERION_THRESHOLD


def test_gate_alone_is_not_the_answer():
    leverage = module_leverage_table().set_index("module")
    assert leverage.loc["gate_specificity", "removal_alone_pct"] < 100.0 * CRITERION_THRESHOLD
    # The strongest single lever is not the gate.
    assert leverage["removal_alone_pct"].idxmax() != "gate_specificity"


def test_joint_potency_can_reach_what_no_single_module_can():
    reached = 100.0 * clearance_at(2.0, 0.35, 24.0)
    assert reached >= 100.0 * CRITERION_THRESHOLD


def test_clearance_is_monotone_in_potency_and_source():
    values = [clearance_at(k, 0.35, 24.0) for k in (1.0, 1.2, 1.5, 2.0)]
    assert values == sorted(values)
    sources = [clearance_at(1.5, s, 24.0) for s in (0.0, 0.2, 0.35, 0.5)]
    assert sources == sorted(sources)


def test_more_time_relaxes_the_potency_requirement():
    short = minimum_potency(0.35, 24.0)
    long = minimum_potency(0.35, 120.0)
    assert short is not None and long is not None
    assert long < short


def test_without_source_reduction_two_years_is_unreachable():
    assert minimum_potency(0.0, 24.0) is None


def test_operating_points_flag_precedent_and_feasibility():
    operating = clearance_operating_points()
    assert operating["source_within_precedent"].all()
    assert operating["admissible"].any()
    unreachable = operating.loc[
        operating["source_reduction"].eq(0.0) & operating["horizon_months"].eq(24.0)
    ].iloc[0]
    assert not unreachable["feasible"]
    assert not unreachable["admissible"]


def test_cheapest_point_prefers_the_shortest_admissible_horizon():
    operating = clearance_operating_points()
    cheapest = cheapest_clearance_point(operating)
    assert cheapest is not None
    assert cheapest["horizon_months"] == 24.0
    assert cheapest["source_reduction"] <= SOURCE_REDUCTION_PRECEDENT_MAX
    assert cheapest["minimum_potency_multiplier"] > 1.0


def test_cheapest_point_returns_none_when_nothing_is_admissible():
    operating = clearance_operating_points()
    operating["admissible"] = False
    assert cheapest_clearance_point(operating) is None


def test_the_reported_operating_point_actually_hits_the_target():
    cheapest = cheapest_clearance_point(clearance_operating_points())
    reached = clearance_at(
        cheapest["minimum_potency_multiplier"],
        cheapest["source_reduction"],
        cheapest["horizon_months"],
    )
    assert reached >= CRITERION_THRESHOLD - 1e-6


# --------------------------------------------------------------------------
# Step 3: the retraction. Nominal arithmetic was luck.
# --------------------------------------------------------------------------


def test_the_nominal_two_factor_set_is_not_robust():
    """Retraction guard: +20.3 nominal, but almost never under uncertainty."""
    two = tuple(
        f
        for f in REMUSCULARISATION_FACTORS
        if f.name in {"endothelial_co_delivery", "epicardial_stromal_co_delivery"}
    )
    result = remuscularisation_robustness(two, draws=3000, chronic_penalty_low=1.0)
    assert result.nominal_gain_pp >= 20.0
    assert result.probability_of_target < 0.05
    assert not result.robust


def test_no_factor_combination_is_robust_at_eighty_percent():
    assert minimum_robust_set(draws=3000) == ()


def test_even_all_four_factors_fall_short():
    result = remuscularisation_robustness(
        REMUSCULARISATION_FACTORS, draws=4000, chronic_penalty_low=1.0
    )
    assert result.probability_of_target < 0.8
    assert result.nominal_gain_pp == pytest.approx(SCAR_RECOVERABLE_LVEF_PP)


def test_a_chronic_haircut_strictly_reduces_the_probability():
    no_haircut = remuscularisation_robustness(
        REMUSCULARISATION_FACTORS, draws=4000, chronic_penalty_low=1.0
    )
    haircut = remuscularisation_robustness(
        REMUSCULARISATION_FACTORS, draws=4000, chronic_penalty_low=0.5
    )
    assert haircut.probability_of_target < no_haircut.probability_of_target


def test_sensitivity_reports_the_best_available_set_when_none_is_robust():
    """Guards against reporting '0 factors' as if no factors were needed."""
    sensitivity = chronic_penalty_sensitivity(penalties=(1.0, 0.5), draws=2000)
    assert not sensitivity["robust_set_exists"].any()
    assert sensitivity["best_available_n_factors"].eq(len(REMUSCULARISATION_FACTORS)).all()
    assert sensitivity["best_probability_of_target"].gt(0.0).all()


def test_sensitivity_probability_falls_as_the_haircut_worsens():
    sensitivity = chronic_penalty_sensitivity(penalties=(1.0, 0.75, 0.5, 0.25), draws=2000)
    values = list(sensitivity["best_probability_of_target"])
    assert values == sorted(values, reverse=True)


def test_robustness_is_reproducible():
    first = remuscularisation_robustness(REMUSCULARISATION_FACTORS, draws=500, seed=42)
    second = remuscularisation_robustness(REMUSCULARISATION_FACTORS, draws=500, seed=42)
    assert first.probability_of_target == second.probability_of_target


def test_robustness_search_covers_every_subset():
    search = robust_remuscularisation_search(draws=300)
    assert len(search) == 2 ** len(REMUSCULARISATION_FACTORS)
    assert search["probability_of_target"].between(0.0, 1.0).all()
    assert search["status"].eq(RESOLUTION_TAG).all()


def test_gains_never_exceed_full_scar_replacement():
    result = remuscularisation_robustness(REMUSCULARISATION_FACTORS, draws=2000)
    assert result.median_gain_pp <= SCAR_RECOVERABLE_LVEF_PP
    assert result.p05_gain_pp <= result.median_gain_pp


def test_robustness_rejects_invalid_inputs():
    with pytest.raises(ValueError):
        remuscularisation_robustness((), draws=0)
    with pytest.raises(ValueError):
        remuscularisation_robustness((), chronic_penalty_low=1.5)


# --------------------------------------------------------------------------
# Report and outputs.
# --------------------------------------------------------------------------


def test_report_states_the_retraction(matrix):
    text = resolution_report(matrix)
    assert DISCLAIMER in text
    assert "imprecise" in text
    assert "arithmetic luck" in text
    assert "None of this is evidence" in text


def test_outputs_are_tagged_and_quarantined(tmp_path, matrix):
    target = tmp_path / "res"
    summary = run_resolution_study(target, matrix=matrix)
    assert summary["status"] == "scenario_only"
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 9
    for path in (target / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].eq(RESOLUTION_TAG).all(), path.name


def test_summary_records_that_remuscularisation_is_unsolved(tmp_path, matrix):
    summary = run_resolution_study(tmp_path / "res", matrix=matrix)
    assert summary["minimum_robust_remuscularisation_set"] == ()
    assert summary["cheapest_clearance_point"] is not None


def test_missing_matrix_is_a_clear_error(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(FileNotFoundError, match="composite"):
        run_resolution_study(tmp_path / "res")


def test_cli_resolve_subcommand(tmp_path, capsys, matrix):
    from ihd_solutions.cli import main

    source = tmp_path / "matrix.csv"
    matrix.to_csv(source, index=False)
    target = tmp_path / "cli-res"
    assert main(["resolve", "--output-dir", str(target), "--composite-matrix", str(source)]) == 0
    out = capsys.readouterr().out
    assert DISCLAIMER in out
    assert "Wrote 9 artefacts" in out
