"""Tests for the best-of-breed composite assembler.

The load-bearing assertions are the anti-double-counting rules: no summed
capabilities, no stacked molecular targets, no proof credit for an untested
combination, and no event credit for coverage additions.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ihd_solutions.catalog import load_catalog, validate_combination
from ihd_solutions.composite import (
    COMPOSITE_TAG,
    CRITERION_NODES,
    CRITERION_ORDER,
    DISCLAIMER,
    DISEASE_NODES,
    INVENTED_NODES,
    INVENTED_TAG,
    TARGET_OVERLAPS,
    Composite,
    Variant,
    _incidence,
    assemble_composite,
    complete_cure_nodes,
    component_frequency,
    composite_report,
    conflicts,
    evaluate_matrix,
    leave_one_out,
    node_of,
    run_composite_study,
    score_criteria,
    with_invented,
)
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.evidence import CreditMode


@pytest.fixture(scope="module")
def strict_composite(population, catalog, anchors, days) -> Composite:
    return assemble_composite(population, catalog, anchors, days)


@pytest.fixture(scope="module")
def matrix_results():
    return evaluate_matrix(population_ids=("stable_ccs", "ischemic_cm_heart_failure"))


# --------------------------------------------------------------------------
# Rule 1 and 2: no double counting.
# --------------------------------------------------------------------------


def test_one_component_per_node(strict_composite, catalog):
    nodes = [node_of(catalog[tid].exclusive_group) for tid in strict_composite.therapy_ids]
    nodes = [node for node in nodes if node is not None]
    assert len(nodes) == len(set(nodes))


def test_same_molecular_target_is_never_stacked(catalog):
    for overlap in TARGET_OVERLAPS:
        members = [tid for tid, therapy in catalog.items() if therapy.exclusive_group in overlap]
        assert len(members) >= 2
        first = catalog[members[0]]
        for other_id in members[1:]:
            reason = conflicts((first.therapy_id,), catalog[other_id], catalog)
            assert reason is not None


def test_pcsk9_editing_and_pcsk9_antibody_cannot_coexist(catalog):
    reason = conflicts(("verve102_base_editing",), catalog["evolocumab"], catalog)
    assert reason is not None
    assert "same molecular target" in reason


def test_conflicts_allows_a_genuinely_different_node(catalog):
    assert conflicts(("statin_high_intensity",), catalog["colchicine"], catalog) is None


def test_conflicts_rejects_a_second_anti_inflammatory(catalog):
    reason = conflicts(("colchicine",), catalog["ziltivekimab"], catalog)
    assert reason is not None
    assert "node" in reason


def test_assembly_is_a_valid_regimen(strict_composite, catalog):
    assert strict_composite.therapy_ids
    validate_combination(catalog, strict_composite.therapy_ids)


def test_every_step_has_positive_marginal_benefit(strict_composite):
    assert strict_composite.steps
    for step in strict_composite.steps:
        assert step.marginal_events_avoided_per_1000 > 0.0
    incidences = [step.incidence_pct for step in strict_composite.steps]
    assert incidences == sorted(incidences, reverse=True)


def test_invented_components_add_no_score(strict_composite, catalog, anchors, population, days):
    plain = score_criteria(strict_composite, catalog, anchors, population, days)
    invented = score_criteria(with_invented(strict_composite), catalog, anchors, population, days)
    for name in CRITERION_ORDER:
        assert plain[name] == invented[name], name


def test_invented_variant_does_not_change_the_event_curve(
    strict_composite, catalog, anchors, population, days
):
    extended = with_invented(strict_composite)
    assert extended.therapy_ids == strict_composite.therapy_ids
    assert extended.invented_components
    assert _incidence(
        extended.therapy_ids, catalog, anchors, population, days, CreditMode.STRICT
    ) == pytest.approx(
        _incidence(
            strict_composite.therapy_ids, catalog, anchors, population, days, CreditMode.STRICT
        )
    )


# --------------------------------------------------------------------------
# Rule 3: an untested combination earns no proof.
# --------------------------------------------------------------------------


def test_composite_proof_is_zero_even_with_proven_components(
    strict_composite, catalog, anchors, population, days
):
    scores = score_criteria(strict_composite, catalog, anchors, population, days)
    assert scores["randomised_proof"] == 0.0
    assert "does not transfer" in scores["proof_reason"]


def test_no_composite_is_ever_declared_a_cure(matrix_results):
    _, _, criteria = matrix_results
    assert not criteria["is_cure"].any()
    assert criteria["criteria_met"].lt(5).all()


def test_burden_and_restoration_use_max_not_sum(
    strict_composite, catalog, anchors, population, days
):
    scores = score_criteria(strict_composite, catalog, anchors, population, days)
    assert scores["burden_reversal"] <= 0.60
    assert scores["myocardial_restoration"] <= 0.50


# --------------------------------------------------------------------------
# Coverage pass: the myocardium must not be silently dropped.
# --------------------------------------------------------------------------


def test_coverage_pass_adds_the_myocardium_node(strict_composite):
    assert "biovat_hf_patch" in strict_composite.coverage_additions
    assert "biovat_hf_patch" in strict_composite.therapy_ids


def test_coverage_additions_do_not_move_events(population, catalog, anchors, days):
    composite = assemble_composite(population, catalog, anchors, days)
    without = tuple(tid for tid in composite.therapy_ids if tid not in composite.coverage_additions)
    full_risk = _incidence(
        composite.therapy_ids, catalog, anchors, population, days, CreditMode.STRICT
    )
    bare_risk = _incidence(without, catalog, anchors, population, days, CreditMode.STRICT)
    assert full_risk == pytest.approx(bare_risk, abs=1e-12)


def test_coverage_pass_is_idempotent(strict_composite, catalog):
    again, added = complete_cure_nodes(strict_composite.therapy_ids, catalog)
    assert again == strict_composite.therapy_ids
    assert added == ()


def test_coverage_pass_respects_conflicts(catalog):
    chosen = ("biovat_hf_patch",)
    completed, added = complete_cure_nodes(chosen, catalog)
    assert "biovat_hf_patch" not in added
    validate_combination(catalog, completed)


def test_criterion_nodes_are_real_nodes():
    for node in CRITERION_NODES:
        assert node in DISEASE_NODES


def test_coverage_addition_is_recorded_in_the_notes(strict_composite):
    assert any("zero event credit" in note for note in strict_composite.notes)


def test_disabling_completion_drops_the_myocardium(population, catalog, anchors, days):
    composite = assemble_composite(population, catalog, anchors, days, complete_nodes=False)
    assert composite.coverage_additions == ()
    assert "biovat_hf_patch" not in composite.therapy_ids


# --------------------------------------------------------------------------
# Mode dependence.
# --------------------------------------------------------------------------


def test_strict_mode_cannot_include_gene_editing(strict_composite, catalog):
    for therapy_id in strict_composite.therapy_ids:
        assert catalog[therapy_id].exclusive_group not in {"pcsk9_gene", "angptl3_gene"}


def test_surrogate_mode_admits_editing_and_scores_driver_control(
    population, catalog, anchors, days
):
    composite = assemble_composite(population, catalog, anchors, days, mode=CreditMode.SURROGATE)
    groups = {catalog[tid].exclusive_group for tid in composite.therapy_ids}
    assert groups & {"pcsk9_gene", "angptl3_gene"}
    scores = score_criteria(
        composite, catalog, anchors, population, days, mode=CreditMode.SURROGATE
    )
    assert scores["durable_driver_control"] == 1.0


def test_relaxed_modes_never_reduce_benefit(population, catalog, anchors, days):
    composite = assemble_composite(population, catalog, anchors, days)
    strict = _incidence(
        composite.therapy_ids, catalog, anchors, population, days, CreditMode.STRICT
    )
    surrogate = _incidence(
        composite.therapy_ids, catalog, anchors, population, days, CreditMode.SURROGATE
    )
    assert surrogate <= strict + 1e-12


def test_matrix_covers_every_mode_for_every_population(matrix_results):
    _, matrix, criteria = matrix_results
    assert len(matrix) == 2 * len(CreditMode)
    assert set(matrix["mode"]) == {mode.value for mode in CreditMode}
    assert len(criteria) == 2 * len(CreditMode) * len(Variant)


# --------------------------------------------------------------------------
# Empty composites, ablation, outputs.
# --------------------------------------------------------------------------


def test_a_population_with_no_credit_yields_an_empty_composite(catalog, anchors, populations, days):
    composite = assemble_composite(
        populations["non_ischemic_cardiomyopathy"], catalog, anchors, days, complete_nodes=False
    )
    assert composite.therapy_ids == ()
    assert composite.label == "(none)"
    assert any("empty composite is the honest result" in note for note in composite.notes)


def test_leave_one_out_ranks_contributions(strict_composite, catalog, anchors, population, days):
    frame = leave_one_out(strict_composite, catalog, anchors, population, days)
    assert len(frame) == len(strict_composite.therapy_ids)
    assert frame["contribution_per_1000"].is_monotonic_decreasing
    assert frame["contribution_per_1000"].ge(-1e-9).all()


def test_leave_one_out_shows_the_coverage_addition_is_event_neutral(
    strict_composite, catalog, anchors, population, days
):
    frame = leave_one_out(strict_composite, catalog, anchors, population, days)
    row = frame.loc[frame["removed"].eq("biovat_hf_patch")].iloc[0]
    assert row["contribution_per_1000"] == pytest.approx(0.0, abs=1e-9)
    assert not row["indispensable"]


def test_unknown_population_is_rejected():
    with pytest.raises(ValueError, match="unknown population"):
        evaluate_matrix(population_ids=("not_a_population",))


def test_frequency_counts_populations_not_rows(matrix_results):
    assembly, _, _ = matrix_results
    frequency = component_frequency(assembly)
    assert frequency["populations_selected"].le(2).all()
    assert frequency["status"].eq(COMPOSITE_TAG).all()


def test_frequency_handles_an_empty_assembly():
    assert component_frequency(pd.DataFrame()).empty


def test_node_map_covers_the_catalogue(catalog):
    mapped = {group for groups in DISEASE_NODES.values() for group in groups}
    for therapy in catalog.values():
        assert therapy.exclusive_group in mapped, therapy.therapy_id


def test_invented_nodes_are_disjoint_from_catalogue_nodes():
    for node in INVENTED_NODES:
        assert node not in DISEASE_NODES or node == "recurrence_tolerance"


def test_report_states_the_rules_and_the_verdict(matrix_results):
    assembly, matrix, criteria = matrix_results
    text = composite_report(assembly, matrix, criteria, component_frequency(assembly))
    assert DISCLAIMER in text
    assert "combines by maximum" in text
    assert "same molecular target" in text
    assert "untrialled combination scores zero" in text
    assert "Composites declared a cure: 0." in text


def test_outputs_are_tagged_and_quarantined(tmp_path):
    target = tmp_path / "comp"
    summary = run_composite_study(target, population_ids=("stable_ccs",))
    assert summary["status"] == "untrialled_composite"
    assert summary["cures_declared"] == 0
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 8
    assert not any(name.startswith("solutions") for name in written)
    allowed = {COMPOSITE_TAG, INVENTED_TAG}
    for path in (target / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        if "status" in frame.columns and not frame.empty:
            assert set(frame["status"]) <= allowed, path.name


def test_risk_normalisation_can_cross_the_bar(population, catalog, anchors, days):
    """The one criterion a best-of-breed composite actually clears."""
    composite = assemble_composite(population, catalog, anchors, days)
    scores = score_criteria(composite, catalog, anchors, population, days)
    assert scores["risk_normalisation"] >= CRITERION_THRESHOLD


def test_cli_composite_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    target = tmp_path / "cli-comp"
    assert main(["composite", "--output-dir", str(target), "--populations", "stable_ccs"]) == 0
    out = capsys.readouterr().out
    assert DISCLAIMER in out
    assert "Wrote 8 artefacts" in out


def test_cli_composite_rejects_an_unknown_population(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(["composite", "--output-dir", str(tmp_path / "bad"), "--populations", "not_real"])
    assert code == 1
    assert "Error:" in capsys.readouterr().err


def test_assembly_rejects_invalid_limits(population):
    with pytest.raises(ValueError):
        assemble_composite(population, max_components=0)
    with pytest.raises(ValueError):
        assemble_composite(population, min_marginal_per_1000=-1.0)


def test_node_of_returns_none_for_an_unmapped_group():
    assert node_of("not_a_group") is None


def test_composite_defaults_are_available():
    assert load_catalog()
    assert CRITERION_ORDER[-1] == "randomised_proof"
