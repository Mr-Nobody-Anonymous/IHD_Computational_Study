"""Regimen enumeration, ranking and Pareto analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ihd_solutions.catalog import TherapyClass
from ihd_solutions.evidence import CreditMode
from ihd_solutions.population import background_regimen_ids
from ihd_solutions.portfolio import (
    EnumerationSpec,
    compare_populations,
    enumerate_regimens,
    evaluate_regimens,
    expected_regimen_count,
    pareto_front,
    regimen_burden,
)


@pytest.mark.parametrize("max_optional", [0, 1, 2, 3])
def test_enumeration_matches_the_independent_combinatorial_count(catalog, max_optional):
    spec = EnumerationSpec(max_optional=max_optional)
    regimens = enumerate_regimens(catalog, spec)
    assert len(regimens) == expected_regimen_count(catalog, spec)


def test_background_is_always_present_and_first(catalog):
    regimens = enumerate_regimens(catalog, EnumerationSpec(max_optional=2))
    assert regimens[0].therapy_ids == tuple(sorted(background_regimen_ids()))
    for regimen in regimens:
        assert set(background_regimen_ids()) <= set(regimen.therapy_ids)


def test_no_regimen_stacks_two_agents_from_one_exclusive_group(catalog):
    regimens = enumerate_regimens(catalog, EnumerationSpec(max_optional=3))
    for regimen in regimens:
        groups = [catalog[therapy_id].exclusive_group for therapy_id in regimen.therapy_ids]
        assert len(groups) == len(set(groups)), regimen.label


def test_regimens_are_unique(catalog):
    regimens = enumerate_regimens(catalog, EnumerationSpec(max_optional=2))
    keys = [regimen.therapy_ids for regimen in regimens]
    assert len(keys) == len(set(keys))


def test_excluded_classes_are_absent(catalog):
    regimens = enumerate_regimens(catalog, EnumerationSpec(max_optional=3))
    for regimen in regimens:
        for therapy_id in regimen.therapy_ids:
            assert catalog[therapy_id].therapy_class is not TherapyClass.REGENERATIVE


def test_candidate_subset_restricts_enumeration(catalog):
    spec = EnumerationSpec(max_optional=2, candidate_ids=("evolocumab", "colchicine"))
    regimens = enumerate_regimens(catalog, spec)
    assert len(regimens) == expected_regimen_count(catalog, spec)
    allowed = set(background_regimen_ids()) | {"evolocumab", "colchicine"}
    for regimen in regimens:
        assert set(regimen.therapy_ids) <= allowed


def test_candidates_sharing_the_background_group_are_dropped(catalog):
    spec = EnumerationSpec(max_optional=1, candidate_ids=("statin_moderate_intensity",))
    regimens = enumerate_regimens(catalog, spec)
    assert len(regimens) == 1


def test_enumeration_validates_arguments(catalog):
    with pytest.raises(ValueError):
        EnumerationSpec(max_optional=-1)
    with pytest.raises(ValueError, match="background therapy not in catalog"):
        enumerate_regimens(catalog, EnumerationSpec(background=("nope",)))
    with pytest.raises(ValueError, match="unknown candidate ids"):
        enumerate_regimens(catalog, EnumerationSpec(candidate_ids=("nope",)))
    with pytest.raises(ValueError, match="exclusive group"):
        enumerate_regimens(catalog, EnumerationSpec(background=("evolocumab", "enlicitide_oral")))


def test_burden_is_additive(catalog, with_background):
    regimen = with_background("evolocumab")
    expected = sum(catalog[tid].burden_score for tid in regimen.therapy_ids)
    assert regimen_burden(regimen, catalog) == pytest.approx(expected)


def test_ranking_is_sorted_and_has_one_row_per_regimen(catalog, population, anchors, short_days):
    regimens = enumerate_regimens(catalog, EnumerationSpec(max_optional=1))
    frame, simulations = evaluate_regimens(
        regimens, catalog, population, anchors, short_days, modes=(CreditMode.STRICT,)
    )
    assert len(frame) == len(regimens)
    assert len(simulations) == len(regimens)
    values = frame["strict_events_avoided_per_1000"].to_numpy()
    assert np.all(np.diff(values) <= 1e-9)
    assert frame["strict_n_credited"].ge(0).all()


def test_all_modes_are_reported_side_by_side(catalog, population, anchors, short_days):
    regimens = enumerate_regimens(catalog, EnumerationSpec(max_optional=1))
    frame, simulations = evaluate_regimens(regimens, catalog, population, anchors, short_days)
    for mode in CreditMode:
        assert f"{mode.value}_events_avoided_per_1000" in frame
    assert len(simulations) == len(regimens) * len(CreditMode)


def test_surrogate_ranking_never_below_strict_for_ldl_only_additions(
    catalog, population, anchors, short_days
):
    regimens = enumerate_regimens(catalog, EnumerationSpec(max_optional=1))
    frame, _ = evaluate_regimens(regimens, catalog, population, anchors, short_days)
    lipid_rows = frame.loc[
        frame["therapy_ids"].str.contains("enlicitide_oral|azd0780_oral|verve102_base_editing")
    ]
    assert not lipid_rows.empty
    assert (
        lipid_rows["surrogate_events_avoided_per_1000"]
        >= lipid_rows["strict_events_avoided_per_1000"]
    ).all()


def test_evaluate_rejects_empty_inputs(catalog, population, anchors, short_days):
    with pytest.raises(ValueError, match="no regimens"):
        evaluate_regimens([], catalog, population, anchors, short_days)
    regimens = enumerate_regimens(catalog, EnumerationSpec(max_optional=0))
    with pytest.raises(ValueError, match="at least one credit mode"):
        evaluate_regimens(regimens, catalog, population, anchors, short_days, modes=())


def test_pareto_front_keeps_only_non_dominated_rows():
    frame = pd.DataFrame(
        {
            "regimen": ["a", "b", "c", "d"],
            "benefit": [10.0, 20.0, 20.0, 5.0],
            "burden_score": [1.0, 2.0, 3.0, 4.0],
        }
    )
    front = pareto_front(frame, "benefit")
    assert front["regimen"].tolist() == ["a", "b"]


def test_pareto_front_handles_ties_and_missing_columns():
    frame = pd.DataFrame({"regimen": ["a", "b"], "benefit": [1.0, 1.0], "burden_score": [1.0, 1.0]})
    assert len(pareto_front(frame, "benefit")) == 2
    with pytest.raises(ValueError, match="missing the benefit"):
        pareto_front(frame, "nope")


def test_compare_populations_exposes_population_specific_evidence(
    catalog, anchors, populations, short_days, with_background
):
    regimens = [with_background(), with_background("colchicine")]
    frame = compare_populations(
        regimens,
        catalog,
        populations,
        anchors,
        short_days,
        population_ids=["stable_ccs", "post_mi_pci", "ischemic_cm_heart_failure"],
    )
    assert len(frame) == 6
    assert set(frame["population_id"]) == {"stable_ccs", "post_mi_pci", "ischemic_cm_heart_failure"}
    assert frame["scenario_risk_pct"].between(0.0, 100.0).all()


def test_compare_populations_rejects_unknown_ids(
    catalog, anchors, populations, short_days, with_background
):
    with pytest.raises(ValueError, match="unknown population"):
        compare_populations(
            [with_background()], catalog, populations, anchors, short_days, population_ids=["nope"]
        )
