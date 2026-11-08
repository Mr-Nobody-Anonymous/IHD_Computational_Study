"""Base-editing scenarios: calibrated pharmacodynamics, refused safety scoring."""

from __future__ import annotations

import numpy as np
import pytest

from ihd_solutions import editing


def test_calibration_reproduces_the_heart2_anchor():
    reduction = editing.ldl_reduction_from_knockdown(editing.HEART2_PCSK9_KNOCKDOWN)
    assert float(reduction) == pytest.approx(editing.HEART2_LDL_REDUCTION, abs=1e-9)


def test_calibration_sharpness_is_positive_and_reusable():
    sharpness = editing.calibrate_sharpness()
    assert sharpness > 0
    value = editing.ldl_reduction_from_knockdown(0.88, sharpness=sharpness)
    assert float(value) == pytest.approx(0.62, abs=1e-9)


def test_reduction_is_monotone_and_saturating():
    knockdown = np.linspace(0.0, 1.0, 101)
    reduction = editing.ldl_reduction_from_knockdown(knockdown)
    assert np.all(np.diff(reduction) > 0)
    assert reduction[0] == pytest.approx(0.0)
    assert reduction[-1] < editing.MAX_LDL_REDUCTION


@pytest.mark.parametrize("value", [-0.01, 1.01, np.nan])
def test_out_of_range_knockdown_is_rejected(value):
    with pytest.raises(ValueError):
        editing.ldl_reduction_from_knockdown(value)


def test_calibration_rejects_impossible_anchor():
    with pytest.raises(ValueError, match="saturation ceiling"):
        editing.calibrate_sharpness(anchor_reduction=0.9, max_reduction=0.75)


def test_permanent_scenario_holds_its_effect():
    scenario = editing.EditingScenario("permanent", edited_hepatocyte_fraction=0.95)
    frame = editing.simulate_editing(scenario, horizon_years=10.0, step_days=60.0)
    assert frame["pcsk9_knockdown_pct"].iloc[-1] == pytest.approx(
        frame["pcsk9_knockdown_pct"].max()
    )
    assert frame["durability_assumption"].unique().tolist() == ["permanent_scenario"]


def test_waning_scenario_loses_effect_over_time():
    scenario = editing.EditingScenario(
        "waning", edited_hepatocyte_fraction=0.95, durability_half_life_years=2.0
    )
    frame = editing.simulate_editing(scenario, horizon_years=10.0, step_days=30.0)
    assert frame["pcsk9_knockdown_pct"].iloc[-1] < frame["pcsk9_knockdown_pct"].max() * 0.5
    assert frame["ldl_c_mg_dl"].iloc[-1] > frame["ldl_c_mg_dl"].min()


def test_higher_editing_efficiency_gives_more_lowering():
    low = editing.simulate_editing(
        editing.EditingScenario("low", edited_hepatocyte_fraction=0.5), horizon_years=2.0
    )
    high = editing.simulate_editing(
        editing.EditingScenario("high", edited_hepatocyte_fraction=0.95), horizon_years=2.0
    )
    assert high["ldl_reduction_pct"].max() > low["ldl_reduction_pct"].max()


def test_no_effect_before_onset():
    scenario = editing.EditingScenario("onset", onset_days=14.0, plateau_days=28.0)
    frame = editing.simulate_editing(scenario, horizon_years=1.0, step_days=1.0)
    early = frame.loc[frame["day"] < 14.0]
    assert (early["pcsk9_knockdown_pct"] == 0.0).all()


def test_safety_is_never_assessed_numerically():
    scenario = editing.EditingScenario("x")
    frame = editing.simulate_editing(scenario, horizon_years=1.0)
    assert not frame["clinical_safety_assessed"].any()
    assessment = editing.off_target_assessment()
    assert assessment["numeric_off_target_score"] is None
    assert assessment["status"] == "requires_experimental_assay"
    assert len(assessment["required_assays"]) >= 5


def test_required_assays_cover_dna_rna_and_bystanders():
    assays = " ".join(editing.required_safety_assays()).lower()
    for keyword in ("dna", "rna", "bystander", "germline", "long-term"):
        assert keyword in assays


def test_bystander_counting_is_exact():
    report = editing.bystander_report("GCAATACAGATCAACGATCG", "TGG", intended_position=6)
    assert report["adenine_positions"] == [3, 4, 6, 8, 10, 13, 14, 17]
    assert report["adenines_in_window"] == [4, 6, 8]
    assert report["bystander_positions"] == [4, 8]
    assert report["n_bystanders"] == 2


def test_bystander_report_without_intended_position_counts_all_window_adenines():
    report = editing.bystander_report("GCAATACAGATCAACGATCG", "TGG")
    assert report["n_bystanders"] == 3
    assert report["intended_position"] is None


def test_adenine_free_window_has_no_bystanders():
    report = editing.bystander_report("GCGTACCGGATCAACGATCG", "AGG", window=(1, 4))
    assert report["adenines_in_window"] == []
    assert report["n_bystanders"] == 0


@pytest.mark.parametrize(
    "spacer,pam",
    [
        ("ACGT", "TGG"),
        ("GCAATACAGATCAACGATCGN", "TGG"),
        ("GCAATACAGATCAACGATCN", "TGG"),
        ("GCAATACAGATCAACGATCG", "TG"),
        ("GCAATACAGATCAACGATCG", "TAA"),
        ("GCAATACAGATCAACGATCG", "TN G"),
    ],
)
def test_guide_validation_rejects_bad_sequences(spacer, pam):
    with pytest.raises(ValueError):
        editing.validate_guide(spacer, pam)


def test_guide_validation_normalises_case_and_whitespace():
    spacer, pam = editing.validate_guide(" gcaatacagatcaacgatcg ", " tgg ")
    assert spacer == "GCAATACAGATCAACGATCG"
    assert pam == "TGG"


def test_intended_position_must_be_an_adenine():
    with pytest.raises(ValueError, match="not an adenine"):
        editing.bystander_report("GCAATACAGATCAACGATCG", "TGG", intended_position=1)


@pytest.mark.parametrize("window", [(0, 5), (5, 4), (1, 21)])
def test_invalid_editing_windows_are_rejected(window):
    with pytest.raises(ValueError, match="editing window"):
        editing.bystander_report("GCAATACAGATCAACGATCG", "TGG", window=window)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"scenario_id": ""},
        {"scenario_id": "x", "edited_hepatocyte_fraction": 1.2},
        {"scenario_id": "x", "knockdown_ceiling": -0.1},
        {"scenario_id": "x", "plateau_days": 0.0},
        {"scenario_id": "x", "onset_days": 100.0},
        {"scenario_id": "x", "durability_half_life_years": 0.0},
    ],
)
def test_scenario_validation(kwargs):
    with pytest.raises(ValueError):
        editing.EditingScenario(**kwargs)


def test_permanent_flag_and_plateau_knockdown():
    scenario = editing.EditingScenario("x", edited_hepatocyte_fraction=0.8, knockdown_ceiling=0.9)
    assert scenario.is_permanent
    assert scenario.plateau_knockdown == pytest.approx(0.72)
    assert not editing.EditingScenario("y", durability_half_life_years=3.0).is_permanent


@pytest.mark.parametrize(
    "kwargs", [{"horizon_years": 0.0}, {"step_days": 0.0}, {"baseline_ldl_mg_dl": 0.0}]
)
def test_simulation_argument_validation(kwargs):
    with pytest.raises(ValueError):
        editing.simulate_editing(editing.EditingScenario("x"), **kwargs)
