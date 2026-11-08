"""Scenario populations: baseline validation and loader behaviour."""

from __future__ import annotations

import math
from dataclasses import replace

import pytest

from ihd_solutions.population import (
    _as_bool,
    background_regimen_ids,
    load_populations,
)


def test_all_populations_load_with_unique_ids(populations):
    assert len(populations) >= 12
    assert len({p.population_id for p in populations.values()}) == len(populations)


def test_key_trial_populations_are_present(populations):
    expected = {
        "secondary_post_mi",
        "post_mi_pci",
        "post_mi_hscrp_high",
        "stable_ccs",
        "ascvd_ckd_inflammation",
        "statin_intolerant_ascvd",
        "t2d_ascvd",
        "ascvd_obesity_no_diabetes",
        "ascvd_high_triglycerides",
        "ascvd_high_lpa",
        "hefh_primary",
        "ischemic_cm_heart_failure",
    }
    assert expected <= set(populations)


def test_log_rate_sigma_reproduces_a_log_symmetric_interval(populations):
    population = populations["stable_ccs"]
    sigma = population.log_rate_sigma
    centre = math.sqrt(population.baseline_rate_ci_low * population.baseline_rate_ci_high)
    assert math.exp(math.log(population.baseline_annual_mace_rate) - 1.959963984540054 * sigma) == (
        pytest.approx(population.baseline_annual_mace_rate**2 / population.baseline_rate_ci_high)
    )
    assert centre > 0


def test_disease_flags_are_parsed_as_booleans(populations):
    assert populations["ascvd_ckd_inflammation"].has_ckd
    assert populations["t2d_ascvd"].has_diabetes
    assert populations["ischemic_cm_heart_failure"].has_heart_failure
    assert not populations["hefh_primary"].has_ascvd


def test_background_regimen_is_statin_plus_antithrombotic(catalog):
    ids = background_regimen_ids()
    assert set(ids) <= set(catalog)
    assert "statin_high_intensity" in ids
    assert "aspirin_secondary_prevention" in ids


@pytest.mark.parametrize(
    "changes",
    [
        {"population_id": ""},
        {"baseline_annual_mace_rate": 0.0},
        {"baseline_annual_mace_rate": 1.0},
        {"baseline_rate_ci_low": 0.5},
        {"baseline_rate_ci_high": 0.001},
        {"ldl_c_mg_dl": 0.0},
        {"hscrp_mg_l": -1.0},
        {"systolic_bp_mmhg": 0.0},
    ],
)
def test_population_validation(populations, changes):
    with pytest.raises(ValueError):
        replace(populations["stable_ccs"], **changes)


def test_boolean_parsing_accepts_common_spellings_and_rejects_nonsense():
    assert _as_bool(True) is True
    assert _as_bool("TRUE") is True
    assert _as_bool("1") is True
    assert _as_bool("no") is False
    assert _as_bool(False) is False
    with pytest.raises(ValueError):
        _as_bool("perhaps")


def test_loader_rejects_missing_and_malformed_files(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_populations(tmp_path / "missing.csv")
    path = tmp_path / "bad.csv"
    path.write_text("population_id,label\na,b\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing columns"):
        load_populations(path)


def test_loader_rejects_duplicate_population_ids(tmp_path, populations):
    import pandas as pd

    frame = pd.DataFrame(
        [
            {
                "population_id": p.population_id,
                "label": p.label,
                "baseline_annual_mace_rate": p.baseline_annual_mace_rate,
                "baseline_rate_ci_low": p.baseline_rate_ci_low,
                "baseline_rate_ci_high": p.baseline_rate_ci_high,
                "ldl_c_mg_dl": p.ldl_c_mg_dl,
                "apob_mg_dl": p.apob_mg_dl,
                "hscrp_mg_l": p.hscrp_mg_l,
                "systolic_bp_mmhg": p.systolic_bp_mmhg,
                "lp_a_nmol_l": p.lp_a_nmol_l,
                "has_ascvd": p.has_ascvd,
                "has_diabetes": p.has_diabetes,
                "has_ckd": p.has_ckd,
                "has_heart_failure": p.has_heart_failure,
                "is_secondary_prevention": p.is_secondary_prevention,
                "notes": p.notes,
            }
            for p in populations.values()
        ]
    )
    duplicated = pd.concat([frame, frame.iloc[[0]]], ignore_index=True)
    path = tmp_path / "dup.csv"
    duplicated.to_csv(path, index=False)
    with pytest.raises(ValueError, match="duplicate population_id"):
        load_populations(path)
