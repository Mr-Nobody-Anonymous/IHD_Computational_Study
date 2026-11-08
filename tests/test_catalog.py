"""Catalogue loading, validation and combination rules."""

from __future__ import annotations

from dataclasses import replace

import pytest

from ihd_solutions.catalog import (
    DosingProfile,
    Regimen,
    TherapyClass,
    catalog_frame,
    load_catalog,
    therapies_in_class,
    validate_combination,
)


def test_catalog_covers_every_solution_class(catalog):
    represented = {therapy.therapy_class for therapy in catalog.values()}
    assert represented == set(TherapyClass)


def test_catalog_has_the_expected_size_and_unique_ids(catalog):
    assert len(catalog) == 47
    assert len({therapy.therapy_id for therapy in catalog.values()}) == len(catalog)


def test_every_therapy_declares_a_source_url(catalog):
    for therapy in catalog.values():
        assert therapy.source_url.startswith("https://"), therapy.therapy_id


def test_pcsk9_modalities_share_one_exclusive_group(catalog):
    groups = {
        catalog[therapy_id].exclusive_group
        for therapy_id in (
            "evolocumab",
            "alirocumab",
            "inclisiran",
            "enlicitide_oral",
            "azd0780_oral",
        )
    }
    assert groups == {"pcsk9_protein"}


def test_base_editing_is_its_own_group_and_permanent_scenario(catalog):
    therapy = catalog["verve102_base_editing"]
    assert therapy.exclusive_group == "pcsk9_gene"
    assert therapy.profile is DosingProfile.ONE_TIME_EDIT
    assert therapy.is_permanent_scenario


def test_statin_intensity_variants_are_mutually_exclusive(catalog):
    with pytest.raises(ValueError, match="exclusive group"):
        validate_combination(catalog, ("statin_high_intensity", "statin_moderate_intensity"))


def test_two_pcsk9_agents_cannot_be_stacked(catalog):
    with pytest.raises(ValueError, match="exclusive group"):
        validate_combination(catalog, ("evolocumab", "enlicitide_oral"))


def test_unknown_therapy_is_rejected(catalog):
    with pytest.raises(ValueError, match="unknown therapy ids"):
        validate_combination(catalog, ("not_a_therapy",))


def test_regimen_sorts_ids_and_builds_a_label(catalog):
    regimen = validate_combination(catalog, ("evolocumab", "ezetimibe"))
    assert regimen.therapy_ids == ("evolocumab", "ezetimibe")
    assert regimen.label == "evolocumab + ezetimibe"
    assert len(regimen) == 2


def test_regimen_rejects_duplicates_and_empty():
    with pytest.raises(ValueError, match="duplicate therapies"):
        Regimen(("a", "a"))
    with pytest.raises(ValueError, match="at least one therapy"):
        Regimen(())


def test_therapies_in_class_is_sorted(catalog):
    lipid = therapies_in_class(catalog, TherapyClass.LIPID_LOWERING)
    assert [therapy.therapy_id for therapy in lipid] == sorted(t.therapy_id for t in lipid)
    assert len(lipid) >= 8


def test_catalog_frame_has_one_row_per_therapy(catalog):
    frame = catalog_frame(catalog)
    assert len(frame) == len(catalog)
    assert set(frame["therapy_id"]) == set(catalog)


@pytest.mark.parametrize(
    "changes",
    [
        {"ldl_pct": 1.0},
        {"ldl_pct": -0.1},
        {"hscrp_pct": 1.2},
        {"mmp_factor": 0.0},
        {"mmp_factor": 1.5},
        {"plateau_days": 0.0},
        {"onset_days": -1.0},
        {"burden_score": 0.0},
        {"trough_fraction": 0.0},
        {"ldl_pct_ci_low": 0.9},
    ],
)
def test_therapy_validation_rejects_implausible_parameters(catalog, changes):
    with pytest.raises(ValueError):
        replace(catalog["evolocumab"], **changes)


def test_interval_dosing_requires_an_interval(catalog):
    with pytest.raises(ValueError, match="interval dosing requires"):
        replace(catalog["evolocumab"], interval_days=0.0)


def test_plateau_must_follow_onset(catalog):
    with pytest.raises(ValueError, match="plateau_days must be at least onset_days"):
        replace(catalog["evolocumab"], onset_days=100.0)


def test_therapy_requires_identifiers(catalog):
    with pytest.raises(ValueError):
        replace(catalog["evolocumab"], therapy_id="")


def test_procedural_and_fasting_flags_are_read_correctly(catalog):
    assert catalog["pci_stable_cad"].is_procedural
    assert catalog["pci_stable_cad"].symptom_benefit_only
    assert catalog["enlicitide_oral"].requires_fasting
    assert not catalog["azd0780_oral"].requires_fasting


def test_load_catalog_rejects_missing_and_malformed_files(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_catalog(tmp_path / "missing.csv")
    path = tmp_path / "bad.csv"
    path.write_text("therapy_id,name\nx,y\n", encoding="utf-8")
    with pytest.raises(ValueError, match="missing columns"):
        load_catalog(path)


def test_boolean_parsing_rejects_nonsense(catalog):
    from ihd_solutions.catalog import _as_bool

    assert "requires_fasting" in catalog_frame(catalog)
    assert _as_bool(True) is True
    assert _as_bool("yes") is True
    assert _as_bool("0") is False
    with pytest.raises(ValueError):
        _as_bool("maybe")
