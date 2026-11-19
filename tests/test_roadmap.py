from __future__ import annotations

from dataclasses import replace

import pandas as pd
import pytest

from ihd_solutions.catalog import load_catalog
from ihd_solutions.evidence import load_anchors
from ihd_solutions.population import load_populations
from ihd_solutions.roadmap import (
    DATA_INVENTORY_COLUMNS,
    ESTIMAND_COLUMNS,
    audit_anchor_links,
    build_endgame_report,
    build_evidence_review_queue,
    load_estimands,
    load_validation_inventory,
    roadmap_input_audit,
    write_endgame_report,
    write_evidence_review_queue,
)


def test_current_estimands_registry_loads_and_marks_unmodelled_acute_mi():
    populations = set(load_populations())
    frame = load_estimands(population_ids=populations)
    assert len(frame) == 9
    acute = frame.loc[frame["estimand_id"] == "acute_mi_reperfusion"].iloc[0]
    assert acute["model_coverage"] == "not_modelled"
    assert acute["population_id"] == "not_registered"
    assert acute["endpoint_class"] == "NOT_DEFINED"


def test_functional_lvef_estimand_is_not_mislabeled_as_hard_outcome():
    frame = load_estimands()
    row = frame.loc[frame["estimand_id"] == "post_mi_remodeling"].iloc[0]
    assert row["endpoint_class"] == "FUNCTIONAL"
    assert row["effect_measure"] == "mean_difference"
    assert row["model_coverage"] == "not_modelled"


def test_external_data_inventory_records_candidates_not_acquired_data():
    frame = load_validation_inventory()
    assert len(frame) == 4
    assert set(frame["dataset_id"]) == {"DREAM_HF", "REVIVED_BCIS2", "BioVAT_HF", "AATS_LB8"}
    assert (frame["access_status"].str.contains("not verified|unverified|must be reviewed")).any()


def test_roadmap_audit_checks_registry_links_but_does_not_claim_validation():
    result = roadmap_input_audit(load_anchors(), set(load_catalog()), set(load_populations()))
    assert result["estimands_registered"] == 9
    assert result["external_dataset_leads"] == 4
    assert result["external_datasets_accessed"] == 0
    assert result["external_validation_completed"] is False
    assert result["source_content_human_verified"] is False


def _base_estimand() -> dict[str, str]:
    return {
        "estimand_id": "test",
        "disease_stage": "stable",
        "population_id": "stable_ccs",
        "comparator": "standard care",
        "intervention_scope": "test intervention",
        "primary_endpoint": "time to MACE",
        "endpoint_class": "MACE",
        "effect_measure": "hazard_ratio",
        "horizon_years": "5",
        "model_coverage": "scenario_supported",
        "status": "planning_only_not_trial_protocol",
        "rationale": "synthetic test row",
    }


def test_estimand_loader_rejects_duplicate_ids(tmp_path):
    frame = pd.DataFrame([_base_estimand(), _base_estimand()])
    path = tmp_path / "estimands.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="duplicate estimand_id"):
        load_estimands(path, population_ids={"stable_ccs"})


def test_estimand_loader_rejects_unregistered_population_without_blocker(tmp_path):
    row = _base_estimand()
    row["population_id"] = "missing_population"
    path = tmp_path / "estimands.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="unknown population_id"):
        load_estimands(path, population_ids={"stable_ccs"})


def test_estimand_loader_only_allows_unfrozen_endpoint_when_explicitly_blocked(tmp_path):
    row = _base_estimand()
    row.update(
        population_id="not_registered",
        endpoint_class="NOT_DEFINED",
        effect_measure="not_set",
        horizon_years="",
        model_coverage="not_modelled",
    )
    path = tmp_path / "estimands.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    assert (
        load_estimands(path, population_ids={"stable_ccs"}).iloc[0]["population_id"]
        == "not_registered"
    )
    row["horizon_years"] = "1"
    pd.DataFrame([row]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="must not invent horizon"):
        load_estimands(path, population_ids={"stable_ccs"})


def test_estimand_loader_rejects_functional_endpoint_with_hazard_ratio(tmp_path):
    row = _base_estimand()
    row.update(endpoint_class="FUNCTIONAL", effect_measure="hazard_ratio")
    path = tmp_path / "estimands.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="mean_difference"):
        load_estimands(path, population_ids={"stable_ccs"})


def test_estimand_loader_rejects_missing_columns(tmp_path):
    path = tmp_path / "estimands.csv"
    pd.DataFrame([{column: "x" for column in ESTIMAND_COLUMNS[:-1]}]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="estimand registry missing columns"):
        load_estimands(path)


def test_inventory_loader_rejects_invalid_source_url(tmp_path):
    frame = load_validation_inventory()
    frame.loc[0, "source_url"] = "not-a-url"
    path = tmp_path / "inventory.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match=r"HTTP\(S\) URL"):
        load_validation_inventory(path)


def test_inventory_loader_rejects_duplicate_dataset_id(tmp_path):
    frame = load_validation_inventory()
    frame.loc[1, "dataset_id"] = frame.loc[0, "dataset_id"]
    path = tmp_path / "inventory.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="duplicate dataset_id"):
        load_validation_inventory(path)


def test_inventory_loader_rejects_missing_columns(tmp_path):
    path = tmp_path / "inventory.csv"
    pd.DataFrame([{column: "x" for column in DATA_INVENTORY_COLUMNS[:-1]}]).to_csv(
        path, index=False
    )
    with pytest.raises(ValueError, match="validation data inventory missing columns"):
        load_validation_inventory(path)


def test_anchor_audit_rejects_bad_source_and_bad_transport_link():
    anchors = load_anchors()
    first_id = next(iter(anchors))
    anchors[first_id] = replace(anchors[first_id], source_url="relative/path")
    with pytest.raises(ValueError, match="invalid_sources"):
        audit_anchor_links(anchors, set(load_catalog()), set(load_populations()))


def test_anchor_audit_rejects_unknown_transport_target():
    anchors = load_anchors()
    first_id = next(iter(anchors))
    anchors[first_id] = replace(anchors[first_id], transportable_to=("not_a_population",))
    with pytest.raises(ValueError, match="unknown_transport_targets"):
        audit_anchor_links(anchors, set(load_catalog()), set(load_populations()))


def test_estimand_loader_rejects_missing_empty_and_blank_identifier(tmp_path):
    with pytest.raises(FileNotFoundError, match="estimand registry not found"):
        load_estimands(tmp_path / "missing.csv")
    empty_path = tmp_path / "empty.csv"
    pd.DataFrame(columns=ESTIMAND_COLUMNS).to_csv(empty_path, index=False)
    with pytest.raises(ValueError, match="estimand registry is empty"):
        load_estimands(empty_path)
    blank_path = tmp_path / "blank.csv"
    row = _base_estimand()
    row["estimand_id"] = " "
    pd.DataFrame([row]).to_csv(blank_path, index=False)
    with pytest.raises(ValueError, match="estimand_id cannot be blank"):
        load_estimands(blank_path, population_ids={"stable_ccs"})


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("comparator", "", "comparator cannot be blank"),
        ("endpoint_class", "made_up", "unsupported endpoint_class"),
        ("effect_measure", "made_up", "unsupported effect_measure"),
        ("model_coverage", "made_up", "unsupported model_coverage"),
        ("horizon_years", "", "horizon_years must be positive"),
        ("horizon_years", "0", "horizon_years must be positive"),
        ("effect_measure", "mean_difference", "hard outcomes require"),
    ],
)
def test_estimand_loader_rejects_invalid_schema_values(tmp_path, field, value, message):
    row = _base_estimand()
    row[field] = value
    if field == "effect_measure" and value == "mean_difference":
        row["endpoint_class"] = "MACE"
    path = tmp_path / "invalid.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    with pytest.raises(ValueError, match=message):
        load_estimands(path, population_ids={"stable_ccs"})


def test_estimand_loader_rejects_undefined_endpoint_without_blocked_status(tmp_path):
    row = _base_estimand()
    row["endpoint_class"] = "NOT_DEFINED"
    path = tmp_path / "invalid.csv"
    pd.DataFrame([row]).to_csv(path, index=False)
    with pytest.raises(ValueError, match="NOT_DEFINED is permitted only"):
        load_estimands(path, population_ids={"stable_ccs"})


def test_validation_inventory_rejects_missing_empty_blank_and_bad_sample_size(tmp_path):
    with pytest.raises(FileNotFoundError, match="validation data inventory not found"):
        load_validation_inventory(tmp_path / "missing.csv")
    empty_path = tmp_path / "empty.csv"
    pd.DataFrame(columns=DATA_INVENTORY_COLUMNS).to_csv(empty_path, index=False)
    with pytest.raises(ValueError, match="validation data inventory is empty"):
        load_validation_inventory(empty_path)

    valid = load_validation_inventory()
    blank_path = tmp_path / "blank.csv"
    valid.loc[0, "dataset_id"] = " "
    valid.to_csv(blank_path, index=False)
    with pytest.raises(ValueError, match="dataset_id cannot be blank"):
        load_validation_inventory(blank_path)

    invalid_path = tmp_path / "bad_n.csv"
    valid = load_validation_inventory()
    valid.loc[0, "reported_sample_size"] = "not-a-number"
    valid.to_csv(invalid_path, index=False)
    with pytest.raises(ValueError, match="reported_sample_size must be an integer"):
        load_validation_inventory(invalid_path)


def test_validation_inventory_rejects_invalid_ipd_value_and_nonpositive_n(tmp_path):
    frame = load_validation_inventory()
    frame.loc[0, "ipd_plan_declared"] = "maybe"
    path = tmp_path / "bad_ipd.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="invalid ipd_plan_declared"):
        load_validation_inventory(path)

    frame = load_validation_inventory()
    frame.loc[0, "reported_sample_size"] = "0"
    path = tmp_path / "bad_n.csv"
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="reported_sample_size must be positive"):
        load_validation_inventory(path)


def test_anchor_audit_rejects_unknown_catalog_and_primary_population_links():
    anchors = load_anchors()
    with pytest.raises(ValueError, match="unknown_therapies"):
        audit_anchor_links(anchors, set(), set(load_populations()))

    with pytest.raises(ValueError, match="unknown_populations"):
        audit_anchor_links(anchors, set(load_catalog()), set())


def test_endgame_report_states_blockers_and_does_not_claim_evidence(tmp_path):
    report = build_endgame_report(load_anchors(), set(load_catalog()), set(load_populations()))
    assert "datasets accessed: 0" in report
    assert "held-out external validation completed: false" in report
    assert "source contents human-verified: false" in report
    assert "acute_mi_reperfusion" in report
    assert "565 actual enrollment" in report
    assert "version 39" in report and "shows 566" in report
    assert "not acquired or approved" in report
    assert "No cure or 80% target is established" in report

    output = tmp_path / "reports" / "endgame.md"
    assert (
        write_endgame_report(output, load_anchors(), set(load_catalog()), set(load_populations()))
        == output
    )
    assert output.read_text(encoding="utf-8") == report


def test_evidence_review_queue_keeps_all_sources_pending(tmp_path):
    anchors = load_anchors()
    queue = build_evidence_review_queue(anchors)
    assert len(queue) == len(anchors) == 48
    assert set(queue["source_review_status"]) == {"pending_human_source_review"}
    assert queue["reviewer"].eq("").all()
    output = tmp_path / "human_review.csv"
    assert write_evidence_review_queue(output, anchors) == output
    loaded = pd.read_csv(output, keep_default_na=False)
    assert len(loaded) == 48
    assert loaded["source_review_status"].eq("pending_human_source_review").all()
