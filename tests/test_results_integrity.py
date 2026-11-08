"""Integrity checks for the committed pipeline outputs.

These tests bind `results/solutions/` to the code that produced it: every file
must match the SHA-256 recorded in the run manifest, and the headline tables must
respect the same evidence rules as the unit tests. If this file fails after a
model change, regenerate the outputs with `python run_solutions.py simulate`.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results" / "solutions"
MANIFEST = RESULTS / "run_manifest.json"

pytestmark = pytest.mark.skipif(
    not MANIFEST.exists(),
    reason="run outputs not generated; execute `python run_solutions.py simulate`",
)


@pytest.fixture(scope="module")
def manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding="utf-8"))


def test_manifest_declares_a_successful_scoped_run(manifest):
    assert manifest["status"] == "passed"
    assert "not clinical prediction" in manifest["scope"].lower()
    assert manifest["config"]["horizon_years"] > 0
    assert set(manifest["config"]["modes"]) == {"strict", "surrogate", "exploratory"}


def test_every_output_matches_its_recorded_hash_and_size(manifest):
    outputs = manifest["outputs"]
    assert len(outputs) >= 20
    for entry in outputs:
        path = RESULTS / entry["path"]
        assert path.is_file(), f"missing output {entry['path']}"
        payload = path.read_bytes()
        assert len(payload) == entry["bytes"], entry["path"]
        assert hashlib.sha256(payload).hexdigest() == entry["sha256"], entry["path"]


def test_no_unlisted_files_are_left_in_the_output_directory(manifest):
    listed = {entry["path"] for entry in manifest["outputs"]} | {"run_manifest.json"}
    present = {str(path.relative_to(RESULTS)) for path in RESULTS.rglob("*") if path.is_file()}
    assert present - listed == set()


def test_ranking_table_respects_the_credit_rules(manifest):
    frame = pd.read_csv(RESULTS / "tables" / "regimen_ranking.csv")
    assert not frame.empty
    assert frame["regimen"].is_unique or not frame.duplicated(["regimen", "population_id"]).any()
    # Surrogate-only mechanisms cannot outrank themselves into strict credit.
    surrogate_only = frame[frame["regimen"].str.contains("enlicitide|azd0780|verve102|inclisiran")]
    assert (
        surrogate_only["surrogate_events_avoided_per_1000"]
        >= surrogate_only["strict_events_avoided_per_1000"]
    ).all()
    assert (frame["strict_events_avoided_per_1000"] >= -1e-9).all()


def test_neutral_and_unproven_mechanisms_are_recorded_as_uncredited_or_neutral():
    frame = pd.read_csv(RESULTS / "tables" / "credit_decisions.csv")
    strict = frame[frame["mode"] == "strict"]
    for therapy in ("enlicitide_oral", "azd0780_oral", "verve102_base_editing", "pelacarsen"):
        rows = strict[strict["therapy_id"] == therapy]
        assert not rows.empty, therapy
        assert not rows["credited"].any(), therapy
    # CIRT is credited as randomized evidence but must stay numerically neutral.
    cirt = strict[strict["anchor_id"] == "CIRT"]
    assert not cirt.empty
    assert cirt["hazard_ratio"].between(0.99, 1.03).all()


def test_endpoint_specific_effects_keep_heart_failure_endpoints_out_of_mace():
    frame = pd.read_csv(RESULTS / "tables" / "endpoint_specific_effects.csv")
    hf = frame[frame["endpoint_class"].str.contains("hf_composite|functional", case=False)]
    assert not hf.empty
    assert not hf["credited_as_mace"].any()


def test_figures_are_real_png_files():
    figures = sorted((RESULTS / "figures").glob("*.png"))
    assert len(figures) >= 3
    for figure in figures:
        assert figure.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert figure.stat().st_size > 10_000


def test_monte_carlo_intervals_bracket_their_medians():
    frame = pd.read_csv(RESULTS / "tables" / "monte_carlo_summary.csv")
    assert not frame.empty
    assert (frame["events_avoided_p5"] <= frame["events_avoided_median"]).all()
    assert (frame["events_avoided_median"] <= frame["events_avoided_p95"]).all()
    assert frame["probability_of_benefit"].between(0, 1).all()


def test_cure_search_outputs_refuse_a_cure_claim():
    scores = pd.read_csv(RESULTS / "tables" / "cure_criteria_scores.csv")
    assert not scores.empty
    assert scores["cure_readiness_index"].max() < 0.8
    assert scores["criteria_blocking"].str.len().gt(0).all()

    search = json.loads((RESULTS / "tables" / "cure_combination_search.json").read_text())
    assert search["satisfied"] is False
    assert search["selected"] == []
    assert set(search["binding_constraints"]) >= {"burden_reversal", "myocardial_restoration"}
    report = (RESULTS / "reports" / "cure_search_report.txt").read_text()
    assert "VERDICT: no catalogued intervention" in report
    assert "No cure is demonstrated" in report


def test_regression_outputs_never_claim_eradication_or_events():
    grid = pd.read_csv(RESULTS / "tables" / "plaque_regression_grid.csv")
    assert not grid.empty
    assert not grid["eradicated"].any()
    assert (grid["mace_credit_from_imaging"] == 0).all()
    assert grid["calcified_volume_increased"].any()
    imaging = pd.read_csv(RESULTS / "tables" / "plaque_imaging_evidence.csv")
    assert (imaging["mace_credit"] == 0).all()
    fit = json.loads((RESULTS / "tables" / "regression_calibration_fit.json").read_text())
    assert 60.0 < fit["regression_threshold_ldl_mg_dl"] < 90.0
    assert fit["caveats"]


def test_lifetime_scenarios_are_labelled_as_genetic_scenarios():
    frame = pd.read_csv(RESULTS / "tables" / "lifetime_exposure_scenarios.csv")
    assert (frame["evidence_class"] == "mendelian_randomisation_scenario").all()
    for _, group in frame.groupby("modality"):
        ordered = group.sort_values("start_age")
        assert ordered["scenario_lifetime_risk_pct"].is_monotonic_increasing


def test_residual_risk_shows_the_remaining_gap():
    frame = pd.read_csv(RESULTS / "tables" / "residual_risk.csv")
    assert not frame.empty
    assert (frame["excess_risk_remaining_fraction"] > 0).all()
    assert (frame["treated_risk_pct"] < frame["baseline_risk_pct"]).all()


# ---------------------------------------------------------------------------
# The invented REVERT-1 hypothesis must stay out of the evidence manifest.
# ---------------------------------------------------------------------------


def test_evidence_manifest_contains_no_hypothesis_artefacts(manifest):
    payload = json.dumps(manifest).lower()
    assert "revert" not in payload
    assert "hypothes" not in payload


def test_hypothesis_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*revert*"))
    assert not list(RESULTS.rglob("*hypothes*"))


def test_committed_hypothesis_tables_are_all_tagged():
    tables = ROOT / "results" / "hypotheses" / "tables"
    if not tables.exists():
        pytest.skip("hypothesis outputs not generated; run `python run_solutions.py hypothesis`")
    paths = sorted(tables.glob("*.csv"))
    assert len(paths) == 8
    for path in paths:
        frame = pd.read_csv(path)
        assert frame["status"].eq("invented_hypothesis_unvalidated").all(), path.name


def test_evidence_manifest_contains_no_discovery_artefacts(manifest):
    payload = json.dumps(manifest).lower()
    assert "discovery" not in payload
    assert "specification" not in payload


def test_discovery_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*discovery*"))
    assert not list(RESULTS.rglob("*specification*"))


def test_committed_discovery_tables_are_all_tagged():
    tables = ROOT / "results" / "discovery" / "tables"
    if not tables.exists():
        pytest.skip("discovery outputs not generated; run `python run_solutions.py discover`")
    paths = sorted(tables.glob("*.csv"))
    assert len(paths) == 8
    for path in paths:
        frame = pd.read_csv(path)
        assert frame["status"].eq("in_silico_specification_not_evidence").all(), path.name


def test_evidence_manifest_contains_no_composite_artefacts(manifest):
    payload = json.dumps(manifest).lower()
    assert "composite" not in payload


def test_composite_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*composite*"))


def test_committed_composite_tables_declare_no_cure():
    tables = ROOT / "results" / "composite" / "tables"
    if not tables.exists():
        pytest.skip("composite outputs not generated; run `python run_solutions.py composite`")
    criteria = pd.read_csv(tables / "composite_criteria.csv")
    assert not criteria["is_cure"].any()
    assert criteria["randomised_proof"].eq(0.0).all()
    assert criteria["criteria_met"].lt(5).all()
    allowed = {"untrialled_composite_scenario", "invented_component_unvalidated"}
    for path in tables.glob("*.csv"):
        frame = pd.read_csv(path)
        if "status" in frame.columns and not frame.empty:
            assert set(frame["status"]) <= allowed, path.name


def test_evidence_manifest_contains_no_gap_artefacts(manifest):
    assert "results/gaps" not in json.dumps(manifest)


def test_gap_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*gap_*"))


def test_committed_gap_tables_are_all_tagged():
    tables = ROOT / "results" / "gaps" / "tables"
    if not tables.exists():
        pytest.skip("gap outputs not generated; run `python run_solutions.py gaps`")
    for path in tables.glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].eq("gap_closure_scenario_not_evidence").all(), path.name


def test_evidence_manifest_contains_no_resolution_artefacts(manifest):
    assert "results/resolution" not in json.dumps(manifest)


def test_resolution_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*resolution*"))


def test_committed_resolution_tables_are_all_tagged():
    tables = ROOT / "results" / "resolution" / "tables"
    if not tables.exists():
        pytest.skip("resolution outputs not generated; run `python run_solutions.py resolve`")
    for path in tables.glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].eq("residual_resolution_scenario_not_evidence").all(), path.name


def test_evidence_manifest_contains_no_atlas_artefacts(manifest):
    assert "results/atlas" not in json.dumps(manifest)


def test_atlas_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*atlas*"))


def test_committed_atlas_tables_are_all_tagged():
    tables = ROOT / "results" / "atlas" / "tables"
    if not tables.exists():
        pytest.skip("atlas outputs not generated; run `python run_solutions.py atlas`")
    for path in tables.glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].str.contains("atlas").all(), path.name


def test_committed_atlas_scorecard_declares_no_cure_and_no_proof():
    path = ROOT / "results" / "atlas" / "tables" / "atlas_scorecard.csv"
    if not path.exists():
        pytest.skip("atlas outputs not generated")
    frame = pd.read_csv(path)
    assert not frame["is_cure"].any()
    assert frame["evidence_randomised_proof"].eq(0.0).all()
    assert frame["hypothesised_randomised_proof"].eq(0.0).all()


def test_evidence_manifest_contains_no_audit_artefacts(manifest):
    assert "results/audit" not in json.dumps(manifest)


def test_audit_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*audit*"))


def test_committed_audit_tables_are_all_tagged():
    tables = ROOT / "results" / "audit" / "tables"
    if not tables.exists():
        pytest.skip("audit outputs not generated; run `python run_solutions.py audit`")
    for path in tables.glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].str.contains("atlas").all(), path.name


def test_committed_audit_records_the_retraction():
    path = ROOT / "results" / "audit" / "tables" / "audit_operating_point.csv"
    if not path.exists():
        pytest.skip("audit outputs not generated")
    frame = pd.read_csv(path)
    assert frame["nominal_is_misleading"].all()
    assert not frame["robust"].any()
    assert frame["probability_of_target"].max() < 0.02


def test_committed_atlas_scorecard_no_longer_claims_the_retracted_pass():
    path = ROOT / "results" / "atlas" / "tables" / "atlas_scorecard.csv"
    if not path.exists():
        pytest.skip("atlas outputs not generated")
    frame = pd.read_csv(path)
    assert (frame["hypothesised_criteria_met"] == frame["evidence_criteria_met"]).all()
    assert frame["hypothesised_burden_reversal"].max() < 0.8


def test_evidence_manifest_contains_no_evidence_audit_artefacts(manifest):
    assert "results/evidence_audit" not in json.dumps(manifest)


def test_evidence_audit_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*evidence_audit*"))


def test_committed_evidence_audit_tables_are_all_tagged():
    tables = ROOT / "results" / "evidence_audit" / "tables"
    if not tables.exists():
        pytest.skip("evidence audit not generated; run `python run_solutions.py audit-evidence`")
    for path in tables.glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].eq("evidence_layer_uncertainty_audit").all(), path.name


def test_committed_evidence_audit_records_the_strict_retraction():
    path = ROOT / "results" / "evidence_audit" / "tables" / "evidence_audit_matrix.csv"
    if not path.exists():
        pytest.skip("evidence audit not generated")
    frame = pd.read_csv(path)
    strict = frame.loc[frame["mode"].eq("strict")]
    assert strict["nominal_criteria_met"].max() == 1
    assert strict["robust_criteria_met"].max() == 0
    assert not frame["is_cure"].any()
    assert (frame["robust_criteria_met"] <= frame["nominal_criteria_met"]).all()


def test_evidence_manifest_contains_no_upgrade_artefacts(manifest):
    assert "results/upgrade" not in json.dumps(manifest)


def test_upgrade_outputs_live_outside_the_evidence_results_tree():
    assert not list(RESULTS.rglob("*upgrade*"))


def test_committed_upgrade_tables_are_all_tagged():
    tables = ROOT / "results" / "upgrade" / "tables"
    if not tables.exists():
        pytest.skip("upgrade outputs not generated; run `python run_solutions.py upgrade`")
    for path in tables.glob("*.csv"):
        frame = pd.read_csv(path)
        assert frame["status"].eq("atlas2_invented_architecture_unvalidated").all(), path.name


def test_committed_upgrade_records_a_solve_with_margin_and_no_exclusion():
    path = ROOT / "results" / "upgrade" / "tables" / "upgrade_headline.csv"
    if not path.exists():
        pytest.skip("upgrade outputs not generated")
    frame = pd.read_csv(path).set_index("construct")
    solved = frame.loc["atlas2_solved"]
    assert solved["robust"]
    assert solved["robust_with_margin"]
    assert solved["eligible_fraction"] == 1.0
    assert not frame.loc["atlas1_baseline", "robust"]
