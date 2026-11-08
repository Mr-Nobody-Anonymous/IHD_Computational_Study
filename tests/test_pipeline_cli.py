"""End-to-end pipeline, manifest integrity and CLI behaviour."""

from __future__ import annotations

import hashlib
import json

import pandas as pd
import pytest

from ihd_solutions.cli import build_parser, main
from ihd_solutions.evidence import CreditMode
from ihd_solutions.pipeline import (
    EXAMPLE_GUIDE_PROVENANCE,
    PipelineConfig,
    endpoint_specific_effects,
    run,
    write_json,
    write_table,
    write_text,
)


@pytest.fixture(scope="module")
def pipeline_run(tmp_path_factory):
    output_dir = tmp_path_factory.mktemp("solutions")
    config = PipelineConfig(
        output_dir=output_dir,
        horizon_years=3.0,
        step_days=60.0,
        max_optional=1,
        populations=("stable_ccs", "post_mi_pci"),
        monte_carlo_draws=25,
        monte_carlo_top_n=3,
        seed=99,
        plaque_months=36,
    )
    summary = run(config)
    return config, summary


@pytest.mark.slow
def test_pipeline_reports_success(pipeline_run):
    config, summary = pipeline_run
    assert summary["status"] == "passed"
    assert summary["regimens"] > 10
    assert summary["populations"] == 2
    assert summary["outputs"] >= 20


@pytest.mark.slow
def test_manifest_hashes_match_every_written_file(pipeline_run):
    config, _ = pipeline_run
    manifest = json.loads((config.output_dir / "run_manifest.json").read_text())
    assert manifest["status"] == "passed"
    assert "not clinical prediction" in manifest["scope"]
    for record in manifest["outputs"]:
        path = config.output_dir / record["path"]
        assert path.is_file(), path
        payload = path.read_bytes()
        assert len(payload) == record["bytes"]
        assert hashlib.sha256(payload).hexdigest() == record["sha256"]


@pytest.mark.slow
def test_expected_tables_exist_and_are_non_empty(pipeline_run):
    config, _ = pipeline_run
    tables = config.output_dir / "tables"
    expected = [
        "therapy_catalog_resolved.csv",
        "evidence_tier_summary.csv",
        "endpoint_specific_effects.csv",
        "regimen_ranking.csv",
        "credit_decisions.csv",
        "pareto_front.csv",
        "lipid_trajectories.csv",
        "risk_curves.csv",
        "monte_carlo_summary.csv",
        "rank_stability.csv",
        "sensitivity_tornado.csv",
        "plaque_solution_grid.csv",
        "inflammation_solution_grid.csv",
        "editing_scenarios.csv",
        "regeneration_scenarios.csv",
    ]
    for name in expected:
        frame = pd.read_csv(tables / name)
        assert not frame.empty, name


@pytest.mark.slow
def test_figures_are_valid_pngs(pipeline_run):
    config, _ = pipeline_run
    figures = list((config.output_dir / "figures").glob("*.png"))
    assert len(figures) == 4
    for path in figures:
        assert path.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        assert path.stat().st_size > 5000


@pytest.mark.slow
def test_report_states_the_validation_boundary(pipeline_run):
    config, _ = pipeline_run
    text = (config.output_dir / "reports" / "solutions_simulation_report.txt").read_text()
    for phrase in (
        "not a patient-level prediction",
        "hsCRP and Lp(a) changes never create event credit",
        "No base-editor safety conclusion",
        "No remuscularisation",
    ):
        assert phrase in text


@pytest.mark.slow
def test_credit_decisions_table_documents_refusals(pipeline_run):
    config, _ = pipeline_run
    frame = pd.read_csv(config.output_dir / "tables" / "credit_decisions.csv")
    denied = frame.loc[~frame["credited"]]
    assert not denied.empty
    assert denied["reason"].str.len().gt(10).all()
    credited = frame.loc[frame["credited"]]
    assert credited["anchor_id"].str.len().gt(0).all()


@pytest.mark.slow
def test_bystander_and_off_target_outputs_refuse_safety_claims(pipeline_run):
    config, _ = pipeline_run
    bystander = json.loads((config.output_dir / "tables" / "bystander_report.json").read_text())
    assert bystander["provenance"] == EXAMPLE_GUIDE_PROVENANCE
    assert "not the VERVE-102 guide" in bystander["warning"]
    off_target = json.loads((config.output_dir / "tables" / "off_target_status.json").read_text())
    assert off_target["numeric_off_target_score"] is None


@pytest.mark.slow
def test_editing_and_regeneration_tables_carry_no_event_credit(pipeline_run):
    config, _ = pipeline_run
    editing_frame = pd.read_csv(config.output_dir / "tables" / "editing_scenarios.csv")
    assert (editing_frame["mace_credit"] == 0.0).all()
    assert not editing_frame["clinical_safety_assessed"].any()
    regeneration = pd.read_csv(config.output_dir / "tables" / "regeneration_scenarios.csv")
    assert (regeneration["mace_credit"] == 0.0).all()


@pytest.mark.slow
def test_endpoint_specific_table_separates_non_mace_benefits(anchors):
    frame = endpoint_specific_effects(anchors)
    assert not frame.empty
    assert not frame["credited_as_mace"].any()
    assert "HF_COMPOSITE" in set(frame["endpoint_class"])
    assert "FUNCTIONAL" in set(frame["endpoint_class"])


def test_pipeline_rejects_unknown_populations_and_bad_config(tmp_path):
    with pytest.raises(ValueError, match="unknown populations"):
        run(PipelineConfig(output_dir=tmp_path, populations=("nope",), make_figures=False))
    with pytest.raises(ValueError):
        PipelineConfig(horizon_years=0.0)
    with pytest.raises(ValueError):
        PipelineConfig(max_optional=-1)
    with pytest.raises(ValueError):
        PipelineConfig(populations=())
    with pytest.raises(ValueError):
        PipelineConfig(modes=())
    with pytest.raises(ValueError):
        PipelineConfig(monte_carlo_draws=0)


def test_writers_are_atomic_and_reject_non_finite_json(tmp_path):
    frame = pd.DataFrame({"a": [1, 2]})
    path = write_table(tmp_path / "t.csv", frame)
    assert path.is_file()
    assert not path.with_suffix(".csv.tmp").exists()
    text_path = write_text(tmp_path / "t.txt", "hello\n")
    assert text_path.read_text() == "hello\n"
    json_path = write_json(tmp_path / "t.json", {"a": 1})
    assert json.loads(json_path.read_text()) == {"a": 1}
    with pytest.raises(ValueError):
        write_json(tmp_path / "bad.json", {"a": float("nan")})


def test_cli_validate_succeeds(capsys):
    assert main(["validate"]) == 0
    assert "Validated" in capsys.readouterr().out


@pytest.mark.parametrize("command", ["therapies", "anchors", "populations"])
def test_cli_listing_commands(command, capsys):
    assert main([command]) == 0
    output = capsys.readouterr().out
    assert len(output.splitlines()) > 5


def test_cli_regimens_reports_the_combinatorial_check(capsys):
    assert main(["regimens", "--max-optional", "1"]) == 0
    assert "combinatorial check" in capsys.readouterr().out


def test_cli_evaluate_text_and_json(capsys):
    assert main(["evaluate", "--therapies", "statin_high_intensity", "evolocumab"]) == 0
    text = capsys.readouterr().out
    assert "Credit decisions" in text
    assert "CREDITED" in text
    assert main(["evaluate", "--therapies", "enlicitide_oral", "--json"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["mode"] == "strict"


def test_cli_evaluate_shows_mode_dependent_credit(capsys):
    args = ["evaluate", "--therapies", "statin_high_intensity", "enlicitide_oral", "--json"]
    assert main(args + ["--mode", "strict"]) == 0
    strict = json.loads(capsys.readouterr().out)
    assert main(args + ["--mode", "surrogate"]) == 0
    surrogate = json.loads(capsys.readouterr().out)
    assert surrogate["events_avoided_per_1000"] > strict["events_avoided_per_1000"]


def test_cli_error_paths_return_non_zero(capsys):
    assert main(["evaluate", "--therapies", "not_a_drug"]) == 1
    assert main(["evaluate", "--therapies", "evolocumab", "enlicitide_oral"]) == 1
    assert main(["evaluate", "--therapies", "evolocumab", "--population", "nope"]) == 1
    assert "Error:" in capsys.readouterr().err


def test_cli_simulate_writes_manifest(tmp_path, capsys):
    assert (
        main(
            [
                "simulate",
                "--output-dir",
                str(tmp_path),
                "--max-optional",
                "0",
                "--populations",
                "stable_ccs",
                "--modes",
                "strict",
                "--draws",
                "10",
                "--no-figures",
            ]
        )
        == 0
    )
    payload = json.loads(capsys.readouterr().out)
    assert payload["status"] == "passed"
    assert (tmp_path / "run_manifest.json").is_file()


def test_parser_requires_a_command():
    parser = build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([])


def test_all_credit_modes_are_accepted_by_the_parser():
    parser = build_parser()
    for mode in CreditMode:
        args = parser.parse_args(["evaluate", "--therapies", "evolocumab", "--mode", mode.value])
        assert args.mode == mode.value
