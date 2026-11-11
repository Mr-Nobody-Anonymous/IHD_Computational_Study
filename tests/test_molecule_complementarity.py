"""Tests for the evidence-aware combination test (not a chemically fused molecule)."""

import pandas as pd
import pytest

from ihd_solutions.molecule_complementarity import (
    DISCLAIMER,
    REGIMENS,
    STATUS,
    complementarity_report,
    complementarity_table,
    write_complementarity_report,
)


@pytest.fixture(scope="module")
def scenarios():
    return complementarity_table()


def test_all_scenarios_are_valid_untrialled_combinations(scenarios):
    assert len(REGIMENS) == 4
    assert len(scenarios) == 8  # four regimens, strict and surrogate credit modes
    assert (scenarios["status"] == STATUS).all()
    assert not scenarios["combination_randomized_proof"].any()


def test_colchicine_adds_complementary_strict_mode_credit(scenarios):
    strict = scenarios.loc[scenarios["credit_mode"] == "strict"].set_index("scenario")
    assert strict.loc["ldl_plus_colchicine", "increment_vs_ldl_backbone_per_1000"] > 0.0
    assert "transported_population" in strict.loc["ldl_plus_colchicine", "evidence_warnings"]


def test_enlicitide_is_not_given_unearned_hard_outcome_credit(scenarios):
    strict = scenarios.loc[scenarios["credit_mode"] == "strict"].set_index("scenario")
    surrogate = scenarios.loc[scenarios["credit_mode"] == "surrogate"].set_index("scenario")
    assert strict.loc["ldl_plus_enlicitide", "increment_vs_ldl_backbone_per_1000"] == pytest.approx(
        0.0, abs=1e-10
    )
    assert "no_anchor" in strict.loc["ldl_plus_enlicitide", "evidence_warnings"]
    assert surrogate.loc["ldl_plus_enlicitide", "increment_vs_ldl_backbone_per_1000"] > 0.0


def test_report_refuses_to_call_the_regimen_a_new_molecule(scenarios):
    text = complementarity_report(scenarios)
    assert "not a new chemically fused molecule" in text.lower()
    assert "combination-specific efficacy remain unknown" in text
    assert DISCLAIMER in text


def test_cli_exposes_the_complementarity_test():
    from ihd_solutions.cli import build_parser

    args = build_parser().parse_args(["molecule-complementarity"])
    assert args.command == "molecule-complementarity"
    assert str(args.output_dir) == "results/molecule_complementarity"


def test_cli_runs_and_writes_the_complementarity_outputs(tmp_path):
    from ihd_solutions.cli import main

    assert main(["molecule-complementarity", "--output-dir", str(tmp_path)]) == 0
    assert (tmp_path / "complementarity_scenarios.csv").exists()
    assert (tmp_path / "complementarity_report.txt").exists()
    assert (tmp_path / "DISCLAIMER.txt").exists()


def test_report_writer_emits_reproducible_outputs(tmp_path):
    frame = write_complementarity_report(tmp_path)
    assert (tmp_path / "complementarity_scenarios.csv").exists()
    assert (tmp_path / "complementarity_report.txt").exists()
    assert (tmp_path / "DISCLAIMER.txt").exists()
    saved = pd.read_csv(tmp_path / "complementarity_scenarios.csv", keep_default_na=False)
    pd.testing.assert_frame_equal(saved, frame, check_dtype=False)
