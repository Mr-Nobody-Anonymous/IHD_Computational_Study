"""Tests for the ATLAS-1 synthesis.

The load-bearing guard is the falsification sweep: at zero bed coupling the
sequencing idea must contribute exactly nothing. If that ever stops being true, the
"innovation" has become an assumption baked into the baseline.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ihd_solutions.atlas import (
    ALL_ARMS,
    ATLAS_INVENTED_TAG,
    ATLAS_POTENCY,
    DISCLAIMER,
    AtlasArms,
    arm_ablation,
    atlas_report,
    bed_coupling_sweep,
    bed_quality,
    clearance_for,
    evaluate_atlas_matrix,
    graft_outcome,
    run_atlas_study,
    score_atlas,
    translation_floor,
)
from ihd_solutions.cure import CRITERION_THRESHOLD
from ihd_solutions.evidence import CreditMode
from ihd_solutions.gap_solver import SCAR_RECOVERABLE_LVEF_PP


@pytest.fixture(scope="module")
def sweep() -> pd.DataFrame:
    return bed_coupling_sweep(draws=3000)


@pytest.fixture(scope="module")
def matrix_result():
    return evaluate_atlas_matrix(
        population_ids=("stable_ccs", "ischemic_cm_heart_failure"), draws=800
    )


# --------------------------------------------------------------------------
# The falsification guard.
# --------------------------------------------------------------------------


def test_zero_coupling_means_sequencing_contributes_nothing(sweep):
    row = sweep.loc[sweep["bed_coupling"].eq(0.0)].iloc[0]
    assert row["sequencing_gain_pp_probability"] == pytest.approx(0.0, abs=1e-9)
    assert row["translation_floor"] == pytest.approx(0.5)


def test_sequencing_gain_grows_with_coupling(sweep):
    gains = list(sweep["sequencing_gain_pp_probability"])
    assert gains == sorted(gains)
    assert gains[-1] > 20.0


def test_even_perfect_coupling_does_not_make_the_graft_robust(sweep):
    """The honest headline: sequencing helps but does not solve remuscularisation."""
    assert not sweep["robust_at_80pct"].any()
    assert sweep["probability_of_target"].max() < 0.8


def test_disabling_sequencing_restores_the_base_floor():
    assert translation_floor(ALL_ARMS.without("bed_sequencing"), 1.0) == pytest.approx(0.5)
    assert translation_floor(ALL_ARMS, 1.0) > 0.5


def test_coupling_rejects_out_of_range():
    with pytest.raises(ValueError):
        translation_floor(ALL_ARMS, 1.5)
    with pytest.raises(ValueError):
        graft_outcome(ALL_ARMS, bed_coupling=-0.1)


# --------------------------------------------------------------------------
# Inherited parts still behave.
# --------------------------------------------------------------------------


def test_the_operating_point_still_delivers_the_clearance_target():
    assert clearance_for(ALL_ARMS) >= CRITERION_THRESHOLD - 1e-3
    assert ATLAS_POTENCY > 1.0


def test_clearance_requires_the_gated_arm():
    assert clearance_for(ALL_ARMS.without("gated_clearance")) == 0.0


def test_source_reduction_is_load_bearing_for_clearance():
    without = clearance_for(ALL_ARMS.without("source_reduction"))
    assert without < clearance_for(ALL_ARMS)
    assert without < CRITERION_THRESHOLD


def test_bed_quality_is_zero_without_clearance():
    assert bed_quality(ALL_ARMS.without("gated_clearance")) == 0.0
    assert bed_quality(ALL_ARMS) == pytest.approx(1.0)


def test_bed_quality_is_bounded():
    assert 0.0 <= bed_quality(ALL_ARMS) <= 1.0


def test_no_remuscularisation_arm_means_no_gain():
    result = graft_outcome(ALL_ARMS.without("remuscularisation"))
    assert result.probability_of_target == 0.0
    assert result.median_gain_pp == 0.0
    assert not result.robust


def test_graft_gains_never_exceed_full_scar_replacement():
    result = graft_outcome(ALL_ARMS, bed_coupling=1.0, draws=2000)
    assert result.median_gain_pp <= SCAR_RECOVERABLE_LVEF_PP
    assert result.nominal_gain_pp <= SCAR_RECOVERABLE_LVEF_PP
    assert result.p05_gain_pp <= result.median_gain_pp


def test_graft_is_reproducible():
    first = graft_outcome(ALL_ARMS, draws=500, seed=7)
    second = graft_outcome(ALL_ARMS, draws=500, seed=7)
    assert first.probability_of_target == second.probability_of_target


def test_backbone_removal_costs_the_dose_escalation_factor():
    with_backbone = graft_outcome(ALL_ARMS, draws=3000)
    without = graft_outcome(ALL_ARMS.without("proven_backbone"), draws=3000)
    assert without.probability_of_target < with_backbone.probability_of_target


# --------------------------------------------------------------------------
# Ablation.
# --------------------------------------------------------------------------


def test_every_arm_is_load_bearing():
    ablation = arm_ablation(draws=2500)
    assert ablation["load_bearing"].all()
    assert len(ablation) == len(ALL_ARMS.names())


def test_editing_is_not_wrongly_labelled_inert():
    """It contributes durability, which clearance and graft mass cannot see."""
    ablation = arm_ablation(draws=2500).set_index("removed_arm")
    row = ablation.loc["durable_editing"]
    assert row["durable_driver_control_lost"]
    assert row["load_bearing"]
    assert row["clearance_loss_pp"] == pytest.approx(0.0, abs=1e-9)


def test_unknown_arm_is_rejected():
    with pytest.raises(ValueError, match="unknown arm"):
        ALL_ARMS.without("not_an_arm")


def test_arm_register_is_complete():
    assert set(ALL_ARMS.names()) == {
        "proven_backbone",
        "durable_editing",
        "gated_clearance",
        "source_reduction",
        "remuscularisation",
        "bed_sequencing",
    }
    assert all(getattr(ALL_ARMS, name) for name in ALL_ARMS.names())
    assert not AtlasArms(bed_sequencing=False).bed_sequencing


# --------------------------------------------------------------------------
# Two scorecards kept apart; no cure claim.
# --------------------------------------------------------------------------


def test_invented_arms_never_raise_the_evidence_scorecard(population, catalog, anchors, days):
    card = score_atlas(population, CreditMode.STRICT, catalog, anchors, days, draws=800)
    assert card.evidence_scores["randomised_proof"] == 0.0
    assert card.hypothesised_scores["randomised_proof"] == 0.0
    assert card.evidence_scores["myocardial_restoration"] <= 0.5
    assert card.hypothesised_criteria_met >= card.evidence_criteria_met


def test_proof_stays_zero_in_both_scorecards(matrix_result):
    _, scores = matrix_result
    assert scores["evidence_randomised_proof"].eq(0.0).all()
    assert scores["hypothesised_randomised_proof"].eq(0.0).all()


def test_no_configuration_is_declared_a_cure(matrix_result):
    _, scores = matrix_result
    assert not scores["is_cure"].any()
    assert scores["hypothesised_criteria_met"].lt(5).all()


def test_the_nominal_burden_pass_is_retracted_not_scored(matrix_result):
    """Pinned retraction: the nominal clearance pass must never count as a criterion.

    See ihd_solutions.audit. Under delivery uncertainty P(clearance target) is ~0,
    so the hypothesised burden score falls back to the evidence-only score.
    """
    _, scores = matrix_result
    assert (scores["hypothesised_burden_reversal"] >= scores["evidence_burden_reversal"]).all()
    assert scores["hypothesised_burden_reversal"].max() < CRITERION_THRESHOLD
    # The retracted number is still recorded, explicitly labelled, for audit trail.
    assert scores["retracted_nominal_burden_reversal"].max() >= CRITERION_THRESHOLD


def test_the_invented_arms_add_no_robust_criterion(matrix_result):
    _, scores = matrix_result
    assert (scores["hypothesised_criteria_met"] == scores["evidence_criteria_met"]).all()


def test_matrix_covers_every_mode(matrix_result):
    events, scores = matrix_result
    assert len(scores) == 2 * len(CreditMode)
    assert set(scores["mode"]) == {m.value for m in CreditMode}
    assert len(events) == len(scores)


def test_unknown_population_is_rejected():
    with pytest.raises(ValueError, match="unknown population"):
        evaluate_atlas_matrix(population_ids=("not_real",), draws=200)


def test_graft_probability_is_reported_alongside_the_median(matrix_result):
    """A median pass must never be mistaken for a robust one."""
    _, scores = matrix_result
    assert "graft_probability_of_target" in scores.columns
    assert scores["graft_probability_of_target"].between(0.0, 1.0).all()
    assert not scores["graft_robust"].any()


# --------------------------------------------------------------------------
# Report and outputs.
# --------------------------------------------------------------------------


def test_report_leads_with_the_disclaimer_and_the_falsifier(matrix_result, sweep):
    _, scores = matrix_result
    text = atlas_report(scores, sweep, arm_ablation(draws=1500))
    assert DISCLAIMER in text
    assert "assumption" in text
    assert "At zero coupling sequencing adds +0.0 pp" in text
    assert "It is not a cure." in text
    assert "WITHDRAWN" in text
    assert "THESE NUMBERS ASSUME PHASE 1 SUCCEEDED" in text


def test_outputs_are_quarantined_and_tagged(tmp_path):
    target = tmp_path / "atlas"
    summary = run_atlas_study(target, population_ids=("stable_ccs",), draws=800)
    assert summary["status"] == "invented_untrialled"
    assert summary["cures_declared"] == 0
    assert (target / "DISCLAIMER.txt").read_text().strip() == DISCLAIMER
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 9
    assert not any(name.startswith("solutions") for name in written)


def test_invented_tables_carry_the_invented_tag(tmp_path):
    target = tmp_path / "atlas"
    run_atlas_study(target, population_ids=("stable_ccs",), draws=800)
    scorecard = pd.read_csv(target / "tables" / "atlas_scorecard.csv")
    assert scorecard["status"].eq(ATLAS_INVENTED_TAG).all()
    sweep = pd.read_csv(target / "tables" / "atlas_bed_coupling_sweep.csv")
    assert sweep["status"].eq(ATLAS_INVENTED_TAG).all()


def test_cli_atlas_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    target = tmp_path / "cli-atlas"
    code = main(
        ["atlas", "--output-dir", str(target), "--populations", "stable_ccs", "--draws", "600"]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert DISCLAIMER in out
    assert "Wrote 9 artefacts" in out


def test_cli_atlas_rejects_a_bad_coupling(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(
        [
            "atlas",
            "--output-dir",
            str(tmp_path / "bad"),
            "--populations",
            "stable_ccs",
            "--bed-coupling",
            "2.0",
            "--draws",
            "200",
        ]
    )
    assert code == 1
    assert "Error:" in capsys.readouterr().err
