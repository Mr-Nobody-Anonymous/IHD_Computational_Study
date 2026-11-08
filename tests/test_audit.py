"""Tests for the ATLAS-1 self-audit.

These pin a retraction. If a future edit makes the nominal clearance operating
point look robust again, or makes the sequencing gain independent of whether the
bed was actually cleared, one of these must fail.
"""

from __future__ import annotations

import pandas as pd
import pytest

from ihd_solutions.atlas import ATLAS_POTENCY
from ihd_solutions.audit import (
    AUDIT_DISCLAIMER,
    asymptotic_clearance_probability,
    audit_report,
    audited_criteria,
    clearance_draws,
    clearance_robustness,
    coupled_graft_outcome,
    coupled_sweep,
    minimum_robust_potency,
    potency_sweep,
    run_audit_study,
    source_threshold_for_perfect_modules,
)
from ihd_solutions.cure import CRITERION_THRESHOLD


@pytest.fixture(scope="module")
def operating_point():
    return clearance_robustness(draws=1200)


@pytest.fixture(scope="module")
def sweep() -> pd.DataFrame:
    return coupled_sweep(draws=1200)


# --------------------------------------------------------------------------
# The retraction itself.
# --------------------------------------------------------------------------


def test_the_nominal_operating_point_passes_but_is_not_robust(operating_point):
    assert operating_point.nominal >= CRITERION_THRESHOLD
    assert not operating_point.robust
    assert operating_point.nominal_is_misleading


def test_probability_at_the_claimed_operating_point_is_negligible(operating_point):
    assert operating_point.probability_of_target < 0.02


def test_the_median_is_nowhere_near_the_nominal(operating_point):
    assert operating_point.median < 0.2
    assert operating_point.nominal - operating_point.median > 0.5


def test_burden_reversal_is_retracted_in_the_corrected_scorecard():
    criteria = audited_criteria(draws=1200).set_index("criterion")
    row = criteria.loc["burden_reversal"]
    assert row["nominal_passes"]
    assert not row["robust_passes"]
    assert row["probability_of_target"] < 0.02


def test_myocardial_restoration_fails_both_ways():
    criteria = audited_criteria(draws=1200).set_index("criterion")
    row = criteria.loc["myocardial_restoration"]
    assert not row["nominal_passes"]
    assert not row["robust_passes"]


# --------------------------------------------------------------------------
# Reachability: what potency would actually be required.
# --------------------------------------------------------------------------


def test_clearance_is_monotone_in_potency_under_common_random_numbers():
    low = clearance_draws(1.5, draws=400)
    high = clearance_draws(5.0, draws=400)
    assert float(high.mean()) > float(low.mean())


def test_probability_rises_with_potency():
    sweep = potency_sweep(potencies=(1.284, 3.0, 20.0), draws=1000)
    probabilities = list(sweep["probability_of_target"])
    assert probabilities == sorted(probabilities)
    assert not sweep["robust_at_80pct"].any()


def test_no_potency_is_robust_while_the_source_term_is_uncertain():
    assert minimum_robust_potency(draws=500, source_uncertain=True, upper=64.0) is None


def test_a_robust_potency_exists_only_if_source_reduction_is_guaranteed():
    required = minimum_robust_potency(draws=500, source_uncertain=False, upper=64.0)
    assert required is not None
    assert required > 5.0 * ATLAS_POTENCY


def test_perfect_modules_still_need_a_minimum_source_reduction():
    threshold = source_threshold_for_perfect_modules()
    assert threshold is not None
    assert 0.0 < threshold < 0.5


def test_the_asymptotic_ceiling_is_below_the_robustness_bar():
    ceiling = asymptotic_clearance_probability()
    threshold = source_threshold_for_perfect_modules()
    assert ceiling == pytest.approx((0.5 - threshold) / 0.5, abs=1e-6)
    assert 0.5 < ceiling < 0.8


def test_zero_source_ceiling_means_the_target_is_unreachable():
    assert asymptotic_clearance_probability(source_reduction=0.0) == 0.0


def test_longer_horizons_do_not_rescue_the_operating_point():
    short = clearance_robustness(months=24.0, draws=800)
    long = clearance_robustness(months=240.0, draws=800)
    assert long.median > short.median
    assert not long.robust
    assert long.probability_of_target < 0.05


# --------------------------------------------------------------------------
# The second error: the arms are coupled, not independent.
# --------------------------------------------------------------------------


def test_zero_coupling_gives_exactly_no_paired_difference(sweep):
    row = sweep.loc[sweep["bed_coupling"].eq(0.0)].iloc[0]
    assert row["overstatement_pp"] == pytest.approx(0.0, abs=1e-9)
    assert row["coupled_probability"] == pytest.approx(row["decoupled_probability"], abs=1e-12)


def test_scoring_the_arms_independently_overstates_the_graft(sweep):
    coupled = sweep.loc[sweep["bed_coupling"].eq(1.0)].iloc[0]
    assert coupled["overstatement_pp"] > 15.0
    assert coupled["coupled_probability"] < coupled["decoupled_probability"]


def test_overstatement_grows_with_coupling(sweep):
    values = list(sweep["overstatement_pp"])
    assert values == sorted(values)


def test_the_bed_is_mostly_not_prepared(sweep):
    assert sweep["median_bed_quality"].max() < 0.15


def test_coupling_the_arms_leaves_the_graft_far_from_robust(sweep):
    assert not sweep["robust_at_80pct"].any()
    assert sweep["coupled_probability"].max() < 0.45


def test_coupled_graft_is_reproducible():
    first = coupled_graft_outcome(draws=600, seed=11)
    second = coupled_graft_outcome(draws=600, seed=11)
    assert first.probability_of_target == second.probability_of_target


def test_coupled_graft_rejects_bad_coupling():
    with pytest.raises(ValueError):
        coupled_graft_outcome(bed_coupling=1.5, draws=200)


def test_clearance_draws_reject_bad_inputs():
    with pytest.raises(ValueError):
        clearance_draws(potency=0.0, draws=100)
    with pytest.raises(ValueError):
        clearance_draws(source_reduction=1.5, draws=100)


# --------------------------------------------------------------------------
# Report and quarantine.
# --------------------------------------------------------------------------


def test_report_states_the_retraction_and_leads_with_the_disclaimer(operating_point, sweep):
    text = audit_report(
        operating_point,
        potency_sweep(potencies=(1.284, 5.0), draws=600),
        sweep,
        audited_criteria(draws=600),
        17.7,
        asymptotic_clearance_probability(),
        source_threshold_for_perfect_modules(),
    )
    assert AUDIT_DISCLAIMER in text
    assert "The nominal pass does not survive" in text
    assert "withdrawn" in text
    assert "2 of 5" in text


def test_outputs_are_quarantined_and_counted(tmp_path):
    target = tmp_path / "audit"
    summary = run_audit_study(target, draws=700)
    assert summary["status"] == "self_audit_retracts_nominal_pass"
    assert summary["burden_reversal_retracted"] is True
    assert summary["cures_declared"] == 0
    assert AUDIT_DISCLAIMER in (target / "DISCLAIMER.txt").read_text()
    written = {p.relative_to(target).as_posix() for p in target.rglob("*") if p.is_file()}
    assert len(written) == summary["outputs"] == 7
    assert not any(name.startswith("solutions") for name in written)


def test_every_table_is_tagged(tmp_path):
    target = tmp_path / "audit"
    run_audit_study(target, draws=700)
    for path in (target / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        assert "status" in frame.columns, path.name
        assert frame["status"].str.contains("atlas").all(), path.name


def test_cli_audit_subcommand(tmp_path, capsys):
    from ihd_solutions.cli import main

    code = main(["audit", "--output-dir", str(tmp_path / "cli-audit"), "--draws", "600"])
    assert code == 0
    out = capsys.readouterr().out
    assert AUDIT_DISCLAIMER in out
    assert "Wrote 7 artefacts" in out
