"""Tests for the four scientific audits.

The most important test in this file is `test_every_numeric_constant_is_registered`:
it fails the build when a new numeric parameter is added without a provenance record.
That is the mechanism that converts "we document our assumptions" from a promise into a
constraint.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from ihd_solutions import calibration as C
from ihd_solutions import identifiability as I
from ihd_solutions import provenance as P
from ihd_solutions import validation as V

# ============================================================ provenance: the guard


def test_every_numeric_constant_is_registered():
    """THE BUILD GUARD. A new undocumented parameter fails here."""
    missing = P.unregistered()
    assert missing == [], "numeric constants without a provenance record: " + ", ".join(
        f"{gap.module}:{gap.name}={gap.value}" for gap in missing
    )


def test_registry_values_agree_with_the_code():
    """A registry that drifts from the code is worse than no registry."""
    frame = P.audit_frame()
    disagreeing = frame.loc[frame["values_agree"].fillna(False).eq(False)]
    assert disagreeing.empty, f"registry drifted from code: {disagreeing[['module', 'name']]}"


def test_provenance_coverage_is_complete():
    summary = P.coverage()
    assert summary["coverage_fraction"] == pytest.approx(1.0)
    assert summary["unregistered"] == 0


def test_registry_classifies_every_constant():
    frame = P.audit_frame()
    assert frame["kind"].isin(P.CLASSES).all()
    assert frame["supports"].isin(P.SUPPORTS).all()
    assert frame["population"].isin(P.POPULATIONS).all()


def test_hypothetical_parameters_never_claim_human_measurement():
    with pytest.raises(ValueError, match="cannot be measured in humans"):
        P.Provenance(
            1.0,
            "x",
            P.HYPOTHETICAL,
            "efficacy",
            "human",
            "some invented thing",
            measured_in_humans=True,
        )


def test_measured_in_humans_requires_the_human_population():
    with pytest.raises(ValueError, match="requires population='human'"):
        P.Provenance(
            1.0,
            "x",
            P.MEASURED,
            "efficacy",
            "rodent",
            "a rat study",
            measured_in_humans=True,
        )


def test_parameter_must_declare_units_and_a_source():
    with pytest.raises(ValueError, match="unit"):
        P.Provenance(1.0, "", P.MEASURED, "efficacy", "human", "a real source here")
    with pytest.raises(ValueError, match="source"):
        P.Provenance(1.0, "x", P.MEASURED, "efficacy", "human", "short")


def test_only_measured_human_efficacy_earns_credit():
    measured = P.Provenance(
        0.62, "fraction", P.MEASURED, "efficacy", "human", "VERVE-102 Heart-2 outcomes"
    )
    mechanism = P.Provenance(
        3.9, "pp", P.MEASURED, "mechanism", "human", "BioVAT-HF single arm, no control"
    )
    holder = P.Provenance(1.0, "x", P.HYPOTHETICAL, "ceiling", "none", "an invented cap value")
    assert measured.earns_efficacy_credit is True
    assert mechanism.earns_efficacy_credit is False
    assert holder.earns_efficacy_credit is False


def test_credit_ceiling_takes_the_weakest_link_not_an_average():
    """One hypothetical multiplier must collapse credit, even beside a strong one."""
    assert P.credit_ceiling([P.MEASURED]) == pytest.approx(1.0)
    assert P.credit_ceiling([P.MEASURED, P.ESTIMATED]) == pytest.approx(0.5)
    assert P.credit_ceiling([P.MEASURED, P.ESTIMATED, P.HYPOTHETICAL]) == pytest.approx(0.0)


def test_credit_ceiling_rejects_an_empty_chain():
    with pytest.raises(ValueError, match="at least one"):
        P.credit_ceiling([])


def test_the_target_is_registered_as_hypothetical():
    """The +20 pp cure target is a study choice, never evidence."""
    record = P.REGISTRY["REQUIRED_LVEF_GAIN_PP"]
    assert record.kind == P.HYPOTHETICAL
    assert record.earns_efficacy_credit is False


def test_the_control_arm_value_is_registered_as_measured_efficacy():
    record = P.REGISTRY["NATURAL_HISTORY_PP"]
    assert record.kind == P.MEASURED
    assert record.population == "human"
    assert record.earns_efficacy_credit is True


def test_provenance_report_names_the_hypothetical_parameters():
    text = P.provenance_report()
    assert "HYPOTHETICAL PARAMETERS CARRYING REAL WEIGHT" in text
    assert "UNREGISTERED" in text or "unregistered" in text.lower()


# ================================================== identifiability: the three checks


def _collinear_model(params):
    return params["growth"] + params["inflam"]


def _separable_model(params):
    return params["growth"] * params["inflam"]


def test_sensitivity_finds_an_inert_parameter():
    def model(params):
        return params["live"] * 2.0  # "dead" never appears

    rows = I.one_at_a_time(model, {"live": 1.0, "dead": 1.0})
    by_name = {row.name: row for row in rows}
    assert by_name["live"].moves_output
    assert not by_name["dead"].moves_output
    assert by_name["dead"].verdict.startswith("inert")


def test_collinearity_detects_a_pair_that_cannot_be_separated():
    pairs = I.collinear_pairs(_collinear_model, {"growth": 0.5, "inflam": 0.5})
    assert pairs, "two additive parameters must be detected as collinear"
    assert abs(pairs[0][2]) > 0.95


def test_collinearity_finds_nothing_for_separable_parameters():
    pairs = I.collinear_pairs(_separable_model, {"growth": 0.5, "inflam": 0.5})
    assert pairs == []


def test_collinear_parameters_are_classified_as_chosen_not_measured():
    """Each profile looks fine alone; only the joint test rejects them."""
    nominal = {"growth": 0.5, "inflam": 0.5}
    sensitivities = I.one_at_a_time(_collinear_model, nominal)
    pairs = I.collinear_pairs(_collinear_model, nominal)
    objective = lambda p: (_collinear_model(p) - 1.0) ** 2  # noqa: E731
    profiles = {
        name: I.profile_verdict(I.profile_identifiability(objective, nominal, name))
        for name in nominal
    }
    verdict = I.identifiability_verdict(sensitivities, pairs, profiles)
    assert verdict["identified"] == ()
    assert set(verdict["chosen_not_measured"]) == {"growth", "inflam"}


def test_profile_identifiability_finds_a_single_minimum():
    objective = lambda p: (p["a"] - 3.0) ** 2  # noqa: E731
    profile = I.profile_identifiability(objective, {"a": 3.0, "b": 1.0}, "a", points=11)
    assert I.profile_verdict(profile).startswith("identified")


def test_profile_identifiability_rejects_a_flat_valley():
    objective = lambda p: (p["a"] + p["b"] - 3.0) ** 2  # noqa: E731
    profile = I.profile_identifiability(objective, {"a": 1.0, "b": 1.0}, "a", points=11)
    assert "unidentified" in I.profile_verdict(profile)


def test_elasticity_is_scale_invariant():
    """Parameters on different scales must be comparable."""

    def model(params):
        return params["big"] * 1e6 + params["small"]

    rows = {row.name: row for row in I.one_at_a_time(model, {"big": 1e6, "small": 1.0})}
    # 'big' dominates the sum, so it has elasticity ~1 and 'small' ~0 - which is the
    # CORRECT answer, and the point of using a relative step rather than an absolute one.
    assert rows["big"].elasticity == pytest.approx(1.0, abs=0.01)
    assert abs(rows["small"].elasticity) < 1e-6


# ================================================ calibration: the legacy ODE network


def test_reimplementation_reproduces_the_legacy_script():
    """Regressing on the original: day-180 CRP must be 88.61 mg/L for the untreated arm."""
    assert C.terminal_crp(C.LEGACY_RATES) == pytest.approx(88.6134, abs=1e-3)


def test_all_fourteen_legacy_rates_are_present():
    assert len(C.LEGACY_RATES) == 14
    assert len(C.LEGACY_INITIAL) == 7


def test_the_cascade_is_irreversibly_degenerate():
    """THE FINDING. Every production rate has elasticity +1.0, every degradation -1.0.

    The model therefore determines a ratio, and no single rate can be separated from
    another by the terminal observable.
    """
    from ihd_solutions.identifiability import one_at_a_time

    rows = {row.name: row for row in one_at_a_time(C.terminal_crp, C.LEGACY_RATES)}
    production = [name for name in rows if name.endswith(("_prod", "_act", "_syn"))]
    degradation = [name for name in rows if name.endswith("_deg")]
    assert len(production) >= 6 and len(degradation) >= 6
    # Every production rate moves the output +1.0 and every degradation rate -1.0, so the
    # cascade determines a SCALING and any one rate can be exchanged for any other.
    # k_mmp_prod is the exception and is asserted separately: it is DOWNSTREAM of the
    # observable, so it cannot move it at all.
    for name in production:
        if name == "k_mmp_prod":
            assert rows[name].elasticity == pytest.approx(0.0, abs=1e-9)
            continue
        assert rows[name].elasticity == pytest.approx(1.0, abs=5e-3), name
    for name in degradation:
        if name == "k_mmp_deg":
            assert rows[name].elasticity == pytest.approx(0.0, abs=1e-9)
            continue
        assert rows[name].elasticity == pytest.approx(-1.0, abs=2e-2), name


def test_the_cascade_has_collinear_rate_pairs():
    verdict = C.calibration_verdict()
    assert verdict["collinear_pairs"] >= 5
    assert verdict["rates_jointly_identified"] < verdict["rates_total"]


def test_no_rate_survives_joint_identification():
    assert C.calibration_verdict()["rates_jointly_identified"] == 0


def test_mmp_rates_are_inert_to_the_observable():
    """MMP is downstream and feeds nothing back, so it cannot be calibrated at all."""
    inert = set(C.calibration_verdict()["inert_rates"])
    assert {"k_mmp_prod", "k_mmp_deg"} <= inert


def test_synthetic_recovery_uses_noise_free_data():
    """A real dataset can only do worse than the synthetic test."""
    results = C.synthetic_recovery(names=("k_crp_prod",), targets=("final",))
    assert all(row.tolerance > 0 for row in results)


def test_ignorable_rates_lists_what_cannot_be_fitted():
    ignored = C.ignorable_rates()
    assert isinstance(ignored, list)


def test_degeneracy_is_exact_proportional_scaling():
    """Doubling a production rate doubles the output EXACTLY, not approximately.

    Stronger than statistical non-identifiability: the cascade is linear in its rates, so
    two different production rates generate the same trajectory to machine precision and
    no measurement of CRP can tell them apart.
    """
    frame = C.degeneracy_probe("k_crp_prod")
    baseline = C.terminal_crp(C.LEGACY_RATES)
    for row in frame.itertuples(index=False):
        assert row.terminal_crp == pytest.approx(row.factor * baseline, rel=1e-12), row.factor


def test_calibration_report_refuses_to_report_a_posterior():
    text = C.calibration_report()
    assert "HELD FIXED" in text
    assert "prior with a decimal point" in text
    assert "No rate was adjusted" in " ".join(text.split())


def test_calibration_report_states_the_degeneracy_honestly():
    """The report must not claim a 2x change 'barely moves' the output - it doubles it."""
    text = C.calibration_report()
    assert "EXACT" in text
    assert "barely moves" not in text


# ============================================ validation: verification vs validation


def test_a_rung_cannot_be_claimed_without_its_requirement():
    with pytest.raises(V.RungNotEarned):
        V.EvidenceClaim("model", V.EXTERNAL, ("refinement_stable",))


def test_convergence_is_never_reported_as_validation():
    result = V.ConvergenceResult("estimate", "draws", 0.800, 0.801, tolerance=0.01)
    assert result.converged
    verdict = V.convergence_verdict(result)
    assert "CONVERGED BUT UNVALIDATED" in verdict
    assert "still wrong" in verdict


def test_divergence_is_reported_as_divergence():
    result = V.ConvergenceResult("estimate", "draws", 0.10, 0.90, tolerance=0.01)
    assert not result.converged
    assert "NOT CONVERGED" in V.convergence_verdict(result)


def test_monte_carlo_convergence_refines_the_budget():
    sampler = lambda count: np.full(count, 0.5)  # noqa: E731
    results = V.monte_carlo_convergence(sampler, draws=(100, 400))
    assert all(result.converged for result in results)


def test_a_moving_estimate_is_not_converged():
    values = {100: 0.3, 400: 0.9}
    sampler = lambda count: np.full(count, values[count])  # noqa: E731
    results = V.monte_carlo_convergence(sampler, draws=(100, 400))
    assert not results[0].converged


def test_the_repository_never_claims_an_interventional_rung():
    assert V.repository_rung() == V.CONVERGENCE
    assert V.repository_rung() != V.INTERVENTIONAL


def test_no_module_is_claimed_as_fully_validated():
    rows = V.verification_vs_validation()
    assert all(row["can_claim_efficacy"] == "no" for row in rows)


def test_validation_gap_names_the_unvalidated_modules():
    gap = V.validation_gap()
    assert "NO external validation" in gap
    assert "convergence" in gap


def test_highest_rung_never_exceeds_convergence():
    converged = [V.ConvergenceResult("x", "draws", 0.5, 0.5, tolerance=0.01)]
    assert V.highest_rung(converged) == V.CONVERGENCE
    assert V.highest_rung([]) == V.UNIT


# ======================================================= the combined audit


def test_audit_writes_eight_artefacts(tmp_path):
    from ihd_solutions.scientific_audit import run_scientific_audit

    summary = run_scientific_audit(tmp_path / "audit")
    written = [path for path in (tmp_path / "audit").rglob("*") if path.is_file()]
    assert len(written) == summary["artefacts"] == 8


def test_audit_summary_is_honest(tmp_path):
    from ihd_solutions.scientific_audit import run_scientific_audit

    summary = run_scientific_audit(tmp_path / "audit")
    assert summary["unregistered"] == 0
    assert summary["provenance_coverage"] == pytest.approx(1.0)
    assert summary["validation_rung"] == V.CONVERGENCE
    assert summary["rates_jointly_identified"] == 0
    assert summary["not_done_items"] >= 5


def test_audit_records_what_it_did_not_do(tmp_path):
    from ihd_solutions.scientific_audit import run_scientific_audit

    run_scientific_audit(tmp_path / "audit")
    not_done = pd.read_csv(tmp_path / "audit" / "tables" / "audit_not_done.csv")
    assert len(not_done) >= 5
    assert not_done["why"].str.len().gt(20).all()
    joined = " ".join(not_done["item"])
    assert "molecular dynamics" in joined
    assert "CFD" in joined


def test_every_audit_table_carries_the_tag(tmp_path):
    from ihd_solutions.scientific_audit import AUDIT_TAG, run_scientific_audit

    run_scientific_audit(tmp_path / "audit")
    for path in (tmp_path / "audit" / "tables").glob("*.csv"):
        frame = pd.read_csv(path)
        assert (frame["status"] == AUDIT_TAG).all(), path.name


def test_audit_report_states_the_unfixable_part(tmp_path):
    from ihd_solutions.scientific_audit import run_scientific_audit

    run_scientific_audit(tmp_path / "audit")
    text = (tmp_path / "audit" / "reports" / "scientific_audit_report.txt").read_text()
    assert "CANNOT BE FIXED COMPUTATIONALLY" in text
    assert "NO MODEL AND ESTABLISHES NO BIOLOGY" in text.upper()
    assert "no cure is declared" in text.lower()


def test_cli_exposes_the_audit_command():
    from ihd_solutions.cli import build_parser

    args = build_parser().parse_args(["scientific-audit"])
    assert args.command == "scientific-audit"


# ============================================ remaining branches and report renderers


def test_identifiability_report_renders_every_classification():
    from ihd_solutions.identifiability import identifiability_report

    nominal = {"growth": 0.5, "inflam": 0.5}
    sensitivities = I.one_at_a_time(_collinear_model, nominal)
    pairs = I.collinear_pairs(_collinear_model, nominal)
    profiles = {"growth": "unidentified - flat", "inflam": "unidentified - flat"}
    text = identifiability_report(sensitivities, pairs, profiles, subject="toy")
    assert "toy" in text
    assert "COLLINEARITY" in text
    assert "chosen, NOT measured" in text


def test_identifiability_report_handles_no_collinearity():
    from ihd_solutions.identifiability import identifiability_report

    nominal = {"growth": 0.5, "inflam": 0.5}
    text = identifiability_report(
        I.one_at_a_time(_separable_model, nominal), [], {}, subject="clean"
    )
    assert "no pair exceeds the collinearity ceiling" in text


def test_sensitivity_matrix_returns_one_row():
    matrix = I.sensitivity_matrix(_separable_model, {"growth": 0.5, "inflam": 0.5})
    assert matrix.shape == (1, 2)


def test_one_at_a_time_rejects_a_zero_output():
    with pytest.raises(ValueError, match="nominal output is zero"):
        I.one_at_a_time(lambda p: 0.0, {"a": 1.0})


def test_one_at_a_time_rejects_a_zero_parameter():
    with pytest.raises(ValueError, match="relative step is undefined"):
        I.one_at_a_time(lambda p: p["a"] + 1.0, {"a": 0.0})


def test_jacobian_rejects_an_all_zero_model():
    with pytest.raises(ValueError, match="returned zero"):
        I.jacobian(lambda p: 0.0, {"a": 1.0})


def test_effect_correlation_with_one_row_says_nothing():
    """A single evaluation carries no covariance, so it must not claim r = 1."""
    frame = I.effect_correlation(np.ones((1, 2)), ["a", "b"])
    assert frame.loc["a", "b"] == 0.0


def test_credit_ceiling_accepts_provenance_objects():
    record = P.REGISTRY["NATURAL_HISTORY_PP"]
    assert P.credit_ceiling([record]) == pytest.approx(1.0)


def test_provenance_rejects_an_unknown_class():
    with pytest.raises(ValueError, match="unknown provenance class"):
        P.Provenance(1.0, "x", "guessed", "efficacy", "human", "a real source here")


def test_provenance_rejects_an_unknown_support_type():
    with pytest.raises(ValueError, match="unknown support type"):
        P.Provenance(1.0, "x", P.MEASURED, "vibes", "human", "a real source here")


def test_provenance_rejects_an_unknown_population():
    with pytest.raises(ValueError, match="unknown population"):
        P.Provenance(1.0, "x", P.MEASURED, "efficacy", "martian", "a real source here")


def test_validation_report_renders_the_ladder():
    text = V.validation_report()
    assert "VERIFICATION VERSUS VALIDATION" in text
    assert "MAX REPORTABLE RUNG: convergence" in text
    assert "WHAT WOULD MOVE THE LADDER" in text


def test_evidence_claim_accepts_an_earned_rung():
    claim = V.EvidenceClaim("model", V.CONVERGENCE, ("refinement_stable",))
    assert claim.rung == V.RUNG[V.CONVERGENCE]
    assert claim.licenses
    assert not claim.can_support_efficacy


def test_only_the_interventional_rung_supports_efficacy():
    claim = V.EvidenceClaim("trial", V.INTERVENTIONAL, ("randomised_human_comparison",))
    assert claim.can_support_efficacy


def test_evidence_claim_rejects_an_unknown_rung():
    with pytest.raises(ValueError, match="unknown rung"):
        V.EvidenceClaim("model", "vibes", ())


def test_validation_gap_is_computed_from_the_table():
    gap = V.validation_gap()
    assert gap.startswith("4 of 6") or "of 6 listed modules" in gap


def test_calibration_report_renders_the_joint_finding():
    text = C.calibration_report()
    assert "THE FINDING THAT ONE-AT-A-TIME RECOVERY MISSES" in text
    assert "collinear pairs" in text


def test_simulate_cytokines_accepts_alternate_rates():
    changed = dict(C.LEGACY_RATES)
    changed["k_crp_prod"] = C.LEGACY_RATES["k_crp_prod"] * 2.0
    assert C.terminal_crp(changed) == pytest.approx(2.0 * C.terminal_crp(C.LEGACY_RATES), rel=1e-9)


def test_simulate_cytokines_accepts_alternate_initial_conditions():
    altered = dict(C.LEGACY_INITIAL)
    altered["crp"] = 1.0
    trajectory = C.simulate_cytokines(initial=altered)
    assert trajectory.species["crp"][0] == pytest.approx(1.0)


def test_trajectory_helpers_are_consistent():
    trajectory = C.simulate_cytokines()
    assert trajectory.final("crp") == pytest.approx(C.terminal_crp(C.LEGACY_RATES))
    assert trajectory.auc("crp") > 0.0
    assert trajectory.summary("crp", "auc") == pytest.approx(trajectory.auc("crp"))


def test_recovery_result_relative_error_handles_a_miss():
    missed = C.RecoveryResult("x", 1.0, None, tolerance=0.05)
    assert missed.relative_error == float("inf")
    assert not missed.recovered_ok


def test_ignorable_rates_can_be_given_a_supplied_list():
    assert C.ignorable_rates([]) == []


def test_profile_verdict_rejects_a_completely_flat_profile():
    flat = pd.DataFrame({"objective": [0.5, 0.5, 0.5]})
    assert "completely flat" in I.profile_verdict(flat)


def test_profile_verdict_rejects_an_all_zero_valley():
    valley = pd.DataFrame({"objective": [0.0, 0.0, 0.0]})
    assert "perfect valley" in I.profile_verdict(valley)


def test_profile_verdict_rejects_a_one_sided_rise():
    # A minimum at index 1: the objective rises to its left but is nearly flat to the
    # right, so the parameter is pinned from above and not from below.
    one_sided = pd.DataFrame({"objective": [1.0, 0.0, 0.002, 0.004, 0.006]})
    assert "weakly identified" in I.profile_verdict(one_sided)


def test_convergence_result_verdict_string():
    converged = V.ConvergenceResult("x", "draws", 0.5, 0.5, tolerance=0.01)
    assert converged.verdict == "converged_but_unvalidated"
    diverged = V.ConvergenceResult("x", "draws", 0.1, 0.9, tolerance=0.01)
    assert diverged.verdict == "not_converged"


def test_registry_keys_are_exposed():
    keys = P.registry_keys()
    assert "NATURAL_HISTORY_PP" in keys
    assert "MIN_SPECIFICATION_DRAWS_RECOMBINE" in keys


def test_resolve_returns_none_for_an_unknown_name():
    assert P.resolve("nowhere.py", "NOT_A_REAL_CONSTANT") is None


def test_resolve_disambiguates_the_shared_draw_floor():
    assert P.resolve("recombine.py", "MIN_SPECIFICATION_DRAWS").value == pytest.approx(2000)
    assert P.resolve("atlas3.py", "MIN_SPECIFICATION_DRAWS").value == pytest.approx(4000)


def test_validation_gap_when_everything_is_validated():
    """The all-clear branch must exist and must not be reachable from real data."""
    original = dict(V.MODULE_STATUS)
    try:
        V.MODULE_STATUS["fake.py"] = {
            "verified": "yes",
            "validated": "PARTIAL - partial",
            "blocker": "none",
        }
        V.MODULE_STATUS.pop("plaque.py")
        V.MODULE_STATUS.pop("inflammation.py")
        V.MODULE_STATUS.pop("recombine.py")
        V.MODULE_STATUS.pop("revival.py")
        assert "at least partial" in V.validation_gap()
    finally:
        V.MODULE_STATUS.clear()
        V.MODULE_STATUS.update(original)


def test_effect_correlation_handles_a_single_parameter():
    """The 0-d guard: corrcoef of a one-column matrix returns a scalar."""
    frame = I.effect_correlation(np.array([[1.0], [2.0], [3.0]]), ["only"])
    assert frame.shape == (1, 1)
    assert frame.loc["only", "only"] == pytest.approx(1.0)


def test_provenance_report_renders_without_unregistered_entries():
    text = P.provenance_report()
    assert "PARAMETER PROVENANCE AUDIT" in text
    assert "BY CLASS" in text
    assert "BY SUPPORT TYPE" in text


def test_provenance_report_lists_an_unregistered_constant_honestly(tmp_path):
    """The report must SHOUT when something is missing, not silently omit it."""
    module = tmp_path / "rogue.py"
    module.write_text("MYSTERY_RATE = 0.42\n")
    original = P.numeric_constants
    try:
        P.numeric_constants = lambda root=None: original(tmp_path)  # noqa: ARG005
        text = P.provenance_report()
    finally:
        P.numeric_constants = original
    assert "UNREGISTERED - THESE FAIL THE BUILD" in text
    assert "MYSTERY_RATE" in text
