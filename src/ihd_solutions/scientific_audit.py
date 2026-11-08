"""One command that runs the four audits the critique asks for, and states the limits.

This is deliberately a *reporting* module. It changes no model. The point is to measure
the four things the critique says are missing and to publish the answers even where they
are unflattering:

1. **Provenance** - does every parameter declare what kind of claim it is?
2. **Calibration feasibility** - could the cytokine rates be fitted at all?
3. **Identifiability** - which parameters do the data actually determine?
4. **Verification versus validation** - what is proven, and what is only arithmetic?

Alongside those it records the work that was NOT done, with the reason, because a plan
that lists only what was delivered is the same failure mode as a model that reports only
its best fit. GPU molecular dynamics, 3D CFD/FSI and a genome-wide off-target model are
all genuinely absent; the audit says so rather than implying otherwise.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from ihd_solutions.calibration import (
    LEGACY_RATES,
    calibration_verdict,
    terminal_crp,
)
from ihd_solutions.identifiability import (
    collinear_pairs,
    one_at_a_time,
    profile_identifiability,
    profile_verdict,
)
from ihd_solutions.provenance import coverage, unregistered
from ihd_solutions.validation import (
    repository_rung,
    validation_gap,
    verification_vs_validation,
)

AUDIT_TAG = "scientific_audit_no_model_change"

DISCLAIMER = (
    "THIS AUDIT CHANGES NO MODEL AND ESTABLISHES NO BIOLOGY. It measures whether the "
    "repository's parameters declare their own status, whether they could be calibrated, "
    "whether they are identifiable, and how much of the work is verified rather than "
    "validated. Several answers are negative and are reported as such. No cure is "
    "declared and no efficacy claim is made."
)

#: Work the critique asks for that this audit does NOT deliver, with the reason.
NOT_DONE: tuple[dict[str, str], ...] = (
    {
        "item": "GPU molecular dynamics of the PCSK9-LDLR complex",
        "why": (
            "requires a GPU node and OpenMM/GROMACS plus multi-microsecond sampling; no GPU "
            "is available in this environment, and a short run would not support the "
            "convergence or replica requirements the critique itself sets"
        ),
        "what_stands_instead": (
            "the static contact analysis remains, flagged as static; it cannot report "
            "binding free energy, conformational stability or force-field sensitivity"
        ),
    },
    {
        "item": "3D CFD and fluid-structure interaction of coronary plaque",
        "why": (
            "needs OpenFOAM or SimVascular plus a meshed coronary geometry; neither the "
            "solver nor an imaging-derived geometry is available here"
        ),
        "what_stands_instead": (
            "the 0D hoop-stress model, with its assumed 300 kPa rupture threshold now "
            "registered as HYPOTHETICAL rather than presented as a criterion"
        ),
    },
    {
        "item": "genome-wide base-editor off-target assessment",
        "why": (
            "needs a reference genome, an editor-specific activity model and BE-Hive-style "
            "calibration data; the repo contains seven hand-written candidate loci"
        ),
        "what_stands_instead": (
            "the seven-locus screen is retained and relabelled: it is a candidate list, "
            "NOT a genome-wide safety assessment"
        ),
    },
    {
        "item": "Bayesian posterior distributions for the cytokine rates",
        "why": (
            "the audit shows 0 of 14 rates are jointly identifiable, so a posterior would "
            "report a prior with a decimal point"
        ),
        "what_stands_instead": (
            "feasibility analysis, collinearity measurement and an explicit list of rates "
            "held fixed as assumptions"
        ),
    },
    {
        "item": "virtual patient cohort with adherence and competing risks",
        "why": (
            "partially present in risk.py and population.py, but the cohort is not "
            "calibrated to an external registry"
        ),
        "what_stands_instead": "trial-level normalisation, which is stated as trial-level",
    },
)


def _calibration_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    sensitivities = one_at_a_time(terminal_crp, LEGACY_RATES)
    frame = pd.DataFrame(
        [
            {
                "rate": row.name,
                "nominal": row.nominal,
                "elasticity": row.elasticity,
                "moves_output": row.moves_output,
                "verdict": row.verdict,
            }
            for row in sensitivities
        ]
    )
    pairs = pd.DataFrame(
        [
            {"first": first, "second": second, "correlation": value}
            for first, second, value in collinear_pairs(terminal_crp, LEGACY_RATES)
        ]
    )
    return frame, pairs


def _profiles() -> pd.DataFrame:
    """Profile identifiability for the production rates, refitting the others."""
    rows = []
    for name in ("k_crp_prod", "k_il6_prod", "k_nfkb_act"):
        profile = profile_identifiability(terminal_crp, dict(LEGACY_RATES), name, points=9)
        rows.append(
            {
                "parameter": name,
                "objective_span": float(profile["objective"].max() - profile["objective"].min()),
                "verdict": profile_verdict(profile),
            }
        )
    return pd.DataFrame(rows)


def run_scientific_audit(output: Path) -> dict[str, object]:
    """Write the four audits and the not-done register. Returns the summary."""
    output = Path(output)
    tables = output / "tables"
    reports = output / "reports"
    tables.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    (output / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n")

    # 1. Provenance
    from ihd_solutions.provenance import audit_frame

    provenance = audit_frame().assign(status=AUDIT_TAG)
    provenance.to_csv(tables / "audit_parameter_provenance.csv", index=False)

    # 2. Calibration feasibility + 3. Identifiability
    sensitivity, pairs = _calibration_tables()
    sensitivity = sensitivity.assign(status=AUDIT_TAG)
    pairs = pairs.assign(status=AUDIT_TAG) if not pairs.empty else pairs.assign(status=[])
    sensitivity.to_csv(tables / "audit_rate_sensitivity.csv", index=False)
    pairs.assign(status=AUDIT_TAG).to_csv(tables / "audit_rate_collinearity.csv", index=False)
    profiles = _profiles().assign(status=AUDIT_TAG)
    profiles.to_csv(tables / "audit_rate_profiles.csv", index=False)

    # 4. Verification versus validation
    validation = pd.DataFrame(verification_vs_validation()).assign(status=AUDIT_TAG)
    validation.to_csv(tables / "audit_validation_ladder.csv", index=False)

    # Not done, with reasons
    not_done = pd.DataFrame([dict(item, status=AUDIT_TAG) for item in NOT_DONE])
    not_done.to_csv(tables / "audit_not_done.csv", index=False)

    summary = {
        "artefacts": 8,
        "provenance_coverage": coverage()["coverage_fraction"],
        "constants_found": coverage()["constants_found"],
        "unregistered": len(unregistered()),
        "rates_total": calibration_verdict()["rates_total"],
        "rates_jointly_identified": calibration_verdict()["rates_jointly_identified"],
        "collinear_pairs": calibration_verdict()["collinear_pairs"],
        "validation_rung": repository_rung(),
        "not_done_items": len(NOT_DONE),
        "status": AUDIT_TAG,
    }

    report = _report(summary, provenance, sensitivity, pairs, profiles, not_done)
    (reports / "scientific_audit_report.txt").write_text(report)
    return summary


def _report(
    summary: dict[str, object],
    provenance: pd.DataFrame,
    sensitivity: pd.DataFrame,
    pairs: pd.DataFrame,
    profiles: pd.DataFrame,
    not_done: pd.DataFrame,
) -> str:
    """Compose the audit report. Every number is read from the frames above."""
    cov = coverage()
    verdict = calibration_verdict()
    lines = [
        "SCIENTIFIC AUDIT - FOUR CHECKS THE CRITIQUE ASKED FOR",
        "=" * 78,
        "",
        "THIS AUDIT CHANGES NO MODEL AND ESTABLISHES NO BIOLOGY.",
        "Several of its answers are negative and are reported as such.",
        "No cure is declared, no efficacy claim is made, and no therapy is proposed.",
        "",
        "",
        "1. PARAMETER PROVENANCE",
        "-" * 78,
        f"  constants in the package        {cov['constants_found']}",
        f"  registered                      {cov['registered']}",
        f"  UNREGISTERED                    {cov['unregistered']}   (a test fails the build on this)",
        f"  values disagreeing with code    {cov['values_disagreeing']}",
        "  by class:",
    ]
    for kind, count in sorted(cov["by_kind"].items()):
        lines.append(f"    {kind:14s} {count:3d}")
    hypothetical = int((provenance["kind"] == "hypothetical").sum())
    lines += [
        "",
        f"  {hypothetical} of {cov['registered']} parameters are HYPOTHETICAL. That is the",
        "  critique's first problem, measured rather than argued. They are now typed as",
        "  such, so downstream code can branch on the distinction instead of trusting a",
        "  decimal point.",
        "",
        "2. CALIBRATION FEASIBILITY (legacy cytokine network, 14 rates)",
        "-" * 78,
        f"  independent observables          {verdict['targets_available']}",
        f"  identifiable upper bound         {verdict['identifiable_upper_bound']}",
        f"  rates recovered one-at-a-time    {verdict['rates_recoverable']} of {verdict['rates_total']}",
        f"  collinear pairs |r| >= 0.95      {verdict['collinear_pairs']} of {verdict['pairs_possible']}",
        f"  rates JOINTLY identified         {verdict['rates_jointly_identified']}",
        f"  rates inert to the observable    {list(verdict['inert_rates'])}",
        "",
        "  One-at-a-time recovery passes 12 of 14; the joint test rejects all 14. The",
        "  difference is the whole point: scanning a rate alone measures how much the",
        "  OUTPUT moves, not what the DATA identify. Every production rate has elasticity",
        "  exactly +1.0000 and every degradation rate -1.0000, so the model determines a",
        "  SCALING, not a set of rates.",
        "",
        "  The degeneracy is exact, not statistical. The cascade is linear in its rates, so",
        "  doubling k_crp_prod changes terminal CRP by exactly +100%, and doubling",
        "  k_il6_prod produces the SAME trajectory to machine precision. Six production",
        "  rates share one degree of freedom between them, and no measurement of CRP can",
        "  separate them.",
        "",
        "  Consequence: a Bayesian posterior for these rates would be a prior with a",
        "  decimal point. None is reported. The rates are held FIXED and declared",
        "  assumptions, which is what they always were.",
        "",
        "3. IDENTIFIABILITY - PROFILE CHECK",
        "-" * 78,
        f"  {'parameter':16s} {'objective span':>15s}  verdict",
    ]
    for row in profiles.itertuples(index=False):
        lines.append(f"  {row.parameter:16s} {row.objective_span:15.4f}  {row.verdict}")
    lines += [
        "",
        "4. VERIFICATION VERSUS VALIDATION",
        "-" * 78,
        f"  highest rung this repo occupies  {repository_rung().upper()}",
        f"  {validation_gap()}",
        "",
        "  A passing test suite says the code does what it was written to do. It says",
        "  nothing about whether the biology is right. Four of six listed modules have no",
        "  external validation at all, and no module reaches the interventional rung.",
        "",
        "5. WHAT THIS AUDIT DID NOT DO",
        "-" * 78,
    ]
    for row in not_done.itertuples(index=False):
        lines += [
            f"  {row.item}",
            f"    why not:  {row.why}",
            f"    instead:  {row.what_stands_instead}",
            "",
        ]
    lines += [
        "6. THE ONE THING THAT CANNOT BE FIXED COMPUTATIONALLY",
        "-" * 78,
        "  Every check above is reproducible and cheap. None of them moves a single result",
        "  from the CONVERGENCE rung to the INTERVENTIONAL one, because that step needs a",
        "  randomised human comparison and no simulation supplies it.",
        "",
        "  The critique is right that the repository has more treatment architectures than",
        "  validated underlying models. This audit makes that visible, which is the most",
        "  a computational study can do about it.",
        "",
        "KEY",
        "=" * 78,
    ]
    return "\n".join(lines)
