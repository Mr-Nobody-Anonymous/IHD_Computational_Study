"""Command line interface for the solution simulator."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

from ihd_solutions.catalog import catalog_frame, load_catalog, validate_combination
from ihd_solutions.evidence import CreditMode, load_anchors
from ihd_solutions.pipeline import DEFAULT_OUTPUT_DIR, PipelineConfig, run
from ihd_solutions.population import load_populations
from ihd_solutions.portfolio import EnumerationSpec, enumerate_regimens, expected_regimen_count
from ihd_solutions.risk import simulate_events
from ihd_solutions.units import time_grid


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ihd-solutions",
        description=(
            "Simulate and compare candidate solutions for ischaemic heart disease. "
            "Scenario modelling only; not clinical prediction."
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("therapies", help="list the therapy catalogue")
    subparsers.add_parser("anchors", help="list the trial-anchor registry")
    subparsers.add_parser("populations", help="list scenario populations")

    regimens = subparsers.add_parser("regimens", help="enumerate admissible regimens")
    regimens.add_argument("--max-optional", type=int, default=2)

    evaluate = subparsers.add_parser("evaluate", help="evaluate a single regimen")
    evaluate.add_argument("--therapies", nargs="+", required=True)
    evaluate.add_argument("--population", default="stable_ccs")
    evaluate.add_argument("--mode", choices=[mode.value for mode in CreditMode], default="strict")
    evaluate.add_argument("--horizon-years", type=float, default=5.0)
    evaluate.add_argument("--json", action="store_true", help="emit JSON instead of text")

    simulate = subparsers.add_parser("simulate", help="run the full simulation suite")
    simulate.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    simulate.add_argument("--horizon-years", type=float, default=5.0)
    simulate.add_argument("--step-days", type=float, default=30.0)
    simulate.add_argument("--max-optional", type=int, default=2)
    simulate.add_argument("--populations", nargs="+", default=None)
    simulate.add_argument("--modes", nargs="+", choices=[m.value for m in CreditMode], default=None)
    simulate.add_argument("--draws", type=int, default=120)
    simulate.add_argument("--seed", type=int, default=20261106)
    simulate.add_argument("--no-figures", action="store_true")

    hypothesis = subparsers.add_parser(
        "hypothesis",
        help="run the quarantined REVERT-1 invented-therapy study (NOT evidence)",
    )
    hypothesis.add_argument("--output-dir", type=Path, default=Path("results/hypotheses"))
    hypothesis.add_argument("--draws", type=int, default=400)
    hypothesis.add_argument("--seed", type=int, default=20261107)
    hypothesis.add_argument("--months", type=float, default=24.0)

    discover = subparsers.add_parser(
        "discover",
        help="search for a cure SPECIFICATION in simulation only (never a cure)",
    )
    discover.add_argument("--output-dir", type=Path, default=Path("results/discovery"))
    discover.add_argument("--draws", type=int, default=4000)
    discover.add_argument("--seed", type=int, default=20261107)
    discover.add_argument("--baseline-risk", type=float, default=0.1393)
    discover.add_argument("--reference-risk", type=float, default=0.0198)

    composite = subparsers.add_parser(
        "composite",
        help="assemble the best part per disease node and test across all populations/modes",
    )
    composite.add_argument("--output-dir", type=Path, default=Path("results/composite"))
    composite.add_argument("--populations", nargs="+", default=None)
    composite.add_argument("--horizon-years", type=float, default=5.0)
    composite.add_argument("--max-components", type=int, default=8)

    gaps = subparsers.add_parser(
        "gaps", help="attempt to close the three capability gaps (scenario, not evidence)"
    )
    gaps.add_argument("--output-dir", type=Path, default=Path("results/gaps"))
    gaps.add_argument(
        "--composite-matrix",
        type=Path,
        default=Path("results/composite/tables/composite_test_matrix.csv"),
    )

    resolve = subparsers.add_parser(
        "resolve", help="resolve the three residual gaps step by step (scenario, not evidence)"
    )
    resolve.add_argument("--output-dir", type=Path, default=Path("results/resolution"))
    resolve.add_argument(
        "--composite-matrix",
        type=Path,
        default=Path("results/composite/tables/composite_test_matrix.csv"),
    )

    atlas = subparsers.add_parser(
        "atlas", help="ATLAS-1 synthesis: best part of every construct, sequenced (invented)"
    )
    atlas.add_argument("--output-dir", type=Path, default=Path("results/atlas"))
    atlas.add_argument("--populations", nargs="+", default=None)
    atlas.add_argument("--bed-coupling", type=float, default=0.5)
    atlas.add_argument("--draws", type=int, default=3000)

    audit = subparsers.add_parser(
        "audit", help="self-audit ATLAS-1: does its operating point survive uncertainty?"
    )
    audit.add_argument("--output-dir", type=Path, default=Path("results/audit"))
    audit.add_argument("--draws", type=int, default=2500)

    evidence_audit = subparsers.add_parser(
        "audit-evidence", help="audit the evidence layer's criterion scores under trial uncertainty"
    )
    evidence_audit.add_argument("--output-dir", type=Path, default=Path("results/evidence_audit"))
    evidence_audit.add_argument("--populations", nargs="+", default=None)
    evidence_audit.add_argument("--draws", type=int, default=300)

    upgrade = subparsers.add_parser(
        "upgrade", help="ATLAS-2: upgrade the parts the audits broke (invented)"
    )
    upgrade.add_argument("--output-dir", type=Path, default=Path("results/upgrade"))
    upgrade.add_argument("--draws", type=int, default=800)

    atlas3 = subparsers.add_parser(
        "atlas3", help="ATLAS-3: staged retry + feedback for the graft arm (invented)"
    )
    atlas3.add_argument("--output-dir", type=Path, default=Path("results/atlas3"))
    atlas3.add_argument("--draws", type=int, default=2500)

    atlas4 = subparsers.add_parser(
        "atlas4", help="ATLAS-4: precedent-corrected graft model + exhaustive design search"
    )
    atlas4.add_argument("--output-dir", type=Path, default=Path("results/atlas4"))
    atlas4.add_argument("--draws", type=int, default=3000)

    solutions = subparsers.add_parser(
        "gap-solutions", help="search + design solutions for the six open science gaps"
    )
    solutions.add_argument("--output-dir", type=Path, default=Path("results/gap_solutions"))
    solutions.add_argument("--draws", type=int, default=1200)

    recombine = subparsers.add_parser(
        "recombine", help="compose the load-bearing mechanism of each arm and test overlaps"
    )
    recombine.add_argument("--output-dir", type=Path, default=Path("results/recombine"))
    recombine.add_argument("--draws", type=int, default=2500)

    revival = subparsers.add_parser(
        "revival",
        help="charge the composition for its missing control arm and re-search",
    )
    revival.add_argument("--output-dir", type=Path, default=Path("results/revival"))
    revival.add_argument("--draws", type=int, default=2500)

    audit_sci = subparsers.add_parser(
        "scientific-audit",
        help="provenance, calibration feasibility, identifiability, verification vs validation",
    )
    audit_sci.add_argument("--output-dir", type=Path, default=Path("results/scientific_audit"))

    unanswered = subparsers.add_parser(
        "unanswered",
        help="what this study cannot answer, with the cost of each gap quantified",
    )
    unanswered.add_argument("--output-dir", type=Path, default=Path("results/unanswered"))

    frontier = subparsers.add_parser(
        "frontier",
        help="what must hold before a simulation result is worth a lab's time",
    )
    frontier.add_argument("--output-dir", type=Path, default=Path("results/frontier"))

    interventions = subparsers.add_parser(
        "interventions",
        help="what the literature says about each transfer-gate item, and what happened here",
    )
    interventions.add_argument("--output-dir", type=Path, default=Path("results/interventions"))

    subparsers.add_parser("validate", help="validate all input registries and exit")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return _dispatch(args)
    except (ValueError, FileNotFoundError, AssertionError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1


def _dispatch(args: argparse.Namespace) -> int:
    if args.command == "interventions":
        from ihd_solutions.interventions import write_intervention_report

        summary = write_intervention_report(str(args.output_dir))
        for key, value in summary.items():
            print(f"{key}: {value}")
        print((args.output_dir / "intervention_report.txt").read_text())
        return 0

    if args.command == "frontier":
        from ihd_solutions.simulation_frontier import write_frontier_report

        summary = write_frontier_report(args.output_dir)
        for key, value in summary.items():
            print(f"{key}: {value}")
        print((args.output_dir / "frontier_report.txt").read_text())
        return 0

    if args.command == "unanswered":
        from ihd_solutions.gaps import write_gap_report

        summary = write_gap_report(args.output_dir)
        for key, value in summary.items():
            print(f"{key}: {value}")
        print((args.output_dir / "gap_report.txt").read_text())
        return 0

    if args.command == "scientific-audit":
        from ihd_solutions.scientific_audit import run_scientific_audit

        summary = run_scientific_audit(args.output_dir)
        for key, value in summary.items():
            print(f"{key}: {value}")
        print((args.output_dir / "reports" / "scientific_audit_report.txt").read_text())
        return 0

    if args.command == "revival":
        from ihd_solutions.revival import run_revival_study

        summary = run_revival_study(args.output_dir, draws=args.draws)
        for key, value in summary.items():
            print(f"{key}: {value}")
        print((args.output_dir / "reports" / "revival_report.txt").read_text())
        return 0

    if args.command == "recombine":
        from ihd_solutions.recombine import run_recombine_study

        summary = run_recombine_study(args.output_dir, draws=args.draws)
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "recombine_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "gap-solutions":
        from ihd_solutions.gap_solutions import run_gap_solutions_study

        summary = run_gap_solutions_study(args.output_dir, draws=args.draws)
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "gap_solutions_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "atlas4":
        from ihd_solutions.atlas4 import run_atlas4_study

        summary = run_atlas4_study(args.output_dir, draws=args.draws)
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "atlas4_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "atlas3":
        from ihd_solutions.atlas3 import run_atlas3_study

        summary = run_atlas3_study(args.output_dir, draws=args.draws)
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "atlas3_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "upgrade":
        from ihd_solutions.upgrade import run_upgrade_study

        summary = run_upgrade_study(args.output_dir, draws=args.draws)
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "upgrade_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "audit-evidence":
        from ihd_solutions.evidence_audit import run_evidence_audit_study

        summary = run_evidence_audit_study(
            args.output_dir,
            population_ids=tuple(args.populations) if args.populations else None,
            draws=args.draws,
        )
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "evidence_audit_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "audit":
        from ihd_solutions.audit import run_audit_study

        summary = run_audit_study(args.output_dir, draws=args.draws)
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "audit_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "atlas":
        from ihd_solutions.atlas import run_atlas_study

        summary = run_atlas_study(
            args.output_dir,
            population_ids=tuple(args.populations) if args.populations else None,
            bed_coupling=args.bed_coupling,
            draws=args.draws,
        )
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "atlas_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "resolve":
        from ihd_solutions.resolution import run_resolution_study

        summary = run_resolution_study(args.output_dir, matrix=pd.read_csv(args.composite_matrix))
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "resolution_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "gaps":
        from ihd_solutions.gap_solver import run_gap_study

        summary = run_gap_study(args.output_dir, matrix=pd.read_csv(args.composite_matrix))
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "gap_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "composite":
        from ihd_solutions.composite import run_composite_study

        summary = run_composite_study(
            args.output_dir,
            population_ids=tuple(args.populations) if args.populations else None,
            horizon_years=args.horizon_years,
            max_components=args.max_components,
        )
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "composite_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "discover":
        from ihd_solutions.discovery import run_discovery_study

        summary = run_discovery_study(
            args.output_dir,
            draws=args.draws,
            seed=args.seed,
            baseline_risk=args.baseline_risk,
            reference_risk=args.reference_risk,
        )
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "discovery_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "hypothesis":
        from ihd_solutions.hypothesis import run_hypothesis_study

        summary = run_hypothesis_study(
            args.output_dir, draws=args.draws, seed=args.seed, months=args.months
        )
        print(summary["disclaimer"])
        print()
        print((args.output_dir / "reports" / "revert1_design_report.txt").read_text())
        print(f"Wrote {summary['outputs']} artefacts to {summary['output_dir']}")
        return 0
    if args.command == "therapies":
        frame = catalog_frame(load_catalog())
        print(frame.to_string(index=False))
        return 0
    if args.command == "anchors":
        anchors = load_anchors()
        frame = pd.DataFrame(
            [
                {
                    "anchor_id": anchor.anchor_id,
                    "therapy_id": anchor.therapy_id,
                    "population_id": anchor.population_id,
                    "endpoint_class": anchor.endpoint_class.value,
                    "effect_unit": anchor.effect_unit.value,
                    "hr": anchor.hr,
                    "neutral_ci": anchor.is_neutral,
                    "evidence_tier": anchor.evidence_tier.value,
                }
                for anchor in sorted(anchors.values(), key=lambda a: a.anchor_id)
            ]
        )
        print(frame.to_string(index=False))
        return 0
    if args.command == "populations":
        populations = load_populations()
        frame = pd.DataFrame(
            [
                {
                    "population_id": population.population_id,
                    "label": population.label,
                    "annual_mace_rate": population.baseline_annual_mace_rate,
                    "ldl_c_mg_dl": population.ldl_c_mg_dl,
                    "hscrp_mg_l": population.hscrp_mg_l,
                }
                for population in sorted(populations.values(), key=lambda p: p.population_id)
            ]
        )
        print(frame.to_string(index=False))
        return 0
    if args.command == "regimens":
        catalog = load_catalog()
        spec = EnumerationSpec(max_optional=args.max_optional)
        regimens = enumerate_regimens(catalog, spec)
        expected = expected_regimen_count(catalog, spec)
        if len(regimens) != expected:
            raise AssertionError(f"enumerated {len(regimens)} but expected {expected}")
        for regimen in regimens:
            print(regimen.label)
        print(f"\n{len(regimens)} admissible regimens (combinatorial check: {expected})")
        return 0
    if args.command == "evaluate":
        catalog = load_catalog()
        anchors = load_anchors()
        populations = load_populations()
        if args.population not in populations:
            raise ValueError(f"unknown population: {args.population}")
        regimen = validate_combination(catalog, tuple(args.therapies))
        simulation = simulate_events(
            regimen,
            catalog,
            populations[args.population],
            anchors,
            time_grid(args.horizon_years),
            mode=CreditMode(args.mode),
        )
        row = simulation.summary_row()
        if args.json:
            print(json.dumps(row, indent=2, sort_keys=True, default=str))
        else:
            for key, value in row.items():
                print(f"{key}: {value}")
            print("\nCredit decisions:")
            for decision in simulation.decisions:
                state = "CREDITED" if decision.credited else "denied"
                print(f"  {decision.therapy_id}: {state} - {decision.reason}")
        return 0
    if args.command == "simulate":
        from ihd_solutions.pipeline import DEFAULT_POPULATIONS

        config = PipelineConfig(
            output_dir=args.output_dir,
            horizon_years=args.horizon_years,
            step_days=args.step_days,
            max_optional=args.max_optional,
            populations=tuple(args.populations) if args.populations else DEFAULT_POPULATIONS,
            modes=tuple(CreditMode(m) for m in args.modes)
            if args.modes
            else (CreditMode.STRICT, CreditMode.SURROGATE, CreditMode.EXPLORATORY),
            monte_carlo_draws=args.draws,
            seed=args.seed,
            make_figures=not args.no_figures,
        )
        summary = run(config)
        print(json.dumps(summary, indent=2, sort_keys=True))
        return 0
    if args.command == "validate":
        catalog = load_catalog()
        anchors = load_anchors()
        populations = load_populations()
        missing = sorted({anchor.therapy_id for anchor in anchors.values()} - set(catalog))
        if missing:
            raise ValueError(f"anchors reference unknown therapies: {missing}")
        unknown_populations = sorted(
            {anchor.population_id for anchor in anchors.values()} - set(populations)
        )
        if unknown_populations:
            raise ValueError(f"anchors reference unknown populations: {unknown_populations}")
        print(
            f"Validated {len(catalog)} therapies, {len(anchors)} anchors, "
            f"{len(populations)} populations."
        )
        return 0
    raise ValueError(f"unsupported command: {args.command}")
