"""Run every solution simulation and write auditable tables, figures and a manifest."""

from __future__ import annotations

import hashlib
import json
import platform
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ihd_solutions import cure, editing, inflammation, regeneration, regression
from ihd_solutions.catalog import (
    Regimen,
    Therapy,
    TherapyClass,
    catalog_frame,
    load_catalog,
    validate_combination,
)
from ihd_solutions.evidence import (
    CreditMode,
    EndpointClass,
    TrialAnchor,
    load_anchors,
    summarise_tiers,
)
from ihd_solutions.lipids import trajectory_frame
from ihd_solutions.plaque import simulate_plaque
from ihd_solutions.population import Population, background_regimen_ids, load_populations
from ihd_solutions.portfolio import (
    EnumerationSpec,
    enumerate_regimens,
    evaluate_regimens,
    expected_regimen_count,
    pareto_front,
)
from ihd_solutions.risk import decisions_frame, simulate_events
from ihd_solutions.uncertainty import monte_carlo_regimen, rank_stability, tornado
from ihd_solutions.units import DAYS_PER_YEAR, require_positive, time_grid

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "solutions"

#: Explicitly synthetic guide used only to demonstrate sequence-level bystander counting.
EXAMPLE_SPACER = "GCAATACAGATCAACGATCG"
EXAMPLE_PAM = "TGG"
EXAMPLE_GUIDE_PROVENANCE = "synthetic_example_not_a_clinical_guide"

DEFAULT_POPULATIONS = (
    "secondary_post_mi",
    "post_mi_pci",
    "stable_ccs",
    "ascvd_ckd_inflammation",
    "hefh_primary",
    "ischemic_cm_heart_failure",
)


@dataclass
class PipelineConfig:
    output_dir: Path = DEFAULT_OUTPUT_DIR
    horizon_years: float = 5.0
    step_days: float = 30.0
    max_optional: int = 2
    populations: tuple[str, ...] = DEFAULT_POPULATIONS
    modes: tuple[CreditMode, ...] = (
        CreditMode.STRICT,
        CreditMode.SURROGATE,
        CreditMode.EXPLORATORY,
    )
    monte_carlo_draws: int = 120
    monte_carlo_top_n: int = 6
    seed: int = 20261106
    make_figures: bool = True
    plaque_months: int = 60
    extra_therapies: tuple[str, ...] = field(default=())

    def __post_init__(self) -> None:
        self.output_dir = Path(self.output_dir)
        require_positive("horizon_years", self.horizon_years)
        require_positive("step_days", self.step_days)
        require_positive("monte_carlo_draws", self.monte_carlo_draws)
        require_positive("monte_carlo_top_n", self.monte_carlo_top_n)
        require_positive("plaque_months", self.plaque_months)
        if self.max_optional < 0:
            raise ValueError("max_optional must be non-negative")
        if not self.populations:
            raise ValueError("at least one population is required")
        if not self.modes:
            raise ValueError("at least one credit mode is required")


def write_table(path: Path, frame: pd.DataFrame) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)
    return path


def write_json(path: Path, payload: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False, default=str) + "\n",
        encoding="utf-8",
    )
    temporary.replace(path)
    return path


def write_text(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)
    return path


def endpoint_specific_effects(anchors: dict[str, TrialAnchor]) -> pd.DataFrame:
    """Real benefits measured on endpoints that are not MACE composites."""
    rows = []
    for anchor in sorted(anchors.values(), key=lambda a: a.anchor_id):
        if anchor.endpoint_class is EndpointClass.MACE:
            continue
        rows.append(
            {
                "anchor_id": anchor.anchor_id,
                "therapy_id": anchor.therapy_id,
                "population_id": anchor.population_id,
                "endpoint": anchor.endpoint,
                "endpoint_class": anchor.endpoint_class.value,
                "hazard_ratio": anchor.hr if anchor.hr is not None else np.nan,
                "ci_low": anchor.ci_low if anchor.ci_low is not None else np.nan,
                "ci_high": anchor.ci_high if anchor.ci_high is not None else np.nan,
                "neutral_ci": anchor.is_neutral,
                "evidence_tier": anchor.evidence_tier.value,
                "credited_as_mace": False,
                "source_url": anchor.source_url,
                "notes": anchor.notes,
            }
        )
    return pd.DataFrame(rows)


def _figure_regimens(frame: pd.DataFrame, output_dir: Path, column: str) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    top = frame.head(12).iloc[::-1]
    figure, axis = plt.subplots(figsize=(11, 7), dpi=160)
    try:
        axis.barh(top["regimen"], top[column], color="#2563eb", alpha=0.85)
        axis.set_xlabel("Scenario events avoided per 1000 over the horizon")
        axis.set_title(
            "Ranked solution scenarios (evidence-credited only)\n"
            "Scenario modelling, not a patient prediction"
        )
        axis.grid(True, axis="x", alpha=0.3)
        figure.tight_layout()
        path = output_dir / "figures" / "regimen_ranking.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path)
        return path
    finally:
        plt.close(figure)


def _figure_trajectories(frame: pd.DataFrame, output_dir: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(13, 5), dpi=160)
    try:
        for label, group in frame.groupby("regimen", sort=False):
            axes[0].plot(group["year"], group["ldl_c_mg_dl"], label=label, linewidth=1.6)
            axes[1].plot(
                group["year"], group["cumulative_ldl_mmol_years"], label=label, linewidth=1.6
            )
        axes[0].axhline(55, color="crimson", linestyle="--", linewidth=1, label="55 mg/dL goal")
        axes[0].set_ylabel("LDL-C (mg/dL)")
        axes[0].set_xlabel("Year")
        axes[0].set_title("LDL-C trajectories under combined therapy")
        axes[1].set_ylabel("Cumulative LDL-C exposure (mmol/L·years)")
        axes[1].set_xlabel("Year")
        axes[1].set_title("Cumulative exposure drives the scenario benefit")
        for axis in axes:
            axis.grid(True, alpha=0.3)
        axes[0].legend(fontsize=7, loc="upper right")
        figure.tight_layout()
        path = output_dir / "figures" / "lipid_trajectories.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path)
        return path
    finally:
        plt.close(figure)


def _figure_uncertainty(summary: pd.DataFrame, output_dir: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    data = summary.head(10).iloc[::-1]
    median = data["events_avoided_median"].to_numpy()
    low = data["events_avoided_p5"].to_numpy()
    high = data["events_avoided_p95"].to_numpy()
    figure, axis = plt.subplots(figsize=(11, 6), dpi=160)
    try:
        positions = np.arange(len(data))
        axis.errorbar(
            median,
            positions,
            xerr=[median - low, high - median],
            fmt="o",
            color="#059669",
            ecolor="#6b7280",
            capsize=4,
        )
        axis.set_yticks(positions, data["regimen"], fontsize=8)
        axis.axvline(0, color="crimson", linestyle="--", linewidth=1)
        axis.set_xlabel("Events avoided per 1000 (median and 90% interval)")
        axis.set_title("Parameter uncertainty across efficacy, adherence and baseline risk")
        axis.grid(True, axis="x", alpha=0.3)
        figure.tight_layout()
        path = output_dir / "figures" / "uncertainty_intervals.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path)
        return path
    finally:
        plt.close(figure)


def _figure_tornado(frame: pd.DataFrame, output_dir: Path) -> Path:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(10, 5), dpi=160)
    try:
        data = frame.iloc[::-1]
        base = float(data["base_output"].iloc[0])
        axis.barh(
            data["parameter"],
            data["low_output"] - base,
            color="#f59e0b",
            alpha=0.9,
            label="lower input",
        )
        axis.barh(
            data["parameter"],
            data["high_output"] - base,
            color="#2563eb",
            alpha=0.7,
            label="higher input",
        )
        axis.axvline(0, color="black", linewidth=1)
        axis.set_xlabel("Change in events avoided per 1000 versus base case")
        axis.set_title("One-at-a-time sensitivity")
        axis.grid(True, axis="x", alpha=0.3)
        axis.legend(fontsize=8)
        figure.tight_layout()
        path = output_dir / "figures" / "sensitivity_tornado.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path)
        return path
    finally:
        plt.close(figure)


def run(config: PipelineConfig | None = None) -> dict[str, Any]:
    """Execute the full solution-simulation suite and return a run summary."""
    config = config or PipelineConfig()
    catalog = load_catalog()
    anchors = load_anchors()
    populations = load_populations()
    unknown = [pid for pid in config.populations if pid not in populations]
    if unknown:
        raise ValueError(f"unknown populations: {unknown}")
    days = time_grid(config.horizon_years, config.step_days)
    spec = EnumerationSpec(max_optional=config.max_optional)
    regimens = enumerate_regimens(catalog, spec)
    if len(regimens) != expected_regimen_count(catalog, spec):
        raise AssertionError("regimen enumeration disagrees with the combinatorial count")

    tables = config.output_dir / "tables"
    reports = config.output_dir / "reports"
    outputs: list[Path] = [
        write_table(tables / "therapy_catalog_resolved.csv", catalog_frame(catalog)),
        write_table(tables / "evidence_tier_summary.csv", summarise_tiers(anchors)),
        write_table(tables / "endpoint_specific_effects.csv", endpoint_specific_effects(anchors)),
    ]

    ranking_frames: list[pd.DataFrame] = []
    all_simulations = []
    for population_id in config.populations:
        population = populations[population_id]
        frame, simulations = evaluate_regimens(
            regimens, catalog, population, anchors, days, config.modes
        )
        ranking_frames.append(frame)
        all_simulations.extend(simulations)
    ranking = pd.concat(ranking_frames, ignore_index=True)
    outputs.append(write_table(tables / "regimen_ranking.csv", ranking))
    outputs.append(write_table(tables / "credit_decisions.csv", decisions_frame(all_simulations)))

    primary_mode = config.modes[0]
    primary_population = config.populations[0]
    column = f"{primary_mode.value}_events_avoided_per_1000"
    primary_frame = ranking.loc[ranking["population_id"].eq(primary_population)].reset_index(
        drop=True
    )
    outputs.append(write_table(tables / "pareto_front.csv", pareto_front(primary_frame, column)))

    top_labels = primary_frame.head(5)["regimen"].tolist()
    label_to_regimen = {regimen.label: regimen for regimen in regimens}
    top_regimens = [label_to_regimen[label] for label in top_labels]
    population = populations[primary_population]
    trajectory_rows = []
    curve_rows = []
    for regimen in top_regimens:
        simulation = simulate_events(regimen, catalog, population, anchors, days, mode=primary_mode)
        trajectory_rows.append(
            trajectory_frame(simulation.biomarkers, regimen.label, primary_population)
        )
        curve_rows.append(
            pd.DataFrame(
                {
                    "regimen": regimen.label,
                    "population_id": primary_population,
                    "mode": primary_mode.value,
                    "year": simulation.years,
                    "relative_risk": simulation.relative_risk,
                    "scenario_cumulative_incidence_pct": 100.0 * simulation.cumulative_incidence,
                    "baseline_cumulative_incidence_pct": 100.0 * simulation.baseline_incidence,
                }
            )
        )
    trajectories = pd.concat(trajectory_rows, ignore_index=True)
    outputs.append(write_table(tables / "lipid_trajectories.csv", trajectories))
    outputs.append(
        write_table(tables / "risk_curves.csv", pd.concat(curve_rows, ignore_index=True))
    )

    mc_results = [
        monte_carlo_regimen(
            regimen,
            catalog,
            population,
            anchors,
            days,
            mode=primary_mode,
            draws=config.monte_carlo_draws,
            seed=config.seed + index,
        )
        for index, regimen in enumerate(top_regimens[: config.monte_carlo_top_n])
    ]
    mc_summary = pd.DataFrame([result.summary() for result in mc_results]).sort_values(
        "events_avoided_median", ascending=False, kind="stable"
    )
    outputs.append(write_table(tables / "monte_carlo_summary.csv", mc_summary))
    outputs.append(write_table(tables / "rank_stability.csv", rank_stability(mc_results)))

    reference = (
        top_regimens[0] if top_regimens else validate_combination(catalog, background_regimen_ids())
    )

    def adherence_model(taking_fraction: float) -> float:
        from ihd_solutions.pharmacology import Adherence

        adherence = {
            therapy_id: Adherence(taking_fraction=taking_fraction)
            for therapy_id in reference.therapy_ids
        }
        return simulate_events(
            reference, catalog, population, anchors, days, primary_mode, adherence=adherence
        ).events_avoided_per_1000

    def horizon_model(years: float) -> float:
        grid = time_grid(years, config.step_days)
        return simulate_events(
            reference, catalog, population, anchors, grid, primary_mode
        ).events_avoided_per_1000

    def baseline_model(rate: float) -> float:
        from dataclasses import replace

        altered = replace(
            population,
            baseline_annual_mace_rate=rate,
            baseline_rate_ci_low=min(population.baseline_rate_ci_low, rate),
            baseline_rate_ci_high=max(population.baseline_rate_ci_high, rate),
        )
        return simulate_events(
            reference, catalog, altered, anchors, days, primary_mode
        ).events_avoided_per_1000

    sensitivity = tornado(adherence_model, {"adherence_fraction": (1.0, 0.5, 1.0)})
    sensitivity = pd.concat(
        [
            sensitivity,
            tornado(horizon_model, {"horizon_years": (config.horizon_years, 1.0, 10.0)}),
            tornado(
                baseline_model,
                {
                    "baseline_annual_rate": (
                        population.baseline_annual_mace_rate,
                        population.baseline_rate_ci_low,
                        population.baseline_rate_ci_high,
                    )
                },
            ),
        ],
        ignore_index=True,
    ).sort_values("span", ascending=False, kind="stable")
    outputs.append(write_table(tables / "sensitivity_tornado.csv", sensitivity))

    outputs.append(
        write_table(
            tables / "plaque_solution_grid.csv",
            _plaque_grid(catalog, populations, primary_population, config),
        )
    )
    outputs.append(
        write_table(tables / "inflammation_solution_grid.csv", _inflammation_grid(catalog))
    )
    outputs.append(write_table(tables / "editing_scenarios.csv", _editing_grid()))
    outputs.append(write_json(tables / "bystander_report.json", _bystander_payload()))
    outputs.append(write_json(tables / "off_target_status.json", editing.off_target_assessment()))
    outputs.append(
        write_table(
            tables / "regeneration_scenarios.csv",
            regeneration.simulate_regeneration(
                draws=max(200, config.monte_carlo_draws), seed=config.seed
            ),
        )
    )
    outputs.append(
        write_json(
            tables / "regeneration_surveillance.json", regeneration.surveillance_requirements()
        )
    )

    # ---- cure search: criteria scoring, combination cover, lifetime exposure ----
    cure_assessments = cure.assess_catalog(catalog, anchors, population, days)
    outputs.append(write_table(tables / "cure_criteria_scores.csv", cure_assessments))
    cover = cure.search_cure_combination(cure_assessments)
    outputs.append(write_json(tables / "cure_combination_search.json", _cover_payload(cover)))
    relaxed = cure.search_cure_combination(cure_assessments, threshold=0.5)
    outputs.append(
        write_json(tables / "cure_combination_search_relaxed.json", _cover_payload(relaxed))
    )
    lifetime = cure.lifetime_grid()
    outputs.append(write_table(tables / "lifetime_exposure_scenarios.csv", lifetime))
    outputs.append(
        write_table(tables / "regression_calibration.csv", regression.calibration_frame())
    )
    outputs.append(
        write_json(tables / "regression_calibration_fit.json", regression.calibration_summary())
    )
    outputs.append(
        write_table(
            tables / "plaque_imaging_evidence.csv", regression.imaging_anchor_frame(anchors)
        )
    )
    regression_addons = tuple(
        sorted(
            therapy.therapy_id
            for therapy in catalog.values()
            if therapy.ldl_pct > 0
            and therapy.therapy_id not in background_regimen_ids()
            and therapy.exclusive_group
            not in {catalog[tid].exclusive_group for tid in background_regimen_ids()}
        )
    )
    regression_grid = regression.regression_grid(
        catalog,
        populations[primary_population],
        background_regimen_ids(),
        regression_addons,
        years=2.0,
    )
    outputs.append(write_table(tables / "plaque_regression_grid.csv", regression_grid))

    residuals = []
    for population_id in config.populations:
        frame = ranking.loc[ranking["population_id"].eq(population_id)]
        if frame.empty:
            continue
        best_label = frame.sort_values(column, ascending=False).iloc[0]["regimen"]
        best_regimen = label_to_regimen[best_label]
        sim = simulate_events(
            best_regimen,
            catalog,
            populations[population_id],
            anchors,
            days,
            mode=primary_mode,
        )
        residuals.append(cure.residual_risk(sim, populations[population_id]))
    outputs.append(
        write_table(
            tables / "residual_risk.csv",
            pd.DataFrame([item.__dict__ for item in residuals]),
        )
    )
    outputs.append(
        write_text(
            reports / "cure_search_report.txt",
            cure.cure_search_report(cure_assessments, cover, lifetime, residuals),
        )
    )

    if config.make_figures:
        outputs.extend(
            [
                _figure_regimens(primary_frame, config.output_dir, column),
                _figure_trajectories(trajectories, config.output_dir),
                _figure_uncertainty(mc_summary, config.output_dir),
                _figure_tornado(sensitivity, config.output_dir),
            ]
        )

    report = _report_text(
        config, catalog, anchors, regimens, ranking, primary_frame, mc_summary, column
    )
    outputs.append(write_text(reports / "solutions_simulation_report.txt", report))
    manifest = _manifest(config, outputs, regimens, ranking)
    write_json(config.output_dir / "run_manifest.json", manifest)
    return {
        "status": "passed",
        "regimens": len(regimens),
        "populations": len(config.populations),
        "modes": [mode.value for mode in config.modes],
        "outputs": len(outputs),
        "output_dir": str(config.output_dir),
    }


def _plaque_grid(
    catalog: dict[str, Therapy],
    populations: dict[str, Population],
    population_id: str,
    config: PipelineConfig,
) -> pd.DataFrame:
    population = populations[population_id]
    rows = []
    extras = [
        "methotrexate_low_dose",  # inactive comparator: no LDL and no MMP effect
        None,
        *sorted(
            therapy.therapy_id
            for therapy in catalog.values()
            if therapy.ldl_pct > 0 or therapy.mmp_factor < 1.0
        ),
    ]
    for extra in extras:
        if extra == "methotrexate_low_dose":
            ids = (extra,)
            label = "inactive_reference"
        else:
            ids = background_regimen_ids() + ((extra,) if extra else ())
            label = "background_only" if extra is None else extra
        try:
            regimen = validate_combination(catalog, ids)
        except ValueError:
            continue
        outcome = simulate_plaque(regimen, catalog, population, months=config.plaque_months)
        rows.append(
            {
                "regimen": regimen.label,
                "population_id": population_id,
                "added_therapy": label,
                "first_threshold_crossing_month": (
                    outcome.first_crossing_month if outcome.first_crossing_month is not None else -1
                ),
                "crossed_within_horizon": outcome.crossed,
                "intact_months": outcome.intact_months,
                "min_cap_thickness_um": outcome.min_cap_um,
                "max_hoop_stress_kpa": outcome.max_stress_kpa,
                "final_core_area_mm2": outcome.final_core_area_mm2,
                "interpretation": "thin-wall proxy thresholds; not a rupture prediction",
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["crossed_within_horizon", "intact_months"], ascending=[True, False], kind="stable"
    )


def _inflammation_grid(catalog: dict[str, Therapy]) -> pd.DataFrame:
    rows = []
    interventions = [inflammation.UNTREATED] + [
        inflammation.intervention_for(therapy)
        for therapy in sorted(catalog.values(), key=lambda t: t.therapy_id)
        if therapy.therapy_class
        in {TherapyClass.ANTI_INFLAMMATORY, TherapyClass.LIPID_LOWERING, TherapyClass.METABOLIC}
    ]
    for intervention in interventions:
        frame = inflammation.simulate_inflammation(intervention, days=180.0, step_days=0.5)
        final = frame.iloc[-1]
        equilibrium = inflammation.steady_state(intervention=intervention)
        rows.append(
            {
                "therapy_id": intervention.therapy_id,
                "hscrp_day180_mg_l": float(final["hsCRP"]),
                "il6_index_day180": float(final["IL6"]),
                "mmp_index_day180": float(final["MMP"]),
                "analytic_equilibrium_hscrp": float(equilibrium[4]),
                "distance_to_equilibrium": float(abs(final["hsCRP"] - equilibrium[4])),
                "event_credit": 0.0,
                "interpretation": "biomarker scenario only; hsCRP is not a validated surrogate",
            }
        )
    return pd.DataFrame(rows)


def _editing_grid() -> pd.DataFrame:
    scenarios = [
        editing.EditingScenario("permanent_high_efficiency", edited_hepatocyte_fraction=0.95),
        editing.EditingScenario("permanent_moderate_efficiency", edited_hepatocyte_fraction=0.70),
        editing.EditingScenario(
            "waning_5_year_half_life",
            edited_hepatocyte_fraction=0.95,
            durability_half_life_years=5.0,
        ),
        editing.EditingScenario(
            "narrow_window_editor",
            editor_version="ABE9",
            edited_hepatocyte_fraction=0.85,
        ),
    ]
    frames = [
        editing.simulate_editing(scenario, horizon_years=10.0, step_days=60.0)
        for scenario in scenarios
    ]
    combined = pd.concat(frames, ignore_index=True)
    rows = []
    for scenario_id, group in combined.groupby("scenario_id", sort=True):
        years = group["year"].to_numpy()
        reduction = group["ldl_reduction_pct"].to_numpy() / 100.0
        mean_reduction = float(np.trapezoid(reduction, years) / (years[-1] - years[0]))
        rows.append(
            {
                "scenario_id": scenario_id,
                "editor_version": group["editor_version"].iloc[0],
                "durability_assumption": group["durability_assumption"].iloc[0],
                "plateau_pcsk9_knockdown_pct": float(group["pcsk9_knockdown_pct"].max()),
                "plateau_ldl_reduction_pct": float(group["ldl_reduction_pct"].max()),
                "ten_year_mean_ldl_reduction_pct": 100.0 * mean_reduction,
                "final_ldl_c_mg_dl": float(group["ldl_c_mg_dl"].iloc[-1]),
                "clinical_safety_assessed": False,
                "mace_credit": 0.0,
            }
        )
    return pd.DataFrame(rows)


def _cover_payload(cover: cure.CoverResult) -> dict[str, Any]:
    return {
        "threshold": cover.threshold,
        "satisfied": cover.satisfied,
        "selected": list(cover.selected),
        "coverage": cover.coverage,
        "missing": list(cover.missing),
        "ceiling": cover.ceiling,
        "binding_constraints": list(cover.binding_constraints),
        "notes": list(cover.notes),
        "interpretation": (
            "Criteria combine by maximum, never by summation. A binding ceiling means no "
            "combination of catalogued interventions can satisfy that criterion."
        ),
    }


def _bystander_payload() -> dict[str, Any]:
    report = editing.bystander_report(EXAMPLE_SPACER, EXAMPLE_PAM, intended_position=6)
    report["provenance"] = EXAMPLE_GUIDE_PROVENANCE
    report["warning"] = (
        "This spacer is a synthetic teaching example. It is not the VERVE-102 guide and must "
        "not be used for any experimental or clinical purpose."
    )
    return report


def _report_text(
    config: PipelineConfig,
    catalog: dict[str, Therapy],
    anchors: dict[str, TrialAnchor],
    regimens: list[Regimen],
    ranking: pd.DataFrame,
    primary_frame: pd.DataFrame,
    mc_summary: pd.DataFrame,
    column: str,
) -> str:
    background_label = validate_combination(catalog, background_regimen_ids()).label
    lines = [
        "Ischaemic heart disease: solution simulation report",
        "=" * 68,
        "SCOPE: scenario modelling of candidate interventions. Not a clinical trial,",
        "not a patient-level prediction, not a treatment recommendation, not a cure.",
        "",
        f"Therapies in catalogue: {len(catalog)}",
        f"Trial anchors: {len(anchors)}",
        f"Admissible regimens per population: {len(regimens)}",
        f"Populations: {', '.join(config.populations)}",
        f"Evidence modes: {', '.join(mode.value for mode in config.modes)}",
        f"Horizon: {config.horizon_years} years; grid step {config.step_days} days",
        f"Monte Carlo draws per regimen: {config.monte_carlo_draws}; seed {config.seed}",
        "",
        "CREDITING RULES",
        "- strict: only randomised hard-endpoint anchors matched or declared transportable",
        "- surrogate: additionally transports the LDL/pressure relationships to unproven agents",
        "- exploratory: additionally admits observational and non-MACE endpoint evidence",
        "- hsCRP and Lp(a) changes never create event credit in any mode",
        "- heart-failure, functional and mortality endpoints are never relabelled as MACE",
        "- neutral published intervals are held at HR = 1.0 (RiskParameters.neutral_ci_policy)",
        "",
        "BEST REGIMEN PER POPULATION (headline mode)",
    ]
    for population_id, group in ranking.groupby("population_id", sort=False):
        ordered = group.sort_values(column, ascending=False)
        best = ordered.iloc[0]
        background = ordered[ordered["regimen"] == background_label]
        reference = f"{background.iloc[0][column]:.1f}" if not background.empty else "n/a"
        lines.append(
            f"  {population_id}: background {reference} -> {best[column]:.1f} "
            f"events avoided/1000 with {best['regimen']}"
        )
    lines.extend(
        [
            "",
            f"TOP REGIMENS ({config.populations[0]}, {config.modes[0].value} mode)",
        ]
    )
    for row in primary_frame.head(10).itertuples(index=False):
        mapping = row._asdict()
        lines.append(
            f"  {mapping['regimen']}: {mapping[column]:.1f} events avoided/1000; "
            f"burden {mapping['burden_score']:.0f}"
        )
    lines.extend(["", "UNCERTAINTY (median and 90% interval, events avoided per 1000)"])
    for row in mc_summary.itertuples(index=False):
        lines.append(
            f"  {row.regimen}: {row.events_avoided_median:.1f} "
            f"[{row.events_avoided_p5:.1f}, {row.events_avoided_p95:.1f}]; "
            f"P(benefit)={row.probability_of_benefit:.2f}"
        )
    lines.extend(
        [
            "",
            "WHAT THIS SUITE DOES NOT ESTABLISH",
            "- No plaque-rupture prediction for an individual artery",
            "- No base-editor safety conclusion; dedicated off-target assays are required",
            "- No remuscularisation or mortality claim for regenerative cell sheets",
            "- No claim that biomarker improvement implies fewer clinical events",
            f"- Total simulated years per curve: {config.horizon_years}; "
            f"days per year convention: {DAYS_PER_YEAR}",
        ]
    )
    return "\n".join(lines) + "\n"


def _manifest(
    config: PipelineConfig,
    outputs: list[Path],
    regimens: list[Regimen],
    ranking: pd.DataFrame,
) -> dict[str, Any]:
    import matplotlib
    import scipy

    records = []
    for path in outputs:
        payload = path.read_bytes()
        records.append(
            {
                "path": str(path.relative_to(config.output_dir)),
                "bytes": len(payload),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        )
    return {
        "status": "passed",
        "scope": "scenario simulation; not clinical prediction",
        "config": {
            "horizon_years": config.horizon_years,
            "step_days": config.step_days,
            "max_optional": config.max_optional,
            "populations": list(config.populations),
            "modes": [mode.value for mode in config.modes],
            "monte_carlo_draws": config.monte_carlo_draws,
            "seed": config.seed,
        },
        "counts": {
            "regimens_per_population": len(regimens),
            "ranking_rows": int(len(ranking)),
            "outputs": len(records),
        },
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "outputs": records,
    }
