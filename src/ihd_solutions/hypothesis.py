"""REVERT-1: an invented, explicitly hypothetical therapeutic architecture.

NOT A THERAPY. NOT EFFICACY EVIDENCE. Nothing in this module is a clinical result,
and nothing here can enter the evidence-graded simulator: hypotheses live outside
`therapy_catalog.csv` and `trial_anchors.csv` entirely, write to their own output
directory, and carry a disclaimer on every row.

The design hypothesis
---------------------
Every published attempt at *reversing* (not merely slowing) atherosclerosis
attacks one node. The literature says the nodes are coupled, and that attacking
one alone can backfire:

* Engulfing an apoptotic cell loads the macrophage with cholesterol. Without
  ABCA1/ABCG1 efflux capacity the efferocyte dies, enlarging the necrotic core -
  so boosting efferocytosis alone can make lesions worse.
* LXR agonism raises efflux but causes hepatic lipogenesis unless it is confined
  to lesional macrophages.
* Dissolving crystalline cholesterol generates oxysterols, the endogenous LXR
  ligand, which raises both ABCA1/ABCG1 *and* MerTK.

REVERT-1 proposes one lesion-gated chassis carrying the whole loop, so that the
plaque's own crystalline cholesterol becomes the signal that drives its removal.
The novelty is the *coupling and the gate*, not the individual parts: each module
has a published precedent, all of them preclinical.

This module models the coupling to ask a falsifiable question: does co-delivery
produce clearance that no single module achieves, and what would the effect sizes
have to be to satisfy the cure criteria in :mod:`ihd_solutions.cure`?
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

import numpy as np
import pandas as pd

from ihd_solutions.units import require_fraction, require_nonnegative, require_positive

DISCLAIMER = (
    "INVENTED HYPOTHESIS - not a therapy, not a trial result, not efficacy evidence. "
    "All parameters are assumptions with stated ranges. Excluded from all "
    "evidence-graded outputs."
)

#: Marker written into every hypothesis output row.
HYPOTHESIS_TAG = "invented_hypothesis_unvalidated"


class PrecedentLevel(str, Enum):
    """How much of a module has actually been demonstrated, and in what."""

    HUMAN_MECHANISM = "human_mechanism"
    MOUSE_GENETIC = "mouse_genetic"
    MOUSE_PHARMACOLOGIC = "mouse_pharmacologic"
    IN_VITRO = "in_vitro"
    ENGINEERING_ONLY = "engineering_only"


@dataclass(frozen=True)
class ModuleSpec:
    """One payload module and the published precedent it is borrowed from."""

    module_id: str
    target: str
    role: str
    precedent: PrecedentLevel
    precedent_reference: str
    source_url: str
    known_failure_mode: str

    def __post_init__(self) -> None:
        if not self.module_id or not self.source_url.startswith("https://"):
            raise ValueError("module_id and an https source_url are required")


#: The six payload modules plus the two systemic arms. Every precedent is real;
#: the combination is the invention.
REVERT1_MODULES: tuple[ModuleSpec, ...] = (
    ModuleSpec(
        "M0_gate",
        "Lesional macrophage AND oxidised-phospholipid dual recognition",
        "Confines payload release to plaque macrophages (the AND gate)",
        PrecedentLevel.MOUSE_PHARMACOLOGIC,
        "Macrophage-targeted nanoparticles delivering an LXR agonist avoided hepatic "
        "lipogenesis while raising ABCA1/ABCG1",
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC7530139/",
        "Off-target release: systemic CD47 blockade causes anaemia; systemic LXR "
        "agonism causes hepatic steatosis",
    ),
    ModuleSpec(
        "M1_solubilise",
        "Crystalline cholesterol",
        "Dissolves crystals and generates the endogenous LXR ligand that fuels M2/M3",
        PrecedentLevel.MOUSE_PHARMACOLOGIC,
        "Cyclodextrin dissolved cholesterol crystals, raised oxysterols and drove "
        "LXR-dependent regression of established plaques",
        "https://www.semanticscholar.org/paper/Cyclodextrin-promotes-atherosclerosis-regression-Zimmer-Grebe/4c3bf0237601aa3dc9756240b1aa5c1570868e06",
        "A separate Apoe-null study failed to reproduce lesion regression; "
        "ototoxicity limits systemic cyclodextrin dosing",
    ),
    ModuleSpec(
        "M2_efflux",
        "LXR-beta / ABCA1 / ABCG1",
        "Exports the liberated cholesterol so the efferocyte survives its own meal",
        PrecedentLevel.MOUSE_GENETIC,
        "Abca1/Abcg1-deficient efferocytes die after engulfing apoptotic cells; "
        "efflux preserves viability during efferocytosis",
        "https://www.ahajournals.org/doi/10.1161/circresaha.110.217281",
        "Hepatic LXR activation causes lipogenesis and hypertriglyceridaemia",
    ),
    ModuleSpec(
        "M3_efferocytosis",
        "MerTK cleavage / CD47-SIRPalpha",
        "Restores clearance of the apoptotic cells that form the necrotic core",
        PrecedentLevel.MOUSE_GENETIC,
        "Cleavage-resistant MerTK gave smaller necrotic cores and thicker caps; "
        "CD47 blockade restored phagocytosis and regressed established lesions",
        "https://www.nature.com/articles/nrcardio.2016.125",
        "Systemic CD47 blockade destroys red cells; MerTK signalling is profibrotic "
        "in liver, so the gate must hold",
    ),
    ModuleSpec(
        "M4_influx",
        "SR-A / CD36 scavenger uptake",
        "Stops re-loading the cell while it is being unloaded",
        PrecedentLevel.MOUSE_PHARMACOLOGIC,
        "Two-pronged nanoparticle (SR-A siRNA plus LXR ligand) reduced plaque free "
        "and esterified cholesterol more than either payload alone",
        "https://www.sciencedirect.com/science/article/abs/pii/S0142961220305792",
        "Scavenger receptors also serve host defence; chronic knockdown may impair "
        "bacterial clearance",
    ),
    ModuleSpec(
        "M5_resolution",
        "NLRP3 / ADAM17",
        "Breaks the inflammation-to-MerTK-cleavage feed-forward loop",
        PrecedentLevel.MOUSE_GENETIC,
        "ADAM17-mediated MerTK cleavage impairs efferocytosis and resolution; "
        "inflammasome activity drives the cleavage",
        "https://onlinelibrary.wiley.com/doi/full/10.1002/mco2.70764",
        "Infection risk with sustained inflammasome blockade (see ZEUS serious "
        "infections in the human IL-6 data)",
    ),
    ModuleSpec(
        "M6_mineral",
        "Microcalcification nucleation",
        "Stops new cap micro-deposits while existing mineral matures to stable macro",
        PrecedentLevel.HUMAN_MECHANISM,
        "Microcalcification tracks instability while dense macrocalcification is "
        "stabilising; statins reduce the microcalcified fraction as density rises",
        "https://doi.org/10.3390/cells15161496",
        "Removing mature calcium could destabilise the cap - the design must NOT "
        "attempt total decalcification",
    ),
    ModuleSpec(
        "S1_driver",
        "Hepatic PCSK9 / ANGPTL3 / LPA",
        "One-time editing shuts off the substrate supply for life",
        PrecedentLevel.HUMAN_MECHANISM,
        "VERVE-102 gave durable 62% LDL-C reduction to 18 months; CTX310 gave 53% "
        "LDL-C reduction at one year",
        "https://www.nejm.org/doi/full/10.1056/NEJMoa2601283",
        "No outcomes data; base-editor off-targets need dedicated assays",
    ),
    ModuleSpec(
        "S2_tolerance",
        "ApoB / oxLDL-specific regulatory T cells",
        "Prevents recurrence without chronic immunosuppression",
        PrecedentLevel.MOUSE_PHARMACOLOGIC,
        "CAR-Tregs against oxidised LDL reduced foam-cell formation and limited "
        "plaque development in mice",
        "https://academic.oup.com/cardiovascres/advance-article/doi/10.1093/cvr/cvag020/8471361",
        "Foxp3 instability could convert tolerance into pathogenic effector cells",
    ),
)


def module_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "module_id": spec.module_id,
                "target": spec.target,
                "role": spec.role,
                "precedent_level": spec.precedent.value,
                "precedent": spec.precedent_reference,
                "known_failure_mode": spec.known_failure_mode,
                "source_url": spec.source_url,
                "status": HYPOTHESIS_TAG,
            }
            for spec in REVERT1_MODULES
        ]
    )


@dataclass(frozen=True)
class ClearanceParameters:
    """Assumed module strengths. Every value is a hypothesis, not a measurement."""

    gate_specificity: float = 0.90
    crystal_solubilisation: float = 0.60
    efflux_gain: float = 0.80
    efferocytosis_gain: float = 0.70
    influx_block: float = 0.50
    resolution: float = 0.60
    nucleation_block: float = 0.60
    plasma_driver_reduction: float = 0.70

    def __post_init__(self) -> None:
        for name in (
            "gate_specificity",
            "crystal_solubilisation",
            "efflux_gain",
            "efferocytosis_gain",
            "influx_block",
            "resolution",
            "nucleation_block",
            "plasma_driver_reduction",
        ):
            require_fraction(name, getattr(self, name), allow_one=True)

    def delivered(self) -> ClearanceParameters:
        """Scale payload strength by how much of it actually reaches the lesion.

        The gate is not a safety detail - it is what makes the efficacious dose
        reachable at all, because the systemic dose is capped by anaemia and
        hepatic steatosis.
        """
        factor = self.gate_specificity
        return ClearanceParameters(
            gate_specificity=1.0,
            crystal_solubilisation=self.crystal_solubilisation * factor,
            efflux_gain=self.efflux_gain * factor,
            efferocytosis_gain=self.efferocytosis_gain * factor,
            influx_block=self.influx_block * factor,
            resolution=self.resolution * factor,
            nucleation_block=self.nucleation_block * factor,
            plasma_driver_reduction=self.plasma_driver_reduction,
        )


#: Named arms for the ablation experiment. The point of the design is that the
#: whole is supposed to beat every part; this is where that claim can fail.
ARMS: dict[str, tuple[str, ...]] = {
    "untreated": (),
    "driver_only": ("plasma_driver_reduction",),
    "efflux_only": ("crystal_solubilisation", "efflux_gain"),
    "efferocytosis_only": ("efferocytosis_gain",),
    "efflux_plus_efferocytosis": ("crystal_solubilisation", "efflux_gain", "efferocytosis_gain"),
    "two_pronged_published": ("crystal_solubilisation", "efflux_gain", "influx_block"),
    "revert1_full": (
        "plasma_driver_reduction",
        "crystal_solubilisation",
        "efflux_gain",
        "efferocytosis_gain",
        "influx_block",
        "resolution",
        "nucleation_block",
    ),
    "revert1_ungated": (
        "plasma_driver_reduction",
        "crystal_solubilisation",
        "efflux_gain",
        "efferocytosis_gain",
        "influx_block",
        "resolution",
        "nucleation_block",
    ),
}

_ZERO_STRENGTHS = (
    "crystal_solubilisation",
    "efflux_gain",
    "efferocytosis_gain",
    "influx_block",
    "resolution",
    "nucleation_block",
    "plasma_driver_reduction",
)


def arm_parameters(
    arm: str, base: ClearanceParameters = ClearanceParameters()
) -> ClearanceParameters:
    """Build one ablation arm by zeroing every module it does not contain."""
    if arm not in ARMS:
        raise ValueError(f"unknown arm {arm!r}; choose from {sorted(ARMS)}")
    active = ARMS[arm]
    values = {name: (getattr(base, name) if name in active else 0.0) for name in _ZERO_STRENGTHS}
    gate = 0.25 if arm == "revert1_ungated" else base.gate_specificity
    return ClearanceParameters(gate_specificity=gate, **values)


@dataclass(frozen=True)
class ClearanceResult:
    """Trajectory of the rupture-prone substrate under one arm."""

    arm: str
    months: np.ndarray
    free_cholesterol: np.ndarray
    necrotic_debris: np.ndarray
    microcalcification: np.ndarray
    macrophage_viability: np.ndarray
    notes: tuple[str, ...] = field(default=(DISCLAIMER,))

    @property
    def substrate(self) -> np.ndarray:
        """Rupture-prone burden: lipid + uncleared debris + micro-mineral.

        Mature macrocalcification is deliberately excluded: it is stabilising, so
        counting it as disease to be removed would be a design error.
        """
        return self.free_cholesterol + self.necrotic_debris + self.microcalcification

    @property
    def substrate_removal_fraction(self) -> float:
        start = float(self.substrate[0])
        return 0.0 if start <= 0 else 1.0 - float(self.substrate[-1]) / start

    @property
    def necrotic_core_change_fraction(self) -> float:
        start = float(self.necrotic_debris[0])
        return 0.0 if start <= 0 else float(self.necrotic_debris[-1]) / start - 1.0

    @property
    def microcalcification_removal_fraction(self) -> float:
        start = float(self.microcalcification[0])
        return 0.0 if start <= 0 else 1.0 - float(self.microcalcification[-1]) / start

    @property
    def phagocyte_depletion(self) -> bool:
        """True when a quarter or more of the functional phagocyte pool is lost."""
        return bool(np.min(self.macrophage_viability) < 0.75)

    def necrotic_core_change_at(self, month: float) -> float:
        """Core change at an interim timepoint - short studies can be misled."""
        require_nonnegative("month", month)
        start = float(self.necrotic_debris[0])
        if start <= 0:
            return 0.0
        index = int(np.argmin(np.abs(self.months - month)))
        return float(self.necrotic_debris[index]) / start - 1.0


@dataclass(frozen=True)
class LesionState:
    """Starting lesion, in arbitrary-but-consistent burden units."""

    free_cholesterol: float = 60.0
    necrotic_debris: float = 40.0
    microcalcification: float = 8.0
    macrophage_viability: float = 1.0

    def __post_init__(self) -> None:
        for name in ("free_cholesterol", "necrotic_debris", "microcalcification"):
            require_positive(name, getattr(self, name))
        require_fraction("macrophage_viability", self.macrophage_viability, allow_one=True)
        require_positive("macrophage_viability", self.macrophage_viability)


def simulate_clearance(
    params: ClearanceParameters,
    arm: str = "custom",
    months: float = 24.0,
    step_months: float = 0.25,
    lesion: LesionState = LesionState(),
    source_reduction: float = 0.0,
) -> ClearanceResult:
    """Integrate the coupled clearance loop with explicit Euler steps.

    The coupling that matters: engulfing debris adds cholesterol to the cell;
    cholesterol above what efflux can export kills the phagocyte; dead phagocytes
    become debris. Boosting uptake without export therefore feeds the necrotic
    core instead of clearing it. That behaviour is emergent here, not asserted.
    """
    require_positive("months", months)
    require_positive("step_months", step_months)
    require_fraction("source_reduction", source_reduction, allow_one=True)
    if step_months > 0.5:
        raise ValueError("step_months must be <= 0.5 for integrator stability")

    effective = params.delivered()
    grid = np.arange(0.0, months + step_months * 0.5, step_months, dtype=float)
    free = np.empty(grid.size)
    debris = np.empty(grid.size)
    micro = np.empty(grid.size)
    viable = np.empty(grid.size)
    free[0], debris[0], micro[0], viable[0] = (
        lesion.free_cholesterol,
        lesion.necrotic_debris,
        lesion.microcalcification,
        lesion.macrophage_viability,
    )

    # Rate constants, per month. Assumed, and swept in the Monte Carlo.
    # Baseline is calibrated so an untreated lesion slowly progresses with a
    # viable phagocyte pool: defective clearance is the disease, not mass death.
    influx_rate = 3.0
    debris_production = 1.5
    export_capacity = 5.0
    uptake_capacity = 1.2  # defective efferocytosis at baseline
    cholesterol_per_corpse = 2.5
    load_tolerance = 70.0  # cholesterol a phagocyte can hold without dying
    overload_death = 0.02
    recovery_rate = 0.15
    phagocyte_pool = 12.0
    nucleation_rate = 0.020
    maturation_rate = 0.030

    for index in range(grid.size - 1):
        f, d, c, v = free[index], debris[index], micro[index], viable[index]

        # M1: crystals dissolved, producing the oxysterol signal that fuels M2/M3.
        # Solubilisation is not itself export - it unblocks export and yields a
        # small oxysterol exit flux.
        solubilised = effective.crystal_solubilisation * 0.25 * f
        signal = min(1.0, solubilised / 8.0)
        # M2: efflux capacity rises with the endogenous LXR signal.
        capacity = export_capacity * (1.0 + 2.0 * effective.efflux_gain * signal)
        efflux = capacity * f / (f + 20.0)
        # M4 + S1: less cholesterol arriving in the first place.
        influx = (
            influx_rate * (1.0 - effective.influx_block) * (1.0 - effective.plasma_driver_reduction)
        )
        # M3: uptake scales with how many phagocytes are still alive. M5 protects
        # the clearance machinery from inflammatory cleavage.
        uptake = (
            uptake_capacity
            * (1.0 + 2.0 * effective.efferocytosis_gain)
            * (1.0 + 0.5 * effective.resolution)
            * v
            * d
            / (d + 15.0)
        )

        # The coupling. Cholesterol above what the cell can hold or export kills
        # it; the dead cell deposits both its debris and its cholesterol, which
        # loads the survivors. Uptake without export therefore feeds the core.
        overload = max(0.0, f - load_tolerance)
        death_rate = overload_death * overload / 20.0 / (1.0 + effective.resolution)
        death_flux = death_rate * v * phagocyte_pool

        df = influx + cholesterol_per_corpse * (uptake + death_flux) - efflux - solubilised * 0.15
        # An optional source-reduction module (M7, ER-stress/CHOP) lowers the rate
        # at which new apoptotic cells appear. Defaults to 0 so every existing
        # artefact and ceiling is unchanged; see ihd_solutions.gap_solver.
        dd = debris_production * (1.0 - source_reduction) + death_flux - uptake
        dv = recovery_rate * (1.0 - v) - death_rate * v
        dc = nucleation_rate * d * (1.0 - effective.nucleation_block) - maturation_rate * c

        free[index + 1] = max(0.0, f + step_months * df)
        debris[index + 1] = max(0.0, d + step_months * dd)
        micro[index + 1] = max(0.0, c + step_months * dc)
        viable[index + 1] = float(min(1.0, max(0.0, v + step_months * dv)))

    notes = [
        DISCLAIMER,
        "coupled loop: uptake loads cholesterol, overload kills phagocytes, dead "
        "phagocytes become debris",
        "rupture substrate excludes mature macrocalcification, which is stabilising",
        "rate constants are assumed and swept in the Monte Carlo, not measured",
    ]
    if arm == "revert1_ungated":
        notes.append("ungated arm assumes only 25% of payload reaches the lesion")
    return ClearanceResult(arm, grid, free, debris, micro, viable, tuple(notes))


def ablation_study(
    base: ClearanceParameters = ClearanceParameters(),
    months: float = 24.0,
    lesion: LesionState = LesionState(),
) -> pd.DataFrame:
    """Run every arm so the co-delivery claim can actually fail."""
    rows = []
    for arm in ARMS:
        result = simulate_clearance(arm_parameters(arm, base), arm, months, lesion=lesion)
        rows.append(
            {
                "arm": arm,
                "modules": ";".join(ARMS[arm]) or "none",
                "substrate_removal_pct": 100.0 * result.substrate_removal_fraction,
                "necrotic_core_change_pct": 100.0 * result.necrotic_core_change_fraction,
                "necrotic_core_change_12m_pct": 100.0 * result.necrotic_core_change_at(12.0),
                "microcalcification_removal_pct": 100.0
                * result.microcalcification_removal_fraction,
                "final_macrophage_viability": float(result.macrophage_viability[-1]),
                "phagocyte_depletion": result.phagocyte_depletion,
                "status": HYPOTHESIS_TAG,
            }
        )
    return (
        pd.DataFrame(rows)
        .sort_values("substrate_removal_pct", ascending=False)
        .reset_index(drop=True)
    )


@dataclass(frozen=True)
class TargetProductProfile:
    """What any invented cure would have to achieve, back-calculated."""

    required_substrate_removal: float
    required_hazard_ratio: float
    required_lvef_gain_pp: float
    best_current_substrate_removal: float
    best_current_lvef_gain_pp: float
    notes: tuple[str, ...] = field(default=(DISCLAIMER,))

    @property
    def substrate_gap_factor(self) -> float:
        if self.best_current_substrate_removal <= 0:
            return math.inf
        return self.required_substrate_removal / self.best_current_substrate_removal

    @property
    def myocardial_gap_factor(self) -> float:
        if self.best_current_lvef_gain_pp <= 0:
            return math.inf
        return self.required_lvef_gain_pp / self.best_current_lvef_gain_pp


def required_hazard_ratio(
    baseline_risk: float,
    reference_risk: float,
    allowed_excess_remaining: float = 0.2,
) -> float:
    """Hazard ratio needed to leave at most ``allowed_excess_remaining`` excess risk.

    Constant-hazard arithmetic over the same horizon, matching the event model:
    the target risk is the reference plus the permitted share of the excess.
    """
    require_fraction("baseline_risk", baseline_risk)
    require_fraction("reference_risk", reference_risk)
    require_fraction("allowed_excess_remaining", allowed_excess_remaining, allow_one=True)
    if reference_risk >= baseline_risk:
        raise ValueError("reference risk must be below baseline risk")
    target = reference_risk + allowed_excess_remaining * (baseline_risk - reference_risk)
    return float(math.log1p(-target) / math.log1p(-baseline_risk))


def target_product_profile(
    baseline_risk: float = 0.1393,
    reference_risk: float = 0.0198,
    criterion_threshold: float = 0.8,
    best_current_substrate_removal: float = 0.431,
    best_current_lvef_gain_pp: float = 3.9,
    normal_lvef_pp: float = 50.0,
    advanced_hf_lvef_pp: float = 25.0,
) -> TargetProductProfile:
    """Turn "invent a cure" into a specification with numbers.

    Defaults are taken from this repository's own simulated results: 13.9%
    five-year risk in chronic coronary syndrome against a 2.0% low-risk reference,
    43.1% lipid-pool removal for the best regimen, and BioVAT-HF's +3.9 LVEF points.
    """
    require_fraction("criterion_threshold", criterion_threshold, allow_one=True)
    require_positive("best_current_substrate_removal", best_current_substrate_removal)
    require_positive("best_current_lvef_gain_pp", best_current_lvef_gain_pp)
    if advanced_hf_lvef_pp >= normal_lvef_pp:
        raise ValueError("advanced HF ejection fraction must be below normal")
    required_hr = required_hazard_ratio(
        baseline_risk, reference_risk, allowed_excess_remaining=1.0 - criterion_threshold
    )
    return TargetProductProfile(
        required_substrate_removal=float(criterion_threshold),
        required_hazard_ratio=required_hr,
        required_lvef_gain_pp=float(normal_lvef_pp - advanced_hf_lvef_pp),
        best_current_substrate_removal=float(best_current_substrate_removal),
        best_current_lvef_gain_pp=float(best_current_lvef_gain_pp),
        notes=(
            DISCLAIMER,
            "substrate target is a proposed refinement: remove the rupture-prone "
            "fraction, not the stabilising macrocalcification",
            "randomised proof cannot be back-calculated; it requires a trial",
        ),
    )


def monte_carlo_hypothesis(
    draws: int = 400,
    seed: int = 20261107,
    months: float = 24.0,
    required_removal: float = 0.8,
) -> pd.DataFrame:
    """Sweep the assumed module strengths over their full plausible range.

    Uniform [0, 1] on every module is deliberately generous: it asks how often the
    architecture could hit the target even if every module worked as well as it
    conceivably could. A low probability here is informative.
    """
    require_positive("draws", draws)
    require_fraction("required_removal", required_removal, allow_one=True)
    rng = np.random.default_rng(seed)
    names = (
        "gate_specificity",
        "crystal_solubilisation",
        "efflux_gain",
        "efferocytosis_gain",
        "influx_block",
        "resolution",
        "nucleation_block",
        "plasma_driver_reduction",
    )
    samples = rng.uniform(0.0, 1.0, size=(int(draws), len(names)))
    rows = []
    for row in samples:
        params = ClearanceParameters(**dict(zip(names, row, strict=True)))
        result = simulate_clearance(params, "revert1_full", months=months, step_months=0.5)
        rows.append(
            {
                **dict(zip(names, row, strict=True)),
                "substrate_removal": result.substrate_removal_fraction,
                "necrotic_core_change": result.necrotic_core_change_fraction,
                "phagocyte_depletion": result.phagocyte_depletion,
                "meets_substrate_target": result.substrate_removal_fraction >= required_removal,
                "status": HYPOTHESIS_TAG,
            }
        )
    return pd.DataFrame(rows)


def best_case_horizon_sweep(
    horizons: tuple[float, ...] = (12.0, 24.0, 36.0, 60.0, 120.0, 240.0),
    base: ClearanceParameters = ClearanceParameters(),
    required_removal: float = 0.8,
) -> pd.DataFrame:
    """Central assumptions against a deliberately impossible best case.

    The best case sets every module to 1.0 - perfect delivery, perfect payload.
    It bounds what the architecture could ever do, so a target above that bound
    is unreachable by tuning and needs a different design.
    """
    require_fraction("required_removal", required_removal, allow_one=True)
    perfect = ClearanceParameters(gate_specificity=1.0, **{n: 1.0 for n in _ZERO_STRENGTHS})
    rows = []
    for horizon in horizons:
        require_positive("horizon", horizon)
        for label, params in (("central_assumptions", base), ("perfect_modules", perfect)):
            result = simulate_clearance(params, label, months=horizon, step_months=0.5)
            rows.append(
                {
                    "horizon_months": horizon,
                    "scenario": label,
                    "substrate_removal_pct": 100.0 * result.substrate_removal_fraction,
                    "necrotic_core_change_pct": 100.0 * result.necrotic_core_change_fraction,
                    "meets_substrate_target": (
                        result.substrate_removal_fraction >= required_removal
                    ),
                    "status": HYPOTHESIS_TAG,
                }
            )
    return pd.DataFrame(rows)


def architecture_ceiling(
    required_removal: float = 0.8, max_months: float = 240.0
) -> dict[str, float | bool | None]:
    """Asymptotic clearance ceiling and the earliest time the target is met."""
    require_fraction("required_removal", required_removal, allow_one=True)
    require_positive("max_months", max_months)
    perfect = ClearanceParameters(gate_specificity=1.0, **{n: 1.0 for n in _ZERO_STRENGTHS})
    trajectory = simulate_clearance(perfect, "perfect_modules", months=max_months, step_months=0.5)
    fractions = 1.0 - trajectory.substrate / float(trajectory.substrate[0])
    hits = np.flatnonzero(fractions >= required_removal)
    return {
        "ceiling_removal_fraction": float(np.max(fractions)),
        "asymptotic_removal_fraction": float(fractions[-1]),
        "months_to_target": (None if hits.size == 0 else float(trajectory.months[int(hits[0])])),
        "target_reachable": bool(hits.size > 0),
        "eradication": bool(fractions[-1] >= 0.999),
        "required_removal": float(required_removal),
    }


def sensitivity_summary(draws: pd.DataFrame) -> pd.DataFrame:
    """Rank-order correlation of each module with substrate removal."""
    from ihd_solutions.uncertainty import spearman_sensitivity

    names = [
        name
        for name in draws.columns
        if name
        not in {
            "substrate_removal",
            "necrotic_core_change",
            "phagocyte_depletion",
            "meets_substrate_target",
            "status",
        }
    ]
    inputs = {name: draws[name].to_numpy() for name in names}
    coefficients = spearman_sensitivity(inputs, draws["substrate_removal"].to_numpy())
    frame = coefficients[["parameter", "spearman_rho", "p_value"]].rename(
        columns={"parameter": "module"}
    )
    frame["status"] = HYPOTHESIS_TAG
    return frame.sort_values("spearman_rho", key=np.abs, ascending=False).reset_index(drop=True)


@dataclass(frozen=True)
class KillCriterion:
    """An observation that would falsify the design, with the cheapest test first."""

    order: int
    hypothesis_tested: str
    falsifying_observation: str
    cheapest_experiment: str
    readout: str


KILL_CRITERIA: tuple[KillCriterion, ...] = (
    KillCriterion(
        1,
        "Coupling: uptake without export is counterproductive",
        "Boosting efferocytosis alone shrinks the necrotic core as well as co-delivery does",
        "Apoe-null or Ldlr-null mice, four arms: anti-CD47 alone, macrophage-targeted "
        "LXR agonist alone, both, control; 12 weeks on Western diet",
        "Necrotic core area and free/macrophage-associated TUNEL ratio on aortic root "
        "sections; the design is only worth building if the combination beats both "
        "single arms on core area",
    ),
    KillCriterion(
        2,
        "Fuel: dissolved crystal cholesterol really does drive the loop",
        "Crystal solubilisation adds no oxysterol signal or no ABCA1/MerTK induction "
        "in human plaque macrophages",
        "Ex vivo human carotid endarterectomy explants, cyclodextrin versus vehicle, 48 hours",
        "Oxysterol mass spectrometry plus ABCA1/ABCG1/MerTK expression; no induction "
        "means the self-fuelling claim is dead",
    ),
    KillCriterion(
        3,
        "Gate: lesion-restricted release is achievable at a useful dose",
        "Dual-recognition particles cannot exceed ~25% lesional delivery, or systemic "
        "exposure still causes anaemia or hepatic steatosis",
        "Biodistribution and toxicology in a large animal with labelled particles; "
        "haemoglobin, reticulocytes, liver triglyceride",
        "Lesion-to-liver ratio and dose-limiting toxicity; without the gate the "
        "efficacious dose is unreachable",
    ),
    KillCriterion(
        4,
        "Mineral: suppressing nucleation does not destabilise the cap",
        "Blocking microcalcification nucleation thins the cap or raises rupture "
        "features instead of shifting mineral to stable macro",
        "Same animals, OCT/NaF-PET plus histology for cap thickness and micro-deposit "
        "size distribution",
        "Cap thickness and micro versus macro mineral fraction; a thinner cap kills "
        "module M6 outright",
    ),
    KillCriterion(
        5,
        "Durability: tolerance prevents recurrence after withdrawal",
        "Lesions regrow at control rates once the chassis is withdrawn despite CAR-Treg transfer",
        "Withdrawal-and-rechallenge design in the mouse regression model",
        "Substrate burden 12 weeks after withdrawal; recurrence means this is a "
        "chronic therapy, not a cure",
    ),
)


def kill_criteria_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "order": item.order,
                "hypothesis_tested": item.hypothesis_tested,
                "falsifying_observation": item.falsifying_observation,
                "cheapest_experiment": item.cheapest_experiment,
                "readout": item.readout,
                "status": HYPOTHESIS_TAG,
            }
            for item in KILL_CRITERIA
        ]
    )


def design_report(
    ablation: pd.DataFrame,
    profile: TargetProductProfile,
    draws: pd.DataFrame,
    sensitivity: pd.DataFrame,
    ceiling: dict | None = None,
) -> str:
    """The design write-up, with the numbers that could kill it."""
    full = ablation.loc[ablation["arm"].eq("revert1_full")].iloc[0]
    eff_only = ablation.loc[ablation["arm"].eq("efflux_only")].iloc[0]
    effero_only = ablation.loc[ablation["arm"].eq("efferocytosis_only")].iloc[0]
    published = ablation.loc[ablation["arm"].eq("two_pronged_published")].iloc[0]
    probability = float(draws["meets_substrate_target"].mean())
    lines = [
        "REVERT-1: an invented architecture for clearing the rupture-prone substrate",
        "=" * 74,
        DISCLAIMER,
        "",
        "THE IDEA",
        "Do not add another lipid-lowering drug. Close a loop inside the lesion: dissolve",
        "the plaque's own cholesterol crystals, let the resulting oxysterols drive both",
        "cholesterol export and corpse clearance, keep the phagocyte alive while it eats,",
        "and gate the whole payload to the lesion so the dose is reachable.",
        "",
        "WHY COUPLING IS THE CLAIM",
        "Engulfing a dead cell loads the macrophage with cholesterol. Without export",
        "capacity the efferocyte dies and becomes necrotic core itself. So efferocytosis",
        "alone can enlarge the lesion, and export alone cannot clear existing debris.",
        "",
        "SIMULATED ABLATION (24 months, assumed parameters)",
    ]
    for row in ablation.itertuples(index=False):
        flag = " [PHAGOCYTE DEPLETION]" if row.phagocyte_depletion else ""
        lines.append(
            f"  {row.arm:28s} substrate removed {row.substrate_removal_pct:6.1f}%  "
            f"necrotic core {row.necrotic_core_change_pct:+7.1f}%{flag}"
        )
    lines.extend(
        [
            "",
            f"  full architecture vs best published-style two-pronged arm: "
            f"{full.substrate_removal_pct:.1f}% vs {published.substrate_removal_pct:.1f}%",
            f"  efferocytosis alone: necrotic core {effero_only.necrotic_core_change_pct:+.1f}%"
            f" (export alone {eff_only.necrotic_core_change_pct:+.1f}%)",
            "",
            "WHAT A CURE WOULD REQUIRE (back-calculated, not hoped for)",
            f"  rupture-substrate removal: {profile.required_substrate_removal:.0%} "
            f"(best current {profile.best_current_substrate_removal:.1%}; "
            f"gap x{profile.substrate_gap_factor:.1f})",
            f"  hazard ratio needed to normalise risk: {profile.required_hazard_ratio:.3f}",
            f"  LVEF gain for the scar arm: +{profile.required_lvef_gain_pp:.0f} points "
            f"(best current +{profile.best_current_lvef_gain_pp:.1f}; "
            f"gap x{profile.myocardial_gap_factor:.1f})",
            "",
            "MONTE CARLO OVER THE FULL ASSUMED RANGE",
            f"  draws: {len(draws)}; P(substrate removal >= "
            f"{profile.required_substrate_removal:.0%}) = {probability:.1%}",
            f"  median substrate removal: {draws['substrate_removal'].median():.1%}",
            f"  draws with phagocyte depletion: {draws['phagocyte_depletion'].mean():.1%}",
            "  modules most associated with clearance:",
        ]
    )
    for row in sensitivity.head(4).itertuples(index=False):
        lines.append(f"    {row.module}: rho {row.spearman_rho:+.2f}")
    if ceiling is not None:
        months_to = ceiling["months_to_target"]
        lines.extend(
            [
                "",
                "THE CEILING, WITH EVERY MODULE SET TO PERFECT",
                f"  asymptotic substrate removal: {ceiling['asymptotic_removal_fraction']:.1%}",
                f"  earliest time the {ceiling['required_removal']:.0%} target is met: "
                + (
                    "never, at any horizon tested"
                    if months_to is None
                    else f"{months_to:.0f} months"
                ),
                f"  full eradication: {ceiling['eradication']}",
                "  Perfect modules are not achievable. This is an upper bound on the",
                "  architecture, so the honest reading is that REVERT-1 could approach",
                "  the burden-reversal bar only with near-perfect delivery over years,",
                "  and it would still leave residual substrate behind.",
            ]
        )
    lines.extend(
        [
            "",
            "A NOTE ON CALCIUM",
            "This design deliberately does NOT try to decalcify the artery. Dense",
            "macrocalcification is stabilising; microcalcification in the cap is what",
            "raises local stress. Module M6 suppresses new nucleation and lets existing",
            "deposits mature. Total-calcium removal would be a design error.",
            "",
            "KILL CRITERIA, CHEAPEST FIRST",
        ]
    )
    for item in KILL_CRITERIA:
        lines.append(
            f"  {item.order}. {item.hypothesis_tested}\n"
            f"     falsified if: {item.falsifying_observation}\n"
            f"     cheapest test: {item.cheapest_experiment}"
        )
    lines.extend(
        [
            "",
            "HONEST STATUS",
            "- Every module has a published precedent; all precedents are preclinical or",
            "  mechanistic. The combination, the gate and the coupling claim are invented.",
            "- The simulation is a toy model with assumed rate constants. It cannot",
            "  establish efficacy, and it is excluded from all evidence-graded outputs.",
            "- Even if every module worked, the scar arm is a separate problem: this",
            "  architecture does not regenerate myocardium.",
            "- Nothing here is a treatment recommendation or a patient-level prediction.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_hypothesis_study(
    output_dir,
    draws: int = 400,
    seed: int = 20261107,
    months: float = 24.0,
) -> dict:
    """Write the hypothesis artefacts to their own directory, never the evidence one."""
    from pathlib import Path

    root = Path(output_dir)
    (root / "tables").mkdir(parents=True, exist_ok=True)
    (root / "reports").mkdir(parents=True, exist_ok=True)

    ablation = ablation_study(months=months)
    profile = target_product_profile()
    samples = monte_carlo_hypothesis(draws=draws, seed=seed, months=months)
    sensitivity = sensitivity_summary(samples)
    sweep = best_case_horizon_sweep()
    ceiling = architecture_ceiling()

    ablation.to_csv(root / "tables" / "revert1_ablation.csv", index=False)
    module_frame().to_csv(root / "tables" / "revert1_modules.csv", index=False)
    samples.to_csv(root / "tables" / "revert1_monte_carlo.csv", index=False)
    sensitivity.to_csv(root / "tables" / "revert1_sensitivity.csv", index=False)
    kill_criteria_frame().to_csv(root / "tables" / "revert1_kill_criteria.csv", index=False)
    sweep.to_csv(root / "tables" / "revert1_horizon_sweep.csv", index=False)
    pd.DataFrame([{**ceiling, "status": HYPOTHESIS_TAG}]).to_csv(
        root / "tables" / "revert1_ceiling.csv", index=False
    )
    pd.DataFrame(
        [
            {
                "required_substrate_removal": profile.required_substrate_removal,
                "required_hazard_ratio": profile.required_hazard_ratio,
                "required_lvef_gain_pp": profile.required_lvef_gain_pp,
                "best_current_substrate_removal": profile.best_current_substrate_removal,
                "best_current_lvef_gain_pp": profile.best_current_lvef_gain_pp,
                "substrate_gap_factor": profile.substrate_gap_factor,
                "myocardial_gap_factor": profile.myocardial_gap_factor,
                "status": HYPOTHESIS_TAG,
            }
        ]
    ).to_csv(root / "tables" / "cure_target_product_profile.csv", index=False)
    (root / "reports" / "revert1_design_report.txt").write_text(
        design_report(ablation, profile, samples, sensitivity, ceiling), encoding="utf-8"
    )
    (root / "DISCLAIMER.txt").write_text(DISCLAIMER + "\n", encoding="utf-8")
    return {
        "status": "hypothesis_only",
        "disclaimer": DISCLAIMER,
        "output_dir": str(root),
        "arms": len(ablation),
        "draws": len(samples),
        "probability_of_meeting_substrate_target": float(samples["meets_substrate_target"].mean()),
        "probability_of_phagocyte_depletion": float(samples["phagocyte_depletion"].mean()),
        "architecture_ceiling": ceiling,
        "outputs": 10,
    }
