"""Separate software verification from biological validation, and prove nothing else.

The critique names the failure exactly: *"separate software verification from biological
validation."* A test suite that passes says the code does what it was written to do. It
says nothing about whether the biology is right. This module keeps those two claims apart
so they cannot be confused, and gives every result a rung on an explicit evidence ladder.

Five rungs, in increasing strength. Each one is a *different kind of claim*, and a result
may only occupy the top rung it has actually climbed:

``unit``            the code computes what it says. No biology implied.
``convergence``     the output is stable under mesh/step/sample refinement. Still no
                    biology: a converged wrong answer is still wrong.
``internal``        calibrated or fitted to data that was also used to build it - in
                    sample. Circular by construction.
``external``        predicts data held out from the fitting. The first rung with real
                    predictive content.
``interventional``  predicts the effect of an intervention in humans, established by a
                    randomised comparison. The only rung that can support an efficacy
                    claim.

Nothing in this repository reaches ``interventional``. Several things reach
``convergence``. Tests enforce that ordering, and enforce that climbing a rung is never
inferred from a number being large.

Two things this module deliberately refuses to do:

1. Call a passing test suite validation. `verification_vs_validation()` reports both and
   will show, for any module, that verification is complete and validation is not.
2. Treat a converged simulation as predictive. `convergence_verdict()` passes the same
   standard to a wrong-but-stable model as to a right one, and says so.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# ----------------------------------------------------------------------- the ladder

UNIT = "unit"
CONVERGENCE = "convergence"
INTERNAL = "internal"
EXTERNAL = "external"
INTERVENTIONAL = "interventional"

LADDER = (UNIT, CONVERGENCE, INTERNAL, EXTERNAL, INTERVENTIONAL)
RUNG = {name: index for index, name in enumerate(LADDER)}

#: What each rung licenses a result to claim.
LICENSES = {
    UNIT: "the code does what it was written to do - and nothing about biology",
    CONVERGENCE: "the number is stable under refinement - a converged wrong answer is wrong",
    INTERNAL: "the model fits the data it was built from - circular by construction",
    EXTERNAL: "the model predicts data withheld from fitting - the first predictive rung",
    INTERVENTIONAL: "a randomised human comparison supports an efficacy claim",
}

#: The weakest input each rung must have. A rung cannot be claimed without these.
REQUIRES = {
    UNIT: (),
    CONVERGENCE: ("refinement_stable",),
    INTERNAL: ("in_sample_fit",),
    EXTERNAL: ("held_out_fit",),
    INTERVENTIONAL: ("randomised_human_comparison",),
}


class RungNotEarned(Exception):
    """Raised when a result claims a rung it has not climbed."""


@dataclass(frozen=True)
class EvidenceClaim:
    """A claim about a result, with the rung it is allowed to occupy."""

    subject: str
    claimed: str
    satisfied: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.claimed not in RUNG:
            raise ValueError(f"unknown rung {self.claimed!r}")
        missing = [item for item in REQUIRES[self.claimed] if item not in self.satisfied]
        if missing:
            raise RungNotEarned(f"{self.subject}: claims {self.claimed!r} without {missing}")

    @property
    def rung(self) -> int:
        return RUNG[self.claimed]

    @property
    def licenses(self) -> str:
        return LICENSES[self.claimed]

    @property
    def can_support_efficacy(self) -> bool:
        return self.claimed == INTERVENTIONAL


# --------------------------------------------------------------- convergence testing


@dataclass(frozen=True)
class ConvergenceResult:
    """The outcome of refining a discretisation and watching the answer."""

    subject: str
    parameter: str  # what was refined: step size, mesh size, draw count
    coarse: float
    fine: float
    tolerance: float

    @property
    def relative_change(self) -> float:
        scale = max(abs(self.fine), 1e-12)
        return abs(self.coarse - self.fine) / scale

    @property
    def converged(self) -> bool:
        return self.relative_change <= self.tolerance

    @property
    def verdict(self) -> str:
        if not self.converged:
            return "not_converged"
        return "converged_but_unvalidated"


def convergence_verdict(result: ConvergenceResult) -> str:
    """Always returns a string containing 'unvalidated' when it converges.

    The point is rhetorical precision, encoded: convergence is never allowed to be
    reported as validation. A model can be perfectly stable and perfectly wrong.
    """
    if not result.converged:
        return (
            f"NOT CONVERGED: {result.subject} changes by "
            f"{100 * result.relative_change:.2f}% on refining {result.parameter} "
            f"(tolerance {100 * result.tolerance:.2f}%). No claim may be made until it "
            f"stabilises."
        )
    return (
        f"CONVERGED BUT UNVALIDATED: {result.subject} is stable to "
        f"{100 * result.relative_change:.2f}% on refining {result.parameter}. This says "
        f"the answer does not depend on the discretisation. It says nothing about "
        f"whether the answer is right, because a converged wrong model is still wrong."
    )


def monte_carlo_convergence(
    sampler,
    draws: tuple[int, ...] = (500, 2000, 8000),
    tolerance: float = 0.01,
    subject: str = "estimate",
) -> list[ConvergenceResult]:
    """Refine the sampling budget and check the estimate settles.

    This is what `MIN_SPECIFICATION_DRAWS` guards against in the ATLAS series: a
    specification selected on too few draws moves when the budget changes, and a
    specification that moves is not a specification.
    """
    estimates = [float(np.mean(np.asarray(sampler(count), dtype=float))) for count in draws]
    results = []
    for index in range(len(draws) - 1):
        results.append(
            ConvergenceResult(
                subject=subject,
                parameter=f"draws {draws[index]} -> {draws[index + 1]}",
                coarse=estimates[index],
                fine=estimates[index + 1],
                tolerance=tolerance,
            )
        )
    return results


def highest_rung(results: list[ConvergenceResult]) -> str:
    """The best rung a set of convergence results can reach: never above convergence."""
    if not results or not all(result.converged for result in results):
        return UNIT
    return CONVERGENCE


# ------------------------------------------------- verification vs validation report

#: What verification covers, per module, and what validation would require.
MODULE_STATUS: dict[str, dict[str, str]] = {
    "plaque.py": {
        "verified": "unit tests, monotonicity, saturation bounds",
        "validated": "NO - calibrated to a regression grid, not to longitudinal imaging",
        "blocker": "no longitudinal plaque imaging dataset is ingested",
    },
    "inflammation.py": {
        "verified": "unit tests, conservation, positivity",
        "validated": "NO - reaction rates are asserted, not fitted",
        "blocker": "no calibration step exists; see the inflammation work order",
    },
    "editing.py": {
        "verified": "unit tests, floors, monotonicity",
        "validated": "PARTIAL - HEART-2 reductions are measured human values",
        "blocker": "no genome-wide off-target model; 7 candidate loci only",
    },
    "regression.py": {
        "verified": "unit tests, grid bounds",
        "validated": "PARTIAL - calibrated on trial-level rates",
        "blocker": "calibration and evaluation share the same grid",
    },
    "recombine.py": {
        "verified": "86 tests, 100% coverage, budget-stable specification",
        "validated": "NO - every mechanism is invented or borrowed from animal work",
        "blocker": "no randomised human test of the composition exists",
    },
    "revival.py": {
        "verified": "42 tests, 100% coverage, monotone correction, zero rung exact",
        "validated": "NO - correction rests on one external trial arm",
        "blocker": "the control arm is borrowed, not randomised alongside the graft",
    },
}

#: A result may never be reported above the rung its evidence supports.
MAX_REPORTABLE_RUNG = CONVERGENCE


def verification_vs_validation() -> list[dict[str, str]]:
    """One row per module: what is proven and what is not."""
    rows = []
    for module, status in sorted(MODULE_STATUS.items()):
        rows.append(
            {
                "module": module,
                "verified": status["verified"],
                "validated": status["validated"],
                "blocker": status["blocker"],
                "max_reportable_rung": MAX_REPORTABLE_RUNG,
                "can_claim_efficacy": "no",
            }
        )
    return rows


def repository_rung() -> str:
    """The highest rung anything in this repository may occupy."""
    return MAX_REPORTABLE_RUNG


def validation_gap() -> str:
    """A single sentence, computed rather than asserted."""
    rows = verification_vs_validation()
    unvalidated = [row["module"] for row in rows if row["validated"].startswith("NO")]
    if not unvalidated:
        return "every listed module has at least partial external validation"
    return (
        f"{len(unvalidated)} of {len(rows)} listed modules have NO external validation "
        f"({', '.join(sorted(unvalidated))}). The highest rung this repository can "
        f"occupy is {repository_rung()}, which licenses arithmetic stability and "
        f"nothing about biology."
    )


def validation_report() -> str:
    """The audit the critique asked for, as text."""
    lines = [
        "VERIFICATION VERSUS VALIDATION",
        "=" * 78,
        "",
        "Two different questions, routinely confused:",
        "  verification - does the code do what it was written to do?",
        "  validation    - is what it was written to do correct about the world?",
        "",
        "THE LADDER (a result may occupy only the rung it has climbed)",
        "-" * 78,
    ]
    for rung in LADDER:
        lines.append(f"  {rung:15s} {LICENSES[rung]}")
    lines += [
        "",
        "PER MODULE",
        "-" * 78,
    ]
    for row in verification_vs_validation():
        lines += [
            f"  {row['module']}",
            f"    verified:  {row['verified']}",
            f"    validated: {row['validated']}",
            f"    blocker:   {row['blocker']}",
        ]
    lines += [
        "",
        f"  MAX REPORTABLE RUNG: {repository_rung()}",
        "",
        "FINDING",
        "-" * 78,
        f"  {validation_gap()}",
        "",
        "  This is the finding that matters most, and it is not a bug. The repository",
        "  verifies far more than it validates. That is normal for a simulation study,",
        "  but it means every biological number here is a specification of what would",
        "  have to be true, not a measurement of what is.",
        "",
        "WHAT WOULD MOVE THE LADDER",
        "-" * 78,
        "  INTERNAL -> EXTERNAL requires data held out from fitting. There is none:",
        "  every calibration and every evaluation draw on the same trial-level grid.",
        "  EXTERNAL -> INTERVENTIONAL requires a randomised human comparison. Only the",
        "  Phase 3 that BioVAT's investigators have announced could supply it, and no",
        "  amount of simulation can.",
    ]
    return "\n".join(lines)
