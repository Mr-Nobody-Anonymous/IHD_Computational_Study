# Scientific audit — answering the critique, and what it does not fix

**THIS AUDIT CHANGES NO MODEL AND ESTABLISHES NO BIOLOGY.** It measures four things the
critique says are missing and reports the answers, including the negative ones. No cure is
declared and no efficacy claim is made.

```bash
.venv/bin/python run_solutions.py scientific-audit
```

Artefacts: `results/scientific_audit/` — **8 files**, all tagged
`scientific_audit_no_model_change`.

---

## 1. What the critique got right

Most of it. Specifically, and these are now measured rather than argued:

| Critique | Verified? | Evidence |
|---|---|---|
| Parameters are manually invented and drive the outcome | **YES** | **32 of 92** are hypothetical (27 of 73 before the unanswered-question audit added its own). The plaque model assigns MMP activity 1.0 / 0.70 / 0.45 / 0.15 to four regimens with no source at all |
| The inflammatory rates are uncalibrated assumptions | **YES, and worse** | 0 of 14 rates are jointly identifiable. See §3 |
| The CRISPR screen cannot establish genome-wide safety | **YES** | It is 7 hand-written candidate loci, not a genome-wide search |
| ATLAS assumes additive combination | **YES** | Full query in `index 6`; `revival.py` already withdrew ATLAS-5 for a different, larger reason |
| Verification is not validation | **YES** | No module reaches the INTERVENTIONAL rung; 4 of 6 have no external validation |
| The original scripts have no tests | **YES** | `src/0*.py` are outside the test suite and outside the package |

## 2. What the critique got wrong, or describes a state that has moved

I want to be precise here rather than agreeable, because accepting a wrong diagnosis leads
to fixing the wrong thing.

**a) The legacy scripts are not the study.** They are 5 standalone scripts at the repo root
of `src/`. The analysis lives in `src/ihd_solutions/` — 34 modules, 1102 tests, 98.68%
coverage. The critique reviews the former and attributes its faults to the latter.

**b) ATLAS-5 has already been withdrawn — for a bigger defect than the one named.** The
critique's Problem 6 says ATLAS-5 assumes unproven additivity. That was addressed two turns
ago: 64 compositions enumerated with constraints held on, all 15 pairwise interactions
measured, **14.08 pp of double-counted benefit refused**. Then `revival.py` found the larger
problem — the whole design rests on an **uncontrolled** observation — and took it from
97.08% to **57.70%** at 20,000 draws, concluding `robust_set_exists = False`. There is no ATLAS-5 to fix.

**c) "Too many invented biological parameters" is right about the legacy scripts and
partly wrong about the package.** The package's parameters are almost all *measured human
values with citations* (BioVAT, VERVE-102, CATO, REVIVED, CTT). The problem is not that
they are undocumented — it is that **27 are assumptions dressed as parameters**. That is a
different and narrower fault, and §3 measures it.

**d) The GPU plan cannot deliver what it promises, by the critique's own standard.** It
says itself that simulations "cannot independently establish that a new treatment works in
a living human." Agreed. But the plan then lists System F as producing a "virtual
cardiovascular trial platform" and acceptance criteria such as "robustness across patient
heterogeneity." No GPU changes the rung. **You cannot get from CONVERGENCE to
INTERVENTIONAL with more compute**, and the plan does not say what would.

## 3. What was actually built, and what it found

### 3.1 A parameter registry that fails the build

`provenance.py`. Every numeric module constant must declare:

- its **unit**, **source** (≥8 chars), and **population**
- a **class**: `measured` / `estimated` / `hypothetical`
- a **support type**: `efficacy` / `mechanism` / `arithmetic` / `design` / `ceiling`

A missing record raises, and `test_every_numeric_constant_is_registered` fails the build.
This is the difference between documentation and a constraint.

**Result: 92 of 92 registered, 0 unregistered, 0 values disagreeing with the code.**

| class | count |
|---|---|
| measured | 41 |
| estimated | 19 |
| **hypothetical** | **32** |

Two guards are enforced structurally:

- A hypothetical parameter **may not** claim `measured_in_humans=True`.
- `measured_in_humans` requires `population="human"`.

And one that matters for the whole study: **a parameter can be measured and still not
support the claim being made.** `BIOVAT_LVEF_GAIN_3MO_PP` is `measured`, population
`human`, but supports `mechanism`, not `efficacy` — because it comes from a single-arm
study with no control. `earns_efficacy_credit` is `False` for it. Only **12** of 92
parameters earn efficacy credit.

`credit_ceiling()` takes the **minimum** across an input chain, not an average, so one
hypothetical multiplier collapses credit to zero however many real measurements sit beside
it.

### 3.2 Calibration feasibility — the finding that matters

`calibration.py` re-implements the legacy 14-rate cytokine network and verifies it
reproduces the original to 3 decimal places (day-180 CRP **88.613** mg/L).

Then it asks whether the rates could be fitted. The answer is decisive:

| measure | value |
|---|---|
| independent observables | **2** |
| identifiable upper bound | 2 |
| rates "recovered" one-at-a-time | 12 of 14 |
| collinear rate pairs (\|r\| ≥ 0.95) | **7 of 91** |
| rates **jointly** identified | **0 of 14** |

**Every production rate has elasticity exactly +1.0000; every degradation rate −0.99.**
The cascade is *linear in its rates*, so:

> Doubling `k_crp_prod` changes terminal CRP by **exactly +100%**. Doubling `k_il6_prod`
> produces the **same trajectory to machine precision**.

Six production rates share one degree of freedom. This is stronger than statistical
non-identifiability — it is an exact symmetry, so **no measurement of CRP can ever separate
them**, no matter how precise.

The one-at-a-time test "recovers" 12 of 14 because it measures how much the *output* moves,
not what the *data* identify. The joint test rejects all 14.

**Consequence: no Bayesian posterior is reported.** A posterior for these rates would be a
prior with a decimal point. All 14 are held **fixed** and declared assumptions — which is
what they always were.

### 3.3 Identifiability tooling

`identifiability.py` gives the three checks: one-at-a-time elasticity, effect-collinearity
via a sampled Jacobian, and profile identifiability (sweep one parameter, re-fit the rest).
A parameter is `identified` only if it passes **all three**. Anything else is
`chosen_not_measured`.

### 3.4 Verification is separated from validation

`validation.py` defines a five-rung ladder:

| rung | licenses |
|---|---|
| `unit` | the code does what it was written to do — nothing about biology |
| `convergence` | the number is stable under refinement — a converged wrong answer is wrong |
| `internal` | the model fits the data it was built from — circular |
| `external` | predicts data withheld from fitting — the first predictive rung |
| `interventional` | a randomised human comparison supports an efficacy claim |

Claiming a rung without its requirement raises `RungNotEarned`. `convergence_verdict()`
always returns a string containing **"CONVERGED BUT UNVALIDATED"** — convergence is never
reported as validation.

**Highest rung this repository occupies: `convergence`.** 4 of 6 listed modules have no
external validation at all.

### 3.5 The legacy scripts are now labelled, in the file and at runtime

The five root-level `src/0*.py` scripts were the target of the critique. Rather than
delete them (they are the reproduction target for `calibration.py`) they now:

- carry a **NON-USABLE** banner as the first line of their docstring,
- state explicitly **what they can and cannot support**,
- emit a `UserWarning` on import, so nobody reads a number without seeing the label,
- and are covered by `tests/test_legacy_scripts_are_labelled.py`, which fails the build if
  the banner is removed, if the PLAQUE script stops admitting its invented MMP
  activities, if the CRISPR script starts claiming a genome-wide screen, or if any
  package module starts importing them.

The plaque script is the worst offender and now says so itself: cap thickness hardcoded
at 180 um, a 300 kPa rupture threshold now registered as hypothetical, and four regimen
MMP activities (1.00 / 0.70 / 0.45 / 0.15) with **no source at all**.

## 4. What is still missing, stated plainly

`audit_not_done.csv` records the work **not** delivered and why. It is the most important
table in the audit:

| item | why not |
|---|---|
| GPU MD of PCSK9–LDLR | needs a GPU node, OpenMM/GROMACS, multi-µs sampling; no GPU here, and a short run would fail the critique's own replica requirement |
| 3D CFD / FSI of plaque | needs OpenFOAM or SimVascular **plus** a meshed coronary geometry; neither is present |
| genome-wide base-editor off-target model | needs a reference genome and BE-Hive-style calibration data |
| Bayesian posteriors for the cytokine rates | 0 of 14 rates are identifiable, so a posterior would be fiction |
| Virtual cohort with adherence and competing risks | partially present; not calibrated to an external registry |

## 5. The honest answer to "what should I do next"

The critique's priority order — fix parameters, build validation, *then* integrate — is
right. But its first item is already done in the only form that can be done
computationally: the parameters are **typed**, and 27 are now marked as assumptions. That
is a ceiling on what any downstream model can claim, and no amount of GPU time raises it.

**The critique proposes a six-month GPU programme to make "computational evidence so
rigorous that any future laboratory experiments have a clear scientific justification."**
The second half of that sentence is already true, and it did not require a GPU:

- 20 pp LVEF needs **5.13×** the demonstrated effect at 3 months, or 2.90× at 17 months.
- Maturation supplies 1.77× **as an uncontrolled observation** — REVIVED's therapy arm
  gives **+3.4 pp free**.
- The graft needs a randomised comparison against guideline therapy with an **outcome**
  endpoint, because STICH shows LVEF change is not associated with death here.

That *is* a clear scientific justification. What it is not is new evidence, and a GPU does
not make it new evidence.

**One thing I would add to the critique:** its acceptance criteria include "the system must
be allowed to conclude that a treatment has no detectable advantage." Good — but it should
also require the system to say **which parameters it cannot estimate**, because a model
that reports a posterior for an unidentifiable rate has already failed that rule silently.
That is what §3.2 exists to prevent.

---

## 6. Code map

| symbol | purpose |
|---|---|
| `provenance.Provenance`, `REGISTRY` | 73 typed parameters |
| `provenance.unregistered()`, `audit_frame()` | the build guard |
| `provenance.credit_ceiling()` | weakest-link credit across an input chain |
| `calibration.simulate_cytokines()` | faithful re-implementation of `src/04` |
| `calibration.synthetic_recovery()` | plant-and-recover feasibility test |
| `calibration.joint_identifiability()` | the exact-degeneracy measurement |
| `calibration.calibration_verdict()` | 0 of 14 identifiable |
| `identifiability.one_at_a_time()`, `collinearity()`, `profile_identifiability()` | the three checks |
| `validation.EvidenceClaim`, `LADDER` | rung enforcement |
| `validation.convergence_verdict()` | never returns "validated" |
| `scientific_audit.run_scientific_audit()` | 8 artefacts |

Tests: `tests/test_scientific_audit.py` (**79**) + `tests/test_legacy_scripts_are_labelled.py`
(**20**), 35 s combined. All four new modules are at **100% line coverage**. The build
guards are `test_every_numeric_constant_is_registered` and
`test_no_package_module_imports_a_legacy_script`.
