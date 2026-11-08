# What the literature says about each gate item — and what happened when it was applied

**The eight transfer-gate items are not invented here.** Every one corresponds to an
established method with a name, a citation, and usually a software package. This document
records what each method is, what the literature says about whether such problems are
solvable, and — the part that matters — **what happened when the method was applied to this
repository.**

```bash
.venv/bin/python run_solutions.py interventions   # -> results/interventions/ (6 artefacts)
```

**8 items, 8 instruments applied, 0 resolved, 8 still failed.**

---

## The headline: one of my own conclusions was wrong

I reported that the inflammatory model **cannot rank the four mechanisms** it must choose
between, because they differ by **0.000252 mg/L** in terminal CRP.

That number is correct. **The conclusion drawn from it was too strong.**

The day-180 terminal values of CRP are read at a horizon that is a *design choice*, on a
species the cascade *does not measure*, and at a point where the trajectory *has not
converged*. Applying optimal experimental design ([Atkinson & Fedorov 1975](https://doi.org/10.2307/2335289);
[Box & Hill 1967](https://doi.org/10.1080/00401706.1967.10490441)) to the **same model**, over
time and over reachable observables:

| reading | separation | patients per arm |
|---|---|---|
| terminal CRP, day 180 | 0.000252 mg/L | **66,000,000,000** |
| CRP time course, day 22 | **9.954 mg/L** | **6.9** |

**A 39,507-fold improvement from a better design, not from more compute.**

So the honest statement is not "the model cannot rank them". It is:

> The model ranks them — **through a time course of CRP or IL-6, and only through that.**
> The terminal value carries almost none of the information. What the repository lacks is a
> time course, not a ranking capability.

One detail that matters: the mechanism *nodes* (NLRP3, caspase-1, IL-1β) separate by more in
absolute terms, but only because **they are the perturbed quantity themselves** — that is
circular, and the table flags it. CRP is downstream of every knob, which is why its
separation is evidence about the chain rather than a restatement of the perturbation.

---

## The eight instruments

| # | gate item | established instrument | primary source |
|---|---|---|---|
| 1 | a pre-specified endpoint, effect size, population | **estimand framework** — five attributes before any estimator | ICH E9(R1), Step 4, Nov 2019 |
| 2 | the model can rank the mechanisms | **T-optimality / Bayesian entropy design** | Atkinson & Fedorov 1975; Box & Hill 1967 |
| 3 | the model distinguishes the interventions | **structural identifiability** | DAISY (Bellu 2007); SIAN (Hong 2019); STRIKE-GOLDD (Villaverde 2016) |
| 4 | parameters determined by data, not chosen | **practical identifiability via profile likelihood** | Raue et al. 2009, *Bioinformatics* 25:1923 |
| 5 | at least one external validation | **TRIPOD+AI / REMARK** | Collins et al. 2024, *BMJ* 385:e078378; McShane 2005 |
| 6 | predictions are falsifiable | **surrogate endpoint validation** | Prentice 1989; Buyse & Molenberghs 1998 |
| 7 | uncertainty reaches the hypothetical parameters | **sloppy models** | Gutenkunst et al. 2007, *PLoS Comput Biol* 3:e189 |
| 8 | a stated comparator | **computational reproducibility** | Tiwari et al. 2021 (51% of 455 models reproduced directly) |

---

## Item by item: what happened here

### 1. The estimand — **still fails**, and it is a decision, not data

The template is now filled in as far as it honestly can be:

| attribute | state |
|---|---|
| **population** | **NOT SPECIFIED** — 19 populations registered, none declared the target |
| **treatment** | **STATED** — the architecture vs guideline therapy alone |
| **variable** | **PARTIAL** — LVEF (surrogate) and a composite outcome are both in use and are not interchangeable |
| **intercurrent events** | **NOT SPECIFIED** — death, transplant, immunosuppression withdrawal, repeat dosing |
| **summary** | **NOT SPECIFIED** — mean difference, risk difference and hazard ratio all appear |

Under ICH E9(R1), **an unstated intercurrent-event strategy means the treatment effect is
not defined.** That is why this stays item one even though the framework for fixing it is
completely settled. All three missing decisions are free to make.

### 2. Ranking — **corrected, still fails on the evidence held**

The instrument exists and works: the design is resolvable in **7 patients per arm**. The item
still fails because **no such time course exists in `data/`** — the ranking is computed on
nobody's data. The fix is a dataset, not a model.

### 3. Mechanism distinguishability — **still fails, structurally**

LDL lowering and NF-kB inhibition are identical to machine precision, because `oxldl_factor`
multiplies `k_nfkb_act` and appears nowhere else. **No identifiability tool can separate two
parameters that are the same parameter** — DAISY, SIAN and StructuralIdentifiability.jl all
work on the model you give them. The fix is a new ODE term, not a new tool.

### 4. Parameter identification — **still fails**, and the literature says it is normal

Profile likelihood gives the correct answer here: an **infinite** confidence interval. That is
not a failure of the method. The structural result (0 of 14 rates jointly identifiable)
stands.

### 5. External validation — **still fails**, unambiguously

TRIPOD+AI is explicit: **a random split is not external validation**, and reproducing your
own legacy script to three decimals does not count. 0 of 7 named validation targets exist.

### 6. Falsifiability — **still fails**, but the number now exists

The surrogate literature gives the threshold: treatment effects on LVEF explain at most part
of the effect on mortality, with **R² 0.59 and STE +6 pp**. The repository's own corrected
estimate is **below that threshold**. So the repository can now say precisely why its
prediction cannot license an efficacy claim — which is progress, not a pass.

### 7. Uncertainty propagation — **still fails**, and sloppiness is the strongest defence available

Gutenkunst et al. show that in **every** systems-biology model examined, sensitivities span
many decades; the response is **to validate predictions rather than to pin parameters**. That
is a genuine defence of this repository's position — *and it requires the second half*,
prediction validation, which is item 5 and has not been done. Doing the first half and not
the second is the worst of both: unidentifiable parameters *and* unvalidated predictions.

### 8. The comparator — **still fails**, but it is now named

Guideline medical therapy, whose **unaided** effect on LVEF is **+3.4 pp** in REVIVED. Every
gain this repository reports is a difference from baseline, not from that. Subtracting it
does not require a new experiment — it requires editing a number.

---

## What the search settles, and what it does not

**Settled by the literature:** the methods. Prespecification, optimal design, structural and
practical identifiability, sloppy-model analysis, external validation and surrogate
validation are all mature, with published checklists and named software.

**Not settled by any literature:** this study's **data**, its **population**, its **endpoint
choice**, or its **intercurrent-event policy**.

That is the whole answer to "search whether these were solved". The *techniques* were solved
years ago and are free to use. The *evidence* was never solved by anyone, and cannot be.

---

## Code map

| symbol | purpose |
|---|---|
| `interventions.gate_matrix()` | the 8 items joined to the live gate, so a verdict cannot drift |
| `interventions.INSTRUMENTS` | method, citation and honest verdict for each |
| `interventions.t_optimal_design()` | which observable, read when, best separates the mechanisms |
| `interventions.discrimination_power()` | that design converted into patients per arm |
| `interventions.terminal_value_was_the_wrong_coordinate()` | **the correction to this repository's own claim** |
| `interventions.estimand_specification()` | the ICH E9(R1) template, filled in as far as it can be |
| `interventions.write_intervention_report()` | 6 artefacts |

A test asserts the register agrees with `simulation_frontier.transfer_gate()` — so a verdict
here can never drift from what the code computes.

Tests: `tests/test_interventions.py` — **31 tests**, `interventions.py` at **100%**.
