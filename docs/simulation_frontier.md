# The simulation-first frontier

**The question.** *"Find the solution in simulation only, then transfer to the lab."* What
would have to exist for that plan to work — and what exists here?

```bash
.venv/bin/python run_solutions.py frontier     # -> results/frontier/ (7 artefacts)
```

**Transfer gate: 0 of 8 checks pass. Transfer authorised: NO.**

No model is changed, nothing is tuned, and no claim is upgraded.

---

## The result that decides the plan

The inflammatory model contains six mechanisms a laboratory could target. It orders them
into **coarse classes** — upstream versus downstream is visible. But it **cannot rank the
four downstream mechanisms at all:**

| mechanism (halved) | terminal CRP |
|---|---|
| untreated | 88.6134 |
| LDL lowering *or* NF-kB inhibition | 44.3508 |
| NLRP3 inflammasome | 44.3070 |
| caspase-1 | 44.3067 |
| IL-1β | 44.3067 |
| IL-6 | 44.3067 |

The four distinct anti-inflammatory targets differ by **0.000252 mg/L** — one part in five
million. At a between-patient CV of 30% and an assay CV of 10%, separating them would need
**66 billion patients per arm.**

This is not a resolution problem more compute fixes. It is structural: each parameter is a
multiplicative factor on a linear chain, so the chain product is invariant to *which* factor
you adjust.

### And one comparison is a parameter against itself

**LDL lowering and NF-kB inhibition are the same perturbation to machine precision** (max
difference over every species and every time point: **exactly 0.000e+00**). `oxldl_factor`
multiplies `k_nfkb_act` and appears nowhere else. Anyone reading this model as "LDL drug
versus anti-inflammatory drug" is comparing a parameter with itself.

### Why this blocks the plan specifically

Choosing **which** experiment to run requires knowing which mechanism matters — which is
what the experiment is for. Selecting the target by fitting the unfitted is circular. That
is the reason a simulation cannot terminate in "the solution", however much compute it gets.

---

## The transfer gate

Eight preconditions for "simulation done, now transfer". **Zero pass.**

| # | requirement | state |
|---|---|---|
| 1 | **a pre-specified endpoint, effect size and population defining "the solution"** | **FAIL** — never specified. This is the single most important missing item |
| 2 | the model can rank the mechanisms it is asked to choose between | **FAIL** — 0.000252 mg/L apart |
| 3 | the model distinguishes the interventions it claims to compare | **FAIL** — LDL ≡ NF-kB, difference 0.0 |
| 4 | fitted parameters are determined by data rather than chosen | **FAIL** — 0 of 14 rates identifiable; 32 of 97 parameters hypothetical |
| 5 | at least one external validation | **FAIL** — 0 of 7 targets present in `data/` |
| 6 | predictions are falsifiable | **FAIL** — every test asserts the code does what it was written to do; none asserts a quantity external data could contradict |
| 7 | uncertainty is propagated | **FAIL** — the interval covers seeds and scenarios, **not** the 32 hypothetical parameters, because they have no distribution to propagate |
| 8 | the transfer has a stated comparator | **FAIL** — REVIVED's guideline-therapy arm moves LVEF +3.4 pp unaided; the architecture has never been compared with it here |

A gate that cannot pass is not a gate, so a test forces it to pass when the checks are
satisfied. It is a real checklist, not a rhetorical device.

---

## What is missing, in order

1. **Define "the solution."** An endpoint, a threshold, a population. Without it, any result
   can be narrated as success afterwards — this is not a technical gap, it is the whole gap.
2. **One external validation set.** No GPU, no patients. 0 of 7 targets exist in `data/`.
3. **A second independent observable** in the cytokine chain. Fourteen rates with one
   observable is a one-dimensional axis; the model needs independent species to break it.
4. **Uncertainty distributions for the 32 hypothetical parameters** — or the honest admission
   that the output interval excludes them.
5. **A falsifiable quantity per module** — a number external data could contradict.
6. **A comparator**, so the intervention can be separated from background therapy.

None of these is a GPU workload. Items 1, 5 and 6 are free.

---

## What simulation CAN hand over

Being unable to name a winner is not being unable to contribute. The correct output of a
simulation-first programme is an **experiment specification**:

| # | measurement | resolves | needs |
|---|---|---|---|
| 1 | serial IL-6 and CRP time course under a defined stimulus | whether the chain structure is right at all | ≥6 time points, inter-individual SD reported |
| 2 | randomised dose-response to an NLRP3/IL-1β inhibitor | mechanism attribution, which the model cannot supply | ≥3 doses, placebo, randomised |
| 3 | paired serial cap thickness and LDL in the same patients | whether cap thinning follows LDL as assumed | ≥2 time points, ≥200 patients |
| 4 | randomised cell-therapy arm against guideline therapy | the entire efficacy question | ≥26 per arm for LVEF; an outcome endpoint for survival |

Plus five deliverables a laboratory cannot produce: a falsifiable prediction with stated
resolution; the **minimum detectable difference** for any proposed trial; which combinations
are redundant because they share a pathway; arithmetic limits on what a trial of a given
size could show; and an ordered list of experiments by information per cost.

---

## The honest bottom line

> **CORRECTED — see [`literature_register.md`](literature_register.md).** The 0.000252 mg/L
> reading above uses the day-180 terminal value of an *unobserved* species at a
> *design-chosen* horizon. Over a CRP time course the same four mechanisms separate by
> **9.954 mg/L around day 22** — **39,507×** more, resolvable with **7 patients per arm**.
> The model **can** rank them through a time course; what the repository lacks is the time
> course, not the capability. The claim below is therefore narrowed: a model cannot choose
> an experiment it has no data for.

**Simulation cannot find "the solution" — not because there is not enough compute, but
because a model cannot choose the experiment that would tell you whether the model is right.**

What it can do is narrower and genuinely valuable: rule things out, bound what a trial could
show, and specify the measurement that would discriminate. That is worth a laboratory's
time. "The solution" is not a thing this or any simulation can hand over, and a programme
that promises it will end by narrating whatever the first experiment returns as success —
which is precisely why gate item 1 comes first.

Tests: `tests/test_simulation_frontier.py` — **32 tests**, `simulation_frontier.py` at
**100%** line coverage.
