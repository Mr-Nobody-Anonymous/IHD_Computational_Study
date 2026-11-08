# Resolving the three residuals, step by step

> **STEP 2 PARTIALLY RETRACTED — see [`atlas_self_audit.md`](atlas_self_audit.md).**
> The clearance operating point below (1.284× potency, source reduction 0.50,
> 24 months, 80.01% removal) is a **nominal** result at central module strengths.
> Under delivery uncertainty its median is 1.83% and P(clearance ≥ 80%) is **0.00%**.
> No potency makes it robust while source reduction is uncertain (ceiling 55.4%);
> with source reduction guaranteed, robustness needs ~17.7× potency, not 1.284×.
> Clearance is therefore **not** a solved specification. Steps 1 and 3 are
> unaffected: the evidence-maturity finding concerns real trial anchors, and the
> remuscularisation retraction already used distributional scoring.


**Status: scenario arithmetic.** Each step converts a vague shortfall into a costed
requirement. None of it is efficacy evidence.

```bash
python run_solutions.py resolve      # -> results/resolution/ (9 artefacts)
```

Code: `src/ihd_solutions/resolution.py` · Tests: `tests/test_resolution.py` (35 tests)

| Residual | Before | After this turn |
|---|---|---|
| Hazard ratio | "closed in 1, rest need 16–32%" | **Solved as a procurement problem** — a named list of outcomes trials |
| Clearance | "≥0.223 source reduction, delivery and time binding" | **Solved as an engineering spec** — one potency number per operating point |
| Remuscularisation | "closable on paper, two factors" | **RETRACTED — not robust.** P(+20) = 0.06%, not ~100% |

Two of three moved forward. The third went backwards, which is the more important result.

---

## Step 1 — Hazard ratio: it's an evidence-maturity gap, not a biology gap

Recomputing the required HR per baseline and asking which credit mode first closes it:

| First closing mode | Populations |
|---|---|
| strict | **2** (one is the reference population) |
| surrogate | **6** |
| exploratory | **11** |

**All 19 close by exploratory mode.** Nothing new needs to be invented — the mechanisms
already exist. What's missing is outcomes-tier proof. Ranked by populations unlocked:

| Agent | Unlocks | Current tier | Needs |
|---|---|---|---|
| **verve102_base_editing** | **17** | RANDOMIZED_SURROGATE | RANDOMIZED_OUTCOMES |
| smoking_cessation | 15 | OBSERVATIONAL | RANDOMIZED_OUTCOMES |
| canakinumab_tet2_guided | 11 | OBSERVATIONAL | RANDOMIZED_OUTCOMES |
| retatrutide | 6 | RANDOMIZED_SURROGATE | RANDOMIZED_OUTCOMES |
| obicetrapib | 4 | RANDOMIZED_SURROGATE | RANDOMIZED_OUTCOMES |
| diet_ldl_lowering | 2 | RANDOMIZED_SURROGATE | RANDOMIZED_OUTCOMES |

**A correctness fix worth noting:** my first version of this table also listed aspirin,
bempedoic acid and icosapent ethyl as "needing outcomes trials" — absurd, since they
already have them. Those agents enter relaxed composites only because the *greedy
ordering shifts* when others gain credit. They're now labelled `composition_reordering`
and excluded from the lever list, with a test asserting aspirin can never appear there.

**Verdict: solved as a procurement problem.** One PCSK9-editing outcomes trial unlocks
17 of 19 populations. That is the single highest-leverage item in this entire study, and
it requires no new biology whatsoever.

---

## Step 2 — Clearance: a correction, then an operating point

### The correction

I previously wrote "near-perfect delivery." That was imprecise. Potentiating each module
**alone** to its maximum, at the precedented source reduction of 0.35 over 24 months:

| Module pushed to max | Removal reached |
|---|---|
| crystal_solubilisation | **63.3%** |
| efferocytosis_gain | 62.3% |
| gate_specificity | 62.3% |
| resolution | 60.0% |
| *all modules perfect* | *82.7%* |

**No single lever reaches 80%, and the gate isn't even the strongest one.** The
requirement is joint performance across all modules simultaneously. So the solver
searches one interpretable knob: a uniform potency multiplier on every central strength.

### The operating points

| Source reduction | 24 mo | 36 mo | 60 mo | 120 mo |
|---|---|---|---|---|
| 0.00 | **unreachable** | 1.627× | 1.403× | 1.332× |
| 0.20 | **unreachable** | 1.302× | 1.156× | 1.101× |
| **0.35** | **1.426×** | 1.152× | 1.041× | 1.000× |
| **0.50** | **1.284×** | 1.058× | **1.000×** | 1.000× |

Every entry is the minimum potency multiplier that reaches 80% clearance. Both source
values are inside the CHOP-deletion precedent bound (≤0.50).

**Cheapest admissible point: source reduction 0.50, 24 months, every module at 1.28×
its central strength.** Or relax the timeline: at 60 months with source 0.50, central
strengths suffice (1.00×).

**Verdict: solved as an engineering spec.** The trade-off is explicit — buy time or buy
potency. Without any source reduction, 24 months is unreachable at *any* potency, which
reconfirms that the source module is the necessary ingredient.

---

## Step 3 — Remuscularisation: retracted

Last turn I reported that two co-delivery factors give ×5.2 against a ×5.13 requirement,
so the gap "closes on paper." I flagged the margin as nil. Testing that properly destroys
the result.

Replace each point estimate with **Uniform(1, published effect)** — the factor delivers
anywhere from nothing to its full published benefit — and add a chronic-injury haircut:

| Factor set | Nominal | Median | **P(reach +20)** |
|---|---|---|---|
| **endothelial + epicardial/stromal** (last turn's answer) | +20.3 | +10.1 | **0.06%** |
| all four factors | +25.0 | +20.2 | **50.9%** |

And with the chronic-injury translation haircut, using all four factors:

| Haircut floor | Best P(+20) | Robust at 80%? |
|---|---|---|
| 1.00 (none) | 51.1% | no |
| 0.75 | 38.6% | no |
| 0.50 | 26.0% | no |
| 0.25 | 17.8% | no |

**No combination of the four precedented factors is robust at 80%, under any haircut.**

Why the nominal answer misled: multiplying point estimates assumes *every* factor
delivers its full published effect simultaneously. The expected product of the uniform
draws is far below the product of the maxima. A ×5.2 nominal against ×5.13 required isn't
a narrow pass — it's a coin-flip at best, and worse once translation is discounted.

**Verdict: not solved. The earlier "closable on paper" is withdrawn.** A test now pins
this retraction so it cannot silently regress.

The honest requirement: remuscularisation needs either a **larger single effect** than
any published factor, or **more independent factors** than the four available — and it
needs the chronic-injury efficacy question answered, because that assumption alone moves
the probability from 51% to 26%.

---

## What changed, and what to do about it

| Residual | Verdict | Next action |
|---|---|---|
| Hazard ratio | **Solved — procurement** | One PCSK9-editing outcomes trial unlocks 17/19 populations |
| Clearance | **Solved — engineering spec** | Hit 1.28× central potency across all modules + 0.50 source reduction, or take 60 months |
| Remuscularisation | **Not solved — retracted** | Answer the chronic-injury efficacy question; find a larger or additional factor |

**The ranking has inverted.** Remuscularisation was previously the gap with a paper
solution and the hazard ratio was the one already closed. Now the hazard ratio has a
concrete procurable path, clearance has a numbered specification, and remuscularisation
is the genuinely open problem — which is consistent with it having the largest gap
(×5.13) all along. The paper solution was the artefact, not the progress.

**Method note.** Both corrections this turn came from the same move: replacing a point
estimate with a distribution or a per-population recomputation. The ×5.2 pass and the
"aspirin needs a trial" error were both consequences of treating a single number as if
it were robust. That is the failure mode this repository exists to catch.

---

## Reproduce

```bash
python run_solutions.py composite     # prerequisite
python run_solutions.py resolve
python -m pytest tests/test_resolution.py -q
```

| Artefact | Contents |
|---|---|
| `step1_population_closure.csv` | First closing mode and residual per population |
| `step1_evidence_upgrade_levers.csv` | Agents ranked by populations unlocked, reordering separated |
| `step2_module_leverage.csv` | One-at-a-time potentiation — none suffices |
| `step2_operating_points.csv` | Minimum potency per (source, horizon) |
| `step3_robustness_search.csv` | All 16 subsets with P(target), not nominal gain |
| `step3_chronic_penalty_sensitivity.csv` | Best available set at each haircut floor |
| `resolution_summary.csv` | Cheapest clearance point, robust set (empty), closure counts |
