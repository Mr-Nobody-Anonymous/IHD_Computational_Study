# Mechanism recombination — ATLAS-5

**THIS IS A COMPOSITION OF UNVALIDATED MECHANISMS, NOT A THERAPY AND NOT A RESULT.**
Every part below is borrowed either from an architecture this repository invented, or
from another group's animal or early-phase human work. Randomised proof for the
composition is **zero**, and no simulation can change that. No cure is declared.

```bash
.venv/bin/python run_solutions.py recombine --draws 2500
```

Artefacts: `results/recombine/` — **13 files**, every table carrying
`status = recombine_mechanism_composition_unvalidated`.

---

## 1. What was asked, and what this does

The previous turns built four architectures and then searched for solutions to six open
science gaps. Each was ablated, so we already knew which parts carried the result. The
question here was: **take the successful part of each and build a better one.**

The naive answer is to switch everything on. That answer is wrong, and §6 shows by how
much. The correct answer requires three things that this module implements:

1. **Measure interactions** so two mechanisms doing the same job are not both credited.
2. **Protect the constraints** so the search cannot inflate its score by deleting the
   honesty that produced it.
3. **Charge every precedent axis separately** so a breach on one cannot hide behind
   another.

---

## 2. Corrections applied in code this turn

These are code changes, not documentation notes.

### 2.1 The ATLAS-3 precedent limit was wrong

`MAX_PRECEDENTED_ATTEMPTS` was hard-coded to `2`. CATO (NCT06145035) and HOPE-2 both give
**four** doses in humans. It now reads `max_precedented_doses()` from the registry:

| | before | after |
|---|---|---|
| `MAX_PRECEDENTED_ATTEMPTS` | 2 | **4** (from the registry) |
| `ATLAS3.beyond_precedent` | `True` | **`False`** |

ATLAS-3's three attempts were **never** beyond precedent.

`precedented_option()` had a second, subtler bug: it returned *the limit itself*. That was
harmless when the limit was 2, but after the correction it would have compared three
attempts against a strictly better four and called the difference a cost. It now returns
the **cheapest robust** design inside the limit:

| design | P(target) | adverse risk | robust |
|---|---|---|---|
| **2 attempts — cheapest robust** | **82.50%** | **7.42%** | yes |
| 3 attempts (ATLAS-3) | 93.1% | 8.38% | yes |

The third attempt buys **+10.6 pp** of success for **+0.95 pp** of procedural risk.

### 2.2 The dose-escalation factor's provenance was stale

`REMUSCULARISATION_FACTORS["dose_escalation"]` cited an engineering projection (450M-cell
constructs, ~1400M cells estimated for a 300 g heart). BioVAT-HF (NEJM 2026,
doi:10.1056/NEJMoa2513525) has since **given** 20 EHM units = **800M iPSC-CMs** to **16 of
20 patients**. The provenance is now `human_dose_precedent_escalation_unvalidated`: the
dose is a human precedent, the **1.75× escalation above it** is still unvalidated.
`results/gaps/` regenerated; the stale figure was corrected in `FULL_REPORT.md`,
`gap_closure_attempt.md` and the chemical index.

### 2.3 Four bugs found while building the recombination engine

| bug | effect | fix |
|---|---|---|
| Interaction verdict compared a **probability** difference against **percentage-point** cutoffs | a −14.15 pp overlap was labelled `"additive"` | `INTERACTION_TOLERANCE = 0.005` in probability units |
| Admissibility folded in the **clearance** arm's standing potency breach | every graft composition inadmissible; empty set for the wrong reason | `graft_precedent_breaches` scoped to the arm being searched |
| Response-guided screen tested `gain < required`, which **is** the stop-on-success rule | mechanism was exactly inert by construction | screen on population eligibility, per Vrtovec's baseline-severity signal |
| "Price of admissibility" conflated with "price of cheapness" | report claimed the strongest design needed an unprecedented modification — **it needed none** | three designs compared separately; branches noise-aware |

The last one is the same failure mode as the gap-5 accounting error: **prose asserting a
shape the numbers do not have.** The verdict block is now computed by
`composition_verdict()` and every claim in it is a table lookup.

---

## 3. The mechanism catalogue

Eight parts, decomposed out of ATLAS-2/3/4 and the gap solutions. Six are searchable
benefits; **two are constraints**.

| mechanism | origin | kind | cost axis |
|---|---|---|---|
| `retry_redundancy` | ATLAS-3 / ATLAS-4 | benefit | attempts + procedural risk |
| `cumulative_operator` | ATLAS-4 (corrected from ATLAS-3) | benefit | mechanism assumption |
| `feedback` | ATLAS-3 / ATLAS-4 | benefit | requires reporter |
| `viability_reporter` | gap 3 (hNIS, imaged in patients) | benefit | **genetic modification** |
| `response_guided_selection` | gap 2 (Vrtovec responder analysis) | benefit | screening burden |
| `maturation` | gap 1 (human EHM series) | benefit | time + conditioning |
| `headroom_gating` | ATLAS-4 (human ceiling effect) | **constraint** | none |
| `human_tempering` | ATLAS-4 (10 of 30 responded) | **constraint** | none |

**Constraints lower the reported probability.** Removing `human_tempering` is precisely
what drove ATLAS-4 to a spurious 100%. A search permitted to switch them off would select
against its own honesty, so `exhaustive_subsets()` holds them on in all 64 compositions
and a test asserts `constraints_intact` for every row.

---

## 4. What each part actually carries

Two directions, because one is not enough. Leave-one-out is measured against a
near-saturated base where a part can look inert simply because no headroom remains.

**Leave-one-out** (from everything-on, 2500 draws):

| mechanism | contribution | load-bearing |
|---|---|---|
| `maturation` | **+25.12 pp** | yes |
| `retry_redundancy` | +3.20 pp | yes |
| `cumulative_operator` | +0.08 pp | yes |
| `feedback` | 0.00 pp | **no** |
| `viability_reporter` | 0.00 pp | **no** |
| `headroom_gating` | 0.00 pp | no (masked — see below) |
| `human_tempering` | **−0.24 pp** | yes (constraint) |
| `response_guided_selection` | **−3.52 pp** | yes (negative) |

**Standalone**, add-one-in from an empty benefit set:

| mechanism | gain | risk-adjusted gain | graft breach |
|---|---|---|---|
| `maturation` | **+40.04 pp** | +38.04 pp | none |
| `retry_redundancy` | +20.52 pp | +17.22 pp | none |
| `response_guided_selection` | +0.08 pp | +0.08 pp | none |
| `cumulative_operator` | 0.00 pp | 0.00 pp | none |
| `viability_reporter` | 0.00 pp | 0.00 pp | genetic modification |
| `feedback` | 0.00 pp | 0.00 pp | none |

Three readings:

- **`feedback` and `viability_reporter` are inert.** ATLAS-4 already measured feedback at
  only +0.40 pp once dosing is cumulative. Here it is zero. The reporter cassette exists
  *only* to enable feedback, so both go — and with them the **only genetic modification**
  in the design. That is the single largest cost saving in the whole exercise.
- **`cumulative_operator` is zero alone and valuable in company.** With one dose,
  cumulative and best-of are identical by definition. Its value appears only as an
  interaction (§5).
- **`headroom_gating` reads zero because maturation saturates the first dose**, leaving no
  retries for it to gate. That is masking, not inertness, and it is why the constraint is
  protected structurally rather than kept on merit.

### 4.1 Response-guided selection is a trial-design lever, not a therapy

This is the correction most likely to be misread elsewhere. Vrtovec's responders had
lower pre-dose LVEF (34.9% vs 40.0%), so screening on baseline severity genuinely
concentrates responders — up to **2.34×** in the gap-2 sweep.

But **screening excludes people.** Simulated over the whole population, using the
*central* swept assumption (SD 4 pp, cutoff 37.5 pp → 58.2% responder rate among the
42.5% who pass) rather than the flattering maximum:

- baseline: 33.3% of all patients respond to a repeat dose
- screened: 42.5% pass × 58.2% respond = **24.7% of all patients**

Fewer people benefit. The population effect is **−3.52 pp**, and the mechanism is
correctly excluded. Enrichment shrinks the trial needed to *detect* an effect; it does not
increase the number of patients who *get* one. Reporting only the screened-in subset
would have turned that into a fake therapy gain, so `composite_draws()` simulates
everyone and lets screen failures keep their first-dose gain.

---

## 5. Interactions — measured, so nothing is double-counted

All 15 benefit pairs, on paired draws. Interaction = (gain of both) − (gain of first) −
(gain of second).

| pair | both | sum of parts | interaction | verdict |
|---|---|---|---|---|
| `retry_redundancy` + `maturation` | +47.12 | +61.20 | **−14.08 pp** | **overlapping** |
| `retry_redundancy` + `response_guided_selection` | +14.72 | +23.04 | **−8.32 pp** | **overlapping** |
| `response_guided_selection` + `maturation` | +41.28 | +42.96 | −1.68 pp | overlapping |
| `retry_redundancy` + `cumulative_operator` | +31.12 | +20.64 | **+10.48 pp** | **synergistic** |
| all remaining 11 pairs | — | — | 0.00 pp | additive |

**The overlap is real and it was the thing to look for.** Maturation and retry redundancy
both push the graft toward the same recoverable-scar ceiling. Crediting +40.04 and +20.52
additively would claim **61.20 pp**; the measured joint gain is **47.12 pp**. **14.08 pp of
that sum was the same biology counted twice.** (Both figures come from the interaction
table at one draw count; the standalone table above runs at another, so its +40.04 and
+20.52 are not the addends used here.)

The synergy is equally real: the cumulative operator is worth nothing alone and +10.48 pp
in the presence of retries, because adding increments is only different from replacing
them when there is more than one.

---

## 6. Naive stacking is worse than selecting

All 64 subsets enumerated, constraints held on, draw count floored at
`MIN_SPECIFICATION_DRAWS = 2000` so **the specification cannot move with the sampling
budget** (verified identical at 300, 2000 and 3600 draws).

| design | mechanisms | P(target) | risk-adjusted |
|---|---|---|---|
| everything on | all 6 | 94.40% | 89.48% |
| **selected — ATLAS-5** | **3** | **97.92%** | **92.53%** |
| constraints only | 0 | 51.68% | 49.10% |

**Switching everything on is 3.52 pp worse than selecting**, because
`response_guided_selection` has a negative population effect and drags the stack down.
This is the concrete reason the search was worth running: **recombination is a selection
problem, not an accumulation problem.**

Search outcome: 64 compositions, 33 robust with margin, **0 saturation artefacts**, 16
admissible after every filter.

---

## 7. ATLAS-5

Named only after the search ran, and defined by what survived it.

**`retry_redundancy` + `cumulative_operator` + `maturation`**, four doses available,
both constraints on. No feedback loop. **No reporter cassette. No genetic modification.**

At 2500 draws:

| quantity | value |
|---|---|
| P(LVEF gain ≥ 20 pp) | **97.92%** |
| risk-adjusted | **92.53%** |
| mean doses **actually used** | **1.10** |
| cumulative adverse risk | 5.50% |
| saturated fraction | 79.5% |

At 4000 draws: **96.63%**, risk-adjusted 91.31%. In the 2000-draw search grid: 97.15%,
risk-adjusted 91.75%. The *specification* is fixed; only the estimate moves.

**The mean dose count is 1.10, not 4.** Maturation makes the first dose succeed often
enough that the retry machinery is rarely invoked. The composite is cheap in practice and
not only on paper — which is also why `retry_redundancy` contributes only +3.20 pp in
leave-one-out despite being worth +20.52 pp alone.

### 7.1 Comparison, computed rather than asserted

| comparison | value |
|---|---|
| strongest composition anywhere | 97.20% (5 mechanisms, includes the reporter) |
| strongest **admissible** | **97.15%** (ATLAS-5, 3 mechanisms) |
| best **single** mechanism | 90.60% (maturation alone) |
| gain over best single mechanism | **+6.55 pp** |
| ATLAS-4 best admissible (committed reference) | 89.60% |
| **gain over ATLAS-4** | **+7.55 pp** |
| price of admissibility | 0.05 pp — **within sampling noise** |
| price of cheapness | 0.00 pp |

So the composite genuinely beats both the best single part and the previous best
architecture, and the gain survives charging procedural risk. The reporter buys 0.05 pp —
noise — for a genetic modification with no cardiac precedent. **It is dropped.**

The **frugal alternative** is maturation alone: 90.60%, risk-adjusted 86.07%, one
mechanism, one dose. It is reported alongside, not hidden.

### 7.2 The bill — every axis charged separately

| axis | value | precedented limit | beyond precedent |
|---|---|---|---|
| doses | 4 | 4 | no |
| clearance potency cap | 2.40× | 1.284× | **YES** |
| clearance redundant paths | 3 | 3 | no |
| genetic modifications | **0** | 0 | **no** |

The graft arm is now **fully inside precedent on every axis it controls**. The remaining
breach is the clearance arm's potency requirement, carried forward from the gap-5 sweep
and **not forgiven** by anything the graft arm achieved.

Both arms together (at 1250 draws): graft **96.40%**, clearance **90.40%**, product
**87.15%**. The product assumes
independence; measured bed coupling says the arms are not independent, so the joint is
labelled optimistic. **Randomised proof: 0.00%.**

---

## 8. Does the selection survive its own assumptions?

The largest assumption in ATLAS-5 is a **transfer**: the maturation ratio (1.77×) is
measured in epicardial EHM patches, while the model's enhancement-factor stack is
assembled from other groups' animal work. Multiplying one by the other is unvalidated. So
it is swept, not asserted — 36 scenarios over multiplier × translation penalty × dose
count:

| maturation multiplier | floor 0.4 | floor 0.7 | floor 1.0 |
|---|---|---|---|
| 1.00× (off) | 37.33% | 51.67% | 66.20% |
| 1.25× | 75.27% | 87.33% | 92.93% |
| 1.50× | 75.27% | 87.33% | 92.93% |
| 1.77× | 75.27% | 87.33% | 92.93% |

**18 of 36 scenarios are robust with margin (50%).** ATLAS-5 is *not* robust to a severe
translation penalty on the animal factors combined with no maturation benefit. That is
reported as a fraction, not hidden behind the point estimate.

### 8.1 A declared blind spot in the model

The three lower rows are **identical**. Increments clamp at the recoverable-scar ceiling,
so above **1.25×** the model cannot distinguish a larger multiplier from a smaller one —
it cannot resolve the full 1.77× at all. `maturation_resolution_limit()` computes this
plateau and the report prints it as a **MODEL BLIND SPOT**.

The consequence: the maturation credit is *what this model can resolve*, not the full
human ratio. A model that reported the full 1.77× as realised benefit would be claiming
resolution it does not have.

---

## 9. What this does not establish

- **No randomised proof.** Every mechanism is borrowed from an invented architecture or
  from animal / early-phase human work. The composition has never been tested in anything.
- **No new evidence.** Composing mechanisms rearranges assumptions; it does not
  manufacture data. The evidence credit for every invented part remains zero.
- **The 50% sensitivity figure is the honest headline** for robustness, not the 97.92%
  point estimate.
- **The clearance arm still breaches precedent** on potency, by 1.87×.
- **The saturation fraction is 79.5%.** Most simulated patients hit the modelled ceiling,
  which is a property of the ceiling and the factor stack, not a demonstrated clinical
  result.
- **P = 100% never appears** and would be rejected as an artefact if it did. Zero
  saturation artefacts across all 64 compositions.

**What recombination bought:** a design that is stronger than any single arm (+7.55 pp
over ATLAS-4), *simpler* than the architecture it replaces (3 mechanisms, no genetic
modification, 1.10 procedures on average), and explicit about the 14.08 pp of
double-counting it refused to claim. It is a sharper specification of what to test. It is
not a cure.

---

## 10. Code map

| symbol | purpose |
|---|---|
| `Mechanism`, `MECHANISMS` | catalogue with origin, kind and cost axis |
| `Composite`, `COMPOSITE`, `ATLAS5` | points in mechanism space + precedent predicates |
| `composite_draws()` | paired simulation, whole population, constraints switchable |
| `evaluate_composite()` | probability, risk, saturation, Monte Carlo error |
| `mechanism_marginals()` | leave-one-out |
| `mechanism_standalone()` | add-one-in |
| `pairwise_interactions()` | all 15 pairs, overlap/synergy verdicts |
| `exhaustive_subsets()` | all 64 compositions, budget-floored |
| `select_composite()` | recommended (max risk-adjusted) + cheapest admissible |
| `composition_verdict()` | computed comparison; separates the two prices |
| `stacking_contrast()` | everything-on vs selected vs constraints-only |
| `maturation_sensitivity()`, `maturation_resolution_limit()` | assumption sweep + blind spot |
| `cost_ledger()` | every axis charged separately |
| `run_recombine_study()` | 13 artefacts |

Tests: `tests/test_recombine.py` — **102 tests**, `recombine.py` at **100% coverage**,
including regressions for all four engine bugs, the population-vs-screened-in accounting,
budget-independent specification, the noise-aware price comparison, every alternative
report branch, and the artefact count.
