# What the best design was missing — and what correcting it costs

**THIS IS A CORRECTION, NOT A THERAPY. It subtracts a documented control-arm effect
from an uncontrolled observation. It removes attribution; it does not add evidence.
Randomised proof for every design in this study remains 0.00%. No cure is declared.**

```bash
.venv/bin/python run_solutions.py revival --draws 2500
```

Artefacts: `results/revival/` — **12 files**. The core tables use
`status = recombination_natural_history_corrected_unvalidated`; the counterfactual uplift
table is separately tagged `hypothetical_uniform_uplift_sensitivity_unvalidated`. The
new evidence-gate table only credits a candidate-specific uplift when all prespecified
evidence requirements are met; otherwise the credit is zero.

---

## 1. The task, and why it could not be done the way it was asked

The request was to search the literature for what ATLAS-5 is missing, add that part to
it, and show the result is better. I did the search. **The thing that is missing cannot
be added by recombination, because it is not a mechanism — it is a control arm.**

Every part of ATLAS-5 rests, ultimately, on one uncontrolled observation: BioVAT-HF
reported **+3.9 LVEF points at three months** in 16 patients, *with no control group*.
The trial's own authors say so: "the lack of a control group, low patient number, short
follow-up, and use of surrogate end points will have to be addressed in follow-up
investigations."

Recombination rearranges assumptions. It cannot manufacture a comparator. So instead of
adding a mechanism that does not exist, this module **charges the missing control arm**
and asks whether anything survives.

**It does not.**

---

## 2. What the deep search found

Five deficiencies, each documented, each with a direction it must move the estimate.

| # | Deficiency | Evidence | Charged as |
|---|---|---|---|
| 1 | **No contemporaneous control** | REVIVED-BCIS2 OMT arm gained **+3.4 pp at 6 months unaided** (n=353). BioVAT's +3.9 (90% CI 0.9–6.8) is **not distinguishable** from it | attribution |
| 2 | **A second uncontrolled patch trial was null** | AATS LB8: n=8, LVEF 28.7 → 33.0% at 52 weeks, **p = 0.156** | replication |
| 3 | **The endpoint does not track outcomes here** | STICH, 10 years, n=601: LVEF rose **only with viability, irrespective of treatment assignment**, and **no association with subsequent death** | endpoint validity |
| 4 | **Viability does not identify who benefits** | REVIVED: viability testing did not select responders; extent of **non-viable** myocardium predicted prognosis | mechanism access |
| 5 | **Transmurality, not viability, decides recovery** | no LGE → **91%** recover; >75% transmural → **8%**. Our population has established scar | population definition |

Sources: Perera et al. *NEJM* 2022 (REVIVED-BCIS2); Perera et al. *JAMA Cardiol* 2023
(viability substudy); Panza et al. *NEJM* 2019 (STICH viability, 10-year); Kim et al.
*Circulation* 2000 and *AJR* 2005 (segmental recovery by transmurality); AATS late-
breaking LB8 (allogeneic iPSC-CM patches in ischaemic cardiomyopathy).

**Deficiency 1 is the one that matters.** The comparator is not a placebo — it is
guideline-directed medical therapy, which those patients were already receiving, and
which returns function on its own.

---

## 3. The result

### 3.1 The ladder

| rung | step | P(target) | of original |
|---|---|---|---|
| 0 | as published (no control group) | **97.08%** | 100.0% |
| 1 | − unaided gain, upper range (+3.4 pp) | **57.70%** | 59.4% |
| 2 | − unaided gain, lower range (+1.1 pp) | **73.81%** | 76.0% |
| 3 | second uncontrolled patch trial | n/a (null) | — |

**Attribution loss: 39.37 pp.** Roughly two-fifths of the headline was what the
guideline therapy arm achieved for free.

### 3.2 The verdict and the expanded search

The re-search now varies both mechanism composition and dose count: **64 compositions ×
1–4 precedented doses = 256 composition-dose pairs**. Every pair uses the same Monte
Carlo budget, the two honesty constraints remain on, and selection first prioritises
clearing the 80% target threshold, then corrected risk-adjusted probability and parsimony.
Risk adjustment is
computed per simulated patient as `I(target) × (1 − per-dose risk)^doses_used`, assuming
independent per-dose risk and treating an adverse event as disqualifying; it no longer
applies an average-dose risk to every simulated patient. The full study uses one
shared Monte Carlo budget (20,000 draws) for the ladder, sensitivity analysis, and search.

At **20,000 draws per pair**, the selected simulated design is `retry_redundancy +
maturation` at four permitted doses. I label this search result **ATLAS-5F** below for
tracking only; that name does not imply validation or a real therapy:

| quantity | simulated value |
|---|---:|
| P(LVEF gain ≥ 20 pp), after +3.4 pp comparator correction | **57.81%** |
| risk-adjusted probability | **54.62%** |
| mean doses used | **1.12** |
| any composition clears the model's 80% robustness threshold? | **No** |

ATLAS-5's three-mechanism composition is effectively tied: at the same draw count it
scores about **57.70%** on the corrected target probability and **54.53%** risk-adjusted.
The alternative's risk-adjusted edge is only **0.09 percentage points**, far below the
approximately **0.35 pp Monte Carlo standard error** of a 57–58% target probability at
20,000 draws. At 2,500 draws the ranking can reverse or tie.
So the simulation does **not** establish that the simpler design is biologically superior;
it says the extra mechanism does not earn a reliable modeled advantage. Both are far below
the 80% threshold, and **no composition is robust after correction**.

The old 64-composition, fixed-four-dose result did not test the full design space. Also,
its raw and corrected calculations used unequal draw counts when callers requested fewer
than 2,000 draws. The revised search now floors both sides to the same budget and includes
the precedented dose-count choices. Its better result is a more complete and fair search,
not evidence that the intervention became more effective.

ATLAS-4's **89.60%** is an uncorrected result; comparing it directly with this corrected
search is not an apples-to-apples performance comparison. ATLAS-4 needs the same
correction before a fair numerical comparison.

The correct report is therefore still **`robust_set_exists = False`**, not a percentage
dressed up as a solution.

### 3.3 Counterfactual: the size of an unsupported uplift

I added a transparent what-if calculation—not a proposed molecule or supported mechanism.
If a hypothetical extra LVEF effect were added uniformly to each simulated patient, the
20,000-draw run requires:

| Counterfactual | Uniform extra LVEF gain | P(target) | Risk-adjusted P |
|---|---:|---:|---:|
| Current corrected model | 0.00 pp | 57.81% | 54.62% |
| Just enough to hit 80% target probability | **3.37 pp** | **80.00%** | 75.57% |
| Just enough to hit 80% risk-adjusted probability | **4.34 pp** | 84.69% | **80.00%** |

These shifts are requirements of this simulation's distribution, **not effects known to be
produced by a drug**. They cannot be assigned to colchicine, enlicitide, or another
molecule without evidence. The generated table is
`results/revival/tables/revival_hypothetical_uplift.csv` and every row is explicitly tagged
as an unvalidated hypothetical.

### 3.4 Evidence gate: don't credit what has not been shown

The counterfactual answers how large an extra effect would need to be; it does not say
that a candidate produces it. I added a fail-closed evidence gate so a proposed uplift
cannot enter the corrected score merely because a user supplies a number. To receive any
credit, a record must identify a randomized concurrent control, the same modeled
population, the same candidate-specific intervention, a prespecified LVEF endpoint, a
traceable publication or trial-registry identifier, and a **positive lower 95% confidence
bound** for the between-arm LVEF effect, plus an independently validated patient-level
response model with effect draws aligned to the simulated cohort. A mean effect or its
confidence bound is **not** applied uniformly to every simulated patient. Source truth
and response-model validation still require human verification; code cannot validate a
citation or create trial data.

No qualifying candidate-specific trial record is present in this run. Therefore the gate
credits **0.00 pp**, and the result remains **57.81% target / 54.62% risk-adjusted**;
neither 80% bar is met. It also records the remaining 3.37 / 4.34 pp thresholds. This is
the best defensible model-side fix: prevent unsupported effects from being treated as
capability, and make the missing evidence requirements executable and testable. The
The literature search did find a relevant but non-qualifying signal: DREAM-HF was a
randomized, double-blind, sham-controlled phase 3 trial of mesenchymal precursor cells
(MPCs), with 565 randomized HFrEF participants (about 56% ischemic). The paper reports a
small 12-month LVEF signal among echo-evaluable participants (P=0.041; stronger in the
hsCRP-high subgroup, P=0.008), but the primary endpoint was neutral (HR 1.17, 95% CI
0.81–1.69; P=0.41), and its authors state primary and secondary endpoints were negative.
MPCs are not the iPSC cardiac patch construct represented by ATLAS-5, and the publication
does not provide the validated patient-level response draws required by this model's
gate. It is a real lead to investigate, **not a defensible 3.37-point credit**. Source:
Perin et al., *JACC* 2023, doi:10.1016/j.jacc.2022.11.061 (`NCT02032004`).

The actual missing capability—demonstrated, candidate-specific LVEF benefit plus a
validated response distribution in the target population—requires evidence, not another
mechanism or code edit. See `results/revival/tables/revival_evidence_gated_uplift.csv`.

### 3.5 Sensitivity: the comparator is swept, not chosen

Estimates below use 20,000 Monte Carlo draws.

| unaided gain (pp) | P(target) | median gain | robust |
|---|---|---|---|
| 0.00 | 97.08% | 25.00 | yes |
| 1.10 | 73.81% | 23.36 | **no** |
| 2.25 | 65.98% | 22.21 | **no** |
| 3.40 | 57.70% | 21.06 | **no** |

The design is robust at **one** comparator value — zero — which is the assumption being
tested, not an independent case. At every documented value it fails.

I also re-searched **all 256 composition-dose pairs separately at each nonzero documented
comparator**, so this does not assume the same mechanism mix is best at every point:

| unaided gain (pp) | best design | P(target) | risk-adjusted | either 80% bar met? |
|---:|---|---:|---:|---|
| 1.10 | ATLAS-5F, four doses | **73.91%** | **69.82%** | **No** |
| 2.25 | ATLAS-5F, four doses | **66.09%** | **62.43%** | **No** |
| 3.40 | ATLAS-5F, four doses | **57.81%** | **54.62%** | **No** |

Thus even the most favorable documented comparator and the best composition/dose search
remain below 80%; selecting or composing within the available model cannot reach the bar.
The machine-readable frontier is `results/revival/tables/revival_comparator_frontier.csv`.

---

## 4. The second deficiency, which correcting the first exposes

Subtract the unaided gain and the model still claims a **median +20.6 pp** of LVEF.

That is **3.0× the largest LVEF gain any remuscularisation approach has measured in a
human** (+6.9 pp, BioVAT-HF at a mean 17 months).

So the surrogate threshold comparison — corrected gain clears the 8 pp threshold that
predicts a mortality effect — is **crossed by construction, not by prediction**. It is
labelled `comparison_usable = False` rather than reported as a pass. A model operating
outside the measured range cannot be validated against a threshold inside it.

This is the deepest finding: **ATLAS-5's factor stack is not an extrapolation from
human data, it is a different order of magnitude from human data.** Gap 1 already said
this in its own units — closing +20 pp needs 5.13× the demonstrated effect at three
months, or 2.90× at seventeen months. The composition was credited with reaching it
anyway, because the factors multiply.

---

## 5. The temptation that is *not* a gain

Restricting to patients **without transmural scar** raises the conditional probability
of recovery to 91%.

It is claimed as **zero population-level gain**, and a test enforces that.

| transmurality | conditional recovery | population gain claimed |
|---|---|---|
| no LGE | 91% | **0%** |
| LGE 1–50% | 44% | **0%** |
| LGE 51–75% | 43% | **0%** |
| LGE >75% | 8% | **0%** |

Excluding patients does not make anyone better. This is the same trap the recombination
search already had to reject for response-guided selection (−3.52 pp), and the same one
the ATLAS-4 admissibility filter had to reject for rodent optimism.

---

## 6. What this means for every earlier result

| result | status after correction |
|---|---|
| ATLAS-5, 97.08% at 20k draws | **withdrawn** — 57.70% after +3.4 pp correction, not robust |
| ATLAS-5 vs ATLAS-4 | **not a fair comparison** — ATLAS-4 is uncorrected; apply same correction first |
| ATLAS-3 @ 4 doses, 96.92% | **same defect** — every part rests on the same uncontrolled base |
| ATLAS-4, 89.60% | **same defect** — as is any LVEF-optimising design here |
| Gap 1's 5.13× / 2.90× arithmetic | **stands** — it is a requirement, not a claim |
| Gap 1's `1.77×` maturation multiplier | **stands as an observation, not as a credit** |
| Gap 6 PCSK9 outcomes trial | **unaffected** — randomised, measured durability, hard endpoints |
| Legacy docking/screening | **non-usable, unchanged** |

The distinction that matters: a *requirement* derived from human data is not damaged by
a missing control arm. A *credit* for an effect is.

---

## 7. Was the correction itself honest?

Four bugs were caught in my own correction before shipping, and two of them were
optimistic in the wrong direction:

1. **The comparator was floored at zero.** Truncating `N(3.4, 5)` at zero rents the mean
   up to **4.16 pp**, so the model would have subtracted *more* than the literature
   documents and called it conservative. Patients in REVIVED both improved and worsened;
   the comparator must be signed. The physical bound belongs on the outcome.
2. **A zero comparator was drawn as noise.** `N(0, 5)` leaves the mean where it was but
   changes the variance, so the "no correction" rung disagreed with the published number
   for no reason. Zero now means exactly zero, and a test asserts the rung reproduces
   `evaluate_composite` to 1e-12.
3. **A saturated mean was nearly read as a surrogate pass.** Corrected in §4.
4. **The ledger was not restamped** with this module's tag.

Guards enforced by tests, not by prose:

- A correction may **only ever lower** an estimate. Every comparator is tested.
- More correction never yields more probability (monotone).
- The zero-comparator rung **is** the published figure, exactly.
- No re-search may manufacture a robust design.
- Stratification claims **zero** population-level gain.
- The corrected median must not silently exceed the largest human observation.
- Every table carries the tag; the artefact count is pinned at 9.

---

## 8. The honest closing statement

**There was nothing to add.** The request assumed a mechanism was missing and that
finding it would make ATLAS-5 better. The search found the opposite: the design was
missing evidence, and evidence cannot be composed. Charging that deficiency correctly
does not produce a better version of ATLAS-5 — it removes ATLAS-5.

What survives is narrower and more useful:

- The **composition logic** holds. Re-searching under correction still selects a
  coherent set, just a smaller and weaker one — and it is the frugal design, not the
  elaborate one, which is consistent with maturation being the only part carrying real
  weight (42% of the effect; drop it and the design loses 17.24 pp).
- The **requirement arithmetic** stands. +20 pp LVEF needs 5.13× the demonstrated
  three-month effect, or 2.90× at seventeen months.
- The **gap is now precisely specified**: a randomised comparison of the graft against
  guideline therapy, with an outcome endpoint, at a dose that can plausibly move LVEF by
  more than the ~3.4 points guideline therapy moves it unaided.

That is what BioVAT's investigators have already announced for Phase 3. It is also the
only thing in this entire study that could ever settle the question — and it is not
something a simulation can do.

---

## 9. Code map

| symbol | purpose |
|---|---|
| `Correction`, `CORRECTIONS`, `correction_table()` | the five documented deficiencies |
| `NATURAL_HISTORY_PP`, `NATURAL_HISTORY_RANGE` | +3.4 pp documented, swept 1.1–3.4 |
| `LVEF_CHANGE_SD_PP` | reuses the gap-4 family, not a new invented number |
| `LARGEST_HUMAN_OBSERVATION_PP` | 6.9 pp, the ceiling of the observed range |
| `natural_history_draws()` | signed per-patient comparator; zero means zero |
| `revival_draws()`, `evaluate_revival()` | paired corrected simulation |
| `correction_ladder()` | rung-by-rung attribution loss |
| `natural_history_sensitivity()` | comparator swept, not chosen |
| `surrogate_translation()` | STE comparison, with unusable-comparison labelling |
| `transmurality_stratification()` | the temptation, claimed at zero |
| `revival_search()` | re-run all 64 compositions across precedented dose counts at a chosen comparator |
| `comparator_frontier()` | repeat all 256 composition-dose pairs at each documented nonzero comparator |
| `adjudicate_lvef_uplift_evidence()` | fail-closed gate; requires controlled evidence plus a validated patient-level effect model |
| `evidence_gated_uplift_table()` | records the currently supportable uplift and unmet evidence requirements |
| `revival_verdict()`, `revival_report()` | computed verdict, no asserted prose |
| `run_revival_study()` | 12 artefacts |

Tests: `tests/test_revival.py` — **52 tests**. The full suite is run with coverage and
all four own-correction bugs remain pinned as regressions alongside the evidence-gate
fail-closed and qualifying-record paths.
