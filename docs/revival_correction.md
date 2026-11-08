# What the best design was missing — and what correcting it costs

**THIS IS A CORRECTION, NOT A THERAPY. It subtracts a documented control-arm effect
from an uncontrolled observation. It removes attribution; it does not add evidence.
Randomised proof for every design in this study remains 0.00%. No cure is declared.**

```bash
.venv/bin/python run_solutions.py revival --draws 2500
```

Artefacts: `results/revival/` — **9 files**, all carrying
`status = recombination_natural_history_corrected_unvalidated`.

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
| 0 | as published (no control group) | **97.92%** | 100.0% |
| 1 | − unaided gain, upper range (+3.4 pp) | **57.56%** | 58.8% |
| 2 | − unaided gain, lower range (+1.1 pp) | **73.32%** | 74.9% |
| 3 | second uncontrolled patch trial | n/a (null) | — |

**Attribution loss: 40.36 pp.** Roughly two-fifths of the headline was what the
guideline therapy arm achieved for free.

### 3.2 The verdict

| question | answer |
|---|---|
| Is ATLAS-5 still robust? | **NO** |
| Does re-searching the 64 compositions still select it? | **NO** — it selects `retry_redundancy\|maturation` |
| Is *any* composition robust after correction? | **NO** |
| Does it still beat ATLAS-4's 89.60%? | **NO** |

**ATLAS-5 was an artefact of the missing control arm.** The best corrected admissible
design is **55.93%** (`retry_redundancy + maturation`) — the *frugal alternative* that
the recombination search had already reported alongside, which now outranks the design
it lost to.

The correct report is therefore the best available design with
**`robust_set_exists = False`**, not a percentage dressed up as a solution.

### 3.3 Sensitivity: the comparator is swept, not chosen

| unaided gain (pp) | P(target) | median gain | robust |
|---|---|---|---|
| 0.00 | 96.40% | 25.00 | yes |
| 1.10 | 71.92% | 22.93 | **no** |
| 2.25 | 63.92% | 21.78 | **no** |
| 3.40 | 55.12% | 20.63 | **no** |

The design is robust at **one** comparator value — zero — which is the assumption being
tested, not an independent case. At every documented value it fails.

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
| ATLAS-5, 97.92% | **withdrawn** — 57.56%, not robust |
| ATLAS-5 beats ATLAS-4 by +7.55 pp | **withdrawn** — it does not beat it |
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
| `revival_search()` | re-run all 64 compositions under correction |
| `revival_verdict()`, `revival_report()` | computed verdict, no asserted prose |
| `run_revival_study()` | 9 artefacts |

Tests: `tests/test_revival.py` — **42 tests**, `revival.py` at **100% coverage**, all
four own-correction bugs pinned as regressions plus every alternative report branch.
