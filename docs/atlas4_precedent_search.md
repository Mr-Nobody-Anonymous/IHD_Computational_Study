# ATLAS-4 — searching what was missing, and finding most of it already done

> **ATLAS-4 IS A SIMULATION OF UNVALIDATED ARCHITECTURES.** Repeat cell dosing is
> precedented and feasible in humans, but its efficacy is **not established**: the one
> randomised blinded repeat-dose trial with a functional primary endpoint was **null**.
> The cumulative-dosing mechanism is anchored to rodent and porcine work, and the cells
> involved **do not regenerate myocardium**. Nothing here is a therapy, a result, or a
> recommendation. No cure is declared.

Reproduce:

```bash
python run_solutions.py atlas4 --draws 4000
```

Outputs: 11 artefacts in `results/atlas4/` (9 tables, 1 report, `DISCLAIMER.txt`).
Tests: `tests/test_atlas4.py` (59). Registry: `src/ihd_solutions/precedent.py` (12 sources).

---

## 1. The question, and the honest answer up front

The request was to search for what is missing, check whether someone has already done
it, add the feature, test every possibility, and keep going until all criteria read
100%.

Three of those four are done here. The fourth cannot be, and §7 shows why with measured
numbers rather than assertion. Short version: **a configuration that reads exactly
100% does exist in this model, and it is an artefact of assuming patients respond like
rats.** The model finds it, flags it, and excludes it.

## 2. What the literature search found: most of the "missing step" is already done

ATLAS-3 claimed the missing step was a remuscularisation modality that can be
**repeated and measured**. The search found that repeat dosing is not missing at all —
it is in trial right now, in our exact population.

| What | Status | Source |
|---|---|---|
| Repeat dosing, 4 doses, ischaemic cardiomyopathy | **In trial now** — CATO, n=60, 1:1:1 placebo / 1 dose / 4 doses, 100M UC-MSCs IV every 2 months, primary ΔLVEF by cardiac MRI at 12 months. 47 enrolled as of 2025-11-30. | [NCT06145035](https://clinicaltrials.gov/study/NCT06145035) |
| Repeat sequential dosing, humans | **Done** — HOPE-2, four IV CDC doses every 3 months, randomised double-blind (Duchenne, not IHD) | [Lancet 2022](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00012-5/fulltext) |
| Cumulative benefit of splitting a dose | **Done** — 3×12M CPCs gave +13.0±0.8 pp LVEF vs +4.7±0.8 pp for one 36M dose (p<0.01). The two extra injections *could not* be replaced by tripling the single dose. | [Circ Res 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4992406/) |
| Same result in a large animal | **Done** — 3 IV MSC doses in pigs with chronic ischaemic cardiomyopathy beat 1 dose across echo, MRI and haemodynamics | [Cardiovasc Res 2024](https://academic.oup.com/cardiovascres/article/120/15/1939/7737527) |
| Measuring engraftment well enough to predict function | **Done in animals** — day-1 PET reporter-gene signal vs week-2 function: **ρ=0.76** (echo), **ρ=0.80** (MRI) | [Circ Cardiovasc Imaging 2012](https://www.ahajournals.org/doi/10.1161/circimaging.111.969329) |
| **Response-guided retry timing** | **Proposed in print in 2019, still NOT TESTED** | [Circ Res 2019](https://www.ahajournals.org/doi/10.1161/CIRCRESAHA.118.314653) |
| A viability-specific engraftment signal in humans | **Not done** — the animal signal needs a reporter gene the cells must be engineered to carry | — |

**So the genuinely open item is narrower than ATLAS-3 claimed.** Not "repeat the dose" —
that is already randomising patients. The open item is: *measure engraftment in a human
well enough to decide when to repeat.*

### The uncomfortable finding

The one randomised, blinded repeat-dose trial with a functional primary endpoint
(Vrtovec/REMEDIUM, n=60, two transendocardial CD34+ doses 6 months apart) was **null
versus single dosing**. Only **10 of 30** patients responded to the second dose, and
responders had *lower* pre-dose LVEF than non-responders (34.9% vs 40.0%, p=0.02) —
a **ceiling effect**. Two smaller human studies (Yao, Gu) were positive, but one was not
placebo-controlled.

Repeat dosing is precedented and feasible. **Its efficacy in humans is not established.**

## 3. Three corrections to ATLAS-3

### 3.1 The precedent limit was wrong

ATLAS-3 used `MAX_PRECEDENTED_ATTEMPTS = 2` and labelled its three-attempt design
"BEYOND PRECEDENT". CATO and HOPE-2 both give **four** doses. The limit is 4, so
**ATLAS-3's three attempts were never beyond precedent at all.** ATLAS-3 is left in the
tree as the historical record with a correction banner; ATLAS-4 supersedes it.

### 3.2 The retry operator was wrong

ATLAS-3 modelled retries as **best-of-k** — keep the best outcome. Every preclinical
study instead reports **cumulative** increments: +3.4, +4.0, +5.5 pp per successive dose.
Best-of-k cannot represent the field's central finding, because it makes splitting a
fixed total dose irrelevant, when splitting is precisely what works.

Both operators are now implemented and contrasted on identical draws:

| Doses | Best-of-k | Cumulative | Difference |
|---|---|---|---|
| 1 | 49.60% | 49.60% | +0.00 pp |
| 2 | 60.52% | 64.40% | +3.88 pp |
| 3 | 67.85% | 75.12% | +7.27 pp |
| **4** | 74.12% | **82.97%** | **+8.85 pp** |
| 6 | 83.70% | 93.03% | +9.33 pp |

Per-dose retention when splitting is anchored to the rat arithmetic:
`(13.0/3)/4.7 = 0.922` — a third-dose delivers 92% of a full dose's effect. Duration of
exposure matters more than intensity.

### 3.3 The mechanism label was wrong

These cells **do not engraft, do not become cardiomyocytes, and do not regenerate
myocardium**. The pig study states outright that "myocardial regeneration is not
plausible"; the rat studies found new cardiomyocyte formation negligible in all groups.
The benefit is paracrine, anti-inflammatory and anti-fibrotic. Modelling their LVEF gain
as scar replacement overstates the claim, and `SCAR_RECOVERABLE_LVEF_PP` should be read
as a functional ceiling, not as regenerated tissue.

## 4. A consequence: the independence requirement was operator-specific

ATLAS-3's headline specification was "retry outcomes must be independent to at least
≈0.71". That turns out to be a property of the **wrong operator**:

| Mode | Gain at correlation 0 | Gain at correlation 1 | **Retained** |
|---|---|---|---|
| Best-of-k | +24.52 pp | +8.77 pp | **35.8%** |
| Cumulative | +33.38 pp | +30.25 pp | **90.6%** |

Best-of-k loses most of its benefit when doses are correlated, because correlated
attempts are close to one attempt. Cumulative dosing keeps almost all of it: a patient
whose biology is consistently favourable accumulates consistently favourable
increments. It stays robust even at correlation 1.

**ATLAS-3's most emphasised unmeasured quantity mattered less than it claimed.**

## 5. The measurement falsifier

Feedback is only worth what the signal actually measures:

| Modality | Quality | Viability-specific | Gain vs no feedback |
|---|---|---|---|
| Reporter-gene PET | 0.78 | yes | +0.40 pp |
| Iron-oxide MRI | 0.00 | no | **+0.00 pp (exact)** |
| Acute retention only | 0.00 | no | **+0.00 pp (exact)** |

Iron-oxide MRI localises cells but its signal persists in macrophages after the cells
die ([J Nucl Med 2009](https://jnm.snmjournals.org/content/50/7/1088)), and acute
retention does not correlate with engraftment
([J Transl Med 2017](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-017-1157-0)).
A loop built on either is blind, and gains exactly nothing. That is a published
finding, not a modelling convenience.

Note also how *small* the feedback gain is under the correct operator (+0.40 pp). Once
dosing is cumulative and headroom-gated, feedback has little left to add — another
ATLAS-3 emphasis that does not survive the correction.

## 6. Every possibility tested

The full Cartesian grid: 2 operators × 5 dose counts × 3 correlations × 3 learning
levels × 2 headroom settings × 2 responder assumptions = **360 configurations**, all
scored, none pre-filtered.

| | Count |
|---|---|
| Configurations tested | **360** |
| Robust (P ≥ 0.8) | 176 |
| Robust with 2-SE margin | 160 |
| Within precedent | 288 |
| **Saturation artefacts (P = 100%)** | **48** |
| Admissible after all filters | **8** |

**Best admissible point:** cumulative, 4 doses, correlation 0.50, learning 1.00,
human-tempered repeat response — **P 89.60%**, risk-adjusted 82.65%.

### A filter bug I had to fix

My first admissibility filter admitted any robust, precedented, non-saturated point. The
"best" point it returned was the **rodent-optimistic** one at P=99.83% — sitting just
under certainty, and resting on an assumption the human data contradicts. That is how a
grid search launders an assumption into a result. The filter now also requires that the
repeat-response assumption be **no more optimistic than the human observation**
(10/30). Admissible points fell from 32 to 8.

## 7. Can every criterion reach 100%? Measured, not asserted

Catalogue ceilings, computed from the real assessments in the reference population:

| Criterion | Ceiling | At 100%? | Meets 0.8? |
|---|---|---|---|
| durable_driver_control | **1.000** | **YES** | YES |
| burden_reversal | 0.600 | no | no |
| myocardial_restoration | 0.500 | no | no |
| risk_normalisation | 0.061 | no | no |
| randomised_proof | **1.000** | **YES** | YES |

**Two of the five already read 1.000.** They were never the binding constraint. What
blocks the other three:

- **burden_reversal** — capped at 0.6 by the partial-regression rule for any single
  agent. ATLAS-2 reaches the burden target in simulation only via an invented 2.80×
  potency cap at 60 months.
- **myocardial_restoration** — capped at 0.5 without human remuscularisation evidence.
  ATLAS-4 raises the simulated *probability*; it cannot raise the *score*, because the
  score is an evidence cap, not a probability.
- **risk_normalisation** — capped at 0.4 per single agent; the measured catalogue
  maximum in this population is 0.061. No single therapy returns event risk to a
  low-risk reference. That is a property of the effect sizes, not a search failure.
- **randomised_proof** — reaches 1.0 for outcomes-tier catalogued agents, but proof
  attaches to **the exact regimen randomised**. It does not transfer to an untested
  combination, so every invented architecture here scores 0.

**Raising these to 1.0 means deleting the caps.** That would manufacture the number
without producing the evidence, which is the one thing this study has consistently
refused to do.

### And the 100% that does exist

| Design | P(target) | Saturated fraction | Flagged |
|---|---|---|---|
| ATLAS-4, human-tempered (10/30 respond to a repeat) | **82.97%** | 52.2% | clean |
| Same design, rodent optimism (everyone responds) | **100.00%** | 55.8% | **saturation artefact** |

The certain version stays at **99.98%** even with a harsh translation penalty
(`translation_floor = 0.4`). So the near-certainty is driven **almost entirely by the
responder assumption**, not by the biology. It is what you get if you assume every patient
responds to every repeat dose the way a rat does. The human data says a third do.

## 8. Result

| | Value |
|---|---|
| Single dose | 49.60% |
| **ATLAS-4, human-tempered, 4 doses** | **82.97%** |
| Median gain | 25.00 pp (at the functional ceiling) |
| Robust with 2-SE margin | **yes** |
| Beyond precedent | **no** — four doses has human precedent |
| Cumulative adverse risk | 7.98% |
| Risk-adjusted success | 76.35% |

Translation sensitivity (repeat response × translation floor) spans **66.35% to
100.00%**. The honest headline is the human-tempered, unpenalised end: 82.97%. Apply a
0.7 translation floor and it falls to **76.40%** — **no longer robust**. That is reported,
not hidden.

## 9. What to do next

The highest-value experiment is no longer the one ATLAS-3 proposed. It is:

1. **Read out CATO.** It is already randomising patients to four doses against one and
   placebo, with ΔLVEF by MRI at 12 months. It tests the operator question directly and
   needs no new biology. Follow-up completes March 2027.
2. **Test response-guided retry timing.** Proposed in print in 2019, never tested. The
   Vrtovec responder analysis says the signal to guide it is simply pre-dose LVEF —
   available today, no reporter gene required.
3. **Build a viability-specific human engraftment signal.** The animal version needs an
   engineered reporter gene; iron MRI and retention imaging are demonstrably blind.

And unchanged from every prior turn: **a PCSK9-editing outcomes trial** remains the
highest-leverage real-world action for the evidence gap, unlocking 17 of 19 populations.

**No cure is declared. Randomised proof for every architecture here remains zero.**

---

## See also

- [`atlas3_graft.md`](atlas3_graft.md) — the superseded best-of-k version, with its correction banner
- [`atlas2_upgrade.md`](atlas2_upgrade.md) — the clearance arm
- [`FULL_REPORT.md`](FULL_REPORT.md) — complete inventory and the honesty ledger
- [`scientific_methodology_and_limitations.md`](scientific_methodology_and_limitations.md)
