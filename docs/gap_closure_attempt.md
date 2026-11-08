# Trying to close the three gaps

**Status: gap-closure scenario.** Closing a gap in a model is a reason to run an
experiment, not a reason to believe the result. Every precedent below is animal,
construct-level or mechanistic unless explicitly stated.

```bash
python run_solutions.py gaps      # -> results/gaps/ (9 artefacts)
```

Code: `src/ihd_solutions/gap_solver.py` · Tests: `tests/test_gap_solver.py` (39 tests)

| Capability | Needed | Have | Gap | Verdict |
|---|---|---|---|---|
| Remuscularisation | +20 LVEF pts | +3.9 | ×5.13 | **closable — 2 factors, both animal-level** |
| Substrate clearance | 80% | 43.1% | ×1.86 | **closable — needs a new module, plus delivery or time** |
| Stack hazard ratio | ≤0.298 | 0.360 | ×1.21 | **already closed** (in 1 of 19 populations) |

I did not tune parameters until the numbers worked. For each gap I asked what *kind* of
obstacle it was, and in two cases that changed the answer.

---

## Gap 3 — hazard ratio: already closed, and the target was misapplied

The composite from the previous turn achieves **HR 0.2818** in stable coronary disease
against a required 0.298. Closed, with no new biology at all.

But the 0.298 figure was derived from a **13.93% baseline**. Applying it to every
population is a category error: a higher baseline needs a *lower* HR to leave the same
share of excess risk. Recomputed per population:

| Population | Required HR | Achieved | Status |
|---|---|---|---|
| **stable_ccs** | ≤0.298 | **0.282** | **CLOSED** |
| secondary_post_mi | ≤0.257 | 0.309 | needs 16.4% more |
| hefh_primary | ≤0.355 | 0.464 | needs 23.1% more |
| t2d_ascvd | ≤0.267 | 0.351 | needs 23.1% more |
| post_mi_hscrp_high | ≤0.248 | 0.368 | needs 31.7% more |
| statin_intolerant_ascvd | ≤0.267 | 0.398 | needs 32.0% more |

So: **closed in 1 of 19 populations** in strict mode (plus `lifetime_exposure`, which is
already the reference population and isn't a gap). The remaining populations need a
further 16–32% risk reduction on top of the composite. This gap is now a *quantified
remainder*, not a vague shortfall.

---

## Gap 2 — substrate clearance: the plateau was structural

### Diagnosis

With **every module perfect**, the REVERT-1 loop still stalls at 83.5%. I derived the
steady state in closed form and it matched the simulation to three decimals:

- Debris: uptake `5.4·D/(D+15)` balances the source `1.5/month` → **D\* = 5.77**
- Cholesterol: efflux + oxysterol exit balances `2.5 × source` → **F\* = 12.06**
- Total 17.83 of 108 → **83.49%** removal (simulation: 83.49%)

That agreement is the diagnosis. **The ceiling is set by a constant apoptosis source
that none of the six modules touches.** M5 (resolution) protects uptake and reduces the
death rate — it does not reduce the rate at which *new* apoptotic cells appear. No
amount of payload strength fixes a source term.

### The fix: a seventh module aimed at the source

ER-stress/CHOP-mediated macrophage apoptosis is the dominant source in advanced lesions,
and it is causally tractable: **CHOP deletion in Apoe−/− and Ldlr−/− mice reduced lesional
apoptosis ~35% and plaque necrosis ~50%**
([Thorp et al.](https://www.researchgate.net/publication/24402929_Reduced_Apoptosis_and_Plaque_Necrosis_in_Advanced_Atherosclerotic_Lesions_of_Apoe--_and_Ldlr--_Mice_Lacking_CHOP),
[Tabas review](https://www.ahajournals.org/doi/10.1161/circresaha.110.224766)).
CHOP expression correlates with apoptosis in **human** vulnerable plaques. Chemical
chaperones (4-PBA, TUDCA) reduce macrophage apoptosis and stabilise plaques in animals.

Crucially the precedent is **bounded**: CHOP deletion does *not* abolish apoptosis
(IRE1-JNK and Bcl-2 pathways remain), so I capped the module at **0.50**.

### Result — perfect modules

| Source reduction | 24 mo | 36 mo | 60 mo | Asymptote | Within precedent |
|---|---|---|---|---|---|
| 0.00 | 74.9% | 80.4% | 82.3% | 83.5% | yes |
| 0.20 | 79.5% | 84.2% | 85.8% | 87.0% | yes |
| **0.35** | **82.7%** | 86.8% | 88.3% | **89.5%** | **yes** |
| 0.50 | 85.8% | 89.4% | 90.8% | 92.0% | yes |
| 0.65 | 88.7% | 91.9% | 93.2% | 94.4% | **no** |

Minimum for 80% at 24 months: **0.223** — comfortably inside the 0.35 demonstrated
apoptosis reduction. The analytic ceiling tracks the simulation at every strength (tested).

### The honest caveat — central modules

| Source reduction | 24 mo | 60 mo | Asymptote | Meets 80% |
|---|---|---|---|---|
| 0.00 | 43.8% | 62.4% | 66.8% | no |
| 0.35 | 56.5% | 77.8% | **80.4%** | yes, asymptotically |
| 0.50 | 62.0% | 83.0% | 85.2% | yes |

With realistic rather than perfect module strengths, closing the gap in **two years**
needs source reduction **0.989** — far beyond precedent. It closes at **0.338** only over
a much longer horizon.

**So the real finding: delivery efficiency and time are jointly binding.** That is
consistent with the earlier sensitivity result where the gate dominated everything. The
new module is necessary but not sufficient; you also need near-perfect delivery or patience.

One more result worth stating plainly: **eradication requires abolishing the apoptotic
source entirely** (`ceiling(1.0) = 100%`). Nothing short of that fully clears the lesion,
and total suppression has no precedent at all.

---

## Gap 1 — remuscularisation: decomposes, but every factor is animal-level

### The calibration that makes this tractable

BioVAT's +3.9 points is **15.6% of the 25-point deficit**. Primate grafts occupy
**10–17% of scar area** — an independent measurement that lands in the same place. So
modelling LVEF gain as proportional to *fraction of scar remuscularised* is calibrated
twice, and reaching +20 points means replacing **80% of the scar**.

### Precedented multipliers

| Factor | Multiplier | Precedent | Species |
|---|---|---|---|
| Dose escalation | ×1.75 | 20 EHM units = 800M iPSC-CMs in 16 of 20 patients (BioVAT-HF, NEJM 2026); escalation above it unvalidated | human dose precedent |
| **Endothelial co-delivery** | **×2.0** | EC+CM raised primate EF **+8.9 pts vs +1.9** for CM alone; preformed microvessels gave 10× larger grafts | **primate** |
| **Epicardial/stromal co-delivery** | **×2.6** | Doubled CM proliferation, **2.6-fold greater graft size**, improved systolic function | rodent |
| Proliferation induction | ×2.0 | CCND2-overexpressing CMs replaced **>50% of scar** at 6 months; CHIR/FGF patches gave 4× engraftment | rodent |

Sources: [EC+CM in primates](https://pmc.ncbi.nlm.nih.gov/articles/PMC10683868/),
[dose-dependent EHT patches](https://www.ahajournals.org/doi/10.1161/CIRCULATIONAHA.120.047904),
[transplantation review](https://www.nature.com/articles/s44385-025-00048-4).

### The saturation rule

All four factors act on **one limiting quantity** — engrafted contractile mass — so their
product cannot exceed the scar available to replace. Multiplying them without a cap would
produce absurd gains. With the cap:

| Factors | Scar replaced | LVEF gain | Meets +20 |
|---|---|---|---|
| none | 15.6% | +3.9 | no |
| any single factor | ≤40.6% | ≤+10.1 | no |
| **endothelial + epicardial/stromal** | **81.1%** | **+20.3** | **YES** |
| epicardial/stromal + proliferation | 81.1% | +20.3 | YES |
| any three or more | 100% (saturated) | +25.0 | yes |

**Minimum set: two co-delivery changes.** Notably this needs *neither* dose escalation
*nor* genetic modification — it is a change to the **cellular composition of the graft**,
which is the cheapest of the four levers. ×5.2 available against a ×5.13 requirement:
it just clears, with essentially no margin.

### Why I do not believe this yet

- **Every multiplier is rodent, primate or construct-level.** None is human efficacy.
- **Chronic injury is the wrong model.** Most primate work treats **two weeks** after
  infarction. An advanced heart-failure patient has a years-old, stiff, poorly vascularised
  scar. Efficacy there is explicitly unproven.
- **Arrhythmia.** Cell injection produced graft-associated ventricular arrhythmias in
  primates; epicardial patches in 20 macaques showed none, attributed to mechanical rather
  than electrical coupling. That distinction is load-bearing and not fully settled.
- **The margin is nil.** ×5.2 against ×5.13 means any factor underdelivering by 2% breaks it.
- **Immunosuppression** is required, and is not modelled as a cost here.

---

## What this establishes, and what it doesn't

**Established:**

1. The hazard-ratio gap was already closed, and the target needed recomputing per
   population — the remainder is now 16–32% in the populations still open.
2. The clearance plateau had a **structural cause** that payload tuning could never fix.
   Adding a source-reduction module is a design change with bounded animal precedent, and
   it works — but only alongside near-perfect delivery or a long horizon.
3. The remuscularisation gap **decomposes into two precedented co-delivery changes**, and
   the cheapest path is graft composition rather than more cells or gene editing.

**Not established:** that any of it works. Two of three verdicts are
`closable_within_precedent`, which means *the arithmetic permits it using published
animal multipliers*. Zero human efficacy evidence was created here. And even if all three
gaps closed, `randomised_proof` remains untouched — so this is still not a cure.

**The most useful single result:** the clearance ceiling was never a payload problem. It
was a source problem, found by deriving the steady state and checking it against the
simulation. That is the kind of thing simulation is genuinely good for.

---

## Reproduce

```bash
python run_solutions.py composite     # prerequisite: the hazard matrix
python run_solutions.py gaps
python -m pytest tests/test_gap_solver.py -q
```

| Artefact | Contents |
|---|---|
| `gap_verdicts.csv` | Three verdicts with mechanism and residual requirement |
| `gap_hazard_by_population.csv` | Required vs achieved HR, per population and mode |
| `gap_clearance_source_reduction.csv` | Perfect modules, 6 strengths, analytic vs simulated |
| `gap_clearance_central_modules.csv` | The same under realistic strengths |
| `gap_remuscularisation_search.csv` | All 16 factor subsets, with saturation |
| `gap_remuscularisation_factors.csv` | Multipliers, species, precedent, URL |
| `gap_detail.csv` | Minimum strengths and the joint delivery/time constraint |

The source-reduction term defaults to zero in `simulate_clearance`, so **every existing
hypothesis and discovery artefact is unchanged** — a test asserts this explicitly.
