# ATLAS-1 self-audit: the operating point does not survive uncertainty

> **HORIZON CORRECTION — see [`atlas2_upgrade.md`](atlas2_upgrade.md).** The "55.4%
> ceiling at any potency" in section 3 is correct **at 24 months only**. Perfect
> modules need source reduction 0.2230 at 24 months but **none** at 36 months or
> beyond, because the structural ceiling rises above the target on its own. Time and
> source control are substitutes; this audit fixed time while sweeping potency.
> ATLAS-2 exploits that and reaches **86.12%** clearance probability across the whole
> population, so "no potency makes clearance robust" is withdrawn as a general claim.


**This document retracts a claim made in [`atlas_synthesis.md`](atlas_synthesis.md)
and, upstream of it, in [`residual_gap_resolution.md`](residual_gap_resolution.md).**
It is a self-audit of an invented construct, not evidence. Module strengths and the
distributions placed over them are assumptions.

```bash
python run_solutions.py audit      # -> results/audit/ (7 artefacts)
```

---

## 1. Why I ran this

Earlier in this study a remuscularisation "pass" was retracted. A two-factor graft
combination had a nominal gain of +20.3 LVEF points, which looked like a success
until the factors were treated as uncertain: the probability of actually reaching
+20 points was **0.06%**. The lesson recorded was that a nominal product is not an
expected product, because E[∏ U(1, m)] ≪ ∏ m.

I then applied that lesson to the graft arm of ATLAS-1 — and **failed to apply it to
the clearance arm**. The resolution step solved for a potency multiplier of 1.284×
that reaches 80.01% substrate removal at 24 months, and ATLAS-1 scored that as a
burden-reversal pass at 0.8001. That number is a point estimate evaluated at central
module strengths. It had never been given the test that sank the graft claim.

It has now been given that test. It fails.

## 2. Result at the claimed operating point

Each of the eight clearance modules delivers somewhere between nothing and its
assumed central strength, drawn uniformly — the same convention already used for
the graft factors, so both arms are judged by one rule. Source reduction is drawn
on [0, 0.50], because it is no better validated than the modules are.

| Quantity | Value |
|---|---|
| Nominal clearance (the number I reported) | **80.01%** |
| Median clearance under uncertainty | **1.83%** |
| 5th–95th percentile | −38.3% to 41.7% |
| P(clearance ≥ 80%) | **0.00%** |
| Robust at 80%? | **no** |

The nominal pass does not survive. **Burden reversal is retracted.** The negative
5th percentile is not an error: with weak delivery the lesion grows, which is the
same behaviour that made efferocytosis-only therapy harmful in REVERT-1.

## 3. What potency would actually be required

| Potency | Median clearance | P(≥80%) | Robust |
|---|---|---|---|
| **1.284×** (claimed) | 1.83% | **0.00%** | no |
| 2× | 29.8% | 0.36% | no |
| 3× | 50.6% | 4.12% | no |
| 5× | 67.5% | 14.2% | no |
| 10× | 76.8% | 30.4% | no |
| 20× | 78.8% | 41.5% | no |
| 50× | 79.9% | 49.2% | no |
| → ∞ | — | **55.4%** | no |

There is a hard ceiling, and it is structural rather than numerical. Even with
*perfect* modules, reaching 80% at 24 months requires source reduction of at least
**0.2230**. If source reduction is uniform on [0, 0.50], the probability of clearing
that bar is (0.50 − 0.2230)/0.50 = **55.4%**. No potency can beat it, because
potency does not touch the source term.

If source reduction were *guaranteed* at the precedent ceiling of 0.50, robustness
would require roughly **17.7× potency** — not 1.284×. That is more than thirteen
times the escalation the resolution step claimed, and it is an assumption about
eight modules performing simultaneously, not a specification anyone can build.

Longer horizons do not rescue it either: at 240 months the median is 20.0% and
P(≥80%) is 1.1%.

## 4. A second error: the two arms are not independent

This one is subtler and matters more for ATLAS-1's actual claim.

`bed_quality` is derived from achieved clearance. `translation_floor` is derived
from bed quality. So when I scored the sequencing idea — clear the bed, then graft
into it — I evaluated the graft against a bed that was *nominally* 80% cleared.
That silently assumed phase 1 works.

Scored jointly, drawing phase 1 and phase 2 together:

| Bed coupling | Median bed quality | Coupled P(+20) | Decoupled P(+20) | Overstated by |
|---|---|---|---|---|
| 0.00 | 0.023 | 27.3% | 27.3% | **+0.00 pp** |
| 0.25 | 0.023 | 28.1% | 33.7% | +5.60 pp |
| 0.50 | 0.023 | 28.8% | 39.2% | +10.40 pp |
| 0.75 | 0.023 | 29.3% | 45.1% | +15.80 pp |
| 1.00 | 0.023 | 30.3% | 50.6% | **+20.28 pp** |

The comparison is paired — identical factor draws and identical penalty quantiles
pushed through the uncertain floor and the nominal floor — so at zero coupling the
difference is exactly 0.00 pp by construction, and the contrast at higher coupling
is the effect of bed uncertainty alone, with no Monte Carlo noise in it.

Median bed quality is **0.023**, not 1.0. Phase 1 usually fails to prepare the bed,
so there is no prepared bed to graft into. **Sequencing cannot prepare a bed that
phase 1 fails to clear.** The idea is sound in structure and close to worthless in
practice until clearance itself is robust. Its real contribution is a few
percentage points, not the +24 pp I reported.

## 5. Corrected scorecard

| Criterion | Nominal score | Nominal verdict | P(target) | Robust |
|---|---|---|---|---|
| Burden reversal | 0.800 | pass | **0.00%** | **no** |
| Myocardial restoration | 0.609 | fail | 28.8% | no |

`atlas.score_atlas` has been changed to score both invented arms as probabilities
rather than nominal point estimates. The consequence is stark and correct:

- Best criteria met on evidence alone: **2 of 5**.
- Best criteria met if the hypothesis held: **2 of 5**.

**The invented arms add no criterion that survives distributional scoring.** The
claim that ATLAS-1 reaches 3 of 5 where every predecessor reached 2 is withdrawn.
The retracted nominal values are retained in the scorecard as explicitly labelled
`retracted_nominal_*` columns so the audit trail is visible rather than erased.

## 6. What still stands, and what the best candidate is now

Withdrawn:
- Clearance is a solved engineering specification at 1.284× potency.
- ATLAS-1 advances the criterion count to 3 of 5.
- Sequencing is worth up to +24 pp of graft success probability.

Still standing:
- The **structural** clearance ceiling of 83.5% with an unaddressed source term,
  and source reduction as the only lever that moves it. That result is about the
  model's asymptote, not about delivery, and is unaffected.
- **Coupling beats potency** in REVERT-1: efferocytosis without matched efflux is
  worse than nothing (−102 pp), and the gate is the highest-correlated parameter
  (ρ = +0.555).
- **Nominal products are not expected products** — now demonstrated twice, the
  second time against my own candidate.
- The **evidence-maturity finding**: strict credit closes 2 populations, surrogate
  6, exploratory 11; one PCSK9-editing outcomes trial unlocks 17 of 19 with no new
  biology. Nothing in this audit touches that, because it concerns real trial
  anchors rather than invented modules.
- Every ATLAS-1 arm is load-bearing under ablation, and its new idea carries an
  explicit falsifier.

**Best candidate: still ATLAS-1, but on construction quality only.** It is the
best-assembled construct in the study — it reuses only parts that survived their
own ablations, and it states its own falsifier. It is *not* better than its
predecessors on any robust criterion: all of them sit at 2 of 5. There is no
candidate in this study that is a cure, and after this audit there is no candidate
that is even meaningfully ahead.

**Corrected ranking of open problems:**

1. **Clearance under delivery uncertainty** — newly promoted to first. It was
   mistaken for solved; it needs either a guaranteed source reduction plus ~18×
   joint module potency, or a different mechanism entirely.
2. **Remuscularisation** — 28.8% jointly scored, still a coin-flip at best.
3. **Randomised proof** — zero by construction, and no simulation can supply it.

The highest-leverage real-world action is unchanged, and it is the one item in this
study that needs no new biology and no assumed module: **the PCSK9-editing outcomes
trial.**

## 7. Method note

The uniform-on-[0, central] delivery rule is a choice, not a measurement. A more
optimistic rule (say, Beta concentrated near the central strength) would raise the
probabilities. Two things make the conclusion robust to that choice rather than an
artefact of pessimism:

- The **asymptotic ceiling of 55.4%** depends only on the source term and the
  perfect-module threshold, not on the module delivery distribution at all. Even
  infinitely potent modules cannot make clearance robust while source reduction is
  uncertain.
- The **paired coupled-vs-decoupled contrast** is exact at zero coupling and grows
  monotonically, so the direction of the second correction does not depend on the
  distribution either — only its magnitude does.

---

## Reproducibility

- 712 tests, 98.54% coverage on `ihd_solutions`; `ruff` clean, 54 files formatted.
- 25 audit tests pin the retraction, the unreachable-robustness result, the exact
  zero-coupling paired contrast, and the corrected scorecard.
- 6 integrity guards keep audit artefacts out of the evidence manifest and assert
  that the committed ATLAS scorecard no longer claims the retracted pass.
- 33 evidence artefacts SHA-256 verified, 0 mismatches.
- Monte Carlo uses common random numbers across potency sweeps and paired draws
  for the coupled contrast; seeded and reproducible.
