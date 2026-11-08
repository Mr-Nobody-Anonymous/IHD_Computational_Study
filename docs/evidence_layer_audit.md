# Evidence-layer audit: the inherited scores were partly nominal too

**Scenario simulation, not clinical prediction.** This audit propagates the
confidence intervals reported by the source trials, plus adherence and
baseline-rate uncertainty. It does not establish that any mechanism transports to
a new population, and it declares no cures.

```bash
python run_solutions.py audit-evidence     # -> results/evidence_audit/ (6 artefacts)
```

---

## 1. Why this audit existed to be run

[`atlas_self_audit.md`](atlas_self_audit.md) caught a nominal pass in an *invented*
construct. The obvious follow-up question is whether the layer that **every**
construct inherits has the same defect. It does.

`composite.score_criteria` computes `risk_normalisation` from a single
point-estimate simulation:

```python
normalisation = max(0.0, 1.0 - (treated - reference_risk) / (baseline - reference_risk))
```

Meanwhile `uncertainty.monte_carlo_regimen` already propagates LDL efficacy
intervals, hazard-ratio intervals, baseline-rate intervals and adherence. The
machinery to make this criterion distributional was already in the repository. It
was simply never connected to the criterion. That is the whole defect.

The fix was to record per-draw baseline and treated risk in `MonteCarloResult` and
expose `risk_normalisation_draws()`, using exactly the same formula so the nominal
score and the distribution are the same quantity.

## 2. Only one criterion can meaningfully be audited

Rescoring all five under confidence intervals would be meaningless, because four
of them are not effect estimates:

| Criterion | Kind | Audited? | Why |
|---|---|---|---|
| Durable driver control | Structural flag | no | Mechanism presence, not an effect size. **Presence is not evidence of lifetime durability.** |
| Burden reversal | Capped policy credit | no | Ceiling 0.60, below the 0.80 bar regardless of uncertainty |
| Myocardial restoration | Capped policy credit | no | Ceiling 0.50, below the bar regardless |
| **Risk normalisation** | **Effect estimate** | **yes** | Derived from trial hazard ratios and efficacy intervals |
| Randomised proof | Structural zero | no | An untrialled combination earns zero regardless |

This scoping matters. It would have been easy to "audit everything" and produce a
sweeping retraction that was mostly an artefact of applying interval arithmetic to
policy constants.

## 3. Result across all 19 populations × 3 credit modes

| Credit mode | Nominal risk passes | Robust risk passes | Retracted | Best criteria (nominal → robust) |
|---|---|---|---|---|
| Strict | 1 | **0** | 1 | **1 → 0** |
| Surrogate | 7 | **2** | 4 | 2 → 2 |
| Exploratory | 19 | **18** | 1 | 2 → 2 |

The effect is strongly mode-dependent, and that is the interesting part.
Exploratory credit survives almost intact because its effects are large relative
to their intervals. Strict credit does not survive at all.

### The six passes that do not survive

| Population | Mode | Nominal | Median | P(≥0.8) |
|---|---|---|---|---|
| stable_ccs | **strict** | 0.819 | 0.749 | **11.7%** |
| any_ascvd | surrogate | 0.806 | 0.752 | 5.0% |
| post_mi_hscrp_high | surrogate | 0.804 | 0.754 | 5.3% |
| ascvd_obesity_no_diabetes | surrogate | 0.819 | 0.764 | 13.0% |
| secondary_post_mi | surrogate | 0.856 | 0.792 | 43.7% |
| lifetime_exposure | exploratory | 1.000 | 0.500 | 50.0% |

**The `stable_ccs` strict row is the consequential one.** Under strict credit, risk
normalisation was the *only* criterion any population met — the composite's
headline "strict 1/5". At 11.7% probability it does not survive, so **the best
strict result in this study falls to zero of five.**

### One case is too close to call

`statin_intolerant_ascvd` (surrogate) lands at **78.3%**, within two Monte Carlo
standard errors of the 80% bar. It is reported as **borderline**, not as a
retraction. Calling a 78% estimate a confident failure against an 80% bar would
repeat the overclaiming this audit exists to correct, just in the opposite
direction. More draws would settle it.

## 4. A median is not a probability either

`lifetime_exposure` under exploratory credit has a nominal score of **1.000**. Its
distribution is not concentrated near 1 — it is **bimodal**:

- ~50% of draws at exactly 0
- ~50% of draws at exactly 1
- P(≥0.8) = **50%**

Because the score is clipped to [0, 1], mass piles up at both poles with little in
between. The median then sits at whichever pole happens to hold slightly more than
half the draws, and **it flips between 0.000 and 1.000 as the sample size changes**
— I observed both across runs at 200, 300 and 400 draws.

So the median is not a summary of this distribution at all. The audit reports
`mass_at_floor`, `mass_at_ceiling` and a `median_uninformative` flag rather than
letting a median of 1.000 stand in for a coin flip. This is the same lesson as the
nominal-product trap, arriving through a different door.

## 5. What this changes, and what it does not

**Changed:**
- The composite's strict-mode "1 of 5" is retracted → **0 of 5**.
- Four surrogate-mode risk passes and one exploratory pass are retracted.
- `risk_normalisation` is now reported as a probability wherever this audit runs.

**Not changed:**
- Exploratory credit still supports 2 of 5 in 18 of 19 populations.
- The evidence-maturity finding (strict closes 2 populations, surrogate 6,
  exploratory 11; one PCSK9-editing outcomes trial unlocks 17 of 19). That counts
  *credit-mode closure*, not criterion scores under uncertainty, and is a separate
  quantity.
- The REVERT-1 coupling result, the structural clearance ceiling, and the
  remuscularisation retraction.

**Combined with the ATLAS-1 self-audit:** no construct in this study meets more
than two of five criteria robustly, and under strict credit none meets any. The
audits do not rescue a candidate — they lower the ceiling every candidate inherits.

## 6. The pattern worth naming

Three times now, a number that looked like a pass turned out to be a point
estimate standing in for a distribution:

1. A two-factor graft with nominal +20.3 pp → **0.06%** probability.
2. A clearance operating point at nominal 80.01% → **0.00%** probability.
3. An evidence-derived risk score of 0.819 → **11.7%** probability.

The first was caught in an invented construct, the second in my own best
candidate, the third in the layer I had been treating as the trustworthy baseline.
The lesson is not that any particular model was wrong. It is that **every criterion
score in this study must be a probability, and "is this a point estimate?" is a
question to ask of each one individually** — including the ones that come from
real trials, because a trial's point estimate is also a draw.

The recommendation is unchanged, and it is the one item that does not depend on
any of these scores: **a PCSK9-editing outcomes trial** would replace assumed
durability and surrogate credit with a randomised effect estimate — precisely the
deficiency all three audits expose.

---

## Reproducibility

- 741 tests, 98.50% coverage on `ihd_solutions`; `ruff` clean, 56 files formatted.
- 25 audit tests pin the strict retraction, the bimodality diagnostic, the
  borderline exclusion, the criterion scoping, and the formula match between the
  point estimate and the distribution.
- 4 integrity guards keep these artefacts out of the evidence run manifest and
  assert the committed strict retraction.
- 33 evidence artefacts SHA-256 verified, 0 mismatches.
- 300 draws per population × mode, seeded; borderline cases carry their own
  standard error so the reader can see what more sampling would resolve.
