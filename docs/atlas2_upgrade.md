# ATLAS-2: upgrading the parts the audits broke

**Invented, untrialled construct.** The upgrades are architectural hypotheses with
clinical analogues, not demonstrated effects. Every upgrade is ablated and charged
for. This declares no cures.

```bash
python run_solutions.py upgrade      # -> results/upgrade/ (10 artefacts)
```

---

## 1. Diagnosis before prescription

Three audits had established *that* the clearance arm fails under delivery
uncertainty. They had not established *why*, and "why" is what makes a fix possible
rather than a bigger assumed number. Reading the model gave two structural causes:

**Cause 1 — a series bottleneck.** `ClearanceParameters.delivered()` multiplies
**six of eight** module strengths by the single `gate_specificity` draw:

```python
factor = self.gate_specificity
crystal_solubilisation = self.crystal_solubilisation * factor
efflux_gain            = self.efflux_gain * factor
efferocytosis_gain     = self.efferocytosis_gain * factor
influx_block           = self.influx_block * factor
resolution             = self.resolution * factor
nucleation_block       = self.nucleation_block * factor
```

One uncertain parameter throttles almost the whole system, with no redundancy. A
bad gate draw cannot be compensated by any other module.

**Cause 2 — open-loop delivery.** The model picked a potency multiplier and hoped
it landed. Nobody manages lipids that way. Clinical practice measures and titrates
to a target, which is precisely why achieved LDL is far more reliable than
fixed-dose response.

Both are architecture defects, not potency shortages. So ATLAS-2 changes the
architecture and leaves the potency assumption honest.

## 2. The three upgrades

| | Upgrade | Fixes | Real-world analogue |
|---|---|---|---|
| **U1** | Redundant gate and node paths | series bottleneck | several mechanistically distinct lesion-targeting routes; the node works if any one lands |
| **U2** | Closed-loop titration to target | open-loop delivery | treat-to-target titration, the standard of care in lipidology |
| **U3** | Verified source control (run-in) | uncertain source term | respond-and-continue enrichment — **available, but not needed at the solved horizon** |

Redundancy combines by **max**, not by complement product: if two mechanisms both
partially work, counting only the better one is the conservative choice. Titration
escalates within a hard cap and records the potency actually used, so it cannot
quietly exceed its budget.

## 3. Result

| | P(clearance ≥ 80%) | Median clearance |
|---|---|---|
| ATLAS-1 baseline (upgrades off, 24 mo) | **0.00%** | 1.58% |
| **ATLAS-2 (upgrades on, 60 mo)** | **86.12%** | **81.16%** |

- Robust at 80%: **yes**
- Robust with a 2-standard-error margin: **yes** — lower bound **83.68%**
- Eligible fraction: **100%** — nobody is excluded to buy the result
- Unconditional success: **86.12%**, identical to the conditional number

That last line is the one that matters. An earlier version of this work used a
run-in gate and reached 83% *conditional* on eligibility — but only 46.4%
unconditional, because it excluded 44% of patients. That is not a fix, it is a
choice of whom to report on. The shipped configuration excludes no one.

## 4. What it costs

| Cost | Value |
|---|---|
| Potency cap required | **2.80×** |
| Precedented operating point | 1.284× |
| **Precedent gap** | **2.18×** |
| Distinct targeting routes | 3 (3 is the precedented limit) |
| Horizon | **60 months** instead of 24 |

Bisection puts the *minimum* robust cap at 2.590×, but that lands at 80.27% ± 2.1 pp
— sitting on the bar, which is exactly the borderline case the evidence-layer audit
refused to count. The shipped cap of 2.80× buys a real margin instead.

**The remaining gap is potency and time, not architecture.** Getting every clearance
module to 2.18× the precedented potency simultaneously is a strong assumption, and
it is the specific thing that would have to be demonstrated. The architecture is
now sound; the potency is still a hypothesis.

## 5. Ablation at the solved point

| Upgrade removed | P(target) | Loss | Verdict |
|---|---|---|---|
| **Titration** | 6.50% | **+72.83 pp** | load-bearing |
| Redundant nodes | 37.33% | +42.00 pp | load-bearing |
| Redundant gate | 56.33% | +23.00 pp | load-bearing |
| Verified source | 79.33% | 0.00 pp | **not used at this horizon** |

**Titration is the single largest contributor by a wide margin.** Treating to a
measured target beats adding potency — the same lesson as the AND gate in REVERT-1,
where coupling beat raw strength (ρ = +0.555). Control structure is worth more than
amplitude.

Redundancy saturates usefully inside precedent:

| Gate routes | P(target) | Status |
|---|---|---|
| 1 | 29.25% | precedented |
| 2 | 77.50% | precedented |
| 3 | **80.00%** | precedented |
| 4 | 82.75% | beyond precedent |

Three routes suffice. The fix does not depend on an indefensible number of
mechanisms.

## 6. A correction to the earlier audit

The ATLAS-1 self-audit reported a **55.4% probability ceiling at any potency** while
the source term is uncertain, and concluded that no potency could make clearance
robust. That was correct at 24 months and **wrong as a general claim**. Recomputing
the source requirement by horizon:

| Horizon | Source reduction perfect modules need |
|---|---|
| 24 months | **0.2230** |
| 36 months | **0.0000** |
| 60 months | **0.0000** |
| 120 months | **0.0000** |

At 24 months the structural ceiling with an unaddressed source is 74.9%, below the
80% target, so source reduction is genuinely binding. By 36 months the ceiling has
risen above the target on its own, and the source lever **stops binding entirely**.

So the "impossible at any potency" verdict was an artefact of the 24-month horizon.
Time and source control are substitutes, and the earlier analysis had fixed time
while sweeping potency. That is why U3 is switched off in the solved configuration:
a run-in at 60 months would only exclude patients who would have succeeded.

I verified this the honest way — by checking that a *tighter* gate raises the
conditional number while **lowering** the unconditional one, and that the
unconditional probability is invariant to the gate threshold when the gate is not
binding. A conditioning step cannot create success.

## 7. Does the sequencing idea finally pay off?

ATLAS-1's one new idea — clear the bed, then graft into it — collapsed in the
self-audit because phase 1 usually failed, leaving median bed quality at **0.023**.
Scored against the bed the *upgraded* phase 1 actually produces:

| Bed coupling | Median bed quality | Graft P(+20 pp) | Sequencing gain |
|---|---|---|---|
| 0.00 | 1.000 | 27.9% | **+0.0 pp** |
| 0.50 | 1.000 | 39.0% | **+11.1 pp** |
| 1.00 | 1.000 | 52.0% | **+24.1 pp** |

**The idea was never structurally wrong. It was starved by a failing phase 1.** With
the bed genuinely prepared it delivers essentially the gain originally claimed —
and the zero-coupling falsifier still returns exactly 0.0 pp, so the gain is still
attributable to the assumption rather than baked into the baseline.

The graft still is not robust (52.0% at perfect coupling). So remuscularisation is
now the **sole** binding biology problem, rather than one of two.

## 8. Verdict

**The clearance arm is fixed as a specification**: robust across the whole
population, with statistical margin, no eligibility exclusion, and every active
upgrade load-bearing under ablation. The price is stated: 2.18× beyond-precedent
potency and a 60-month horizon.

**Still not fixed:**
- **Myocardial restoration** — 52.0% at best, not robust. This is now the binding
  biology gap.
- **Randomised proof** — zero. None of this has been trialled, and no simulation
  can supply it.
- **The potency assumption itself** — 2.18× across every module at once is the
  thing to demonstrate, not a result.

**Best candidate is now ATLAS-2**, and unlike ATLAS-1 it earns that on a robust
criterion rather than on construction quality alone. It is still not a cure.

The real-world recommendation is unchanged: the **PCSK9-editing outcomes trial**
needs no new biology and unlocks 17 of 19 populations. ATLAS-2 sharpened what would
have to be true for plaque clearance; it did not make the evidence appear.

---

## Reproducibility

- 783 tests, 98.56% coverage on `ihd_solutions`; `ruff` clean, 58 files formatted.
- 39 upgrade tests pin the solve with margin, the 100% eligibility, the horizon
  correction, the precedent and redundancy limits, the gate-invariance property,
  and the recovered sequencing gain with its zero-coupling falsifier.
- 4 integrity guards keep these artefacts out of the evidence run manifest and
  assert the committed solve has margin and excludes nobody.
- 33 evidence artefacts SHA-256 verified, 0 mismatches.
- Seeded Monte Carlo; the solved cap is found by bisection, not chosen, and the
  shipped cap adds margin above it.
