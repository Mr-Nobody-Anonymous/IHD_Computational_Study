# ATLAS-3 — the graft arm given a variance architecture

> **SUPERSEDED IN PART BY [`atlas4_precedent_search.md`](atlas4_precedent_search.md).**
> A literature search found three errors here, and ATLAS-4 corrects them:
> 1. **The precedent limit was wrong — now corrected in the code, not just here.**
>    This document says two doses is the precedented maximum and calls three "beyond
>    precedent". CATO (NCT06145035) and HOPE-2 both give **four** doses in humans. Three
>    was never beyond precedent. `MAX_PRECEDENTED_ATTEMPTS` now reads
>    `max_precedented_doses()` from the registry (**4**), so `ATLAS3.beyond_precedent`
>    is `False`. `precedented_option()` was also rewritten: it used to return the limit
>    itself, which after the correction would have compared three attempts against a
>    strictly better four and called the difference a cost. It now returns the
>    **cheapest robust** design within the limit — **2 attempts, P = 82.50%, adverse
>    7.42%** — against ATLAS-3's three at **93.1%** for **+0.95 pp** more risk.
> 2. **The retry operator was wrong.** Best-of-k is used below. Every published
>    repeat-dose study reports **cumulative** increments instead, and shows that
>    splitting a fixed total dose beats giving it at once — which best-of-k cannot
>    represent.
> 3. **The mechanism label was wrong.** These cells do not regenerate myocardium; the
>    effect is paracrine. The LVEF ceiling below is functional, not regenerated tissue.
>
> A consequence: the "measure independence above ≈0.71" specification in §4 is a
> property of the wrong operator. Under cumulative dosing ~91% of the gain survives
> full correlation. The numbers below are retained as the historical record.
>
> **SUPERSEDED FURTHER BY [`mechanism_recombination.md`](mechanism_recombination.md).**
> ATLAS-5 keeps this arm's retry redundancy, replaces its operator, drops its feedback
> loop and reporter cassette as measured-inert, and adds graft maturation.

> **ATLAS-3 IS AN INVENTED, UNTRIALLED CONSTRUCT.** Staged repeatable grafting with a
> measurable feedback signal does not exist. Nothing in this document is a therapy, a
> result, or a recommendation. Randomised proof remains **zero**. Every table in
> `results/atlas3/` carries `status = atlas3_invented_architecture_unvalidated`.

Reproduce:

```bash
python run_solutions.py atlas3 --draws 4000
```

Outputs: 10 artefacts in `results/atlas3/` (8 tables, 1 report, `DISCLAIMER.txt`).
Tests: `tests/test_atlas3.py` (58).

---

## 1. Why this study exists

ATLAS-2 closed the clearance gap, and it did so **without new potency**. It won by
attacking variance: redundant nodes, a conditioning gate, and dose titration. The
ablation showed titration alone was worth **+72.8 pp**.

One biology gap was left open: remuscularisation, the **+20 pp LVEF** criterion. It had
never been given the same treatment. It was still modelled as a single, irreversible,
unmeasured shot. ATLAS-3 applies the ATLAS-2 architecture to it and charges the bill.

## 2. The diagnosis: this was never a capability problem

| Quantity | Value |
|---|---|
| Required delivered multiplier | **5.128×** |
| Product of all four published maxima | **18.200×** |
| Headroom over the requirement | **3.55×** |
| P(target), one graft, **perfect** bed (zero chronic-injury penalty) | **49.60%** |

The four factors are the ones already catalogued in
`results/gaps/tables/gap_remuscularisation_factors.csv`: dose escalation 1.75× (human
construct, animal efficacy), endothelial co-delivery 2.0× (NHP), epicardial stromal
support 2.6× (rodent), induced proliferation 2.0× (rodent).

Even granting every one of them at its published maximum, and even granting a
**perfect** recipient bed, a single graft is a coin flip. The reason is the one this
repository has hit repeatedly: **the nominal product is not the expected product.**
Each factor is uncertain across its published range, and the product of four uncertain
factors has a long left tail. `E[∏] ≪ ∏max`.

A single bigger mechanism is the *wrong* response. Replacing the four factors with one
drawn `U(1, 8)` — a hypothetical mechanism 8× better than anything published — yields
**40.9%**, *worse* than the four modest factors combined. Variance, not ceiling, is the
binding constraint.

## 3. What was added, and what it costs

Three levers, each independently ablatable:

1. **Retry redundancy** — up to *k* staged attempts; the patient keeps the best outcome.
   Modelled as a **maximum, never as extra treated area**. Attempts stop the moment the
   target is met, so a first-attempt responder is never exposed to a second procedure.
2. **Independence** — the degree to which attempts fail for *different* reasons. Swept,
   never assumed (see §4).
3. **Feedback** — a measured engraftment signal that guides the retry, gated by
   `measurement_quality`.

### Result (4000 draws, seed 20261107)

| Construct | P(target) | Median gain | Cumulative adverse risk | Risk-adjusted | Robust w/ margin |
|---|---|---|---|---|---|
| Single graft (today) | **49.60%** | 19.92 pp | 5.00% | 47.12% | no |
| **2 attempts — within precedent** | **82.50%** | 25.00 pp | 7.42% | 76.37% | **yes** |
| ATLAS-3, 3 attempts | **93.10%** | 25.00 pp | 8.38% | 85.30% | yes |

Every extra attempt is charged an adverse-event probability (default 5% per attempt,
swept 0–20%), and success is reported **both raw and discounted by that risk**. A
higher success rate bought with more procedures is not free.

### The third attempt buys headroom, not feasibility

**Two attempts already clear the bar**, and repeat dosing has precedent
(`MAX_PRECEDENTED_ATTEMPTS = 2`). This was verified across 12 independent seeds: two
attempts gives P = 0.8216 on average (range 0.8117–0.8300), robust with margin in 11/12.

So **the number of procedures is not the invention, and it is not what breaks
precedent.** The invention is the **closed loop** — a measured engraftment signal that
decides *whether* to repeat. The third attempt is optional headroom: it lifts 82.5% to
93.1% for **+0.95 pp** of additional adverse-event risk. The two-attempt margin is thin
(2.5 pp), so it is reported rather than relied upon.

> **A methodological correction made during this study.** The first implementation
> reported a minimum of **3** attempts. That was an artefact: at small draw counts the
> 2-SE margin is wide enough to reject a two-attempt design that is genuinely above the
> bar. Both specification searches (`minimum_robust_attempts`,
> `maximum_tolerable_correlation`) now floor their draw count at
> `MIN_SPECIFICATION_DRAWS = 4000`. **A specification must not move when the sampling
> budget moves.** A test pins this.

## 4. Falsifier one — redundancy only counts if attempts are independent

The contrast is **paired**: the same run returns both the best-of-*k* gain and the
first attempt's gain, so at correlation 1 the difference is **exactly zero**, not zero
plus Monte Carlo noise.

| Failure correlation | P(target) | Gain vs. first attempt |
|---|---|---|
| 0.00 | 88.48% | +38.88 pp |
| 0.25 | 89.30% | +35.05 pp |
| 0.50 | 81.38% | +25.08 pp |
| 0.75 | 64.95% | +11.10 pp |
| **1.00** | 50.95% | **+0.00 pp (exact)** |

Three attempts that share a failure mode **are one attempt**.

**Maximum tolerable correlation ≈ 0.71** (spread ±0.02 across seeds). Retry outcomes
must be at least this independent for the design to stay robust. That is the single
most important quantity to *measure* in a real study — and a chronic, poorly
vascularised scar is precisely the kind of shared failure mode that would violate it.

An earlier version of this falsifier used independent random streams and read
**+0.1 pp** at correlation 1. That was noise masquerading as a result, and it was fixed
the same way `coupled_graft_outcome` was fixed in the self-audit.

## 5. Falsifier two — feedback only counts if it can learn

| Feedback learning | P(target) | Gain |
|---|---|---|
| 0.00 | 88.48% | **+0.00 pp (exact)** |
| 0.25 | 91.40% | +2.93 pp |
| 0.50 | 93.10% | +4.63 pp |
| 0.75 | 94.40% | +5.93 pp |
| 1.00 | 95.53% | +7.05 pp |

Feedback is worth having but is **not** the main effect. Redundancy is.

> **A second correction.** The first feedback model raised each factor's lower bound.
> That made the product of the *minima* exceed the requirement, forcing **P = 100%** —
> the nominal-product trap in a new disguise. It was replaced with **one bounded
> aggregate uplift** (`FEEDBACK_MAX_UPLIFT = 0.5`), capped at the product of maxima and
> gated by `measurement_quality`. A clinical signal reports *total* engraftment; it does
> not attribute outcome to individual mechanisms. **Any P = 100% in this model is an
> artefact, not a result** — a test asserts none is reachable.

## 6. Lever ablation — every lever is load-bearing

| Removed lever | P(target) | Loss |
|---|---|---|
| Retry redundancy | 49.60% | **−43.50 pp** |
| Independence | 66.50% | −26.60 pp |
| Feedback | 88.48% | −4.63 pp |
| Measurement | 88.48% | −4.63 pp |

Feedback and measurement ablate identically by construction: with no measurement there
is nothing to learn from, so the two paths must coincide. A test pins that equality.

## 7. Both arms together

The joint and VOI calculations run at half the headline draw budget (2000 draws at
`--draws 4000`), so the graft reads 92.45% here against 93.10% in the headline. Both
are the same construct; the difference is Monte Carlo error.

| Arm | P(target) |
|---|---|
| Clearance (ATLAS-2, 60 months) | 83.50% |
| Graft (ATLAS-3) | 92.45% |
| **Joint, assuming independence** | **77.20%** |
| Randomised proof | **0.00%** |

The joint is a product, which **assumes the arms are independent**. They are not fully
independent in reality — the graft lands in the bed that clearance prepared, and the
self-audit measured that coupling as worth up to **+24.1 pp** to the graft. So 77.20%
is an *optimistic* joint, labelled as such. Even taken at face value it is one number
covering **two of five** criteria, with proof still absent.

## 8. Value of information, recomputed

The original `results/discovery/tables/discovery_value_of_information.csv` ranked
substrate removal first (0.1839) — but it was computed **before ATLAS-2 existed** and is
stale. Asked again of the current architectures:

| Capability | P now | VOI gain | Architecture |
|---|---|---|---|
| **substrate_removal** | 83.50% | **+15.25 pp** | ATLAS-2 (redundancy + titration, 60 mo) |
| lvef_gain_pp | 92.45% | +6.30 pp | ATLAS-3 (retry + feedback) |

**The ranking did not flip, but the reason changed.** Substrate removal is still the
higher-value uncertainty — not because it lacks an architecture (it now has one), but
because its architecture delivers a *lower* probability. Resolving the graft uncertainty
would buy less, because ATLAS-3 already pushes it higher.

## 9. What must be demonstrated before any of this is real

1. **Independence of retry outcomes** above ≈0.71. Correlated failures erase the entire
   gain, and a chronic poorly-vascularised scar is exactly a shared failure mode.
2. **A measurable engraftment signal** good enough to guide a retry decision.
3. **Acceptable cumulative procedural risk.** At 10% per-attempt risk the risk-adjusted
   three-attempt success falls to 77.39%; at 20% it falls to 63.43% and the design stops
   being worth its price.
4. **Human evidence.** All four factors are animal- or NHP-derived. The translation
   penalty is swept (`translation_floor`), never assumed away — the headline runs at the
   *optimistic* end of that sweep.

## 10. Verdict

The graft gap closes **in simulation** under staged retry with feedback. This locates
the missing step precisely: a remuscularisation modality that can be **repeated and
measured**. Repeat dosing already has precedent; the missing piece is the feedback
signal that decides whether to repeat. It is **not more potency** — the published
factors already exceed the requirement 3.55-fold.

**ATLAS-3 is a specification of the experiment to run, not a cure and not a therapy.**
The five-criterion cure test still returns zero cures, and no simulation in this
repository can return anything else: `randomised_proof` is unreachable in silico by
definition. The first real-world action remains the one identified by the evidence
lever analysis — a **PCSK9-editing outcomes trial**, which unlocks 17 of 19 populations.

---

## See also

- [`atlas2_upgrade.md`](atlas2_upgrade.md) — the clearance arm, where this architecture was first validated
- [`atlas_self_audit.md`](atlas_self_audit.md) — nominal passes, the coupling overstatement, horizon specificity
- [`evidence_layer_audit.md`](evidence_layer_audit.md) — six evidence retractions
- [`FULL_REPORT.md`](FULL_REPORT.md) — the complete inventory of agents, compounds, experiments and results
- [`scientific_methodology_and_limitations.md`](scientific_methodology_and_limitations.md)
