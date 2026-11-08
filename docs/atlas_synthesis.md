# ATLAS-1: synthesising the best working part of every construct

> **PARTIALLY RETRACTED — read [`atlas_self_audit.md`](atlas_self_audit.md) first.**
> A later self-audit found that the clearance operating point this document relies
> on was a *nominal* pass: under delivery uncertainty P(clearance ≥ 80%) is **0.00%**,
> not the 80.01% reported below. Because bed quality is derived from clearance, the
> sequencing gain was also overstated by ~20 pp. **ATLAS-1 reaches 2 of 5 criteria,
> not 3.** Sections 4–6 below are corrected in place; the original numbers are kept
> visible and labelled so the error is auditable rather than erased.

**Status: invented, untrialled construct.** ATLAS-1 is not a therapy, a clinical
result, or a recommendation. Its invented arms carry zero evidence credit, its
randomised-proof score is zero by construction, and its one genuinely new idea is
an assumption that is swept to zero rather than asserted.

Run it with:

```bash
python run_solutions.py atlas
```

Outputs land in `results/atlas/` (9 artefacts), deliberately outside the
evidence-graded `results/solutions/` tree and absent from its SHA-256 manifest.

---

## 1. Why assemble rather than invent again

Five constructs preceded this one. Each failed, but each failed *differently*, and
the parts that survived their own tests are reusable. ATLAS-1 keeps only the parts
that earned their place under test, and discards the rest.

| Source construct | Part kept | Why it earned its place | Part discarded |
|---|---|---|---|
| Catalogue / cure search | Node-disjoint proven backbone | Real randomised event benefit, deduplicated by mechanism | Single-agent risk normalisation (ceiling 0.096) |
| Composite best-of-breed | Target-overlap dedup; zero proof for untested combinations | Stopped three PCSK9 modalities counting as three mechanisms | Greedy-by-events myocardial blindness |
| REVERT-1 | AND gate; efflux ↔ efferocytosis coupling | Ablation showed gate removal costs more than any addition | Efferocytosis alone (−102.3 pp, actively harmful) |
| Gap solver, module M7 | Source reduction | The only lever that broke the 83.5% structural clearance ceiling | Unbounded source abolition (no precedent exists) |
| Residual-gap resolution | 1.284× potency operating point; distributions, not points | Reached the 80% clearance target for the first time | Nominal-product "passes" (E[∏U] ≪ ∏max) |
| Discovery | Back-calculated capability target | Turned "cure" into five falsifiable numbers | Feasible-volume optimism (demonstrated volume = 0) |

Nothing was kept on the strength of a narrative. Each row above is pinned by a
test that fails if the part stops doing its job.

## 2. The new idea: sequence, do not stack

Every previous construct delivered all of its arms **simultaneously**. That is
also how the remuscularisation evidence was generated badly: the primate work that
shows large grafts treats animals roughly two weeks after infarction, into a
freshly injured bed. A human candidate presents with a mature, lipid-laden,
inflamed scar. The gap between those two settings was previously modelled as a
flat, unexplained "chronic-injury translation floor" of 0.5.

ATLAS-1's hypothesis is that the floor is not a constant of nature but a property
of the bed — and that phase 1 can change it:

- **Phase 1 (months 0–24):** proven backbone + durable editing + gated clearance +
  source reduction. Clearance reaches 80.0%, inflammation is quietened, the lipid
  and necrotic burden falls.
- **Phase 2 (month 24 onward):** graft remuscularisation into the *prepared* bed.

The link between the phases is the single new assumption, `bed_coupling ∈ [0, 1]`.
At `bed_coupling = 0` the construct degenerates exactly to simultaneous delivery,
so the idea can be falsified by inspection.

### The falsification test

| Bed coupling | Translation floor | P(LVEF +20 pp) | Gain vs unsequenced | Robust at 80%? |
|---|---|---|---|---|
| 0.00 | 0.50 | 26.8% | **+0.0 pp** | no |
| 0.25 | 0.62 | 32.4% | +5.5 pp | no |
| 0.50 | 0.75 | 38.7% | +11.9 pp | no |
| 0.75 | 0.88 | 45.2% | +18.4 pp | no |
| 1.00 | 1.00 | 50.9% | **+24.1 pp** | no |

Two readings, both important:

1. **The idea does real work.** At central coupling it lifts the graft from 26.8%
   to 38.7%; at perfect coupling to 50.9%. That is the largest single improvement
   to the weakest arm found anywhere in this study.
2. **The idea does not solve the problem.** Even granting *perfect* bed
   preparation — an assumption no evidence supports — the graft arm reaches
   50.9%, not the 80% required. Sequencing moves remuscularisation from hopeless
   to a coin-flip. It does not make it a capability.

> **Correction.** The table above scores the graft against a *nominally* cleared
> bed, which assumes phase 1 succeeded. It usually does not. Scored jointly, the
> gain at full coupling is a few points rather than +24.1 pp, and the decoupled
> scoring overstated the graft by **20.3 pp**. See the self-audit, section 4.

A test asserts the zero-coupling gain is exactly 0.0. If a future edit ever makes
the unsequenced baseline benefit from sequencing, the assumption has leaked into
the baseline and the test fails.

## 3. Ablation: what is actually load-bearing

Each arm was removed one at a time and the construct re-scored.

| Arm removed | Clearance | Loss | Graft P | Loss | Load-bearing |
|---|---|---|---|---|---|
| Gated clearance | 0.0% | 80.0 pp | 26.0% | 12.6 pp | yes |
| Source reduction | 66.0% | 14.0 pp | 36.0% | 2.7 pp | yes |
| Durable editing | 80.0% | 0.0 pp | 38.6% | 0.0 pp | yes — durability criterion |
| Proven backbone | 80.0% | 0.0 pp | 13.3% | 25.3 pp | yes |
| Remuscularisation | 80.0% | 0.0 pp | 0.0% | 38.6 pp | yes |
| Bed sequencing | 80.0% | 0.0 pp | 26.0% | 12.6 pp | yes |

All six arms carry weight; none is decoration. Two entries deserve comment:

- **The proven backbone matters to the *graft*, not to clearance** (25.3 pp). Its
  contribution is the manufacturing and delivery maturity that makes a
  dose-escalated patch conceivable at all. Removing it collapses the graft arm.
- **Durable editing looks inert on clearance and graft mass, and is not.** It
  carries the durability criterion, which those two metrics cannot see. An
  earlier version of this ablation wrongly labelled it inert; the fix adds an
  explicit durability column and a test pinning it.

## 4. Tested across everything

19 populations × 3 credit modes = **57 configurations**, each scored on all five
cure criteria, with clearance and remuscularisation evaluated alongside. Two
scorecards are kept strictly apart: what the **evidence** supports, and what
would hold **if the hypothesis were true**.

Scores below are **probabilities of reaching each target under delivery
uncertainty**, as corrected by the self-audit. The retracted nominal values are
shown struck through.

| Criterion | Best evidence score | Best hypothesised score | Threshold 0.8 |
|---|---|---|---|
| Durable driver control | 1.000 ✓ | 1.000 ✓ | passes |
| Atherosclerotic burden reversal | 0.600 ✗ | 0.600 ✗ (was ~~0.800 ✓~~) | **blocked** |
| Myocardial restoration | 0.500 ✗ | 0.500 ✗ (was ~~0.708 ✗~~) | **blocked** |
| Risk normalisation | 1.000 ✓ | 1.000 ✓ | passes |
| Randomised proof | 0.000 ✗ | **0.000 ✗** | **blocked in both** |

- Best on evidence alone: **2 of 5**.
- Best if the hypothesis held: **2 of 5** — the invented arms add no criterion that
  survives distributional scoring. ~~3 of 5, the highest of any construct~~ is
  withdrawn.
- Configurations declared a cure: **0**.

The original burden-reversal "pass" landed at 0.8001 only because the potency
multiplier was back-solved to hit that target with no headroom. The self-audit
showed the deeper problem: at that operating point the *median* clearance is 1.83%
and P(≥80%) is 0.00%. No potency makes it robust while the source term is uncertain
(asymptotic ceiling 55.4%).

Myocardial restoration, scored jointly with the bed phase 1 actually produces,
reaches +20 pp with probability 28.8% (median bed quality 0.023, not 1.0).

Proof is zero in both scorecards, and this is not a gap that more simulation can
close. No combination in ATLAS-1 has ever been randomised.

## 5. What has been discovered across the whole study

**Discovery 1 — the binding constraint on paper is evidence, not biology.**
Counting how many populations close as credit is relaxed: **strict 2, surrogate 6,
exploratory 11**, and all 19 close somewhere. One PCSK9-editing outcomes trial
would unlock **17 of 19** populations. It requires no new biology, only the
maturation of an existing randomised programme. It is the highest-leverage item
found in this entire study.

**Discovery 2 — plaque clearance has a structural ceiling, and one lever moves
it** *(the "solved specification" reading of this was retracted; see below)*. Perfect delivery of every clearance module asymptotes at **83.5%**, because
the debris source keeps refilling the lesion. Module potentiation alone cannot
reach 80% at 24 months: the best single module (crystal solubilisation) reaches
63.3%, and even all-perfect modules reach only 82.7%. Reducing the *source* is
what changes the arithmetic. Bounded at 0.35–0.50 by CHOP-deletion precedent,
source reduction plus a **1.284× potency multiplier at 24 months** reaches the
80% target *at central module strengths*. **Retracted as a specification:** under
delivery uncertainty that operating point delivers a median of 1.83% and P(≥80%) of
0.00%. The structural ceiling and the identity of the lever stand; the claim that
the target is thereby solved does not.

**Discovery 3 — efferocytosis alone is harmful, and the gate is the real
invention.** Boosting corpse uptake without matching export capacity kills the
phagocytes and drives clearance to **−102.3 pp** relative to untreated. The AND
gate that couples efflux to efferocytosis is the single highest-correlated
parameter in REVERT-1 (ρ = +0.555). Coupling, not potency, was the missing idea.

**Discovery 4 — nominal products are not expected products.** A two-factor graft
combination with a nominal +20.3 pp gain has a **0.06%** probability of actually
delivering +20 pp, because E[∏ Uniform(1, m)] ≪ ∏ m. This destroyed an earlier
"pass" and forced a retraction. Every claim in ATLAS-1 is distributional.

**Discovery 5 — remuscularisation is the genuinely open problem.** Ranking
inverted over the course of this study. Clearance turned out to be a solvable
specification; muscle replacement did not. All four published enhancement
factors together, at full effect and no chronic-injury penalty, give **50.9%** —
a coin-flip. Under a 0.5 translation floor, **26.5%**. No subset is robust under
any haircut tested.

**Discovery 6 — sequencing is structurally sound and practically contingent.**
~~Preparing the bed before grafting is worth up to +24.1 pp.~~ **Corrected:** that
figure assumed phase 1 clears the bed, which it does with probability ~0. Scored
jointly the gain is a few points. A preparation step is only as valuable as the
preparation is reliable — the arms of a sequenced construct cannot be scored
independently.

**Discovery 7 — the same error recurred in my own candidate.** The nominal-product
trap that sank the graft claim also invalidated the clearance operating point, and
I did not catch it until I audited ATLAS-1 against the standard I had already written
down. Both arms are now scored as probabilities, and the retracted values are
retained as labelled columns rather than deleted.

## 6. Best candidate

**ATLAS-1**, at the operating point: source reduction 0.50, 24-month phase 1,
1.284× clearance potency, graft at month 24, central bed coupling 0.5.

**Corrected basis.** It is the best candidate on *construction quality* only: it
(a) reuses exclusively parts that survived their own ablations, (b) states its own
falsifier and reports the result of applying it, and (c) has had its own headline
claim audited and withdrawn. It is **not** ahead on any robust criterion — it sits
at 2 of 5, exactly where every predecessor sits. The earlier grounds (b) "3 of 5"
and (c) "improves the weakest arm most" are withdrawn.

It is **not a cure**, and the two reasons are independent:

1. **Neither invented arm is robust** — clearance P(≥80%) = 0.00% at the claimed
   operating point, and myocardial restoration 28.8% scored jointly. These are
   biology and delivery gaps.
2. **Randomised proof is zero** — no part of this combination has been trialled
   together. This is an evidence gap, and no amount of modelling closes it.

The honest next step is unchanged by any of this synthesis: run the PCSK9-editing
outcomes trial. It needs no new biology and unlocks most of the population grid.
The remuscularisation problem needs a different kind of result entirely — a
chronic-injury, human-scale graft demonstration — and until that exists, the
correct answer to "is there a cure here?" remains no.

---

## Reproducibility

- 712 tests, 98.54% coverage on `ihd_solutions`; `ruff` clean.
- 34 ATLAS-specific tests, including the zero-coupling falsifier, the
  editing-not-inert guard, and the no-cure / no-proof assertions.
- 47 therapies, 48 anchors, 19 populations validated.
- 33 evidence artefacts SHA-256 verified; ATLAS artefacts excluded from that
  manifest and tagged `atlas_invented_hypothesis_unvalidated` in every row.
- Monte Carlo draws are seeded; the graft evaluation is reproducible.
