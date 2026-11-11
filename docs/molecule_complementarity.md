# Complementary IHD medicines: what can be tested without inventing a hybrid

**This is a model comparison of separate medicines—not a newly designed molecule, a tested fixed-dose combination, or treatment advice.**

## Selection

ATLAS-5/ATLAS-5F is a simulated cardiomyocyte-patch architecture, not a molecule. The strongest outcome-supported pharmacologic strategy in this study is LDL-C lowering. I used high-intensity statin therapy as the model backbone and tested colchicine as the complementary anti-inflammatory component. I also included enlicitide, an oral PCSK9-directed agent, to test whether its LDL-lowering signal receives hard-outcome credit in the current evidence model.

I did **not** chemically fuse these compounds. A fragment from one medicine does not automatically preserve the other medicine's binding, exposure, selectivity, or safety. The repository's earlier de-novo PCSK9 docking results are explicitly marked non-usable, so they are not a sound basis for creating a new molecular structure (`docs/FULL_REPORT.md`, §D.7).

## Model test

The scenarios were evaluated in the repository's `secondary_post_mi` population over the default five-year horizon, in both strict outcome-credit and surrogate-credit modes. The estimate is **modelled events avoided per 1,000**, not an observed trial result. Every combination is flagged as untrialled, and the population transport warning remains visible.

| Scenario | Strict events avoided / 1,000 | Change vs LDL backbone | Surrogate-mode events avoided / 1,000 | Evidence warning |
|---|---:|---:|---:|---|
| Aspirin + high-intensity statin (backbone) | 75.21 | — | 75.21 | transported population |
| Backbone + colchicine | 102.38 | **+27.17** | 102.38 | transported population |
| Backbone + enlicitide | 75.21 | 0.00 | 93.58 | no hard-outcome anchor; transported population |
| Backbone + colchicine + enlicitide | 102.38 | **+27.17** | 116.93 | no hard-outcome anchor; transported population |

### Interpretation

- **Colchicine is the complementary part that improves the model's strict event estimate.** In this scenario it adds 27.17 modelled events avoided per 1,000 over the LDL backbone. That is not proof of synergy or a measured effect for this exact combination. Trial results vary by context—CLEAR SYNERGY was neutral in ACS—and the model flags population transport.
- **Enlicitide adds surrogate LDL-C credit, not hard-outcome credit.** The strict-mode model therefore gives it no incremental event benefit until an outcome anchor is available. Adding colchicine does not change that evidence limitation.
- The model does not calculate drug-drug interactions, combination-specific adverse effects, or the safety of a fused chemical structure. These rows cannot establish a positive benefit-risk balance.
- The model's `transported_population` warning means the evidence is not fully matched to the simulated population. The combination itself has **zero randomized proof** in this project.

Reproduce with:

```bash
python run_solutions.py molecule-complementarity
```

Outputs: `results/molecule_complementarity/complementarity_report.txt` and
`complementarity_scenarios.csv`. The tests in `tests/test_molecule_complementarity.py`
check that colchicine is scored as a separate component, enlicitide does not receive
unsupported hard-outcome credit, and no scenario is labelled as a validated hybrid.
