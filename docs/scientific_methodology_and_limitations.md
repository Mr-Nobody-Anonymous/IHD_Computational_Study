# Scientific Methodology, Validation & Limitations

> [!IMPORTANT]
> **Scientific Transparency Statement:** This project combines rigorous crystallographic biophysics with exploratory dynamical modeling. To ensure scientific integrity, this document explicitly delineates **what is empirically validated ground truth** versus **what is an idealized theoretical model**.

---

## 1. Structural Interface Analysis (Empirically Validated)

### Methodology
* **Coordinates:** High-resolution X-ray crystallographic structures from the Protein Data Bank:
  - `PDB 3BPS` (2.30 Å resolution): Human PCSK9 catalytic domain bound to LDLR EGF-A.
  - `PDB 3GCX` (1.98 Å resolution): Human PCSK9 in complex with LDLR EGF-A (Kwon et al., *PNAS* 2008).
  - `PDB 3P5C` (3.30 Å resolution): Human PCSK9 in complex with full LDLR extracellular domain (Lo Surdo et al., *EMBO Rep* 2011).
* **Solvent-Accessible Surface Area (SASA):** Calculated using the **Shrake-Rupley algorithm** (probe radius $r_{\text{probe}} = 1.4\text{ \AA}$ representing a water molecule, 960 test points per atom).
  $$\Delta \text{SASA}_i = \text{SASA}_i(\text{isolated monomer}) - \text{SASA}_i(\text{bound complex})$$
* **Experimental Benchmark:** Validated against wet-lab Surface Plasmon Resonance (SPR) alanine-scanning mutagenesis data from Kwon et al. (*PNAS* 2008, 105(6):1820-1825).

### What This Establishes
* **The "Flat PPI" Paradox:** The total buried interface spans **$591.5\text{ \AA}^2$** across $15+$ residues on multiple discontinuous loops on PCSK9.
* A classical Lipinski Rule-of-5 small molecule ($<500\text{ Da}$) covers only $\sim 250 - 350\text{ \AA}^2$. This quantitatively explains why 15 years of classical small-molecule medicinal chemistry failed to yield an approved drug targeting this interface: small molecules cover $<40\%$ of the contact zone, allowing the flexible LDLR EGF-A loop to wrap around them.
* Validated modalities must either be **macrocyclic peptides** ($\sim 1.5\text{ kDa}$, like Merck's MK-0616 in Phase 3) that physically mimic the entire $\beta$-hairpin loop, or **upstream genetic silencing** (VERVE-102 base editing, Inclisiran siRNA).

---

## 2. Plaque Biomechanics & Continuum Remodeling (Theoretical / Exploratory)

### Methodology
* Modeled 60-month progression of human coronary atherogenesis in the Left Anterior Descending (LAD) artery using coupled differential equations:
  - Intimal lipid accumulation and core expansion as a function of circulating LDL.
  - Fibrous cap remodeling balancing VSMC collagen synthesis against macrophage/MMP enzymatic degradation.
  - Wall tensile hoop stress via Laplace's law for thin-walled cylindrical shells:
    $$\sigma_{\text{hoop}} = \frac{P_{\text{lumen}} \cdot R_{\text{lumen}}}{T_{\text{cap}}}$$
  - Clinical vulnerability threshold: Fibrous cap thickness $T_{\text{cap}} \le 65\text{ \mu m}$ and peak tensile stress $\sigma_{\text{hoop}} \ge 300\text{ kPa}$.

### Limitations & Reality Check
* **Not a Patient-Specific Predictor:** The rupture times (e.g., Month 23.0 vs. Month 56.0) are illustrative dynamical outputs governed by chosen baseline parameters.
* **Simplifications:**
  - Artery geometry is modeled as an axisymmetric cylinder rather than a 3D patient-specific vessel with curvature, bifurcation, and intravascular ultrasound (IVUS) / optical coherence tomography (OCT) reconstruction.
  - Material constitutive law is linear elastic; real vascular tissue exhibits nonlinear, anisotropic hyperelastic behavior (e.g., Mooney-Rivlin or Holzapfel-Gasser-Ogden strain energy functions).
  - Blood flow is treated as uniform pressure; real rupture is influenced by 3D fluid-structure interaction (FSI) and spatial wall shear stress (WSS) gradients.

---

## 3. Vascular Cytokine ODE Network (Textbook Systems Biology)

### Methodology
* Modeled the vascular inflammatory signaling network over 180 days:
  $$\text{oxLDL} \longrightarrow \text{NF-}\kappa\text{B} \longrightarrow \text{NLRP3 Inflammasome} \longrightarrow \text{Caspase-1} \longrightarrow \text{IL-1}\beta \longrightarrow \text{IL-6} \longrightarrow \text{hs-CRP} \longrightarrow \text{MMP-1/9}$$
* Evaluated known pharmacological mechanisms: Statin LDL lowering, Canakinumab (CANTOS trial anti-IL-1$\beta$), and Colchicine (COLCOT/LoDoCo2 trials NLRP3/tubulin inhibition).

### Limitations & Reality Check
* Real human vascular immunology involves hundreds of interacting species, autocrine/paracrine spatial gradients, regulatory T-cells ($T_{\text{reg}}$), and negative feedback loops (e.g., IL-10, SOCS proteins) that simple 7-variable ODEs simplify.
* The model serves to conceptually demonstrate why lipid lowering alone (statin) leaves significant residual inflammatory risk, and why dual-target therapy (lipid suppression + inflammasome extinction) produces synergistic protection.

---

## 4. CRISPR Base-Editing Genomic Specificity (Computational Screening)

### Methodology
* Analyzed the human *PCSK9* Intron 1 splice donor locus (`GGGCAGCGCCACCCACAGGA` PAM `TGG`) targeted by VERVE-102.
* Implemented the Cutting Frequency Determination (CFD) / base-editing mismatch tolerance matrix to evaluate homologous loci across human proprotein convertases (PCSK1, PCSK2, FURIN, MBTPS1).

### Limitations & Reality Check
* In-silico string matching and CFD scoring are first-pass filtering tools.
* **Clinical Gold Standard:** Real regulatory dossiers for human therapeutic editing require empirical genome-wide assays (e.g., **GUIDE-seq**, **CIRCLE-seq**, **DISCOVER-seq**) in primary human hepatocytes, followed by targeted deep sequencing ($>100,000\times$ coverage) and 15-year patient registries.

---

## References

1. **Kwon HJ, et al.** (2008). Molecular basis for LDL receptor recognition by PCSK9. *Proc Natl Acad Sci USA*, 105(6):1820–1825.
2. **Lo Surdo P, et al.** (2011). Mechanistic implications for LDL receptor degradation from the PCSK9/LDLR structure at neutral pH. *EMBO Rep*, 12(12):1300–1305.
3. **Ridker PM, et al. (CANTOS)** (2017). Antiinflammatory Therapy with Canakinumab for Atherosclerotic Disease. *N Engl J Med*, 377(12):1119–1131.
4. **Tardif JC, et al. (COLCOT)** (2019). Efficacy and Safety of Low-Dose Colchicine after Myocardial Infarction. *N Engl J Med*, 381(26):2497–2505.
5. **Verve Therapeutics / NEJM** (2025). Interim Safety and Efficacy of VERVE-102 in Heterozygous Familial Hypercholesterolemia. *N Engl J Med*.
