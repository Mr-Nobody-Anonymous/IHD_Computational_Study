# Computational Biophysics & Systems Biology of Ischemic Heart Disease (IHD)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-blue.svg)](https://www.python.org/)
[![PDB Structures](https://img.shields.io/badge/PDB-3BPS%20%7C%203GCX%20%7C%203P5C%20%7C%202P4E-brightgreen.svg)](https://www.rcsb.org/)

A reproducible computational biology framework investigating the structural biophysics, vascular biomechanics, inflammatory signaling networks, and curative therapeutic modalities for **Ischemic Heart Disease (IHD)**.

---

## 📌 Scientific Context & The "Flat PPI" Paradox

Ischemic Heart Disease (IHD) remains the leading cause of human mortality worldwide (~9 million deaths annually). Although standard-of-care pharmacotherapies (statins) and revascularization procedures (PCI stents, CABG surgery) prolong event-free survival, they manage symptoms and slow progression without permanently halting atherogenesis or regrowing infarcted myocardium.

This repository provides an empirical, multi-scale computational investigation addressing the core biological bottlenecks:
1. **The Biophysical Paradox of PCSK9:** Why 15+ years of classical medicinal chemistry failed to yield an approved oral small molecule targeting the PCSK9–LDLR interface (quantified via multi-crystal solvent-accessible surface area burial).
2. **Experimental Alanine-Scanning Validation:** Benchmark of computational residue burial against empirical wet-lab Surface Plasmon Resonance (SPR) binding affinities from Kwon et al. (*PNAS* 2008).
3. **Plaque Rupture Biomechanics:** A 60-month finite-difference simulation of fibrous cap remodeling, necrotic core dynamics, and Laplace tensile hoop stress.
4. **Vascular Systems Immunology:** Ordinary differential equation (ODE) modeling of the $\text{oxLDL} \rightarrow \text{NLRP3} \rightarrow \text{IL-1}\beta \rightarrow \text{hs-CRP}$ inflammatory cascade (CANTOS/COLCOT axis).
5. **CRISPR Base-Editing Specificity:** Computational genomic off-target scoring (CFD) of the clinical-stage VERVE-102 ABE8e guide RNA.

---

## 📂 Repository Structure

```
IHD_Computational_Study/
├── data/
│   ├── structures/                  # High-resolution crystallographic coordinate files (.cif)
│   │   ├── 3BPS.cif                 # PCSK9 + LDLR EGF-A (2.30 Å resolution)
│   │   ├── 3GCX.cif                 # PCSK9 + LDLR EGF-A (1.98 Å resolution, Kwon et al. 2008)
│   │   ├── 3P5C.cif                 # PCSK9 + full LDLR extracellular domain (3.30 Å)
│   │   ├── 2P4E.cif                 # High-resolution PCSK9 apo / domain map (1.97 Å)
│   │   ├── 4NMX.cif                 # PCSK9 with engineered inhibitor peptide (1.90 Å)
│   │   └── 4HHB.cif                 # Human hemoglobin tetramer reference (1.74 Å)
│   └── experimental/
│       └── kwon_2008_alanine_scanning.csv  # Peer-reviewed SPR binding data (PNAS 2008)
├── src/
│   ├── 01_multistructure_sasa_benchmark.py # Multi-structure SASA & SPR correlation
│   ├── 02_interface_contact_mapping.py     # Atomic contacts, H-bonds, and docking grid
│   ├── 03_plaque_biomechanics_model.py     # 60-month coronary atherogenesis simulation
│   ├── 04_immune_cytokine_network.py       # 180-day vascular inflammatory ODE network
│   └── 05_crispr_genomic_safety_cfd.py     # Genome-wide CFD off-target screening
├── results/
│   ├── figures/
│   │   ├── experimental_validation_plot.png # Cross-crystal burial & SPR correlation
│   │   ├── plaque_simulation_curves.png     # 4-panel plaque biomechanics trajectories
│   │   └── cytokine_dynamics.png            # 4-panel inflammatory biomarker timecourses
│   ├── tables/
│   │   ├── multi_structure_sasa_variability.csv
│   │   ├── experimental_vs_computational_correlation.csv
│   │   ├── plaque_biomechanics_timeline.csv
│   │   ├── cytokine_network_data.csv
│   │   └── verve102_offtarget_scores.json
│   └── reports/
│       ├── interface_contact_summary.txt
│       ├── rigorous_sasa_report.txt
│       └── verve102_genomic_safety_report.txt
├── docs/
│   ├── scientific_methodology_and_limitations.md # Critical reality check & assumptions
│   ├── IHD_master_solution_blueprint.md         # Comprehensive multi-scale dossier
│   └── IHD_cure_deep_research.md                # Research review on VERVE-102 & iPSCs
├── .gitignore
├── requirements.txt
└── README.md
```

---

## 🔬 Key Scientific Modules & Results

### 1. Cross-Crystallographic Interface Burial & Alanine-Scanning Benchmark
Using the **Shrake-Rupley algorithm** ($r_{\text{probe}} = 1.4\text{ \AA}$ water probe, 960 points/atom), we evaluated per-residue solvent-accessible surface area burial across three independent crystal structures (`3BPS`, `3GCX`, `3P5C`):

```
ResID   ResName   3GCX (1.98Å)   3BPS (2.30Å)   3P5C (3.30Å)    Mean ± SD (Å²)   Experimental ΔΔG (kcal/mol)
─────────────────────────────────────────────────────────────────────────────────────────────────────────────
379     PHE           81.9           77.4           92.3          83.9 ± 6.2     +3.00 (Binding abolished)
369     ILE           65.2           64.9           50.7          60.3 ± 6.8     +1.87 (23-fold loss)
153     SER           57.9           66.2           26.5          50.2 ± 17.1    +0.10 (Neutral effect)
377     THR           41.1           43.0           54.9          46.3 ± 6.1     +1.31 (9-fold loss)
238     ASP           36.7           43.2           31.9          37.3 ± 4.6     +1.13 (6.7-fold loss)
378     CYS           36.2           37.0           29.2          34.1 ± 3.5     N/A (Disulfide bridge)
374     ASP           36.3           38.8           26.9          34.0 ± 5.1     +0.82 (4-fold loss / D374Y: -1.91)
380     VAL           28.3           37.5           27.0          30.9 ± 4.6     +0.71 (3.3-fold loss)
194     ARG           33.1           20.6           26.0          26.6 ± 5.1     +1.75 (19-fold loss)
```

#### Why Classical Small Molecules Failed:
* Total buried surface area (BSA) is **$591.5\text{ \AA}^2$**, with the altered PCSK9 footprint covering **$2,676\text{ \AA}^2$**.
* A Rule-of-5 small molecule ($<500\text{ Da}$) covers only **$\sim 250 - 350\text{ \AA}^2$** ($<40\%$ of the interface).
* The only viable molecular modalities are:
  1. **Synthetic Macrocycles (e.g., Merck's MK-0616):** Large cyclic peptides ($\sim 1.5\text{ kDa}$) designed to cover $>700\text{ \AA}^2$.
  2. **Gene Silencing (VERVE-102 base editing / Inclisiran siRNA):** Bypasses the flat interface by arresting protein synthesis upstream in hepatocytes.

---

### 2. Plaque Biomechanics & Vulnerability Dynamics
Simulates a 60-month continuum finite-difference model of human Left Anterior Descending (LAD) coronary artery atherogenesis under physiological pulsatile pressure ($140\text{ mmHg} = 18.66\text{ kPa}$):
* **Untreated High-Risk (LDL 190 mg/dL):** Fibrous cap thins below the $65\text{ \mu m}$ critical threshold, and tensile hoop stress surpasses the $300\text{ kPa}$ material failure threshold, predicting rupture at **Month 23.0**.
* **High-Intensity Statin (LDL 95 mg/dL):** Prolongs event-free survival, delaying rupture to **Month 56.0**; persistent residual vascular inflammation ultimately causes late cap failure.
* **VERVE-102 Base Editing (LDL 38 mg/dL):** Drives net cholesterol efflux, reversing necrotic lipid core volume from $0.40\text{ mm}^2$ to **$0.26\text{ mm}^2$** and maintaining thick, stable caps ($>156\text{ \mu m}$).
* **Dual Curative Strategy:** Combines gene silencing with inflammasome blockade, achieving complete plaque stabilization ($180\text{ \mu m}$ cap).

---

### 3. Systems Biology of Vascular Inflammation (CANTOS/COLCOT)
An ODE network tracking the inflammatory cascade over 180 days:
* While statins moderately lower biomarkers, **residual inflammatory risk** remains high ($\text{hs-CRP} \approx 37.7\text{ mg/L}$ in inflamed plaques).
* Targeted inflammasome blockade (Colchicine / Canakinumab) combined with radical LDL reduction drives hs-CRP down to **$1.78\text{ mg/L}$** (safely below the clinical $<2.0\text{ mg/L}$ threshold) and suppresses destructive matrix metalloproteinase (MMP) activity to **$20.8\%$**.

---

### 4. CRISPR Base-Editing Genomic Specificity (VERVE-102)
Off-target cutting frequency determination (CFD) analysis for the *PCSK9* Intron 1 splice donor guide (`GGGCAGCGCCACCCACAGGA` PAM `TGG`):
* Mismatch penalties across the 10 bp seed region ensure negligible off-target cleavage in homologous proprotein convertases (**PCSK1**, **PCSK2**, **FURIN**, **MBTPS1**; all CFD scores $<0.086$).
* The use of **Adenine Base Editing (ABE8e-nCas9)** avoids double-strand breaks (DSBs), thereby eliminating chromosomal translocations, chromothripsis, and p53 damage responses associated with classical Cas9 nucleases.

---

## 🚀 Quickstart & Reproducibility

### Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/Mr-Nobody-Anonymous/IHD_Computational_Study.git
cd IHD_Computational_Study
pip install -r requirements.txt
```

### Running the Analysis Suite

Execute the scripts in sequence:

```bash
# 1. Multi-structure SASA calculation & experimental SPR correlation
python src/01_multistructure_sasa_benchmark.py

# 2. Interface contact mapping and docking grid derivation
python src/02_interface_contact_mapping.py

# 3. 60-month coronary plaque biomechanics simulation
python src/03_plaque_biomechanics_model.py

# 4. 180-day vascular inflammatory cytokine ODE network
python src/04_immune_cytokine_network.py

# 5. CRISPR base-editing genomic specificity & CFD scoring
python src/05_crispr_genomic_safety_cfd.py
```

All generated plots are written to `results/figures/`, tables to `results/tables/`, and reports to `results/reports/`.

---

## ⚠️ Limitations & Scientific Transparency

Please read [`docs/scientific_methodology_and_limitations.md`](docs/scientific_methodology_and_limitations.md) for full context:
* **Structural calculations** are based directly on empirical crystallographic ground truth from the RCSB Protein Data Bank.
* **Biomechanics and cytokine models** are idealized continuum and dynamical systems intended to explore mechanistic concepts; they are not patient-specific clinical diagnostic tools.
* **Genomic off-target scores** are in-silico screening estimates; regulatory dossiers require wet-lab assays (GUIDE-seq, CIRCLE-seq) in primary hepatocytes.

---

## 📚 References

1. **Kwon HJ, et al.** (2008). Molecular basis for LDL receptor recognition by PCSK9. *Proc Natl Acad Sci USA*, 105(6):1820–1825.
2. **Lo Surdo P, et al.** (2011). Mechanistic implications for LDL receptor degradation from the PCSK9/LDLR structure at neutral pH. *EMBO Rep*, 12(12):1300–1305.
3. **Ridker PM, et al. (CANTOS)** (2017). Antiinflammatory Therapy with Canakinumab for Atherosclerotic Disease. *N Engl J Med*, 377(12):1119–1131.
4. **Tardif JC, et al. (COLCOT)** (2019). Efficacy and Safety of Low-Dose Colchicine after Myocardial Infarction. *N Engl J Med*, 381(26):2497–2505.
5. **Verve Therapeutics / NEJM** (2025). Interim Safety and Efficacy of VERVE-102 in Heterozygous Familial Hypercholesterolemia. *N Engl J Med*.

---

## 📄 License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
