# PeriScope

A pairing-free blood–brain world model that nominates peripherally druggable
targets in Parkinson's disease.

Peripheral immune dysfunction is increasingly implicated in Parkinson's
disease (PD), yet no framework connects blood-accessible proteins to brain
glial programs, because paired blood–brain samples do not exist in PD.
PeriScope is a counterfactual cross-tissue world model that learns
blood-state–brain-state coupling from unpaired single-cell atlases anchored
by same-cell protein–RNA (CITE-seq) measurements. A potential-form
equilibrium simulator (two VAE-style tissue encoders, per-tissue and coupling
energy terms, unrolled equilibrium re-solving under clamped protein doses) is
trained on pooled PBMC cohorts and substantia-nigra snRNA-seq with a donor-level
locked test split. Pre-registered destruction gates (random-pairing and
coupling-off retraining) show that directional readouts depend on learned
pairing structure rather than shortcuts. Counterfactual protein perturbation
elicits dose-monotonic, cell-type-coherent brain gene programs, which converge
with spatial-transcriptomic and program-level genetic evidence in the human PD
brain. Evidence integration across nine layers nominates SELP/CD62P, CD22 and
CR1/CD35 as peripherally actionable, plasma-measurable targets whose
modulation acts on the brain without requiring blood–brain-barrier
penetration.

## Repository layout

```
src/pdproduct/        core package (simulators/ccwm.py: CCWM equilibrium model;
                      datasets, provenance run manifests, analysis, QC, CLI)
scripts/              pipeline entry points, in run order:
                      preprocess_* / annotate_* / ingest_*   data freezing
                      step1_world_model/ccwm/                training data + training
                      p3_5_seed_stability/                   stable-edge posterior
                      p3_6_gates/                            diagnostic gates 3-6
                      p4_dose_response/                      dose-response curves
                      p4_5_randko/                           random-protein knockout
                      p6_spatial/                            competitive spatial test
config/               frozen combinatorics (pd_combos.yaml), benchmark grid,
                      paths template (copy paths.example.yaml -> paths.yaml)
contracts/            JSON schemas for run manifests and evidence cards
registry/             method-capability registry
```

Analysis scripts reference a frozen asset layout (see `config/paths.example.yaml`);
they are research artifacts of the archived runs and contain frozen absolute
paths of that layout in a few places.

## Installation

Python >= 3.10.

```bash
pip install -r requirements.txt
pip install -e .
```

Training uses PyTorch 2.5.1 (CUDA build in the archived runs; any recent torch
works for inference). All random seeds, null constructions and multiplicity
corrections are fixed inside the scripts.

## Data sources (public)

- CITE-seq protein–RNA bridge: GSE164378 (Hao et al., Cell 2021)
- PD PBMC cohort: GSE223138 (Xiong et al., npj Parkinson's Disease 2024)
- Reference PBMC cohort: Moquin-Beaudry et al., Brain 2025 (Zenodo 14372434)
- snRNA-seq substantia nigra: GSE178265 (Kamath et al., Nature 2022)
- GeoMx spatial molecular imaging: GSE253975 (Ma et al., Nat Commun 2025)
- PD GWAS summary statistics: Leonard 2025, GCST009325, FinnGen R11

## Citation

If you use this code, please cite the accompanying PeriScope manuscript
(Dong, Wang & Meng).
