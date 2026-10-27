# PeriScope

**A Blood–Brain World Model Nominates Peripherally Druggable Targets in Parkinson's Disease**

The blood–brain barrier limits therapeutic exposure in Parkinson's disease. Choosing a peripheral target requires predicting which brain-cell programs it would shift and in which direction. **PeriScope is a potential-form blood–brain world model** for this cross-tissue therapeutic target analysis.

PeriScope learns from unpaired blood and brain single-cell atlases, with paired RNA–protein measurements anchoring the protein input. A learned coupling potential organizes recipient-state updates: the blood state is held fixed, the recipient brain latent state follows the potential gradient, and a decoder reads out its transcriptional response. Clamping a protein input at observed abundance quantiles makes candidate modulations comparable for a specified recipient cell state.

Computational perturbation, human brain expression and blood-protein Mendelian randomization nominate **CD22, CD40 and CD71 (TFRC)**, with **16 recipient-specific programs and 337 distinct genes**. Their modulation directions depend on the recipient state. Direct target mechanisms link the three proteins to **24 clinical-stage molecules**. Human brain cohorts, spatial transcriptomics and an independent blood–brain response evaluation in pigs provide complementary evidence for the model and its programs.

This repository accompanies the current **Nature Machine Intelligence submission**, updated against the analysis files and figure data on 27 October 2026. The manuscript title is reproduced above.

## Run the released analyses

Two checks require only Python's standard library:

```bash
python reproduce/verify_manifest.py
python reproduce/verify_results.py
```

The first verifies the released file hashes. The second reconstructs gene memberships and all 576 MR estimates, checks the global BH correction, and checks the published benchmark and animal-omission summaries.

To recompute the revised numerical analyses on a CPU, use Python 3.10 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r reproduce/requirements-numerical.txt
python reproduce/recompute_nmi.py --output recomputed/nmi --ma
```

This recomputes the seven-method matched comparison, the additional available-seed aggregation, donor dependence, MR directionality and the Ma cohort test from released inputs. `--ma` reruns the 5,000 expression-matched random sets from donor pseudobulk counts; omit it for the faster checks. Each output directory must be new. Original results remain in `results/`.

The primary comparison in Supplementary Table 1 uses **three shared benchmark seeds, 42–44, for every method**, with equal weighting over four blood routes within each recipient. Pooling the two recipient contrasts reproduces Spearman correlations of **0.259 for PeriScope**, 0.211 for conditional MLP, 0.176 for the label/state marginal, 0.165 for conditional ICNN, 0.159 for the label/state mean, 0.107 for OT–Ridge and 0.042 for CellOT. The ten-model PeriScope aggregation is a separate retained analysis; it does not define this matched comparison.

## Trained models

The [NMI submission release](https://github.com/dong-zhongyuan/PeriScope/releases/tag/nmi-submission-20261027) provides all ten authoritative model checkpoints, seeds **42–51**. Download and verify them without a GitHub login:

```bash
python reproduce/download_checkpoints.py --output /absolute/path/to/periscope_models
python -m pip install torch==2.5.1
python reproduce/check_model_loading.py --models /absolute/path/to/periscope_models
```

The downloader verifies both the archive and every checkpoint against the scientific model registry. The second command sequence tests strict checkpoint loading and a small synthetic CPU inference. For the original analysis directory layout, use an output path ending in `pd_product_assets/results/optimization_20261003/models`.

Discovery uses seeds **42–46** and confirmation **47–51**. The repaired 47–49 checkpoints are part of that ten-model set. They are distinct from the **42–44** models used in the three-run human comparison. Gene and protein vocabularies, centering metadata and state labels are provided in `reproduce/inputs/`.

## Repository guide

| Location | Contents |
| --- | --- |
| `src/pdproduct/` | Core model, preprocessing and package interfaces |
| `scripts/` | Data preparation and analysis entry points |
| `_work/rescreen_20261002/` | Current input preparation |
| `_work/optimization_20261003/`, `_work/checkpoint_repair_20261003/` | Training, registered models, dose queries, confirmation and target integration |
| `_work/benchmark_shared_20261003/` | Human comparison, shared representations and readouts |
| `_work/cellot_results_recovery_20261003/` | Archived CellOT recipe, environment and numerical settings |
| `_work/blood_brain_MR_20261003/` | Instrument acquisition, harmonization and MR |
| `_work/pig_external_validation_20261003/` | Cross-species adapter, predictions and animal-level evaluation |
| `_work/spatial_subtype_20261003/` | Subtype reference and spatial analysis |
| `_work/nmi_revision_20261017/` | Third-cohort, spatial, donor-dependence, matched-benchmark and genetic follow-up scripts |
| `results/` | Complete analysis tables, program memberships, prediction arrays and model identities |
| `results/nmi_revision_20261017/` | Current follow-up results and Ma donor pseudobulk inputs |
| `figure_source_data/` | Current numerical figure inputs |
| `figures/` | Preserved plotting sources, R factories, layouts and asset attribution |
| `reproduce/` | Portable numerical runners, checkpoint download, input metadata and integrity checks |

`_work/` contains tracked scientific source, including sibling imports used by the original analyses. Its historical folder names are retained for compatibility.

## Figure numbering

The plotting project retains its original seven figure identifiers. The NMI submission contains six main figures and one supplementary figure:

| Project identifier | NMI placement | Information |
| --- | --- | --- |
| 1 | Main Fig. 1 | Cohorts and the blood–brain world model |
| 2 | Main Fig. 2 | Protein-input prediction and calibration |
| 3 | Main Fig. 3 | Human method comparison and independent response transfer |
| 4 | Main Fig. 4 | Protein queries and recipient-dependent programs |
| 5 | Supplementary Fig. 1 | Confirmation selectivity, MR and seed robustness |
| 6 | Main Fig. 5 | Human brain and spatial evidence |
| 7 | Main Fig. 6 | Target programs, pharmacology and gene overlap |

## Reproduction scope

[REPRODUCIBILITY.md](REPRODUCIBILITY.md) gives exact commands, input requirements, seed definitions and the verification record. The released tables, compact prediction contrasts, donor pseudobulks and checkpoints support the CPU recomputations above. A complete run from downloaded raw matrices through GPU training, spatial deconvolution and final figure assembly requires additional large input resources and environment setup; it has not been executed from a clean public checkout for this release.

The full analysis environment is listed in `requirements.txt`; install the core package with `pip install -e .`. CellOT's original numerical recipe includes its CPU/BLAS configuration. Spatial deconvolution uses RCTD/CSIDE from [spacexr](https://github.com/dmcable/spacexr). Vendored CellOT code retains its BSD-3-Clause license; figure asset attribution is retained under `figures/assets/`.

## Data resources

| Resource | Use |
| --- | --- |
| [GSE164378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378), Hao et al., Cell 2021 | Paired RNA–protein reference |
| [GSE223138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223138), Xiong et al., npj Parkinson's Disease 2024 | PD blood RNA |
| [Zenodo 14372434](https://zenodo.org/records/14372434), Moquin-Beaudry et al., Brain 2025 | Reference blood RNA |
| [GSE178265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178265), Kamath et al., Nature Neuroscience 2022 | Brain model anchor and held-out brain evaluation |
| [GSE184950](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184950), Wang et al., Science Advances 2024 | Brain resource inventory |
| [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783), Smajić et al., Brain 2022 | Independent matched brain-state evidence |
| [GSE253975](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253975), [GSE253462](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253462), Ma et al., Nature Communications 2025 | Visium tissue and subtype-resolved snRNA-seq cohort |
| Sun et al., Nature 2023; Gudjonsson et al., Nature Communications 2022 | UKB plasma and AGES serum cis-pQTLs |
| Smith et al., Nature Neuroscience 2021 | UKB brain imaging GWAS |
| Olney et al., Journal of Neuroinflammation 2024 | Independent pig systemic inflammation study |

## Citation

Zhongyuan Dong, Lianghua Wang and Xuanlin Meng. *PeriScope: A Blood–Brain World Model Nominates Peripherally Druggable Targets in Parkinson's Disease*. Manuscript, 2026.
