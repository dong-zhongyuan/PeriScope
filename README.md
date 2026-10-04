# PeriScope

PeriScope predicts protein-guided blood–brain transcriptional responses from unpaired single-cell atlases. A potential-form model couples blood and brain latent states; paired CITE-seq measurements anchor the protein input, which can be clamped at observed abundance quantiles.

The current analysis nominates **CD22, CD40 and CD71 (TFRC)** through confirmation gene-set testing, matched human brain evidence, target accessibility, direct clinical pharmacology and blood-protein Mendelian randomization (MR). The final set contains 16 source–recipient response programs and 337 distinct genes.

## Current results

- Ten authoritative PeriScope models use seeds 42–51: discovery 42–46 and confirmation 47–51. Repaired checkpoints 47–49 are included in this set and supply the three-run human benchmark. Checkpoint identities are recorded in `results/optimization_20261003/model_registry.json`.
- The human benchmark compares seven methods across eight blood–brain routes. It separates distribution MMD, mean-expression MSE, condition-aware response correlation and fixed-condition blood-only response correlation. An independent pig inflammatory challenge supplies the fifth readout.
- PeriScope has the highest mean condition-aware disease-direction correlation on the four astrocyte routes. Method comparisons retain all published readouts, including microglial and distribution results.
- MR evaluates 36 circulating proteins against 16 MRI phenotypes (576 tests, one global BH family). UKB supplies 35 proteins; AGES supplies the remaining eligible protein. Five pallidal estimates support CD22, CD40 and TFRC.
- Computational perturbation and MR are parallel requirements. The 16 final programs also meet C1–C6. Protein prediction correlation, relative-antigen rank and seed sensitivity are evaluations, not additional nomination cutoffs. Cell-subtype spatial results remain a separate follow-up.
- The pig study comprises six saline and four LPS animals with paired whole-blood and prefrontal-cortex RNA. It tests aggregate cross-tissue response transfer through one-to-one Ensembl orthologs. It is not a target-specific intervention experiment.

## Start with the released results

No GPU or private file is required for this check:

```bash
python reproduce/verify_results.py
```

It independently reconstructs target gene unions and recipient counts, recomputes all 576 Wald estimates and their BH correction, and checks the benchmark and animal-omission summaries against the detailed tables.

To rerun MR from the supplied harmonized exposure and outcome inputs with base R, keeping the released results untouched:

```bash
mkdir -p recomputed/mr
cp results/blood_brain_MR_20261003/combined_exposure_instruments.csv recomputed/mr/
cp results/blood_brain_MR_20261003/brain_sentinel_statistics.csv recomputed/mr/
Rscript reproduce/run_brain_mr.R recomputed/mr
```

## Repository contents

| Location | Contents |
| --- | --- |
| `src/pdproduct/` | Current core simulator, data processing and package interfaces |
| `scripts/` | Preprocessing, model and analysis entry points |
| `_work/rescreen_20261002/` | Revised input preparation and supporting analysis modules |
| `_work/checkpoint_repair_20261003/` | Checkpoint training and repair implementation |
| `_work/optimization_20261003/` | Ten-seed registry, dose queries, confirmation tests, brain evidence, nomination and robustness |
| `_work/benchmark_shared_20261003/` | Shared-representation human benchmark |
| `_work/cellot_results_recovery_20261003/` | Verified original CellOT numerical recipe and expected results |
| `_work/blood_brain_MR_20261003/` | Exposure selection, outcome retrieval, harmonization and MR |
| `_work/pig_external_validation_20261003/` | Cross-species adapter, prediction and animal-level evaluation |
| `_work/spatial_subtype_20261003/` | Author-annotated subtype reference and spatial analysis |
| `results/` | Detailed current results, MR input rows, model identities and pig prediction arrays |
| `figure_source_data/` | Numerical inputs to the current seven-figure project |
| `reproduce/` | Portable result checks and MR recomputation |

`_work` is tracked source code in this release. The original directory names are preserved because the analysis uses sibling imports and explicit execution paths. They do not indicate disposable files.

## Full analysis reproduction

See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for execution order, input requirements and the distinction between result recomputation and model retraining. The release contains analysis code and detailed results, but does not bundle the large human single-cell matrices or all trained model weights. The full training pipeline has not been rerun from a clean public checkout as part of this release.

Python dependencies are listed in `requirements.txt`; install the core package with `pip install -e .`. The recorded CellOT environment and CPU/BLAS settings are part of its numerical recipe. Spatial deconvolution uses RCTD/CSIDE from [spacexr](https://github.com/dmcable/spacexr). Vendored CellOT source retains its BSD-3-Clause license.

## Data resources

| Resource | Use |
| --- | --- |
| [GSE164378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378), Hao et al., Cell 2021 | Paired RNA–protein reference |
| [GSE223138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223138), Xiong et al., npj Parkinson's Disease 2024 | PD blood RNA |
| [Zenodo 14372434](https://zenodo.org/records/14372434), Moquin-Beaudry et al., Brain 2025 | Reference blood RNA |
| [GSE178265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178265), Kamath et al., Nature Neuroscience 2022 | Brain model anchor and held-out brain evaluation |
| [GSE184950](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184950), Wang et al., Science Advances 2024 | Brain resource inventory |
| [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783), Smajić et al., Brain 2022 | Independent matched brain-state evidence |
| [GSE253975](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253975) and [GSE253462](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253462), Ma et al., Nature Communications 2025 | Visium spatial tissue and annotated brain subtype reference |
| Sun et al., Nature 2023; Gudjonsson et al., Nature Communications 2022 | UKB plasma and AGES serum cis-pQTLs |
| Smith et al., Nature Neuroscience 2021 | UKB brain imaging GWAS |
| Olney et al., Journal of Neuroinflammation 2024 | Independent pig systemic inflammation study |

## Citation

Please cite the accompanying PeriScope manuscript by Zhongyuan Dong, Lianghua Wang and Xuanlin Meng. Detailed source identities and file hashes are supplied with the analysis tables and release manifest.
