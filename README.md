# PeriScope

**A Blood–Brain World Model Maps Peripheral Interventions onto Brain-Cell Responses in Parkinson's Disease**

Treating the brain through the blood bypasses the blood–brain barrier, but which blood target moves which brain cell has been unknown. **PeriScope is a potential-form blood–brain world model** that answers this question computationally: set the dose of any of 213 blood-cell surface proteins, and it predicts how astrocyte and microglial gene programs in Parkinson's disease (PD) respond.

PeriScope learns from unpaired blood and brain single-cell atlases, with paired RNA–protein measurements anchoring the protein input. A learned coupling potential organizes recipient-state updates: the blood state is held fixed, the recipient brain latent state follows the potential gradient, and a decoder reads out its transcriptional response. Clamping a protein input at measured abundance quantiles makes every candidate modulation comparable for a specified recipient cell state.

## Highlights

- **Built on 813,363 blood and brain cells, nuclei and spatial spots** from 110 donors, with 161,159 paired RNA–protein cells calibrating the protein input (median antigen-level Spearman 0.63 in held-out donors).
- **Predicts disease-related brain changes in locked donors.** With the same seeds for every method, PeriScope reaches a disease-direction Spearman of 0.250 over eight blood-to-brain routes: 72% above CellOT, 196% above OT–Ridge and 27% above the strongest comparator. In two independent brain cohorts it reaches 3.0× and 2.3× CellOT.
- **Transfers across species.** In a paired pig endotoxemia study its brain-response correlation is 2.8× that of the next-best method.
- **Depends on blood–brain alignment.** Training with shuffled blood–brain pairing removes the program responses.
- **Adds a protein–brain imaging Mendelian randomization (MR) atlas**: 2,941 UK Biobank plasma protein assays × 3,963 brain imaging phenotypes = 11,655,183 protein–phenotype pairs, 1,006 at full-family q < 0.05, with independent replication, colocalization on two protein platforms and reverse MR.
- **Narrows 213 antigens to two antibody-accessible blood targets.** TFRC (CD71) links blood myeloid cells to iron-sensitive imaging of the substantia nigra and pallidum; IFNGR1 (CD119) provides a second route through astrocytes, and its pathway can be blocked by an approved antibody.

## Repository layout

| Location | Contents |
| --- | --- |
| [`pd-v3/`](pd-v3/README.md) | MR atlas, replication, colocalization, reverse MR, drug evidence, target funnel, random-pairing ablation, benchmark evaluation, the code of the five main figures and the summary tables |
| `src/pdproduct/` | World-model package: core model, preprocessing and interfaces |
| `scripts/` | Data preparation and analysis entry points |
| `_work/` | Input preparation, model training, dose queries, program discovery and confirmation, human benchmark, pig evaluation, spatial analysis and target integration |
| `results/` | Analysis tables, program memberships, prediction arrays and model identities |
| `reproduce/` | Input rebuilding, checkpoint download and loading checks, input metadata |
| `config/`, `contracts/`, `registry/` | Run configuration, artifact and manifest schemas, data and method registries |
| `figure_source_data/`, `figures/` | Figure inputs, plotting assets and asset attribution |

## Data and models

- **MR atlas and result tables.** The complete atlas (primary and sensitivity policy, each in 12 csv.gz parts with headers), the colocalization, replication, drug-evidence and funnel tables, and the trained models (ten full models in three archive parts, and five random-pairing models) are deposited at Zenodo, DOI [10.5281/zenodo.24374535](https://doi.org/10.5281/zenodo.24374535). The same summary tables are in [`pd-v3/result_tables/`](pd-v3/result_tables).

Download and verify the models from Zenodo, then test strict loading and a small CPU inference:

```bash
python reproduce/download_checkpoints.py --output /absolute/path/to/periscope_models
python -m pip install torch==2.5.1
python reproduce/check_model_loading.py --models /absolute/path/to/periscope_models
```

Gene and protein vocabularies, centering metadata and state labels are in `reproduce/inputs/`.

## Reproduce

**Inputs.** The six cell-level input matrices are rebuilt from public sources with

```bash
python reproduce/build_inputs.py --output /absolute/path/to/new_build
```

after installing `reproduce/requirements-inputs.txt` and the R packages SeuratObject and Matrix. The runner downloads and checks the source files, applies the bundled annotations and donor partitions, and verifies all 84 output arrays. `--raw-cache /path/to/raw --offline` reuses an existing source cache; `--resume` continues a stopped build.

**Analyses and figures.** The MR atlas, colocalization, target funnel and the five main figures are run from `pd-v3/`; its [README](pd-v3/README.md) lists the scripts, method definitions and figure commands.

**Environment.** The analysis environment is listed in `requirements.txt`; install the package with `pip install -e .`. Spatial deconvolution uses RCTD/C-SIDE from [spacexr](https://github.com/dmcable/spacexr); colocalization uses the R packages `coloc` and `susieR`. Vendored CellOT code keeps its BSD-3-Clause license; figure asset attribution is under `figures/assets/`.

## Data resources

| Resource | Use |
| --- | --- |
| [GSE164378](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378), Hao et al., Cell 2021 | Paired RNA–protein (CITE-seq) reference |
| [GSE223138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223138), Xiong et al., npj Parkinson's Disease 2024 | PD blood RNA |
| [Zenodo 14372434](https://zenodo.org/records/14372434), Moquin-Beaudry et al., Brain 2025 | Reference blood RNA |
| [GSE178265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178265), Kamath et al., Nature Neuroscience 2022 | Brain anchor and locked test donors |
| [GSE157783](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783), Smajić et al., Brain 2022 | Independent brain cohort |
| [GSE184950](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184950), Wang et al., Science Advances 2024 | Independent brain cohort |
| [GSE253975](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253975), [GSE253462](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253462), Ma et al., Nature Communications 2025 | Substantia nigra spatial transcriptomics and its snRNA-seq reference |
| Olney et al., Journal of Neuroinflammation 2024 | Pig endotoxemia: paired blood and prefrontal cortex |
| [GSE245309](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE245309), Boles et al. 2024; [GSE122766](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122766), Shemer et al., Nature Communications 2018; [GSE109834](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109834), Hirbec et al., Glia 2018 | Mouse systemic inflammation |
| Sun et al., Nature 2023; Gudjonsson et al., Nature Communications 2022 | UKB-PPP plasma and AGES serum protein GWAS and cis-pQTLs |
| Smith et al., Nature Neuroscience 2021; Wang et al., Nature Neuroscience 2022 | BIG40 brain imaging GWAS (discovery and replication); QSM and T2\* GWAS |
| 1000 Genomes Project, European panel | Linkage disequilibrium reference |
| ChEMBL, Open Targets, DrugCentral, Guide to Pharmacology, DGIdb, drug labels | Drug evidence |

## Citation

*PeriScope: A Blood–Brain World Model Maps Peripheral Interventions onto Brain-Cell Responses in Parkinson's Disease*. Manuscript, 2026.

Data: PeriScope: plasma protein–brain imaging Mendelian randomization atlas, colocalization, drug evidence and result tables. Zenodo (2026). doi:[10.5281/zenodo.24374535](https://doi.org/10.5281/zenodo.24374535)
