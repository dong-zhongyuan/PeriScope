# pd-v3: protein–brain imaging atlas, colocalization and target funnel

This directory holds the genetic, pharmacological and benchmark analyses of the manuscript and the code of its five main figures. The blood–brain world model, its ten trained models (seeds 42–51) and the model-side result tables are in the parent repository.

The additions are:

- a plasma protein–brain imaging Mendelian randomization (MR) atlas: 2,941 UK Biobank Pharma Proteomics Project (UKB-PPP) assays × 3,963 brain imaging phenotypes, 11,655,183 protein–phenotype pairs;
- replication of the instruments and of the imaging associations in independent samples;
- single-signal and SuSiE colocalization on two protein platforms (AGES and UKB-PPP), and reverse MR;
- drug evidence for 213 blood-cell surface antigens, with the direction of action resolved per record;
- the target funnel that narrows 213 antigens to TFRC (CD71) and IFNGR1 (CD119);
- the human benchmark on eight blood-to-brain routes and two independent brain cohorts;
- a random-pairing ablation on five additional seeds (52–56);
- the code for the five main figures.

## Layout

| Path | Contents |
| --- | --- |
| `pdv3/` | Shared Python package: MR estimation, multiplicity correction, allele alignment, remote range access |
| `tests/` | Unit tests for the package |
| `papers/01_mr_imaging_scan/` | MR atlas, colocalization, reverse MR and replication: acquisition, analysis and aggregation scripts |
| `papers/01_mr_imaging_scan/results/*/scripts/` | Scripts of the SuSiE colocalization and the BIG40 replication |
| `papers/03_world_model_integration/` | Model evaluation, ablation, drug evidence, the target funnel (`build_target_funnel_g8rel.py`) and TFRC's orientation on the clinical assay (`tfrc_orientation.py`, `extract_qsm_tfrc.py`) |
| `papers/03_world_model_integration/results/benchmark_20261119/src/` | Benchmark evaluation on the eight blood–brain routes |
| `papers/03_world_model_integration/results/random_pairing_20261119/training/` | Training of the random-pairing models, seeds 52–56 |
| `papers/03_world_model_integration/results/drug_evidence_expanded_20261119/code/` | Retrieval, parsing and classification of the drug records |
| `papers/figures_v2_20261119/` | Five main figures: one script per panel, `figN/page/FigureN.pdf` is the assembled page |
| `result_tables/` | Summary tables used by the manuscript (see below) |

## Result tables

`result_tables/` contains the summary tables. The complete MR atlas and the five random-pairing checkpoints are deposited at Zenodo, DOI [10.5281/zenodo.24374535](https://doi.org/10.5281/zenodo.24374535), together with the same summary tables. Each of the two atlas tables (primary and sensitivity policy, 11,655,183 rows) is stored as 12 csv.gz parts of about 125 MB; every part carries the header row and can be read on its own.

| Folder | Contents |
| --- | --- |
| `mr_atlas_tables/` | Selected cis instruments, phenotype catalogue, assay and phenotype coverage, instrument replication between the UKB-PPP discovery and replication sets, and the plot tables of the atlas overview |
| `colocalization_reverse_mr_replication/` | Region, exposure and locus tables of the 1,375 colocalization regions (posterior probabilities for two windows and two priors, reverse MR, replication), SuSiE colocalization pairs and signal pairs, BIG40 replication pairs and summaries |
| `target_funnel_and_drug_evidence/` | Antigen gate matrix, funnel counts, final targets with their MR and colocalization pairs and drug matches; 4,739 antigen–drug–route records (2,367 distinct drugs), the hand-curated overlay the funnel reads, and the per-antigen summary; `tfrc_variant_panel.csv`, TFRC's instrument and the clinical-assay leads on the soluble transferrin receptor scale beside nigral QSM and T2\* |
| `world_model_programs_benchmark_ablation/` | Program evidence and frozen gene memberships, benchmark summaries, route-level ablation (seeds 42–44), random pairing (seeds 52–56) against the full model (seeds 42–51), and three figure tables |

## Method definitions

**Benchmark.** Seven methods are compared on eight routes (cDC, cMono, nMono and pDC to Astro and MG-AP) in the locked test donors and in GSE157783 and GSE184950. The readout is the Spearman correlation between predicted and observed PD–control differences over 4,741 genes, averaged over routes with equal weight. Every method is compared on seeds 42–44; PeriScope's ten seeds (42–51) report seed stability per route.

**MR atlas.** Exposures are all 2,941 UKB-PPP assays; 1,954 have a published lead cis-pQTL, and one instrument per assay is chosen by |β/SE| in the discovery set without reading imaging results. Outcomes are 3,935 BIG40 phenotypes and 28 QSM and T2\* phenotypes. Alleles are aligned on GRCh37. The primary policy drops palindromic variants that cannot be oriented; the sensitivity policy infers their orientation from allele frequency. Estimates are Wald ratios with first-order and delta standard errors. The correction family is the full grid of 11,655,183 pairs; non-estimable pairs enter with P = 1. Significance is Benjamini–Hochberg q < 0.05 over the full family.

**Replication.** Instrument effects are compared between the UKB-PPP discovery and replication sets. Imaging associations are tested for direction in the independent BIG40 replication sample; the query list is fixed before replication values are read.

**Colocalization.** Each region pairs an AGES protein GWAS and the corresponding UKB-PPP assay with the imaging phenotype of largest |z| at the instrument. Single-signal colocalization uses windows of ±500 kb and ±1 Mb with p1 = p2 = 1e-4 and p12 = 1e-5 (1e-6 also reported). SuSiE colocalization uses the 1000 Genomes European panel, L = 5 and coverage 0.95; a protein signal qualifies when the MR instrument lies in its 95% credible set or has r² ≥ 0.8 with a variant in it. The unified posterior probability is the larger of the single-signal and the instrument-carrying SuSiE PP.H4, computed separately for each protein platform.

**TFRC orientation.** TFRC's instrument rs3817672 encodes TfR1 Gly142Ser. The UKB-PPP affinity assay reads its Ser142 allele as lowering plasma TFRC (β = −0.524 SD), whereas the clinical soluble transferrin receptor assay (GCST90451684, 45,330 participants, none from UK Biobank) reads it as raising the receptor (β = 0.058 SD), as do transferrin saturation and total iron-binding capacity. Figure 4i therefore scales TFRC's Wald ratios by β(UKB-PPP) / β(clinical) at rs3817672 (−9.11) and Figure 4d reverses TFRC's signed z; P values, q values and colocalization are unchanged. `tfrc_orientation.py` writes `tfrc_variant_panel.csv` from the clinical-assay GWAS cut to TFRC ± 500 kb and the nigral QSM and T2\* regions cut by `extract_qsm_tfrc.py`.

**Target funnel.** After the five model-side criteria (C1–C5), an antigen needs (6) MR support, (7) a unified PP.H4 ≥ 0.8 for the same protein–phenotype pair on either platform, and (8) a clinical-stage or approved drug, not withdrawn, whose action is not opposite to the model's restorative direction, with all direction-resolved programs of the antigen pointing the same way.

**Drug evidence.** Records are compiled from ChEMBL, Open Targets, DrugCentral, Guide to Pharmacology, DGIdb and drug labels, with manual review. Routes are direct, ligand, complex partner and gene silencing; directions are increase, decrease, depletion and unknown. Depleting agents, receptor-mediated delivery vehicles and links resting on text mining alone are excluded. The gate excludes only the opposite direction, so a record whose direction is unknown is retained; records added by hand are listed in `result_tables/target_funnel_and_drug_evidence/curation_manual_20261122.csv`.

## Figures

Each panel is one script under `papers/figures_v2_20261119/figN/`. Panel a of every figure is rendered in Blender and labelled in R; all other panels are drawn in R, except the evidence wheel (`fig5/c_wheel.py`, pycirclize).

```bash
cd papers/figures_v2_20261119
. _system/env.sh                                                  # sets the R and Python interpreters
$FIGV2_RSCRIPT _system/render_panel.R fig2/b_leaderboard.R 2x2    # render one panel at a block size
$FIGV2_PY _system/assemble.py _system/layout_fig2.json            # assemble a page
_system/run_blender.sh fig2/blender 2a                            # re-render a Blender schematic
$FIGV2_PY fig4/brain/render_brain.py --out fig4/brain             # brain sections of Figure 4j
$FIGV2_PY fig5/c_wheel.py                                         # evidence wheel
```

As with the Blender schematics, the rendered rasters are not in this repository: `fig4/brain/render_brain.py` rebuilds the four sections of Figure 4j from the MNI152 template, the Talairach and Harvard–Oxford atlases and the UK Biobank group-ICA maps (`rfMRI_ICA_d100.nii.gz`, `edge_list_d100.txt` and `rfMRI_GoodComponents_d100_v1.txt` from https://www.fmrib.ox.ac.uk/ukbiobank/).

`_system/check_panel.py` checks panel size, fonts and bar charts. The page is a 4 × 5 grid of 41.25 mm units with 3 mm gaps (174 × 218.25 mm); all text is 7 pt Arial with bold axis titles.

## Environment and paths

Python dependencies are listed in `requirements.txt`; colocalization uses the R packages `coloc` and `susieR`. The scripts read inputs and write results under absolute paths of the authors' server (`/public/home/...`). To run them elsewhere, set the roots in `pdv3/paths.py`, `papers/01_mr_imaging_scan/common.py` and `papers/figures_v2_20261119/_system/common.R` and `env.sh`.

## Data sources

- Plasma protein pQTLs: UKB-PPP (Sun et al., Nature 2023; Supplementary Tables 3, 9 and 16); AGES serum proteome summary statistics (Gudjonsson et al., Nature Communications 2022).
- Brain imaging GWAS: BIG40, discovery and replication samples (Smith et al., Nature Neuroscience 2021); QSM and T2\* (Wang et al., Nature Neuroscience 2022).
- Soluble transferrin receptor GWAS on a clinical assay: Allara et al., Communications Biology 2024 (GCST90451684); iron indices: Bell et al., Communications Biology 2021.
- Linkage disequilibrium reference: 1000 Genomes, European panel.
- Drug evidence: ChEMBL, Open Targets, DrugCentral, Guide to Pharmacology, DGIdb and drug labels.
- Single-cell, spatial, CITE-seq, pig and mouse data: as listed in the parent README.
