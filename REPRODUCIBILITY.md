# Reproducing the NMI submission analyses

The release distinguishes recalculating statistics from frozen inputs, running trained models, and rebuilding the entire analysis from raw data. These operations require different inputs.

## 1. Integrity and numerical checks

From the repository root:

```bash
python reproduce/verify_manifest.py
python reproduce/verify_results.py
```

The 33 standard-library numerical checks rebuild target-level gene unions from the program definitions, verify recipient counts and directions, independently calculate all 576 Wald estimates and global BH values, and check the benchmark and pig animal-omission summaries.

`results/optimization_20261003/mr_parallel_gate_20261003/priority_all_programs.csv` contains the 16 priority programs. `all_programs_parallel_gate.csv` retains all 6,849 programs. Current nomination uses `priority_pass`, C1–C6 and `MR_gate`; legacy diagnostic columns are not additional eligibility thresholds. Exact memberships are in `joint_analysis/programs.json`.

## 2. Revised analyses on CPU

Use Python >=3.10 with `reproduce/requirements-numerical.txt`:

```bash
python reproduce/recompute_nmi.py --output recomputed/nmi --ma
```

| Output | Inputs and calculation |
| --- | --- |
| `mgap_matched.csv` | Frozen gene contrasts; all seven methods at benchmark seeds 42–44, four equal-weight routes per recipient; reproduces Supplementary Table 1 |
| `mgap_allmethods.csv` | Additional available-seed calculation: ten PeriScope models and three runs per comparator; preserved separately |
| `benchmark_pooled_all_available.csv` | Pooled contrasts for the additional available-seed calculation, with explicit seed counts |
| `donor_dependence.csv` | Held-out donor metadata and `bloodonly_profiles.npz`; recomputes group RMS, donor RMS and their ratio for eight routes |
| `mr_steiger.csv` | Released instrument/outcome estimates; repeats the original calculation for all 48 target–outcome pairs |
| `ma_cohort_test.csv` | Released donor × subtype pseudobulk counts; reruns the original 5,000 expression-matched null sets, including all 22 tested programs; requires `--ma` |
| `ma_spatial_complement.csv` | Reconstructs the 16-program join and support indicators from the released Ma and spatial tables; this is not a new spatial deconvolution |
| `verification.json` | Dependency versions and comparison of each recomputed table with the archived result |

The compact input `reproduce/inputs/benchmark_contrasts.npz` contains gene-wise PD-minus-control mean contrasts exported from the existing prediction arrays. The original float precision and route/seed identifiers are retained. `export_provenance.json` records hashes of the source prediction files and exported inputs. Recomputing these statistics does not regenerate predictions.

The primary matched comparison produces PeriScope pooled rho **0.2589259866**, rounded to **0.259** in the manuscript. The separate all-available-seed analysis produces **0.2376738206** for PeriScope. They use different seed sets. The released runner checks both, retaining the original analyses rather than silently substituting one for the other.

The original Ma, Steiger and complement calculations are executed with only input/output locations redirected. The benchmark and donor-dependence wrappers apply the same arithmetic to compact public inputs and compare every output column with the original tables. Output folders must be new and cannot be inside `results/`.

## 3. MR from harmonized input tables

This uses base R and writes to a separate directory:

```bash
mkdir -p recomputed/mr
cp results/blood_brain_MR_20261003/combined_exposure_instruments.csv recomputed/mr/
cp results/blood_brain_MR_20261003/brain_sentinel_statistics.csv recomputed/mr/
Rscript reproduce/run_brain_mr.R recomputed/mr
```

The runner performs allele harmonization, UKB-first exposure selection, Wald estimates, uncertainty and global BH correction across all 576 selected tests. It does not require GWAS re-downloading. The complete family is preserved, including nonsignificant estimates.

## 4. Models, vocabulary and seed identity

```bash
python reproduce/download_checkpoints.py --output /absolute/path/to/periscope_models
python -m pip install torch==2.5.1
python reproduce/check_model_loading.py --models /absolute/path/to/periscope_models
```

`reproduce/checkpoint_release.json` specifies the three public release parts, their sizes and SHA-256 hashes, the complete archive hash and expected members. The downloader verifies all ten checkpoint hashes against `results/optimization_20261003/model_registry.json` before extraction. It can also verify a local download using `--archive /path/to/periscope-models-seeds42-51.tar.gz`.

The model-loading check imports the preserved implementation, strictly loads every state dictionary, checks the feature/context dimensions, and runs a synthetic blood-to-brain forward calculation on CPU. It checks software compatibility; it is not a scientific performance evaluation. The archive contains model parameters/configuration, not the complete training matrices or every comparator checkpoint.

| Use | Seeds |
| --- | --- |
| Discovery programs | 42–46 |
| Confirmation tests | 47–51 |
| Authoritative ten-model set | 42–51, including repaired 47–49 |
| Primary human comparison and Supplementary Table 1 | 42–44 for all seven methods |
| Additional ten-model sensitivity and donor dependence | 42–51 |

`current_benchmark.py`, its loader, the matched-comparison script and the current supplementary methods identify 42–44 for the three-run comparison.

The vocabulary is 4,741 genes and 228 protein measurements. `reproduce/inputs/*_preprocessing.npz` provides the gene order, state vocabulary and recorded normalization metadata. Original prepared matrices are approximately 17 GB across `blood.npz`, `brain_train.npz`, `brain_calibration.npz`, `brain_locked_test.npz`, `citeseq.npz` and `bridge_data.npz`; those cell-level matrices are not bundled here.

## 5. Execution order from original inputs

The preserved scripts use the original roots `/public/home/mengxl/dzy/pd_product` and `/public/home/mengxl/dzy/pd_product_assets`. To prepare a separate execution copy:

```bash
python reproduce/configure_paths.py \
  --assets /absolute/path/to/pd_product_assets \
  --output /absolute/path/to/new_execution_copy
```

This copies source and frozen results, rewrites known root paths in the copy, and substitutes the current Python interpreter. It does not download missing data, install R libraries or resolve every historical helper's dependency. Use the following analysis entry points and inspect their actual input arguments before execution.

1. **Inputs:** obtain the resources linked in the README. Use `scripts/preprocess_*` and `_work/rescreen_20261002/prepare_inputs.py`, retaining donor splits, annotations and the fixed feature vocabulary. `results/optimization_20261003/data_checks.json` records the prepared dimensions. Spatial analysis additionally requires the author-annotated reference and Visium data.
2. **Model training:** `_work/optimization_20261003/train_ccwm.py` and `_work/checkpoint_repair_20261003/train_ccwm.py` preserve the training implementation. To use frozen models, install the release asset under `pd_product_assets/results/optimization_20261003/models`. The repaired canonical files already match the registry. A historical `checkpoint_overrides.json` may refer to a separate repair directory; supply the same hash-verified files there or configure those entries to the canonical files in the execution copy.
3. **Dose queries and confirmation:** `_work/optimization_20261003/run_full.py`, `run_joint_analysis.py` and the registered refresh scripts use prepared arrays, models and frozen discovery definitions. They are downstream runners, not raw-data download installers.
4. **Human comparison:** `_work/benchmark_shared_20261003/current_benchmark.py` and `shared_pipeline.py` use benchmark seeds 42–44. Required protocol files are included. Follow `_work/cellot_results_recovery_20261003/run_verified_recipe.py` for the original CellOT recipe, including its recorded environment and CPU/BLAS settings.
5. **MR and nomination:** `_work/blood_brain_MR_20261003/` provides acquisition and analysis code. `_work/optimization_20261003/apply_mr_parallel_gate.py` combines MR with computational perturbation and C1–C6 over the complete program family.
6. **Brain and spatial evidence:** `brain_competitive.py`, `independent_validation.py`, `_work/spatial_subtype_20261003/` and `_work/nmi_revision_20261017/` preserve the cohort, subtype and donor analyses. The latter directory includes Ma donor robustness, recipient-enriched spatial tests and lineage-matched nulls. The CPU runner above covers only its listed calculations.
7. **Pig response transfer:** `_work/pig_external_validation_20261003/prepare_adapter.py`, `predict_external.py` and `evaluate_external.py` implement the one-to-one ortholog adapter and animal-level comparison. Released input tables, frozen predictions and summaries are under `results/pig_external_validation_20261003/`. Fresh predictions also require the human reference arrays and trained models.
8. **Figures:** `figure_source_data/` is the current export. `figures/` preserves plotting source, copied R factories and layouts; see [figures/README.md](figures/README.md) for external assets and execution requirements. Project Figure 5 is NMI Supplementary Fig. 1; project Figures 6 and 7 are NMI main Figures 5 and 6. Project Figures 1–4 retain their main-figure numbers.

## 6. Verification record and remaining full-run work

The October 2026 maintenance compared all 907 previously tracked files with their server counterparts. Fourteen figure-data files required updates. Eighty-one later analysis scripts and result files were added from the current supplementary analysis archive and verified byte-for-byte against the server. Core scientific implementations and existing result files were preserved.

Checks performed for this release include the 33 original numerical checks, current CPU recomputations including the Ma null test, independent base-R MR recalculation, all ten checkpoint hashes and CPU model-loading/inference checks. See `reproduce/validation_20261027.json` for the recorded results and runtime versions.

A full fresh raw-data download, GPU retraining, refitting every comparator, spatial RCTD/CSIDE run and end-to-end final figure rendering have **not** been executed from the public checkout. These are the remaining steps needed before claiming complete end-to-end reproduction. The published checks establish the narrower, explicitly listed recomputations and model compatibility.
