# Reproducing the current PeriScope analysis

## What can be checked immediately

Run `python reproduce/verify_results.py` from the repository root. This standard-library-only check reads released tables without modifying them. It recomputes the 576 MR estimates, uncertainty and global BH values; rebuilds the 16 priority programs' gene unions from the full program definitions; and checks recipient counts, astrocyte benchmark means and animal-omission ranges.

The base-R command in the README reruns harmonization, UKB-first exposure selection and MR from the supplied exposure and brain sentinel tables. `results/blood_brain_MR_20261003/all_brain_MR_estimates.csv` contains the complete selected test family, not just significant results.

`results/optimization_20261003/mr_parallel_gate_20261003/priority_all_programs.csv` supplies the complete 16-program evidence table. `all_programs_parallel_gate.csv` supplies all 6,849 programs, including nonselected programs. Columns retained from earlier diagnostics do not define the current nomination: use `priority_pass`, the C1–C6 fields and `MR_gate` for the current result. `joint_analysis/programs.json` contains exact memberships.

`figure_source_data/` uses the project figure numbering 1–7. In the Research manuscript, project Figure 2 is supplementary Figure S1; project Figures 3–7 become main Figures 2–6. No numerical values change with journal numbering.

## Execution order for rebuilding from original inputs

The preserved analysis scripts use two original roots:

- `/public/home/mengxl/dzy/pd_product`: the repository.
- `/public/home/mengxl/dzy/pd_product_assets`: downloaded data, prepared arrays, checkpoints and analysis outputs.

On another machine, use `python reproduce/configure_paths.py --assets /absolute/path/to/pd_product_assets --output /absolute/path/to/new_execution_copy`. This creates a separate source execution copy and replaces these two root strings; the released checkout remains unchanged. It also substitutes the archived Python executable with the current interpreter. Additional external R/library and raw-data paths must be configured for the local environment before full execution.

1. **Inputs.** Obtain the public resources in the README and prepare donor-separated RNA and CITE-seq inputs using `scripts/preprocess_*` and `_work/rescreen_20261002/prepare_inputs.py`. Preserve donor splits, annotation mappings and the 4,741-gene/228-ADT vocabulary. Source metadata and `results/optimization_20261003/data_checks.json` describe the prepared input dimensions.
2. **Models.** `_work/optimization_20261003/train_ccwm.py` and `_work/checkpoint_repair_20261003/train_ccwm.py` contain the current training implementation. `model_registry.json` identifies the authoritative ten checkpoints, including repaired 47–49. The large weights and prepared human matrices must be supplied separately; the repository is not a weights download.
3. **Dose programs and confirmation.** `_work/optimization_20261003/run_full.py`/`run_joint_analysis.py` route through the registered model outputs and current downstream refresh. Discovery seeds 42–46 define memberships; confirmation seeds 47–51 test them. The refresh scripts expect the prepared assets and frozen discovery definitions to exist; they are not one-command raw-data installers.
4. **Human benchmark.** `_work/benchmark_shared_20261003/current_benchmark.py` and `shared_pipeline.py` implement the current comparison and readouts. Repaired 47–49 are the PeriScope representation seeds; benchmark run and comparator seed labels are separate. CellOT retains its original shared tissue PCA representation. Use `_work/cellot_results_recovery_20261003/run_verified_recipe.py` for its verified original recipe, preserving the recorded CPU/BLAS settings and expected results. Do not substitute new hyperparameters for the restored result.
5. **MR and nomination.** `_work/blood_brain_MR_20261003/` contains acquisition and analysis code. `reproduce/run_brain_mr.R` performs the portable table-level calculation. `_work/optimization_20261003/apply_mr_parallel_gate.py` intersects MR with perturbation and C1–C6 using the complete result family.
6. **Human tissue and spatial follow-up.** `_work/optimization_20261003/brain_competitive.py`, `independent_validation.py` and the subtype-aware `_work/spatial_subtype_20261003/` scripts preserve donor and recipient identity. Spatial support is reported separately from the final MR-plus-perturbation nomination.
7. **Independent response transfer.** `_work/pig_external_validation_20261003/prepare_adapter.py`, `predict_external.py` and `evaluate_external.py` implement the Ensembl ortholog adapter, frozen predictions and animal-level comparison. The pig RNA inputs, released prediction arrays, per-gene effects and resampling results are under `results/pig_external_validation_20261003/`. Regenerating predictions requires the human model weights and reference inputs; the portable verification does not retrain or re-predict.

## Preservation and verification

`release_manifest.json` records SHA-256 hashes for released code and result files. This release preserves current scientific code and numerical outputs; obsolete migration utilities and local dependency caches are excluded. It does not claim that every historical helper is a supported entry point. Use the entry points above and their actual input requirements.

The publication check covered Python source syntax, the portable numerical verification and a fresh base-R MR calculation. Full GPU retraining, fresh external-data downloading, complete RCTD/CSIDE execution and all figure rendering were not rerun for this repository update.
