# Figure source snapshot

`source/`, `R/`, the layout/scene manifests and `assets/` are copied from the current figure project on the analysis server. Scientific values are in the repository's `figure_source_data/`. The source snapshot preserves the complete plotting scripts and copied template factories; no plots were reimplemented during repository maintenance.

Project numbering 1–7 is stable. NMI places project Figure 5 in Supplementary Fig. 1, and project Figures 6–7 in main Figures 5–6. Project Figures 1–4 remain main Figures 1–4.

These scripts originally run inside `pd_product_figure_design/manuscript_review/figure_redesign_20261004/`, with `data/`, `R/`, `source/`, `build/` and the JSON layout files as siblings. `render_current.R` uses `PD_PROJECT_ROOT` for the parent project. `prepare_current.py` reads analysis tables through the original project's `delivery/` tree; it is an export step, not an analysis runner.

The full visual build also uses cohort plot objects from `pd_product_figure_design/panel_grid/objects`, cohort tables and UMAP coordinates from `template_semantic_review`, `full_template_reuse/inputs/all_cohort_composition.csv`, and workflow illustrations under `pd_product_figure_design/workflow_assets`. Those external figure-project inputs are not all included in this code snapshot. Some historical helpers have additional dependencies and are retained for provenance rather than advertised as portable entry points.

R dependencies referenced by the main renderer include ggplot2, grid, cowplot, jsonlite, ggrepel, ggalluvial, ComplexHeatmap and circlize. Python assembly uses PyMuPDF, pandas and NumPy; Blender-specific steps require Blender. Reproducing typography also requires the original fonts and rendering settings. Asset license/source metadata is preserved in `assets/bioicons/used_asset_license.json`.

The current numerical figure exports were synchronized and checked against the server. The complete figure build and visual layout were not rerun for this repository release.
