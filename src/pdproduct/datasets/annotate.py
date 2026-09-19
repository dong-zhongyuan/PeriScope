"""Marker-based cell type annotation for substantia nigra snRNA-seq.

Transparent scoring: mean log-normalized expression of curated markers per
cell; argmax with a minimum-score floor for 'unassigned'. Markers follow
published SN snRNA conventions (Kamath 2022 and related).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

SN_MARKERS: dict[str, list[str]] = {
    "DA_neuron": ["TH", "SLC6A3", "DDC", "KCNJ6"],
    "Neuron_other": ["RBFOX3", "MAP2", "SNAP25"],
    "Astrocyte": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP"],
    "OPC": ["PDGFRA", "VCAM1", "CSPG4"],
    "Microglia": ["AIF1", "CX3CR1", "CSF1R", "P2RY12"],
    "Pericyte": ["PDGFRB", "RGS5", "CSPG4"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5"],
    "T_cell": ["TRAC", "CD3D", "CD3E"],
    "B_cell": ["MS4A1", "CD79A"],
    "Mono_macro": ["LYZ", "CD14", "FCGR3A"],
}

# genes shared by multiple lines must not double-count (CSPG4: OPC & Pericyte)


def score_cell_types(adata, markers: dict[str, list[str]] | None = None, min_score: float = 0.05):
    """Add obs['cell_type'] and per-type scores. Uses .layers['lognorm'] if present else .X.

    Gene lookup covers both var_names (symbol or ENSG depending on source)
    and the gene_symbol column when present.
    """
    markers = markers or SN_MARKERS
    X = adata.layers["lognorm"] if "lognorm" in adata.layers else adata.X
    gene_index = {str(g): i for i, g in enumerate(adata.var_names)}
    if "gene_symbol" in adata.var.columns:
        for i, sym in enumerate(adata.var["gene_symbol"].astype(str)):
            gene_index.setdefault(sym, i)
    present = {ct: [g for g in gs if g in gene_index] for ct, gs in markers.items()}
    missing = {ct: [g for g in gs if g not in gene_index] for ct, gs in markers.items()}
    scores = {}
    for ct, genes in present.items():
        if not genes:
            continue
        cols = [gene_index[g] for g in genes]
        sub = X[:, cols]
        scores[ct] = np.asarray(sub.mean(axis=1)).ravel()
    score_df = pd.DataFrame(scores, index=adata.obs_names)
    for ct in score_df.columns:
        adata.obs[f"score_{ct}"] = score_df[ct]
    best = score_df.idxmax(axis=1)
    floor_ok = score_df.max(axis=1) >= min_score
    adata.obs["cell_type"] = np.where(floor_ok, best, "unassigned")
    return {
        "n_per_type": adata.obs["cell_type"].value_counts().to_dict(),
        "markers_present": {ct: g for ct, g in present.items()},
        "markers_missing": {ct: g for ct, g in missing.items() if g},
    }
