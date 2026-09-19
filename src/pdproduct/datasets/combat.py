"""COMBAT CITE-seq atlas extraction (read-only reference, D9).

Source atlas: 836,148 cells x 37,694 features (RNA + AB_ ADT merged matrix).
Extracts a donor-stratified subsample with:
  X = top-N expressed RNA genes (official normalized X)
  Y = all ADT (AB_*) channels — MEASURED values, unlike the legacy
      secretome npz whose Y is an inferred target (never use as truth)
Result cached as npz in our asset root for benchmark reuse.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


def extract_combat_bridge(
    h5ad_path: Path,
    out_npz: Path,
    *,
    n_genes: int = 2000,
    cells_per_donor: int = 1500,
    seed: int = 42,
) -> dict:
    import anndata as ad
    import scipy.sparse as sp

    adata = ad.read_h5ad(h5ad_path)  # full load; host has 1TB RAM
    is_adt = np.array([str(v).startswith("AB_") for v in adata.var_names])
    rna_idx = np.where(~is_adt)[0]
    adt_idx = np.where(is_adt)[0]

    # top RNA genes by mean expression (computed on the full matrix)
    Xall = adata.X if not sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    rna_means = np.asarray(Xall[:, rna_idx].mean(axis=0)).ravel()
    top = rna_idx[np.argsort(-rna_means)[:n_genes]]
    top = np.sort(top)

    donors = adata.obs["COMBAT_ID"].astype(str).values
    rng = np.random.default_rng(seed)
    keep = []
    for d in np.unique(donors):
        idx = np.where(donors == d)[0]
        take = min(cells_per_donor, len(idx))
        keep.append(rng.choice(idx, size=take, replace=False))
    keep = np.sort(np.concatenate(keep))

    X = Xall[keep][:, top]
    Y = Xall[keep][:, adt_idx]
    if sp.issparse(X):
        X = np.asarray(X.todense())
    if sp.issparse(Y):
        Y = np.asarray(Y.todense())

    ct_str = adata.obs["cell_type"].astype(str).values[keep]
    labels = sorted(set(ct_str.tolist()))
    ct = np.searchsorted(np.array(labels), ct_str)
    donor_out = donors[keep]
    source = adata.obs["Source"].astype(str).values[keep]

    out_npz = Path(out_npz)
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out_npz,
        X=X.astype(np.float32),
        Y=Y.astype(np.float32),
        ct=ct,
        ct_labels=np.array(labels),
        genes=np.array([str(v) for v in adata.var_names[top]]),
        proteins=np.array([str(v) for v in adata.var_names[adt_idx]]),
        donor=donor_out,
        source=source,
    )
    return {
        "n_cells": int(X.shape[0]),
        "n_genes": int(X.shape[1]),
        "n_adt": int(Y.shape[1]),
        "n_donors": int(len(set(donor_out.tolist()))),
        "cell_types": labels,
        "source_counts": {s: int((source == s).sum()) for s in set(source.tolist())},
        "out_npz": str(out_npz),
    }
