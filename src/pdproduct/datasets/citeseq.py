"""Hao 2021 CITE-seq bridge data access (read-only reference, D9/ADR-0001).

Source: <asset_root>/processed/citeseq_hao/bridge_data.npz (in-project rebuild, pd_edges, 2026-09-09)
161,764 cells x ~2000 genes (RNA; HVG2000 + pd-edge protected) + 228 ADT proteins, 8 donors, 6 cell types.
Reference only — never written to; PD-side use requires independent validation.
Paths must be passed explicitly (G0: no hardcoded absolute paths).
"""
from __future__ import annotations

import numpy as np


def load_hao_bridge(path) -> dict:
    """path is required (str or Path); no default absolute path."""
    z = np.load(path, allow_pickle=True)
    return {
        "X": z["X"],                # (n, g) RNA features (z-scored upstream)
        "Y": z["Y"],                # (n, 228) ADT values
        "ct": z["ct"],              # (n,) cell type index
        "ct_labels": [str(s) for s in z["ct_labels"]],
        "genes": [str(s) for s in z["genes"]],
        "proteins": [str(s) for s in z["proteins"]],
        "donor": z["donor"].astype(str),
        "gene_mu": z["gene_mu"],
        "gene_sd": z["gene_sd"],
        "expr_mean_ct": z["expr_mean_ct"],   # (6, g) historical per-ct RNA means
        "adt_mean_ct": z["adt_mean_ct"],     # (6, 228) historical per-ct ADT means
    }


def donor_split(donors: np.ndarray, n_test: int = 3, seed: int = 42) -> dict[str, list[str]]:
    """Deterministic donor split: unique donors sorted, RNG-shuffled once."""
    uniq = sorted(set(donors.tolist()))
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(uniq))
    test_idx = set(order[:n_test].tolist())
    train = [uniq[i] for i in range(len(uniq)) if i not in test_idx]
    test = [uniq[i] for i in sorted(test_idx)]
    return {"train": train, "test": test}
