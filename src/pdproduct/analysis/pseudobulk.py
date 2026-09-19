"""Donor-level pseudobulk differential expression.

Aggregation: raw counts summed per (donor, cell_type); statistics computed
on donors (8 Ctrl vs 10 Disease) — never on cells (pseudo-replication trap,
HANDOFF_SERVER_V1 §13.2). Discovery-grade screening only: association
claims, no causal wording.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def gene_symbols(adata) -> np.ndarray:
    """Preferred display symbols: gene_symbol column when present, else var_names."""
    if "gene_symbol" in adata.var.columns:
        sym = adata.var["gene_symbol"].astype(str).values
        base = np.array([str(v) for v in adata.var_names])
        return np.where(np.isin(sym, ["", "nan", "None"]), base, sym)
    return np.array([str(v) for v in adata.var_names])


def pseudobulk_by_donor(X, gene_names, donor, cell_type) -> pd.DataFrame:
    """Sum counts per (cell_type, donor); return long DataFrame of raw count matrices."""
    records = []
    for ct in sorted(set(cell_type.tolist())):
        m_ct = cell_type == ct
        for d in sorted(set(donor[m_ct].tolist())):
            m = m_ct & (donor == d)
            s = np.asarray(X[m].sum(axis=0)).ravel()
            records.append((ct, d, s))
    genes = np.array(gene_names)
    return records, genes


def de_donor_level(
    records, genes, condition_of, cell_type_name, min_count=10, min_donor_frac=0.5,
    ctrl_label="Ctrl", disease_label="Disease",
) -> pd.DataFrame:
    """Welch t-test on log1p-CPM per gene across donors, BH-adjusted.

    Filters: gene kept if total pseudobulk counts >= min_count in the
    cell type and expressed (>0) in >= min_donor_frac of donors per group.
    """
    group = [(ct, d, s) for ct, d, s in records if ct == cell_type_name]
    donors = [d for _, d, _ in group]
    M = np.stack([s for _, _, s in group])  # donors x genes raw counts
    lib = M.sum(axis=1, keepdims=True)
    cpm = M / np.clip(lib, 1, None) * 1e6
    logcpm = np.log1p(cpm)

    cond = np.array([condition_of[d] for d in donors])
    a = cond == ctrl_label
    b = cond == disease_label
    keep = (M.sum(axis=0) >= min_count)
    keep &= (M[a] > 0).mean(axis=0) >= min_donor_frac
    keep &= (M[b] > 0).mean(axis=0) >= min_donor_frac

    idx = np.where(keep)[0]
    t, p = stats.ttest_ind(logcpm[b][:, idx], logcpm[a][:, idx], axis=0, equal_var=False)
    log2fc = np.log2((cpm[b][:, idx].mean(axis=0) + 1) / (cpm[a][:, idx].mean(axis=0) + 1))

    res = pd.DataFrame(
        {"gene": genes[idx], "log2fc_disease_vs_ctrl": log2fc, "p": p}
    ).sort_values("p")
    # Benjamini-Hochberg
    n = len(res)
    order = res["p"].rank(method="first")
    res["q"] = (res["p"] * n / order).clip(upper=1.0)
    res["n_donors"] = f"{int(a.sum())}v{int(b.sum())}"
    return res.reset_index(drop=True)
