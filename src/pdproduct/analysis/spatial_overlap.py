"""Spatial colocalization first check (ROI-level, donor-stratified).

Question (fixed before looking): do the two-stage replicated glial candidate
gene sets colocalize with their cell-type marker scores across GeoMx ROIs?
Statistic: within-donor Spearman(candidate-set score, marker score), averaged
over donors; null = same-size random gene sets (permutation).
Evidence class: C3 spatial colocalization, region-level only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def roi_lognorm(adata):
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    lib = X.sum(axis=1, keepdims=True)
    return np.log1p(X / np.clip(lib, 1, None) * 1e6)


def set_score(L, gene_index, genes, min_present: int = 4):
    cols = [gene_index[g] for g in genes if g in gene_index]
    missing = [g for g in genes if g not in gene_index]
    if len(cols) < max(min_present, int(0.2 * len(genes))):
        return None, missing
    return L[:, cols].mean(axis=1), missing


def _expression_bins(L, n_bins: int):
    """Quantile bins of per-gene mean ROI lognorm expression (Amendment 2)."""
    mean_expr = L.mean(axis=0)
    edges = np.quantile(mean_expr, np.linspace(0, 1, n_bins + 1)[1:-1])
    return np.digitize(mean_expr, edges)


def colocalization_test(adata, candidate_genes: list[str], marker_genes: list[str],
                        n_perm: int = 200, seed: int = 42,
                        control_library_size: bool = False, top_n: int | None = None,
                        expr_match_bins: int | None = None, min_marker_genes: int = 4,
                        L=None) -> dict:
    """Optional library-size partial-out (removes GeoMx background correlation)
    and top-N candidate subsetting.

    expr_match_bins: if set (e.g. 20), the permutation null draws size-matched
    random gene sets with the same expression-bin composition as the candidate
    set (bins = quantiles of per-gene mean ROI lognorm expression), instead of
    drawing from the whole genome.
    min_marker_genes: minimum number of marker genes that must be present
    (candidate side keeps the default minimum of 4).
    L: optional precomputed roi_lognorm(adata) matrix (identical values; saves
    recomputation when calling the test many times)."""
    if L is None:
        L = roi_lognorm(adata)
    gene_index = {str(g): i for i, g in enumerate(adata.var_names)}
    if top_n is not None:
        candidate_genes = candidate_genes[:top_n]
    cand, miss_c = set_score(L, gene_index, candidate_genes)
    mark, miss_m = set_score(L, gene_index, marker_genes, min_present=min_marker_genes)
    if cand is None or mark is None:
        return {"error": f"insufficient genes present: candidates missing {len(miss_c)}, markers missing {len(miss_m)}"}

    if control_library_size:
        lib = np.log(np.asarray(adata.X.sum(axis=1)).ravel() + 1)
        lib = (lib - lib.mean()) / (lib.std() + 1e-9)

        def resid(y):
            beta = np.dot(lib, y) / (np.dot(lib, lib) + 1e-9)
            return y - beta * lib

        cand, mark = resid(cand), resid(mark)

    donor = adata.obs["donor"].values.astype(str)
    def per_donor_rho(a, b):
        rhos = {}
        for d in set(donor.tolist()):
            m = donor == d
            if m.sum() > 10 and np.std(a[m]) > 0 and np.std(b[m]) > 0:
                rhos[d] = stats.spearmanr(a[m], b[m]).statistic
        return rhos

    obs_rhos = per_donor_rho(cand, mark)

    rng = np.random.default_rng(seed)
    n_cand_cols = len([gene_index[g] for g in candidate_genes if g in gene_index])
    if expr_match_bins is not None:
        bins = _expression_bins(L, expr_match_bins)
        cand_cols = [gene_index[g] for g in candidate_genes if g in gene_index]
        bin_counts = pd.Series(bins[cand_cols]).value_counts().to_dict()
        bin_members = {b: np.where(bins == b)[0] for b in bin_counts}

        def draw_null_cols():
            return np.concatenate([rng.choice(bin_members[b], size=n, replace=False)
                                   for b, n in bin_counts.items()])
    else:
        pool = np.array([i for i in range(L.shape[1])])

        def draw_null_cols():
            return rng.choice(pool, size=n_cand_cols, replace=False)

    null_means = []
    for _ in range(n_perm):
        cols = draw_null_cols()
        rand_score = L[:, cols].mean(axis=1)
        if control_library_size:
            rand_score = resid(rand_score)
        null_rhos = per_donor_rho(rand_score, mark)
        if null_rhos:
            null_means.append(np.mean(list(null_rhos.values())))
    obs_mean = float(np.mean(list(obs_rhos.values())))
    p_emp = float(np.mean([nm >= obs_mean for nm in null_means])) if null_means else None
    return {
        "n_candidate_genes_present": n_cand_cols,
        "candidate_genes_missing": miss_c,
        "n_marker_genes_present": len([g for g in marker_genes if g in gene_index]),
        "per_donor_rho": {k: float(v) for k, v in obs_rhos.items()},
        "mean_rho": obs_mean,
        "n_donors_scored": len(obs_rhos),
        "null_mode": f"expression_matched_{expr_match_bins}bins" if expr_match_bins else "genomewide",
        "null_mean_of_means": float(np.mean(null_means)) if null_means else None,
        "null_p95": float(np.percentile(null_means, 95)) if null_means else None,
        "perm_p": p_emp,
        "n_perm": n_perm,
    }
