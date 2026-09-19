"""Cell-type mean baseline (the mandatory first baseline of the 2+1 scope).

Prediction = per-cell-type mean ADT profile, computed on training donors only.
Two evaluation modes:
  oracle_ct   — cell type from ground truth (upper reference for the family)
  rna_ct      — cell type classified from RNA via cosine to class means
Null control: donor-internal pairing shuffle of Y (must collapse to ~0 rho).
"""
from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


class CellTypeMeanModel:
    def __init__(self) -> None:
        self.adt_means: np.ndarray | None = None      # (k, n_protein)
        self.rna_class_means: np.ndarray | None = None  # (k, n_gene)

    def fit(self, X: np.ndarray, Y: np.ndarray, ct: np.ndarray) -> "CellTypeMeanModel":
        k = int(ct.max()) + 1
        self.adt_means = np.stack([Y[ct == j].mean(axis=0) for j in range(k)])
        norms = np.linalg.norm(X, axis=1, keepdims=True)
        Xn = X / np.clip(norms, 1e-8, None)
        self.rna_class_means = np.stack(
            [(Xn[ct == j] / np.clip(np.linalg.norm(Xn[ct == j], axis=1, keepdims=True), 1e-8, None)).mean(axis=0) for j in range(k)]
        )
        self.rna_class_means /= np.clip(np.linalg.norm(self.rna_class_means, axis=1, keepdims=True), 1e-8, None)
        return self

    def classify(self, X: np.ndarray) -> np.ndarray:
        norms = np.linalg.norm(X, axis=1, keepdims=True)
        sims = (X / np.clip(norms, 1e-8, None)) @ self.rna_class_means.T
        return sims.argmax(axis=1)

    def predict(self, X: np.ndarray, ct: np.ndarray | None = None) -> np.ndarray:
        if ct is None:
            ct = self.classify(X)
        return self.adt_means[ct]


def per_donor_rho(Y_true: np.ndarray, Y_pred: np.ndarray, donor: np.ndarray) -> dict[str, dict]:
    """Per-donor mean Spearman rho across features (cell-level ranking).

    Constant-column filter uses exact range (ptp): float32 columns that are
    literally constant can show std ~1e-10 from rounding noise and would pass
    a std>0 filter, then poison the mean with spearmanr's constant-input NaN.
    """
    out = {}
    for d in sorted(set(donor.tolist())):
        m = donor == d
        rhos = []
        for j in range(Y_true.shape[1]):
            if np.ptp(Y_true[m, j]) > 0 and np.ptp(Y_pred[m, j]) > 0:
                rhos.append(spearmanr(Y_true[m, j], Y_pred[m, j]).statistic)
        out[d] = {
            "mean_rho": float(np.mean(rhos)),
            "n_proteins_scored": len(rhos),
        }
    return out


def pairing_shuffle_null(Y: np.ndarray, donor: np.ndarray, rng: np.random.RandomState) -> np.ndarray:
    """Shuffle cell pairing within each donor; breaks RNA-protein pairing only."""
    Ys = Y.copy()
    for d in set(donor.tolist()):
        idx = np.where(donor == d)[0]
        Ys[idx] = Ys[rng.permutation(idx)]
    return Ys
