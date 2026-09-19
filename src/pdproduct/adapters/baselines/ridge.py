"""Ridge baseline (2nd of the 2+1 first-version scope).

Multi-output ridge regression RNA -> ADT. Alpha selected by grouped
cross-validation over TRAINING donors only (test donors untouched).
"""
from __future__ import annotations

import numpy as np
from sklearn.linear_model import Ridge


class RidgeBaseline:
    def __init__(self, alphas=(1.0, 10.0, 100.0, 1000.0, 3000.0, 10000.0)) -> None:
        self.alphas = tuple(alphas)
        self.model: Ridge | None = None
        self.selected_alpha: float | None = None

    def fit(self, X: np.ndarray, Y: np.ndarray, groups: np.ndarray, cv_folds: int | None = None) -> "RidgeBaseline":
        """Grouped CV over training donors for alpha; refit on all training cells.

        cv_folds: number of donor groups (grouped K-fold); None = leave-one-
        donor-out. Grouped K-fold is preferred for large donor counts.
        """
        best_alpha, best_score = None, -np.inf
        uniq = np.unique(groups)
        if cv_folds is not None and cv_folds < len(uniq):
            rng = np.random.default_rng(0)
            assignment = {d: i % cv_folds for i, d in enumerate(rng.permutation(uniq))}
            folds = [
                np.array([d for d in uniq if assignment[d] == k]) for k in range(cv_folds)
            ]
        else:
            folds = [np.array([held]) for held in uniq]
        for alpha in self.alphas:
            fold_scores = []
            for held in folds:
                tr = ~np.isin(groups, held)
                va = ~tr
                m = Ridge(alpha=alpha).fit(X[tr], Y[tr])
                pred = m.predict(X[va])
                # per-donor mean rho across proteins (same metric as evaluation)
                from scipy.stats import spearmanr

                rhos = [
                    spearmanr(Y[va, j], pred[:, j]).statistic
                    for j in range(Y.shape[1])
                    if np.std(Y[va, j]) > 0 and np.std(pred[:, j]) > 0
                ]
                fold_scores.append(np.mean(rhos))
            score = float(np.mean(fold_scores))
            if score > best_score:
                best_score, best_alpha = score, alpha
        self.selected_alpha = best_alpha
        self.model = Ridge(alpha=best_alpha).fit(X, Y)
        self.cv_score = best_score
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        return self.model.predict(X)
