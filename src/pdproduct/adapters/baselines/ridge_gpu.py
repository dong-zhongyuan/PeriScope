"""GPU (torch) multi-output ridge — normal equations with shared Gram matrix.

XtX and XtY are computed once on device; each alpha only adds aI and solves
via Cholesky. Orders of magnitude faster than per-alpha sklearn refits on
large multi-output problems. Falls back to CPU tensor if CUDA unavailable.
"""
from __future__ import annotations

import numpy as np


class TorchRidge:
    def __init__(self, alphas=(100.0, 1000.0, 3000.0, 10000.0)) -> None:
        self.alphas = tuple(float(a) for a in alphas)
        self.device = None
        self.coef_ = None  # (n_features, n_targets)
        self.intercept_ = None
        self.selected_alpha: float | None = None

    def _solve(self, XtX, Xty, alpha, n):
        import torch

        A = XtX.clone()
        A.diagonal().add_(alpha)
        L = torch.linalg.cholesky(A)
        return torch.cholesky_solve(Xty, L)

    def fit(self, X: np.ndarray, Y: np.ndarray, groups: np.ndarray | None = None,
            cv_folds: int = 10) -> "TorchRidge":
        import torch

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        n, d = X.shape
        Xt = torch.as_tensor(X, dtype=torch.float32, device=self.device).T.contiguous()
        Yt = torch.as_tensor(Y, dtype=torch.float32, device=self.device)

        mu = Xt.mean(dim=1, keepdim=True)  # (d,1) column means of X
        y_mu = Yt.mean(dim=0, keepdim=True)
        Xtc = Xt - mu
        Ytc = Yt - y_mu
        XtX = Xtc @ Xtc.T
        Xty = Xtc @ Ytc

        # grouped CV over donors for alpha selection
        donors = np.unique(groups) if groups is not None else np.arange(n)
        rng = np.random.default_rng(0)
        assign = {dn: i % max(cv_folds, 1) for i, dn in enumerate(rng.permutation(donors))}
        best_alpha, best_score = self.alphas[0], -np.inf
        from scipy.stats import spearmanr

        for alpha in self.alphas:
            scores = []
            for k in range(max(cv_folds, 1)):
                held = [dn for dn in donors if assign[dn] == k]
                if groups is None:
                    tr = np.ones(n, bool); va = ~tr
                else:
                    va = np.isin(groups, held); tr = ~va
                if va.sum() == 0:
                    continue
                Xtr = torch.as_tensor(X[tr], dtype=torch.float32, device=self.device)
                Ytr = torch.as_tensor(Y[tr], dtype=torch.float32, device=self.device)
                Xva = torch.as_tensor(X[va], dtype=torch.float32, device=self.device)
                mtr = Xtr.mean(0, keepdim=True)
                mytr = Ytr.mean(0, keepdim=True)
                XtXf = (Xtr - mtr).T @ (Xtr - mtr)
                Xtyf = (Xtr - mtr).T @ (Ytr - mytr)
                W = self._solve(XtXf, Xtyf, alpha, n)
                pred = (Xva - mtr) @ W + mytr
                pred_np = pred.cpu().numpy()
                Yva = Y[va]
                rhos = [
                    spearmanr(Yva[:, j], pred_np[:, j]).statistic
                    for j in range(min(Yva.shape[1], 200))  # protein subsample for CV speed
                    if np.std(Yva[:, j]) > 0 and np.std(pred_np[:, j]) > 0
                ]
                scores.append(np.mean(rhos))
            score = float(np.mean(scores))
            if score > best_score:
                best_score, best_alpha = score, alpha

        self.selected_alpha = best_alpha
        self.cv_score = best_score
        W = self._solve(XtX, Xty, best_alpha, n)
        self.coef_ = W.cpu().numpy()
        self.intercept_ = (y_mu.cpu().numpy() - mu.cpu().numpy().T @ self.coef_).ravel()
        return self

    def predict(self, X: np.ndarray) -> np.ndarray:
        import torch

        Xt = torch.as_tensor(X, dtype=torch.float32, device=self.device)
        W = torch.as_tensor(self.coef_, dtype=torch.float32, device=self.device)
        b = torch.as_tensor(self.intercept_, dtype=torch.float32, device=self.device)
        return (Xt @ W + b).cpu().numpy()
