"""OOD risk-coverage evaluation + G5A active-sampling mask-and-reveal.

Implements the two standard evaluation pieces frozen in
docs/G2_PREREGISTRATION.md (§2 risk-coverage, §3 active design).
"""
from __future__ import annotations

import numpy as np
from scipy import stats

from ..adapters.baselines.cell_type_mean import per_donor_rho
from ..adapters.baselines.ridge import RidgeBaseline


def risk_coverage_curve(Y_true, pred_a, pred_b, donor, n_bins: int = 10,
                        uncertainty: np.ndarray | None = None) -> dict:
    """Stratify cell-protein predictions by an uncertainty score.

    Default uncertainty = ensemble disagreement |a-b|; pass `uncertainty`
    (same shape as pred_a) for quantile-interval width etc.
    """
    if uncertainty is None:
        uncertainty = np.abs(pred_a - pred_b)
    disagree = uncertainty
    order = np.argsort(disagree, axis=None)  # global quantile bins over cells x proteins
    n = order.size
    bins = []
    for k in range(n_bins):
        idx = order[int(k * n / n_bins): int((k + 1) * n / n_bins)]
        rows, cols = np.unravel_index(idx, disagree.shape)
        keep_rows, keep_cols = rows, cols
        rhos = []
        for d in np.unique(donor):
            m_d = donor == d
            sel = m_d[keep_rows]
            r_sel, c_sel = keep_rows[sel], keep_cols[sel]
            if len(r_sel) < 30:
                continue
            rv = Y_true[r_sel, c_sel]
            pv = pred_a[r_sel, c_sel]
            if np.std(rv) > 0 and np.std(pv) > 0:
                rhos.append(stats.spearmanr(rv, pv).statistic)
        bins.append(float(np.mean(rhos)) if rhos else None)
    return {"bin_rho_low_to_high_disagreement": bins, "n_bins": n_bins}


def active_sampling_mask_reveal(X_tr, Y_tr, ct_tr, donor_tr, X_te, Y_te, donor_te,
                                budgets=(200, 1000, 5000), seed: int = 42,
                                extra_strategies: bool = True) -> dict:
    """Compare random vs diversity vs residual-margin cell selection.

    Ridge trained on the selected labeled subset; evaluated on held-out
    donors with the standard per-donor rho metric.
    """
    from sklearn.cluster import KMeans

    rng = np.random.default_rng(seed)
    results = {}
    for budget in budgets:
        arm = {}
        # random
        idx_rand = rng.choice(len(X_tr), size=budget, replace=False)
        arm["random"] = _fit_eval(X_tr[idx_rand], Y_tr[idx_rand], X_te, Y_te, donor_te)

        # cell-type stratified quota: rare types get guaranteed coverage
        # (biological motivation: protein profiles are type-specific and
        # random oversamples dominant types)
        quotas = {}
        uq, cnt = np.unique(ct_tr, return_counts=True)
        base = budget * cnt / cnt.sum()
        alloc = np.floor(base).astype(int)
        for j in range(len(uq)):
            if base[j] - alloc[j] > 0:
                alloc[j] += 1
        idx_strat = np.concatenate([
            rng.choice(np.where(ct_tr == t)[0], size=min(int(a), (ct_tr == t).sum()), replace=False)
            for t, a in zip(uq, alloc)
        ])
        arm["celltype_stratified"] = _fit_eval(X_tr[idx_strat], Y_tr[idx_strat], X_te, Y_te, donor_te)

        # committee disagreement (cellmean vs ridge pilot, cell-level)
        pilot = rng.choice(len(X_tr), size=min(500, budget), replace=False)
        pm = RidgeBaseline(alphas=(100.0, 1000.0)).fit(X_tr[pilot], Y_tr[pilot], groups=np.arange(len(pilot)) % 3, cv_folds=3)
        from ..adapters.baselines.cell_type_mean import CellTypeMeanModel
        cm = CellTypeMeanModel().fit(X_tr[pilot], Y_tr[pilot], ct_tr[pilot])
        disagree = np.abs(pm.predict(X_tr) - cm.predict(X_tr, ct=ct_tr)).mean(axis=1)
        idx_qbc = np.argsort(-disagree)[: 2 * budget]
        idx_qbc = rng.choice(idx_qbc, size=budget, replace=False)
        arm["qbc_committee"] = _fit_eval(X_tr[idx_qbc], Y_tr[idx_qbc], X_te, Y_te, donor_te)

        if extra_strategies:
            # diversity: PCA-50 then k-means on training features, nearest cells
            # (raw 2168-dim k-means costs ~30+ min at 36 cores; PCA first)
            from sklearn.decomposition import PCA

            Z = PCA(n_components=50, random_state=seed).fit_transform(X_tr)
            km = KMeans(n_clusters=min(budget, 256), random_state=seed, n_init=3).fit(Z[:: max(1, len(Z) // 20000)])
            labels = km.predict(Z)
        idx_div = []
        for c in np.unique(labels):
            members = np.where(labels == c)[0]
            take = max(1, int(round(budget * len(members) / len(X_tr))))
            idx_div.append(rng.choice(members, size=min(take, len(members)), replace=False))
        idx_div = np.concatenate(idx_div)[:budget]
        arm["diversity"] = _fit_eval(X_tr[idx_div], Y_tr[idx_div], X_te, Y_te, donor_te)

        # k-center greedy coreset on PCA-50 (farthest point sampling)
        Zs = Z[rng.choice(len(Z), size=min(20000, len(Z)), replace=False)]
        chosen = [int(rng.integers(len(Zs)))]
        d = np.linalg.norm(Zs - Zs[chosen[0]], axis=1)
        for _ in range(min(budget, 1000) - 1):
            nxt = int(np.argmax(d))
            chosen.append(nxt)
            d = np.minimum(d, np.linalg.norm(Zs - Zs[nxt], axis=1))
        arm["kcenter_corest"] = _fit_eval(X_tr[chosen], Y_tr[chosen], X_te, Y_te, donor_te)

        # residual-margin: pilot on small random set, then pick high-residual cells
        pilot = rng.choice(len(X_tr), size=min(200, budget), replace=False)
        pilot_model = RidgeBaseline(alphas=(100.0, 1000.0)).fit(X_tr[pilot], Y_tr[pilot], groups=donor_tr[pilot], cv_folds=3)
        resid = np.abs(pilot_model.predict(X_tr) - Y_tr).mean(axis=1)
        top = np.argsort(-resid)[: 2 * budget]
        idx_margin = rng.choice(top, size=budget, replace=False)
        arm["residual_margin"] = _fit_eval(X_tr[idx_margin], Y_tr[idx_margin], X_te, Y_te, donor_te)

        results[str(budget)] = arm
    return results


def _fit_eval(Xs, Ys, X_te, Y_te, donor_te):
    # Selected cells carry no donor grouping (the selection is the treatment);
    # fixed-alpha ridge — alpha chosen a priori, no CV on tiny budgets.
    from sklearn.linear_model import Ridge

    model = Ridge(alpha=1000.0).fit(Xs, Ys)
    pred = model.predict(X_te)
    rhos = per_donor_rho(Y_te, pred, donor_te)
    return {d: v["mean_rho"] for d, v in rhos.items()}
