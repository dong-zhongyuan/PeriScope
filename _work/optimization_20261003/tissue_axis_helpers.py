"""Copied bootstrap/rank helpers; only numerical-identifiability and failure accounting repaired."""
import numpy as np
import pandas as pd
N_BOOT=800

def welch_t(M, c):
    """per-gene donor-level Welch t; M (donors x genes), c (0/1)."""
    m1, m0 = M[c == 1], M[c == 0]
    v1 = m1.var(0, ddof=1) / max(len(m1), 2)
    v0 = m0.var(0, ddof=1) / max(len(m0), 2)
    se = np.sqrt(v1 + v0) + 1e-9
    return (m1.mean(0) - m0.mean(0)) / se


def rank_(x):
    return pd.Series(x).rank().values


def wspearman(a, b, w=None):
    if w is None:
        from scipy.stats import spearmanr
        return float(spearmanr(a, b).statistic)
    ra, rb = rank_(a), rank_(b)
    ra = (ra - ra.mean()) / (ra.std() + 1e-12)
    rb = (rb - rb.mean()) / (rb.std() + 1e-12)
    w = w / w.sum()
    return float((w * ra * rb).sum() / np.sqrt((w * ra * ra).sum() * (w * rb * rb).sum()))


def ols_correct(M, c, C):
    if C is None or not len(C):
        return M - np.zeros(M.shape), False
    if not np.isfinite(C).all():raise ValueError('missing donor covariates')
    C=C[:,np.std(C,axis=0)>1e-12]
    D = np.column_stack([np.ones(len(C)), c.astype(float), C])
    if np.linalg.matrix_rank(D)<D.shape[1]:raise ValueError('nonidentified adjusted contrast')
    beta, *_ = np.linalg.lstsq(D, M, rcond=None)
    return M - D @ beta + D[:, :2] @ beta[:2], True   # keep intercept+cond, drop covariate part


def meta_blocks(blocks):
    """fixed-effect meta across cohort blocks: weighted diff + t by sqrt(n1*n0)."""
    ws = np.array([np.sqrt((b["conds"] == 1).sum() * (b["conds"] == 0).sum()) for b in blocks])
    ws = ws / ws.sum()
    diff = sum(w * b["diff"] for w, b in zip(ws, blocks))
    t = sum(w * b["t"] for w, b in zip(ws, blocks))
    return diff, t


def boot_pair(blk_a_list, blk_b_list, cov_a=None, cov_b_list=None, n_boot=N_BOOT, seed=0):
    """donor bootstrap through the full pipeline (resample -> OLS -> t -> meta -> weighted r)."""
    rng = np.random.RandomState(seed)

    def one_pass(pick_a=None, pick_b=None):
        pick_a = pick_a if pick_a is not None else [None] * len(blk_a_list)
        pick_b = pick_b if pick_b is not None else [None] * len(blk_b_list)
        ab = []
        for blk, pick in zip(blk_a_list, pick_a):
            M, c = blk["means"], blk["conds"]
            if pick is None:
                idx = np.arange(len(c))
            else:
                idx = pick
            Mc, cc = M[idx], c[idx]
            Cc = cov_a[idx] if cov_a is not None else None
            Mc2, _ = ols_correct(Mc, cc, Cc)
            ab.append({"diff": cc[cc == 1].mean() - 0 if False else Mc2[cc == 1].mean(0) - Mc2[cc == 0].mean(0),
                       "t": welch_t(Mc2, cc), "conds": cc})
        bb = []
        for blk, pick, Cv in zip(blk_b_list, pick_b, cov_b_list or [None] * len(blk_b_list)):
            M, c = blk["means"], blk["conds"]
            idx = np.arange(len(c)) if pick is None else pick
            Mc, cc = M[idx], c[idx]
            Cc = Cv[idx] if Cv is not None else None
            Mc2, _ = ols_correct(Mc, cc, Cc)
            bb.append({"diff": McCdiff(Mc2, cc), "t": welch_t(Mc2, cc), "conds": cc})
        da, ta = meta_blocks(ab) if len(ab) > 1 else (ab[0]["diff"], ab[0]["t"])
        db, tb = meta_blocks(bb) if len(bb) > 1 else (bb[0]["diff"], bb[0]["t"])
        wgt = np.abs(ta) * np.abs(tb)
        return wspearman(ta, tb, wgt), wspearman(ta, tb)

    def pick_arm(blk):
        c = blk["conds"]
        return np.concatenate([rng.choice(np.where(c == 1)[0], (c == 1).sum(), True),
                               rng.choice(np.where(c == 0)[0], (c == 0).sum(), True)])

    r_w, r_u = one_pass(None, None)
    bw, bu = [], []
    for _ in range(n_boot):
        pa = [pick_arm(b) for b in blk_a_list]
        pb = [pick_arm(b) for b in blk_b_list]
        try:
            a, b2 = one_pass(pa, pb)
            if np.isfinite(a) and np.isfinite(b2):bw.append(a); bu.append(b2)
        except Exception:
            continue
    bw, bu = np.array(bw), np.array(bu)
    if len(bw)<200:return dict(r_weighted=r_w,r_unweighted=r_u,n_boot=len(bw),n_boot_requested=n_boot,status='insufficient_identifiable_bootstraps')
    return {"r_weighted": round(r_w, 4), "r_unweighted": round(r_u, 4),
            "boot_p_pos_w": round(float((bw <= 0).mean()), 4),
            "ci_w": [round(float(np.percentile(bw, 2.5)), 4), round(float(np.percentile(bw, 97.5)), 4)],
            "boot_p_pos_u": round(float((bu <= 0).mean()), 4),
            "ci_u": [round(float(np.percentile(bu, 2.5)), 4), round(float(np.percentile(bu, 97.5)), 4)],
            "n_boot": len(bw),"n_boot_requested":n_boot,"status":"ok"}


def McCdiff(M, c):
    return M[c == 1].mean(0) - M[c == 0].mean(0)


