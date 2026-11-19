"""Parallel numerical adjustments; thresholds/claims belong to the authors.

Missing p values retain their slots in the declared family (temporarily p=1),
then return NA. Invalid p values fail loudly. BH assumes independence/PRDS;
BY tolerates arbitrary dependence given valid marginal p values. Family BH
controls within each declared family, not globally. Fixed-lambda Storey values
are exploratory sensitivities with an estimated pi0, not guaranteed improved
power or dependence-robust inference. Two-stage BH is NOT hierarchical FDR.
"""
import numpy as np
import pandas as pd


def _prepare(pvalues):
    p = np.asarray(pvalues, dtype=float)
    if p.ndim != 1:
        raise ValueError("pvalues must be one-dimensional")
    if np.any(np.isinf(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("pvalues must be in [0,1] or NaN")
    return np.where(np.isnan(p), 1.0, p), np.isnan(p)


def adjust_pvalues(pvalues, method="bh"):
    p, missing = _prepare(pvalues)
    if method not in {"bh", "by"}:
        raise ValueError("method must be bh or by")
    if not len(p):
        return p
    order = np.argsort(p, kind="stable")
    ranks = np.arange(1, len(p) + 1)
    factor = len(p) * (np.sum(1 / ranks) if method == "by" else 1)
    adjusted = np.minimum.accumulate((p[order] * factor / ranks)[::-1])[::-1]
    result = np.empty_like(p)
    result[order] = np.minimum(1, adjusted)
    result[missing] = np.nan
    return result


def storey_qvalues(pvalues, lam=0.5):
    """Conservative +1 fixed-lambda pi0 estimate; returns (qvalues, pi0).

    pi0=min(1, (1+#(p>lambda))/(m*(1-lambda))). No spline/bootstrap tuning.
    Different lambdas must be reported together, not selected by best result.
    """
    if not np.isfinite(lam) or not 0 < lam < 1:
        raise ValueError("lambda must be between 0 and 1")
    p, missing = _prepare(pvalues)
    if not len(p):
        return p, np.nan
    pi0 = min(1.0, (1 + np.sum(p > lam)) / (len(p) * (1 - lam)))
    # Multiply BEFORE clipping: multiplying clipped BH would underadjust p=1.
    order = np.argsort(p, kind="stable")
    raw = pi0 * len(p) * p[order] / np.arange(1, len(p) + 1)
    result = np.empty_like(p)
    result[order] = np.minimum(1, np.minimum.accumulate(raw[::-1])[::-1])
    result[missing] = np.nan
    return result, float(pi0)


def two_stage_bh(pvalues, *, alpha):
    """Alpha-specific two-stage BH sensitivity (not a hierarchical method).

    https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.fdrcorrection_twostage.html
    The returned values cannot be reinterpreted using a different alpha.
    """
    from statsmodels.stats.multitest import fdrcorrection_twostage
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")
    p, missing = _prepare(pvalues)
    if not len(p):
        return {"adjusted": p, "alpha": alpha, "m0": 0, "alpha_stages": []}
    _, adjusted, m0, stages = fdrcorrection_twostage(p, alpha=alpha, method="bh", maxiter=1)
    adjusted[missing] = np.nan
    return {"adjusted": adjusted, "alpha": float(alpha), "m0": int(m0),
            "alpha_stages": [float(a) for a in stages]}


def adjust_table(df, p_col="p", family_cols=("family",), *, storey_lambdas=(0.25, 0.5, 0.75),
                 two_stage_alpha=None):
    """Return all original rows with adjustments; no filtering or significance flag.

    Family columns must define the chosen hypothesis groups explicitly. Pass ()
    to omit group adjustment. Row count and order are preserved, including NA p.
    """
    if p_col not in df:
        raise ValueError(f"missing p-value column: {p_col}")
    out = df.copy()
    p = pd.to_numeric(out[p_col], errors="raise").to_numpy(dtype=float)
    out["q_bh_global"] = adjust_pvalues(p, "bh")
    out["q_by_global"] = adjust_pvalues(p, "by")
    out["family_n_global"] = len(out)
    out["family_n_observed_global"] = np.isfinite(p).sum()
    if family_cols:
        family_cols = list(family_cols)
        if any(c not in out for c in family_cols):
            raise ValueError(f"missing family columns: {family_cols}")
        out["q_bh_family"] = np.nan
        out["family_n"] = 0
        out["family_n_observed"] = 0
        # Group indices are positional, preserving arbitrary/duplicated df indexes.
        groups = out.groupby(family_cols, dropna=False, sort=False).indices
        q, n, observed = np.full(len(out), np.nan), np.zeros(len(out), int), np.zeros(len(out), int)
        for positions in groups.values():
            q[positions] = adjust_pvalues(p[positions])
            n[positions] = len(positions)
            observed[positions] = np.isfinite(p[positions]).sum()
        out["q_bh_family"], out["family_n"], out["family_n_observed"] = q, n, observed
    for lam in storey_lambdas:
        label = format(lam, "g").replace(".", "p")
        q, pi0 = storey_qvalues(p, lam)
        out[f"q_storey_lambda_{label}"] = q
        out[f"pi0_lambda_{label}"] = pi0
    if two_stage_alpha is not None:
        res = two_stage_bh(p, alpha=two_stage_alpha)
        out["q_two_stage_bh"] = res["adjusted"]
        out["two_stage_alpha"] = res["alpha"]
        out["two_stage_m0"] = res["m0"]
    return out
