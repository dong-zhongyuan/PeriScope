"""Synthetic null calibration, distinct from biological validation.

Merely shuffling labels among real nonzero effects preserves their p-value
distribution and is NOT a valid demonstration of a uniform null. This module
generates zero-mean outcome effects with known standard errors instead.
"""
import numpy as np
from scipy.stats import kstest
from .mr import wald


def null_calibration(seed=20261115, n=10000):
    """Return diagnostics under independent Gaussian zero-effect outcomes.

    Exposure effects/SEs and outcome SEs vary, but outcome beta is freshly drawn
    N(0, SE_y^2). Wald first-order z is sign(beta_x)*beta_y/SE_y. No real data
    are used; this checks numerical calibration, not instrument validity or LD.
    KS P itself is random and must not be used as a biological decision gate.
    """
    if not isinstance(n, (int, np.integer)) or n < 10:
        raise ValueError("n must be an integer >=10")
    rng = np.random.default_rng(seed)
    bx = rng.choice([-1, 1], n) * rng.uniform(0.1, 0.8, n)
    sx = rng.uniform(0.005, 0.015, n)
    sy = rng.uniform(0.005, 0.05, n)
    by = rng.normal(0, sy)
    # Exercise the actual estimator, rather than only testing an algebraic copy.
    p = np.asarray([wald(x, ex, y, ey)["p"] for x, ex, y, ey in zip(bx, sx, by, sy)])
    ks = kstest(p, "uniform")
    return {"kind": "synthetic_zero_effect_calibration", "seed": int(seed), "n": int(n),
            "real_data_used": False, "mean_p": float(p.mean()),
            "fraction_p_below_0p05": float(np.mean(p < 0.05)),
            "ks_statistic": float(ks.statistic), "ks_p": float(ks.pvalue),
            "p_quantiles": {str(q): float(np.quantile(p, q)) for q in (0.01, 0.05, 0.5, 0.95, 0.99)},
            "interpretation": "Numerical self-check only; no biological inference or real-data null validation."}
