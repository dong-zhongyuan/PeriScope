"""Transparent result presentation; no hidden filtering or evidence labels."""
import numpy as np
import pandas as pd
from scipy.stats import binomtest


def rank_effects(df, effect_col="beta", *, scale_cols=("effect_unit",)):
    """Rank all rows by absolute effect within stated comparable scales.

    Required scale columns prevent accidentally comparing arbitrary protein
    units. Pass scale_cols=() only when the caller has verified comparability.
    Missing effects remain at the end of each group. No rows are removed.
    """
    if effect_col not in df or any(c not in df for c in scale_cols):
        raise ValueError("effect or scale column missing")
    out = df.copy()
    effect = pd.to_numeric(out[effect_col], errors="raise")
    if np.isinf(effect.to_numpy(dtype=float)).any():
        raise ValueError("infinite effect")
    out["absolute_effect"] = effect.abs()
    if scale_cols:
        out["effect_rank"] = out.groupby(list(scale_cols), dropna=False)["absolute_effect"].rank(method="min", ascending=False)
    else:
        out["effect_rank"] = out.absolute_effect.rank(method="min", ascending=False)
    return out.sort_values([*scale_cols, "absolute_effect"], ascending=[True] * len(scale_cols) + [False], kind="stable", na_position="last")


def direction_summary(df, group_cols=("protein",), effect_col="beta", *, independent=False):
    """Count directions descriptively; correlated MRI outcomes are not coin flips.

    Exact two-sided binomial P is supplied ONLY when independent=True explicitly
    asserts independent signs under the null. Zeros are excluded from that test;
    missing effects, zero effects, and family size remain separately visible.
    The independence flag does not correct bias from selecting a family after
    seeing directions; that design issue must be recorded outside this function.
    """
    if independent is not True and independent is not False:
        raise ValueError("independent must be an explicit boolean")
    if effect_col not in df or any(c not in df for c in group_cols):
        raise ValueError("effect or group column missing")
    groups = df.groupby(list(group_cols), dropna=False, sort=False) if group_cols else [((), df)]
    rows = []
    for key, block in groups:
        key = key if isinstance(key, tuple) else (key,)
        effects = pd.to_numeric(block[effect_col], errors="raise").to_numpy(dtype=float)
        if np.isinf(effects).any():
            raise ValueError("infinite effect")
        pos, neg = int(np.sum(effects > 0)), int(np.sum(effects < 0))
        n = pos + neg
        rows.append({**dict(zip(group_cols, key)), "n_rows": len(block),
                     "n_positive": pos, "n_negative": neg,
                     "n_zero": int(np.sum(effects == 0)),
                     "n_missing": int(np.isnan(effects).sum()),
                     "majority_fraction_nonzero": max(pos, neg) / n if n else np.nan,
                     "sign_p": float(binomtest(pos, n, 0.5).pvalue) if independent and n else np.nan,
                     "sign_test_assumption": "independent_signs_asserted" if independent else "descriptive_only_correlations_unmodelled"})
    return pd.DataFrame(rows)
