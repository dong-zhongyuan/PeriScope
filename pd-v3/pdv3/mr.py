"""Summary-data MR arithmetic; no biological support decisions.

Wald delta errors assume independent exposure and outcome estimates. IVW uses
first-order weights, assumes independent instruments, and does not model sample
overlap, horizontal pleiotropy, selection bias, or LD. State those assumptions
with results; passing ``independent=True`` is an assertion, not an LD check.
"""
import numpy as np
from scipy.stats import chi2, norm

_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _build(value):
    if value is None:
        return None
    v = str(value).strip().lower()
    aliases = {"grch37": "37", "hg19": "37", "b37": "37", "37": "37",
               "grch38": "38", "hg38": "38", "b38": "38", "38": "38"}
    return aliases.get(v)


def harmonize(exposure_ea, exposure_oa, outcome_ea, outcome_oa, beta_y,
              eaf_x=None, eaf_y=None, *, exposure_build=None, outcome_build=None,
              palindrome="drop", eaf_tolerance=0.1, maf_ceiling=0.42):
    """Align SNP outcome beta/EAF to the exposure effect allele.

    The caller must already have verified variant identity. This function never
    liftover-maps coordinates or checks rsID identity. Missing build metadata is
    retained as ``build_status='unverified'``; explicit different builds fail.
    Only biallelic A/C/G/T SNPs are supported. Palindromes are dropped by default.
    With ``palindrome='infer'``, both EAFs must have MAF <= maf_ceiling and exactly
    one orientation must agree within eaf_tolerance. Population EAF differences
    can defeat this inference. All failures are returned as explicit rows.
    """
    if palindrome not in {"drop", "infer"}:
        raise ValueError("palindrome must be 'drop' or 'infer'")
    if not 0 <= eaf_tolerance < 0.5 or not 0 < maf_ceiling < 0.5:
        raise ValueError("invalid frequency inference tolerances")
    out = {"usable": False, "status": None, "beta_y": np.nan,
           "eaf_y": np.nan, "flipped": None, "build_status": "unverified",
           "eaf_difference": np.nan}

    def fail(reason):
        return {**out, "status": reason}

    bxbuild, bybuild = _build(exposure_build), _build(outcome_build)
    if ((exposure_build is not None and bxbuild is None) or
            (outcome_build is not None and bybuild is None)):
        return fail("unknown_build")
    if bxbuild is not None and bybuild is not None:
        if bxbuild != bybuild:
            return fail("build_mismatch")
        out["build_status"] = "matched"
    alleles = [str(a).strip().upper() for a in
               (exposure_ea, exposure_oa, outcome_ea, outcome_oa)]
    if any(a not in {"A", "C", "G", "T"} for a in alleles):
        return fail("unsupported_allele")
    ea, oa, ya, yo = alleles
    if ea == oa or ya == yo:
        return fail("identical_alleles")
    try:
        beta_y = float(beta_y)
    except (TypeError, ValueError):
        return fail("invalid_beta")
    if not np.isfinite(beta_y):
        return fail("invalid_beta")
    frequencies = []
    for f in (eaf_x, eaf_y):
        if f is None:
            frequencies.append(None)
            continue
        try:
            f = float(f)
        except (TypeError, ValueError):
            return fail("invalid_eaf")
        if np.isnan(f):
            frequencies.append(None)
        elif not np.isfinite(f) or not 0 <= f <= 1:
            return fail("invalid_eaf")
        else:
            frequencies.append(f)
    fx, fy = frequencies
    candidates = []
    for a, b, sign, label in ((ya, yo, 1, "direct"), (yo, ya, -1, "swap"),
                              (ya.translate(_COMPLEMENT), yo.translate(_COMPLEMENT), 1, "strand"),
                              (yo.translate(_COMPLEMENT), ya.translate(_COMPLEMENT), -1, "strand_swap")):
        if (ea, oa) == (a, b):
            candidates.append((sign, label))
    if not candidates:
        return fail("allele_mismatch")
    is_palindromic = {ea, oa} in ({"A", "T"}, {"C", "G"})
    if is_palindromic:
        if palindrome == "drop":
            return fail("palindromic_dropped")
        if fx is None or fy is None:
            return fail("palindromic_missing_eaf")
        if max(min(fx, 1 - fx), min(fy, 1 - fy)) > maf_ceiling:
            return fail("palindromic_ambiguous_eaf")
        candidates = [(s, label) for s, label in candidates
                      if abs(fx - (fy if s == 1 else 1 - fy)) <= eaf_tolerance]
        if len(candidates) != 1:
            return fail("palindromic_frequency_mismatch")
    sign, label = candidates[0]
    aligned_eaf = (fy if sign == 1 else 1 - fy) if fy is not None else np.nan
    return {**out, "usable": True,
            "status": "palindromic_" + label if is_palindromic else label,
            "beta_y": sign * beta_y, "eaf_y": aligned_eaf,
            "flipped": sign == -1,
            "eaf_difference": abs(fx - aligned_eaf) if fx is not None else np.nan}


def _inputs(bx, sx, by, sy):
    values = [np.asarray(x, dtype=float) for x in (bx, sx, by, sy)]
    if any(x.shape != values[0].shape for x in values):
        raise ValueError("effect and SE arrays must have identical shapes")
    if any(not np.all(np.isfinite(x)) for x in values):
        raise ValueError("effect and SE values must be finite")
    if np.any(values[0] == 0) or np.any(values[1] <= 0) or np.any(values[3] <= 0):
        raise ValueError("exposure effects must be nonzero and SEs positive")
    return values


def _normal_summary(beta, se):
    z = beta / se
    return {"se": float(se), "p": float(2 * norm.sf(abs(z))),
            "log_p": float(np.log(2) + norm.logsf(abs(z))),
            "ci_lower": float(beta - norm.ppf(0.975) * se),
            "ci_upper": float(beta + norm.ppf(0.975) * se)}


def wald(bx, sx, by, sy):
    """Scalar Wald ratio with first-order and independent-sample delta SEs."""
    bx, sx, by, sy = _inputs(bx, sx, by, sy)
    if bx.ndim != 0:
        raise ValueError("wald expects scalar values")
    beta = by / bx
    se = abs(sy / bx)
    se_delta = np.hypot(se, by * sx / bx**2)
    first = _normal_summary(beta, se)
    delta = {key + "_delta": value for key, value in _normal_summary(beta, se_delta).items()}
    return {"method": "Wald", "beta": float(beta), **first, **delta,
            "F": float((bx / sx)**2), "n_instruments": 1}


def ivw(bx, sx, by, sy, *, independent=False):
    """IVW through-origin, fixed and multiplicative random errors (no shrinkage).

    Requires >=2 independent instruments. Correlated cis variants must be
    pruned using a suitable LD reference before calling, or analysed elsewhere
    with an LD-aware estimator. Marginal F values do not establish independence.
    """
    if independent is not True:
        raise ValueError("independent=True requires documented LD pruning/independence")
    bx, sx, by, sy = _inputs(bx, sx, by, sy)
    if bx.ndim != 1 or len(bx) < 2:
        raise ValueError("ivw needs one-dimensional arrays with >=2 instruments")
    w = 1 / sy**2
    denominator = np.sum(w * bx**2)
    beta = np.sum(w * bx * by) / denominator
    se = np.sqrt(1 / denominator)
    q = np.sum(w * (by - beta * bx)**2)
    df = len(bx) - 1
    dispersion = max(1.0, float(q / df))
    random = {key + "_random": value
              for key, value in _normal_summary(beta, se * np.sqrt(dispersion)).items()}
    return {"method": "IVW_independent", "beta": float(beta),
            **_normal_summary(beta, se), **random,
            "Q": float(q), "Q_p": float(chi2.sf(q, df)), "df": df,
            "dispersion": dispersion, "F_min": float(np.min((bx / sx)**2)),
            "n_instruments": len(bx)}


def variance_explained(beta, eaf, *, trait_variance=None, standardized=False):
    """2*f*(1-f)*beta^2 / trait_variance for one additive continuous trait SNP.

    This is not liability-scale R2 for binary outcomes. A trait variance must be
    supplied, or standardized=True must assert the beta's outcome has variance 1.
    Do not sum these quantities across correlated instruments.
    """
    if trait_variance is None:
        if standardized is not True:
            raise ValueError("provide trait_variance or assert standardized=True")
        trait_variance = 1.0
    elif standardized and trait_variance != 1:
        raise ValueError("standardized trait must have variance 1")
    if not np.all(np.isfinite([beta, eaf, trait_variance])) or not 0 <= eaf <= 1 or trait_variance <= 0:
        raise ValueError("invalid beta, EAF, or trait variance")
    r2 = 2 * eaf * (1 - eaf) * beta**2 / trait_variance
    if r2 > 1:
        raise ValueError("R2 > 1; verify effect scale and trait variance")
    return float(r2)
