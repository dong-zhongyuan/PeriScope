"""Conditional detectability and assumed-covariance intervals, not overlap correction.

Exposure effects are treated as fixed for detectability. Fieller intervals assume
unselected jointly normal association estimates with a supplied covariance.
Neither calculation establishes valid instruments or an independent sample.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm


def _alpha(alpha):
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError('alpha must be in (0,1)')
    return norm.isf(alpha/2)


def conditional_power(effect, se, *, alpha=.05):
    """Two-sided normal-test power under fixed exposure effect and outcome SE."""
    effect, se = np.broadcast_arrays(np.asarray(effect, float), np.asarray(se, float))
    if not np.isfinite(effect).all() or not np.isfinite(se).all() or (se <= 0).any():
        raise ValueError('finite effect and positive finite SE required')
    z = _alpha(alpha)
    noncentrality = np.abs(effect)/se
    return norm.sf(z-noncentrality) + norm.cdf(-z-noncentrality)


def detectable_effect(se, *, target_power=.8, alpha=.05):
    """Minimum |effect| for the specified conditional normal-test power."""
    _alpha(alpha)
    if not np.isfinite(target_power) or not alpha < target_power < 1:
        raise ValueError('target power must lie strictly between alpha and 1')
    se = np.asarray(se, float)
    if not np.isfinite(se).all() or (se <= 0).any():
        raise ValueError('positive finite SE required')
    upper = _alpha(alpha)+norm.ppf(target_power)+3
    while conditional_power(upper, 1., alpha=alpha) < target_power:
        upper *= 2
    factor = brentq(lambda d: float(conditional_power(d, 1., alpha=alpha))-target_power,
                    0., upper, xtol=1e-13)
    return se*factor


def _associations(bx, sx, by, sy, rho):
    bx, sx, by, sy, rho = np.broadcast_arrays(*[np.asarray(x, float) for x in (bx, sx, by, sy, rho)])
    if not all(np.isfinite(x).all() for x in (bx, sx, by, sy, rho)):
        raise ValueError('association inputs must be finite')
    if (bx == 0).any() or (sx <= 0).any() or (sy <= 0).any() or (np.abs(rho) > 1).any():
        raise ValueError('nonzero bx, positive SEs and error correlation in [-1,1] required')
    return bx, sx, by, sy, rho


def covariance_delta_se(bx, sx, by, sy, *, rho=0.):
    """Delta SE at the ratio estimate; rho is an assumed sampling-error correlation."""
    bx, sx, by, sy, rho = _associations(bx, sx, by, sy, rho)
    theta = by/bx
    # Sum-of-squares form avoids cancellation at singular covariance endpoints.
    return np.hypot(sy-theta*rho*sx, theta*sx*np.sqrt(1-rho*rho))/np.abs(bx)


def bounded_fieller_interval(bx, sx, by, sy, *, rho=0., alpha=.05):
    """Finite Fieller interval only; refuse to misrepresent unbounded cases.

    Inverts (by-theta*bx)^2 <= z^2 Var(by-theta*bx). Inputs are on
    harmonized allele scales. It does not estimate rho or remove selection bias.
    """
    bx, sx, by, sy, rho = _associations(bx, sx, by, sy, rho)
    z = _alpha(alpha)
    a = bx*bx-(z*sx)**2
    if (a <= 0).any():
        raise ValueError('Fieller set is not guaranteed finite: bx^2 <= z^2 sx^2')
    centre = (bx*by-z*z*rho*sx*sy)/a
    radius = z*np.sqrt((by*sx-bx*sy*rho)**2+(1-rho*rho)*sy*sy*a)/a
    return centre-radius, centre+radius
