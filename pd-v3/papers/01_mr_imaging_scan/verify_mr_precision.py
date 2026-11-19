"""Independent numerical checks of the 45-trait precision implementation audit."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import chi2, ncx2, norm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.mr_precision import bounded_fieller_interval

OUT = ROOT/'papers/01_mr_imaging_scan/results/precision_overlap_20261116/validation45'
SOURCE = ROOT/'papers/01_mr_imaging_scan/results/expanded_20261116/combined_45'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    audit = json.loads((OUT/'AUDIT.json').read_text())
    manifest = json.loads((OUT/'MANIFEST.json').read_text())
    for mapping in [audit['input_sha256'], audit['code_sha256'], manifest['artifacts_sha256']]:
        for name, expected in mapping.items():
            assert sha(ROOT/name) == expected, name
    checks = []
    for policy in ['drop', 'infer']:
        index = ['phenotype', 'exposure_id']
        original = pd.read_csv(SOURCE/f'full_family_palindrome_{policy}.csv.gz',
                              usecols=index+['status', 'beta_x', 'se_x', 'beta_y_aligned', 'se_y']).set_index(index)
        result = pd.read_csv(OUT/f'pair_precision_{policy}.csv.gz').set_index(index)
        assert original.index.is_unique and result.index.is_unique
        assert set(original.index) == set(result.index) and len(result) == 132345
        original = original.loc[result.index]
        assert np.array_equal(original.status, result.status)
        ok = result.status.eq('estimated')
        assert result.loc[~ok, result.columns.difference(['family', 'status'])].isna().all().all()
        bx, sx, by, sy = [original.loc[ok, c].to_numpy() for c in ['beta_x', 'se_x', 'beta_y_aligned', 'se_y']]
        s = sy/abs(bx)
        power_error = 0.
        for name, alpha in audit['alpha_references'].items():
            for target in [.8, .9]:
                mde = result.loc[ok, f'MDE{int(target*100)}_{name}'].to_numpy()
                p = ncx2.sf(chi2.isf(alpha, 1), 1, (mde/s)**2)
                power_error = max(power_error, float(np.max(np.abs(p-target))))
        assert power_error < 1e-12
        a = bx*bx-norm.isf(.025)**2*sx*sx
        all_lo, all_hi, all_delta = [], [], []
        zero_checks = 0
        for rho in audit['rho_error_grid']:
            q = norm.isf(.025)**2
            b = -2*bx*by+2*q*rho*sx*sy
            c = by*by-q*sy*sy
            discriminant = b*b-4*a*c
            # Independent direct quadratic formula, unlike the stable sum of
            # squares used by the implementation being verified.
            assert np.min(discriminant) >= 0
            lo, hi = (-b-np.sqrt(discriminant))/(2*a), (-b+np.sqrt(discriminant))/(2*a)
            assert np.array_equal((lo <= 0) & (hi >= 0), abs(by/sy) <= norm.isf(.025))
            zero_checks += len(lo)
            all_lo.append(lo); all_hi.append(hi)
            variance = sy*sy/(bx*bx) + by*by*sx*sx/bx**4 - 2*by*rho*sx*sy/bx**3
            assert variance.min() >= -1e-15
            all_delta.append(np.sqrt(np.maximum(variance, 0)))
        expected = dict(Fieller95_rho0_lower=all_lo[4], Fieller95_rho0_upper=all_hi[4],
                        Fieller95_grid_lower_min=np.min(all_lo, axis=0),
                        Fieller95_grid_upper_max=np.max(all_hi, axis=0),
                        delta_SE_grid_min=np.min(all_delta, axis=0), delta_SE_grid_max=np.max(all_delta, axis=0))
        errors = {}
        for name, values in expected.items():
            actual = result.loc[ok, name].to_numpy()
            errors[name] = float(np.max(abs(values-actual)))
            np.testing.assert_allclose(actual, values, rtol=1e-9, atol=1e-12)
        checks.append(dict(policy=policy, planned_rows=len(result), estimated=int(ok.sum()),
                           power_target_max_abs_error=power_error,
                           quadratic_and_covariance_max_abs_errors=errors, null_checks=zero_checks))
    # Known-parameter Gaussian coverage only: no selection, pleiotropy or actual
    # UKB sampling. This must not be presented as correcting the real study.
    rng = np.random.default_rng(20261116)
    rows = []
    n = 200000
    for theta in [-.25, 0., .25]:
        for rho in [-1., -.75, -.5, -.25, 0., .25, .5, .75, 1.]:
            u, v = rng.standard_normal((2, n))
            bx = .25+.01*u
            by = theta*.25+.02*(rho*u+np.sqrt(1-rho*rho)*v)
            lo, hi = bounded_fieller_interval(bx, .01, by, .02, rho=rho)
            coverage = float(np.mean((lo <= theta) & (hi >= theta)))
            assert abs(coverage-.95) < 6*np.sqrt(.95*.05/n)
            rows.append(dict(true_effect=theta, rho_error=rho, n=n, coverage=coverage,
                             selection_performed=False, true_bx=.25, sx=.01, sy=.02))
    pd.DataFrame(rows).to_csv(OUT/'gaussian_coverage_check.csv', index=False)
    result = dict(complete=True, code_sha256=sha(Path(__file__)), mode='validation45_only',
                  manifest_sha256=sha(OUT/'MANIFEST.json'), checks=checks,
                  simulated_scenarios=len(rows), samples_per_scenario=n,
                  simulated_coverage_range=[min(r['coverage'] for r in rows), max(r['coverage'] for r in rows)],
                  six_MC_standard_error_tolerance=True, actual_overlap_or_selection_validated=False)
    (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
