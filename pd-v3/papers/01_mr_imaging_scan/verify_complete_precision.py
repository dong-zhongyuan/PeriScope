"""Verify every precision result and summary without its production functions."""
import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import chi2, ncx2, norm

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'papers/01_mr_imaging_scan/results'
TOTAL = 11655183
ALPHAS = {'nominal': .05, 'full_plan_bonferroni_reference': .05/TOTAL}
RHOS = [-1., -.75, -.5, -.25, 0., .25, .5, .75, 1.]
KEYS = ['exposure_id', 'phenotype', 'family', 'status']
SOURCE_COLUMNS = KEYS + ['beta_x', 'se_x', 'beta_y_aligned', 'se_y', 'beta', 'se', 'se_delta', 'F']


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, name):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    require(a.shape == b.shape and np.array_equal(np.isnan(a), np.isnan(b)), f'Shape/missingness: {name}')
    np.testing.assert_allclose(a, b, rtol=1e-9, atol=1e-12, equal_nan=True, err_msg=name)
    valid = np.isfinite(a) & np.isfinite(b)
    return float(np.max(abs(a[valid]-b[valid]))) if valid.any() else 0.


def completed_audit(audit, mode):
    require(mode in ('full', 'validation45'), 'Unknown precision verification mode')
    ntraits = 3963 if mode == 'full' else 45
    require(audit.get('complete') is True and audit.get('mode') == mode
            and audit.get('original_full_plan_pairs') == TOTAL
            and audit.get('analysed_pairs_per_policy') == 2941*ntraits,
            'Precision result does not cover the completed declared plan')
    require(audit.get('alpha_references') == ALPHAS and audit.get('rho_error_grid') == RHOS,
            'Precision parameters differ from the frozen specification')
    require(set(audit.get('policies', {})) == {'drop', 'infer'}, 'Both precision policies are required')
    for field in ('actual_overlap_estimated', 'overlap_bias_corrected', 'new_MR_discovery'):
        require(audit.get(field) is False, f'Unsupported interpretation: {field}')
    return ntraits


def context(mode):
    folder = BASE/f'precision_overlap_20261116/{mode}'
    path = folder/'AUDIT.json'
    require(path.is_file(), 'Completed precision audit is unavailable')
    audit = json.loads(path.read_text())
    ntraits = completed_audit(audit, mode)
    source = BASE/('full_atlas_20261116/full' if mode == 'full' else 'expanded_20261116/combined_45')
    evidence = {str(path.relative_to(ROOT)): sha(path)}
    def check(path, digest):
        require(sha(path) == digest, f'Input SHA mismatch: {path}')
        evidence[str(path.relative_to(ROOT))] = digest
    if mode == 'full':
        run = json.loads((source/'RUN.json').read_text())
        require(run.get('complete') is True and run.get('mode') == 'full'
                and run.get('status') == 'complete' and run.get('planned_pairs') == TOTAL
                and run.get('assays') == 2941 and run.get('phenotypes') == 3963
                and run.get('completed_traits') == 3963, 'Original full MR is incomplete')
        check(source/'RUN.json', sha(source/'RUN.json'))
        for policy in ('drop', 'infer'):
            check(source/f'full_family_palindrome_{policy}.csv.gz', run['policies'][policy]['output_sha256'])
    manifest = json.loads((folder/'MANIFEST.json').read_text())
    check(folder/'MANIFEST.json', sha(folder/'MANIFEST.json'))
    for mapping in (audit['input_sha256'], audit['code_sha256'], manifest['artifacts_sha256']):
        for name, digest in mapping.items():
            check(ROOT/name, digest)
    for policy in ('drop', 'infer'):
        name = str((source/f'full_family_palindrome_{policy}.csv.gz').relative_to(ROOT))
        require(name in audit['input_sha256'], 'Precision result belongs to a different MR source')
        check(folder/f'pair_precision_{policy}.csv.gz', audit['policies'][policy]['output_sha256'])
    check(source/'phenotypes.csv', audit['metadata_sha256'])
    check(source/'assay_coverage.csv', audit['assay_coverage_sha256'])
    meta = pd.read_csv(source/'phenotypes.csv')
    cov = pd.read_csv(source/'assay_coverage.csv')
    require(len(meta) == ntraits and meta.phenotype.is_unique, 'Phenotype plan differs')
    require(len(cov) == 2941 and cov.exposure_id.is_unique, 'Assay plan differs')
    pd.testing.assert_frame_equal(pd.read_csv(folder/'phenotypes.csv'), meta,
                                  check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    return folder, source, audit, meta, cov, evidence


def source_blocks(path, mode, meta, cov):
    if mode == 'validation45':
        data = pd.read_csv(path, usecols=SOURCE_COLUMNS, float_precision='round_trip')
        require(len(data) == len(meta)*len(cov)
                and not data.duplicated(['phenotype', 'exposure_id']).any(), 'Invalid legacy Cartesian input')
        require(set(data.phenotype) == set(meta.phenotype), 'Legacy phenotype coverage differs')
        for phenotype, group in data.groupby('phenotype', sort=False):
            require(set(group.exposure_id) == set(cov.exposure_id), 'Legacy assay coverage differs')
            yield group.set_index('exposure_id').loc[cov.exposure_id].reset_index()
    else:
        with pd.read_csv(path, usecols=SOURCE_COLUMNS, chunksize=len(cov),
                         float_precision='round_trip') as reader:
            yield from reader


def paired_keys(original, result, trait, assays):
    require(len(original) == len(result) == len(assays), 'Incomplete source/precision block')
    for frame in (original, result):
        require(frame.exposure_id.tolist() == assays, 'Assay identity/order differs')
        require(frame.phenotype.eq(trait['phenotype']).all() and frame.family.eq(trait['family']).all(),
                'Phenotype metadata differs')
    require(original.status.tolist() == result.status.tolist(), 'Precision changed a source status')


def quadratic(bx, sx, by, sy, rho):
    z2 = norm.isf(.025)**2
    a = bx**2-z2*sx**2
    require((a > 0).all(), 'Fieller denominator is not bounded; finite interval forbidden')
    b = -2*bx*by+2*z2*rho*sx*sy
    c = by**2-z2*sy**2
    discriminant = b*b-4*a*c
    require((discriminant >= 0).all(), 'Invalid Fieller discriminant')
    lo = (-b-np.sqrt(discriminant))/(2*a)
    hi = (-b+np.sqrt(discriminant))/(2*a)
    residual = 0.
    for endpoint in (lo, hi):
        lhs = (by-endpoint*bx)**2
        rhs = z2*(sy**2+endpoint**2*sx**2-2*endpoint*rho*sx*sy)
        scale = lhs+z2*(sy**2+endpoint**2*sx**2+abs(2*endpoint*rho*sx*sy))
        error = abs(lhs-rhs)/np.maximum(scale, np.finfo(float).tiny)
        residual = max(residual, float(np.max(error, initial=0.)))
    require(residual < 1e-10, 'Fieller endpoints fail the original inequality boundary')
    variance = sy**2/bx**2+by**2*sx**2/bx**4-2*by*rho*sx*sy/bx**3
    scale = sy**2/bx**2+by**2*sx**2/bx**4+abs(2*by*rho*sx*sy/bx**3)
    require((variance >= -32*np.finfo(float).eps*scale).all(), 'Negative covariance variance')
    return lo, hi, np.sqrt(np.maximum(variance, 0.)), residual


def quantile_check(record, values, prefix):
    expected = np.percentile(values, [5, 50, 95]) if len(values) else [np.nan]*3
    return close([record[f'{prefix}_{suffix}'] for suffix in ('p05', 'median', 'p95')],
                 expected, prefix+' summary')


def verify_block(original, result, policy, mde, covariance):
    # paired_keys has already checked identical scientific row order. Reader
    # chunks can carry different pandas index labels, which are not identities.
    valid = result.status.eq('estimated').to_numpy()
    derived = result.columns.difference(KEYS)
    require(result.loc[~valid, derived].isna().all().all(), 'Unestimated precision slot was filled')
    require(np.isfinite(result.loc[valid, derived].to_numpy(float)).all(), 'Nonfinite estimated precision')
    bx, sx, by, sy = [original.loc[valid, col].to_numpy(float) for col in
                      ('beta_x', 'se_x', 'beta_y_aligned', 'se_y')]
    require(np.isfinite(np.array([bx, sx, by, sy])).all()
            and (bx != 0).all() and (sx > 0).all() and (sy > 0).all(), 'Invalid source estimate')
    actual = result.loc[valid]
    phenotype, family = original.phenotype.iloc[0], original.family.iloc[0]
    n, nplan = len(actual), len(result)
    se = sy/abs(bx)
    errors = {'conditional_se': close(actual.conditional_se, se, 'Conditional SE')}
    close(original.loc[valid, 'se'], se, 'Original conditional SE')
    power_error, summary_error, residual = 0., 0., 0.
    for label, alpha in ALPHAS.items():
        record = mde.loc[(policy, phenotype, label)]
        require(record['family'] == family and record['planned'] == nplan
                and record['estimated'] == n and record['alpha'] == alpha, 'MDE summary metadata differs')
        for target in (.8, .9):
            name = f'MDE{int(target*100)}'
            values = actual[f'{name}_{label}'].to_numpy(float)
            require((values >= 0).all(), 'Negative detectable effect')
            power = ncx2.sf(chi2.isf(alpha, 1), 1, (values/se)**2)
            power_error = max(power_error, float(np.max(abs(power-target), initial=0.)))
            summary_error = max(summary_error, quantile_check(record, values, name))
    require(power_error < 1e-12, 'Detectable effects fail independent power inversion')
    lo0, hi0, delta0, _ = quadratic(bx, sx, by, sy, 0.)
    lows, highs, deltas = [], [], []
    null_in = abs(by/sy) <= norm.isf(.025)
    for rho in RHOS:
        lo, hi, delta, error = quadratic(bx, sx, by, sy, rho)
        residual = max(residual, error)
        require(np.array_equal((lo <= 0) & (hi >= 0), null_in), 'Fieller null identity differs')
        record = covariance.loc[(policy, phenotype, rho)]
        require(record['planned'] == nplan and record['estimated'] == n
                and record['intervals_including_zero'] == int(null_in.sum()), 'Covariance summary count differs')
        summary_error = max(summary_error,
                            quantile_check(record, (hi-lo)/(hi0-lo0), 'Fieller_width_ratio_to_rho0'),
                            quantile_check(record, delta/delta0, 'delta_SE_ratio_to_rho0'))
        lows.append(lo); highs.append(hi); deltas.append(delta)
    expected = dict(Fieller95_rho0_lower=lo0, Fieller95_rho0_upper=hi0,
                    Fieller95_grid_lower_min=np.min(lows, axis=0), Fieller95_grid_upper_max=np.max(highs, axis=0),
                    delta_SE_grid_min=np.min(deltas, axis=0), delta_SE_grid_max=np.max(deltas, axis=0))
    for name, values in expected.items():
        errors[name] = close(actual[name], values, name)
    return dict(estimated=n, power_target_max_abs_error=power_error,
                arithmetic_max_absolute_differences=errors,
                endpoint_boundary_max_scaled_residual=residual, summary_max_abs_difference=summary_error)


def verify(mode, output):
    require(not output.exists(), 'Precision verification already exists; choose a new output')
    verifier_sha = sha(Path(__file__))
    folder, source, audit, meta, cov, evidence = context(mode)
    mde = pd.read_csv(folder/'conditional_MDE_by_phenotype.csv', float_precision='round_trip').set_index(
        ['policy', 'phenotype', 'reference'])
    covariance = pd.read_csv(folder/'covariance_sensitivity_by_phenotype.csv', float_precision='round_trip').set_index(
        ['policy', 'phenotype', 'rho_error'])
    require(mde.index.is_unique and covariance.index.is_unique, 'Duplicate summary groups')
    expected_mde = {(p, t, a) for p in ('drop', 'infer') for t in meta.phenotype for a in ALPHAS}
    expected_cov = {(p, t, rho) for p in ('drop', 'infer') for t in meta.phenotype for rho in RHOS}
    require(set(mde.index) == expected_mde and set(covariance.index) == expected_cov, 'Summary group coverage differs')
    checks = {}
    for policy in ('drop', 'infer'):
        blocks = source_blocks(source/f'full_family_palindrome_{policy}.csv.gz', mode, meta, cov)
        statuses, metrics, errors = Counter(), Counter(), {}
        with pd.read_csv(folder/f'pair_precision_{policy}.csv.gz', chunksize=len(cov),
                         float_precision='round_trip') as reader:
            for trait in meta.to_dict('records'):
                original, result = next(blocks, None), next(reader, None)
                require(original is not None and result is not None, 'Missing source/precision block')
                paired_keys(original, result, trait, cov.exposure_id.tolist())
                item = verify_block(original, result, policy, mde, covariance)
                statuses.update(result.status.value_counts().to_dict())
                for name in ('power_target_max_abs_error', 'endpoint_boundary_max_scaled_residual', 'summary_max_abs_difference'):
                    metrics[name] = max(metrics[name], item[name])
                for name, value in item['arithmetic_max_absolute_differences'].items():
                    errors[name] = max(errors.get(name, 0.), value)
            require(next(blocks, None) is None and next(reader, None) is None, 'Extra source/precision blocks')
        expected = audit['policies'][policy]
        require(dict(statuses) == expected['status_counts'] and sum(statuses.values()) == expected['rows']
                == len(meta)*len(cov) and statuses['estimated'] == expected['estimated'], 'Precision status/row totals differ')
        checks[policy] = dict(rows=sum(statuses.values()), status_counts=dict(statuses),
                              arithmetic_max_absolute_differences=errors, **dict(metrics))
        print(f"Verified precision {policy}: {sum(statuses.values()):,} planned pairs", flush=True)
    for name, digest in evidence.items():
        require(sha(ROOT/name) == digest, f'Input changed during verification: {name}')
    require(sha(Path(__file__)) == verifier_sha, 'Verifier changed during verification')
    report = dict(complete=True, mode=mode, formal_full_precision_verified=mode == 'full',
                  verified_at_utc=datetime.now(timezone.utc).isoformat(), verifier_sha256=verifier_sha,
                  inputs_sha256=evidence, policies=checks, checked_MDE_summary_rows=len(mde),
                  checked_covariance_summary_rows=len(covariance),
                  versions=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
                  production_precision_functions_imported=False, actual_sample_overlap_estimated=False,
                  overlap_bias_corrected=False, computational_verification_only=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['full', 'validation45'], default='full')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    output = args.output or BASE/f'precision_overlap_20261116/{args.mode}/COMPLETE_PRECISION_VERIFICATION.json'
    result = verify(args.mode, output.resolve())
    print(json.dumps(dict(complete=result['complete'], mode=args.mode, output=str(output)), indent=2))
