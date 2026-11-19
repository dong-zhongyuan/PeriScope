"""Audit a completed atlas without importing or rerunning its MR kernels.

The 45-trait mode validates this reader against real completed output. Full
mode refuses partial input, then examines every declared pair in both policies.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import gc
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.special import erfc, log_ndtr
from scipy.stats import norm
import statsmodels
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'papers/01_mr_imaging_scan/results/full_atlas_20261116'
QCOLS = ['q_bh_global', 'q_by_global', 'q_bh_cis_family', 'q_by_cis_family',
         'q_bh_family', 'q_storey_lambda_0p25', 'q_storey_lambda_0p5',
         'q_storey_lambda_0p75']
NUMERIC = ['beta', 'se', 'p', 'log_p', 'ci_lower', 'ci_upper', 'se_delta',
           'p_delta', 'log_p_delta', 'ci_lower_delta', 'ci_upper_delta', 'F']


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked(path, expected, evidence):
    actual = sha(path)
    require(actual == expected, f'SHA mismatch: {path}')
    evidence[str(path.relative_to(ROOT))] = actual
    return path


def close(actual, expected, name, probability=False):
    a, b = np.asarray(actual, float), np.asarray(expected, float)
    require(a.shape == b.shape, f'Shape mismatch: {name}')
    require(np.array_equal(np.isnan(a), np.isnan(b)), f'Missingness mismatch: {name}')
    # Log P is checked separately, including P underflow. This small absolute
    # tolerance only concerns subnormal/tiny probabilities, not ordinary ones.
    np.testing.assert_allclose(a, b, rtol=1e-10,
                               atol=1e-300 if probability else 1e-12,
                               equal_nan=True, err_msg=name)
    finite = np.isfinite(a) & np.isfinite(b)
    return float(np.max(abs(a[finite] - b[finite]))) if finite.any() else 0.


def completed_run(folder, mode):
    require(mode in ('full', 'validation45'), 'Unknown audit mode')
    path = folder / 'RUN.json'
    require(path.is_file(), 'No completed result RUN; partial atlas cannot be audited as complete')
    run = json.loads(path.read_text())
    ntraits = 3963 if mode == 'full' else 45
    require(run.get('complete') is True and run.get('status') == 'complete'
            and run.get('mode') == mode and run.get('assays') == 2941
            and run.get('phenotypes') == ntraits
            and run.get('completed_traits') == ntraits
            and run.get('planned_pairs') == 2941 * ntraits,
            'Completed run does not cover the exact declared atlas')
    require(set(run.get('policies', {})) == {'drop', 'infer'}, 'Both policies are required')
    return run


def read_context(folder, mode):
    run = completed_run(folder, mode)
    evidence = {str((folder / 'RUN.json').relative_to(ROOT)): sha(folder / 'RUN.json')}
    for mapping in (run['input_sha256'], run['source_code_sha256']):
        for name, digest in mapping.items():
            checked(ROOT / name, digest, evidence)
    cov = pd.read_csv(ROOT / 'data/ukb_ppp/extracted/assay_coverage.csv')
    inst = pd.read_csv(ROOT / 'data/ukb_ppp/extracted/cis_primary_selected.csv')
    require(len(cov) == 2941 and cov.exposure_id.is_unique, 'Assay plan differs')
    require(len(inst) == 1954 and inst.exposure_id.is_unique, 'Cis selection differs')
    combined = ROOT / 'data/imaging_combined'
    manifest = json.loads((combined / 'MANIFEST.json').read_text())
    require(manifest.get('complete') is True, 'Original 45-trait input incomplete')
    prior = pd.read_csv(checked(combined / 'phenotypes_45.csv',
                              manifest['phenotype_plan_sha256'], evidence))
    if mode == 'validation45':
        meta = prior
    else:
        data = ROOT / 'data/big40_atlas'
        plan = json.loads((data / 'PLAN.json').read_text())
        progress = json.loads((data / 'extraction_progress.json').read_text())
        require(progress.get('complete') is True and progress.get('completed') == 3935
                and progress.get('requested') == 3935 and not progress.get('errors'),
                'Original full acquisition is incomplete')
        catalog = pd.read_csv(checked(data / 'phenotypes_3935.csv',
                                     plan['plan_sha256'], evidence), dtype={'idp': str})
        catalog['trait'] = catalog.short_name
        catalog['modality'] = catalog.category
        meta = pd.concat([catalog, prior[prior.phenotype.str.startswith('OxfordQSM2022:')]],
                         ignore_index=True)
        require(catalog.is_qc.eq(True).sum() == 16, 'Original QC membership differs')
    saved = pd.read_csv(checked(folder / 'phenotypes.csv', run['phenotypes_sha256'], evidence))
    require(len(meta) == run['phenotypes'] and meta.phenotype.is_unique, 'Phenotype plan differs')
    for col in ('phenotype', 'family', 'trait'):
        require(saved[col].fillna('').tolist() == meta[col].fillna('').tolist(),
                f'Frozen metadata changed: {col}')
    for policy in ('drop', 'infer'):
        checked(folder / f'full_family_palindrome_{policy}.csv.gz',
                run['policies'][policy]['output_sha256'], evidence)
    return run, cov, inst, meta, evidence


def independent_alignment(g, policy):
    if g.empty:
        return 0
    ea, oa, ya, yo = [g[c].astype(str).str.upper().str.strip() for c in
                      ('ea_x', 'oa_x', 'ea_y', 'oa_y')]
    require(all(s.isin(list('ACGT')).all() for s in (ea, oa, ya, yo)),
            'Estimated unsupported allele')
    require((ea != oa).all() and (ya != yo).all(), 'Estimated identical alleles')
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    forward = (ea.eq(ya) & oa.eq(yo)) | (ea.eq(ya.map(complement)) & oa.eq(yo.map(complement)))
    reverse = (ea.eq(yo) & oa.eq(ya)) | (ea.eq(yo.map(complement)) & oa.eq(ya.map(complement)))
    pal = ea.eq(oa.map(complement))
    require((forward | reverse).all(), 'Estimated allele mismatch')
    require(not pal.any() or policy == 'infer', 'Drop policy estimated a palindrome')
    fx, fy = g.eaf_x.to_numpy(float), g.eaf_y.to_numpy(float)
    p = pal.to_numpy()
    if p.any():
        require(np.isfinite(fx[p]).all() and np.isfinite(fy[p]).all(), 'Missing palindrome EAF')
        require(((fx[p] >= 0) & (fx[p] <= 1) & (fy[p] >= 0) & (fy[p] <= 1)).all(), 'Invalid EAF')
        require((np.maximum(np.minimum(fx[p], 1-fx[p]), np.minimum(fy[p], 1-fy[p])) <= .42).all(),
                'Ambiguous palindrome estimated')
    direct = forward.to_numpy() & (~p | (abs(fx-fy) <= .1))
    swapped = reverse.to_numpy() & (~p | (abs(fx-(1-fy)) <= .1))
    require(np.logical_xor(direct, swapped).all(), 'Non-unique allele orientation')
    sign = np.where(direct, 1., -1.)
    close(g.beta_y_aligned, sign*g.beta_y, 'aligned effect')
    close(g.eaf_y_aligned, np.where(direct, fy, 1-fy), 'aligned EAF')
    require(np.array_equal(g.allele_flipped.to_numpy(), swapped), 'Flip indicator mismatch')
    return int(p.sum())


def independent_estimates(g):
    if g.empty:
        return {c: 0. for c in NUMERIC}
    bx, sx, by, sy = [g[c].to_numpy(float) for c in ('beta_x', 'se_x', 'beta_y_aligned', 'se_y')]
    require(np.isfinite(np.array([bx, sx, by, sy])).all() and (bx != 0).all()
            and (sx > 0).all() and (sy > 0).all(), 'Invalid estimated inputs')
    beta = by/bx
    se = sy/abs(bx)
    delta = np.sqrt(sy**2/bx**2 + by**2*sx**2/bx**4)
    values = {'beta': beta, 'F': (bx/sx)**2}
    for suffix, error in (('', se), ('_delta', delta)):
        z = abs(beta/error)
        values.update({f'se{suffix}': error, f'p{suffix}': erfc(z/np.sqrt(2)),
                       f'log_p{suffix}': np.log(2)+log_ndtr(-z),
                       f'ci_lower{suffix}': beta-norm.isf(.025)*error,
                       f'ci_upper{suffix}': beta+norm.isf(.025)*error})
    require(g.method.eq('Wald').all() and g.n_instruments.eq(1).all(), 'Unexpected estimator')
    return {c: close(g[c], values[c], c, c in ('p', 'p_delta')) for c in NUMERIC}


def external_adjust(p, method, pi0=1.):
    p = np.asarray(p, float)
    missing = np.isnan(p)
    # Scale before the external BH routine clips, preserving the fixed-lambda
    # Storey definition. Missing planned slots participate as P=1 before scaling.
    q = multipletests(np.where(missing, 1., p)*pi0, method=method)[1]
    q[missing] = np.nan
    return q


def audit_policy(folder, policy, run, cov, inst, meta):
    n, size = run['planned_pairs'], len(cov)
    cis_one = cov.exposure_id.isin(inst.exposure_id).to_numpy()
    fixed = inst.set_index('exposure_id').reindex(cov.exposure_id)
    cis = np.tile(cis_one, len(meta))
    p = np.empty(n)
    q = {c: np.empty(n) for c in QCOLS}
    statuses, numerical, families = Counter(), {}, []
    palindromes = 0
    expected = run['policies'][policy]
    require(expected['family_n'] == n and expected['cis_family_n'] == int(cis.sum()),
            'Run correction denominators differ')
    family_sizes = meta.groupby('family').size().to_dict()
    denoms, pi0s = [], {}
    source = folder / f'full_family_palindrome_{policy}.csv.gz'
    with pd.read_csv(source, chunksize=size, float_precision='round_trip') as reader:
        for i, trait in enumerate(meta.to_dict('records')):
            try:
                block = next(reader)
            except StopIteration as error:
                raise ValueError('Incomplete Cartesian result') from error
            require(len(block) == size and block.exposure_id.tolist() == cov.exposure_id.tolist(),
                    'Assay count, order or identity differs')
            for col in ('gene', 'assay'):
                require(block[col].fillna('').tolist() == cov[col].fillna('').tolist(), 'Assay annotation differs')
            for col in ('phenotype', 'family', 'trait'):
                require(block[col].eq(trait[col]).all(), f'Phenotype annotation differs: {col}')
            require(block.palindrome_policy.eq(policy).all(), 'Mixed policy')
            require(block.status.notna().all(), 'Missing explicit status')
            require(np.array_equal(block.has_cis_instrument.to_numpy(), cis_one), 'Cis membership differs')
            for col in ('beta_x', 'se_x', 'eaf_x', 'pos37', 'pos38'):
                close(block[col], fixed[col], f'Frozen exposure tool {col}')
            for col in ('variant', 'ea_x', 'oa_x'):
                require(block[col].fillna('').tolist() == fixed[col].fillna('').tolist(),
                        f'Frozen exposure tool differs: {col}')
            valid = block.status.eq('estimated')
            require(block.loc[~valid, NUMERIC+QCOLS].isna().all().all(), 'Unestimated slot contains a result')
            require(np.isfinite(block.loc[valid, NUMERIC+QCOLS].to_numpy()).all(), 'Estimated numerical value is nonfinite')
            require(((block.loc[valid, ['p', 'p_delta', *QCOLS]] >= 0)
                     & (block.loc[valid, ['p', 'p_delta', *QCOLS]] <= 1)).all().all(), 'P/q outside [0,1]')
            require((block.loc[valid, ['log_p', 'log_p_delta']] <= 0).all().all(), 'Positive log P')
            for col, value in {'family_n_global': n, 'cis_family_n': int(cis.sum()),
                               'family_n': size*family_sizes[trait['family']]}.items():
                require(block[col].eq(value).all(), f'Incorrect denominator: {col}')
            observed = {}
            for col in ('family_n_observed_global', 'cis_family_n_observed', 'family_n_observed'):
                require(block[col].nunique() == 1, f'Inconsistent denominator: {col}')
                # Preserve the actual value: truncating to int could hide a
                # corrupted fractional observed denominator.
                observed[col] = float(block[col].iloc[0])
            denoms.append(observed)
            for label in ('0p25', '0p5', '0p75'):
                col = f'pi0_lambda_{label}'
                require(block[col].nunique() == 1, f'Inconsistent pi0: {col}')
                value = float(block[col].iloc[0])
                require(col not in pi0s or value == pi0s[col], 'Pi0 changed within a policy')
                pi0s[col] = value
            palindromes += independent_alignment(block.loc[valid], policy)
            estimated = block.loc[valid]
            require(estimated.build_x.eq(estimated.build_y).all(), 'Estimated genome builds differ')
            close(estimated.pos_x, estimated.pos_y, 'Estimated source coordinates')
            normalize_chr = lambda s: s.astype(str).str.replace('.0', '', regex=False).replace('23', 'X')
            require(normalize_chr(estimated.chrom_x).eq(normalize_chr(estimated.chrom_y)).all(),
                    'Estimated source chromosomes differ')
            for col, error in independent_estimates(block.loc[valid]).items():
                numerical[col] = max(numerical.get(col, 0.), error)
            xchr = block.chrom_y.astype(str).str.replace('.0', '', regex=False).isin(['X', '23'])
            require(block.loc[xchr, 'n_y'].isna().all(), 'X sample size was invented')
            if str(trait['phenotype']).startswith('OxfordQSM2022:'):
                require(block.n_y.isna().all(), 'QSM trait-specific sample size was invented')
            sl = slice(i*size, (i+1)*size)
            p[sl] = block.p
            for col in QCOLS:
                q[col][sl] = block[col]
            statuses.update(block.status.value_counts().to_dict())
        require(next(reader, None) is None, 'Rows beyond the declared Cartesian result')
    require(dict(statuses) == expected['status_counts'] and statuses['estimated'] == expected['estimated'],
            'Completed result status counts disagree')
    observed = int(np.isfinite(p).sum())
    cis_observed = int(np.isfinite(p[cis]).sum())
    require(all(d['family_n_observed_global'] == observed and d['cis_family_n_observed'] == cis_observed
                for d in denoms), 'Observed global/cis denominator differs')
    corrections = {}
    for name, method, mask in [('q_bh_global', 'fdr_bh', slice(None)),
                               ('q_by_global', 'fdr_by', slice(None)),
                               ('q_bh_cis_family', 'fdr_bh', cis),
                               ('q_by_cis_family', 'fdr_by', cis)]:
        corrections[name] = close(q[name][mask], external_adjust(p[mask], method), name, True)
    filled = np.where(np.isnan(p), 1., p)
    for label, lam in (('0p25', .25), ('0p5', .5), ('0p75', .75)):
        pi0 = min(1., (1 + np.count_nonzero(filled > lam))/(n*(1-lam)))
        close([pi0s[f'pi0_lambda_{label}']], [pi0], 'Fixed-lambda pi0')
        name = f'q_storey_lambda_{label}'
        corrections[name] = close(q[name], external_adjust(p, 'fdr_bh', pi0), name, True)
    for family in meta.family.drop_duplicates():
        traits = meta.family.eq(family).to_numpy()
        mask = np.repeat(traits, size)
        actual_observed = int(np.isfinite(p[mask]).sum())
        require(all(denoms[i]['family_n_observed'] == actual_observed for i in np.flatnonzero(traits)),
                'Observed category denominator differs')
        error = close(q['q_bh_family'][mask], external_adjust(p[mask], 'fdr_bh'), f'BH category {family}', True)
        families.append(dict(family=family, planned=int(mask.sum()), observed=actual_observed,
                             q_bh_max_absolute_difference=error))
    return dict(rows=n, status_counts=dict(statuses), independently_oriented_estimates=observed,
                palindrome_estimates=palindromes, arithmetic_max_absolute_differences=numerical,
                correction_max_absolute_differences=corrections, categories=families,
                all_planned_slots_retained=True, all_missing_result_slots_preserved=True)


def policy_pair(drop, infer):
    for col in ('exposure_id', 'phenotype'):
        require(drop[col].tolist() == infer[col].tolist(), 'Policy pair keys differ')
    previous = drop.status.eq('estimated')
    require(infer.loc[previous, 'status'].eq('estimated').all(), 'Infer policy lost an original estimate')
    for col in NUMERIC:
        close(drop.loc[previous, col], infer.loc[previous, col], f'Policy-invariant {col}', col in ('p', 'p_delta'))
    changed = drop.status.ne(infer.status)
    require(drop.loc[changed, 'status'].eq('palindromic_dropped').all(), 'A non-palindrome status changed')
    require(not infer.status.eq('palindromic_dropped').any(), 'Infer policy retained a drop-only status')
    return Counter(zip(drop.status, infer.status))


def policy_relation(folder, ntraits):
    columns = ['exposure_id', 'phenotype', 'status', *NUMERIC]
    transitions = Counter()
    with pd.read_csv(folder/'full_family_palindrome_drop.csv.gz', usecols=columns,
                     chunksize=2941, float_precision='round_trip') as a, \
         pd.read_csv(folder/'full_family_palindrome_infer.csv.gz', usecols=columns,
                     chunksize=2941, float_precision='round_trip') as b:
        for _ in range(ntraits):
            x, y = next(a, None), next(b, None)
            require(x is not None and y is not None and len(x) == len(y) == 2941,
                    'Policy comparison has incomplete blocks')
            transitions.update(policy_pair(x, y))
        require(next(a, None) is None and next(b, None) is None, 'Extra policy rows')
    return [dict(drop=x, infer=y, pairs=n) for (x, y), n in sorted(transitions.items())]


def verify(folder, mode, output):
    require(not output.exists(), f'Audit output already exists: {output}')
    verifier_sha = sha(Path(__file__))
    run, cov, inst, meta, evidence = read_context(folder, mode)
    checks = {}
    for policy in ('drop', 'infer'):
        checks[policy] = audit_policy(folder, policy, run, cov, inst, meta)
        print(f"Verified {policy}: {checks[policy]['rows']:,} planned pairs", flush=True)
        gc.collect()
    transitions = policy_relation(folder, len(meta))
    for name, digest in evidence.items():
        require(sha(ROOT / name) == digest, f'Input changed during audit: {name}')
    require(sha(Path(__file__)) == verifier_sha, 'Verifier changed during audit')
    report = dict(complete=True, mode=mode, formal_full_atlas_verified=mode == 'full',
                  verified_at_utc=datetime.now(timezone.utc).isoformat(),
                  inputs_sha256=evidence, verifier_sha256=verifier_sha, policies=checks,
                  policy_status_transitions=transitions,
                  drop_estimates_preserved_under_inference=True,
                  independent_adjustment_backend='statsmodels.stats.multitest.multipletests',
                  versions=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
                                statsmodels=statsmodels.__version__),
                  production_mr_or_adjustment_kernels_imported=False,
                  scope='Complete Cartesian coverage, retained missing slots, allele orientation, '
                        'Wald/delta arithmetic and all eight correction columns for both policies. '
                        'Verifies completed result arithmetic, not instrument validity, causal effects, '
                        'unmeasured sample overlap, or independent regional replication.')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['full', 'validation45'], default='full')
    parser.add_argument('--results', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    folder = (args.results or BASE / ('full' if args.mode == 'full' else 'validation45_nativefix')).resolve()
    output = (args.output or folder / 'COMPLETE_RESULT_VERIFICATION.json').resolve()
    result = verify(folder, args.mode, output)
    print(json.dumps(dict(complete=result['complete'], mode=args.mode, output=str(output)), indent=2))
