"""Annotate the fixed MR plan with detectability and covariance sensitivity.

validation45 is an implementation audit, never a new discovery scan. Full mode
requires the completed original full-family result and its verified SHA.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.mr_precision import detectable_effect, covariance_delta_se, bounded_fieller_interval

TOTAL = 2941*3963
RHOS = np.linspace(-1, 1, 9)
ALPHAS = {'nominal': .05, 'full_plan_bonferroni_reference': .05/TOTAL}
COLS = 'exposure_id phenotype family status beta_x se_x beta_y_aligned se_y beta se se_delta F'.split()


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        while block := stream.read(8*1024*1024):
            h.update(block)
    return h.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def inputs(mode):
    if mode == 'full':
        base = ROOT/'papers/01_mr_imaging_scan/results/full_atlas_20261116/full'
        run = json.loads((base/'RUN.json').read_text())
        if not run.get('complete') or run.get('planned_pairs') != TOTAL:
            raise ValueError('Full original MR plan is incomplete; no partial sensitivity scan')
        hashes = {p: run['policies'][p]['output_sha256'] for p in ['drop', 'infer']}
        ntraits = 3963
    else:
        base = ROOT/'papers/01_mr_imaging_scan/results/expanded_20261116/combined_45'
        hashes = {p: json.loads((base/f'full_family_palindrome_{p}.csv.gz_provenance.json').read_text())['result_sha256']
                  for p in ['drop', 'infer']}
        ntraits = 45
    paths = {p: base/f'full_family_palindrome_{p}.csv.gz' for p in hashes}
    for policy, path in paths.items():
        if sha(path) != hashes[policy]:
            raise ValueError(f'MR result SHA mismatch: {path}')
    return base, paths, hashes, ntraits


def groups(path, mode):
    if mode == 'validation45':
        yield from pd.read_csv(path, usecols=COLS).groupby('phenotype', sort=False)
    else:
        # The original complete-atlas writer emits one 2941-assay block per trait.
        for data in pd.read_csv(path, usecols=COLS, chunksize=2941):
            if data.phenotype.nunique() != 1:
                raise ValueError('Unexpected complete-atlas row layout')
            yield data.phenotype.iloc[0], data


def quantiles(values, prefix):
    if not len(values):
        return {prefix+'_'+q: None for q in ['p05', 'median', 'p95']}
    return dict(zip([prefix+'_'+q for q in ['p05', 'median', 'p95']],
                    map(float, np.quantile(values, [.05, .5, .95]))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=['validation45', 'full'], default='validation45')
    args = parser.parse_args()
    base, paths, input_hashes, ntraits = inputs(args.mode)
    output = ROOT/f'papers/01_mr_imaging_scan/results/precision_overlap_20261116/{args.mode}'
    output.mkdir(parents=True, exist_ok=True)
    if (output/'AUDIT.json').exists():
        raise ValueError('Completed audit exists; inspect before overwriting')
    metadata = pd.read_csv(base/'phenotypes.csv')
    assays = pd.read_csv(base/'assay_coverage.csv').exposure_id.tolist()
    assert len(assays) == len(set(assays)) == 2941 and metadata.phenotype.nunique() == ntraits
    metadata.to_csv(output/'phenotypes.csv', index=False)
    audit = dict(complete=False, mode=args.mode, original_full_plan_pairs=TOTAL,
                 analysed_pairs_per_policy=2941*ntraits, alpha_references=ALPHAS, rho_error_grid=RHOS.tolist(),
                 started_at_utc=datetime.now(timezone.utc).isoformat(),
                 input_sha256={str(paths[p].relative_to(ROOT)): h for p, h in input_hashes.items()},
                 metadata_sha256=sha(base/'phenotypes.csv'), assay_coverage_sha256=sha(base/'assay_coverage.csv'),
                 code_sha256={str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), ROOT/'pdv3/mr_precision.py',
                     ROOT/'papers/01_mr_imaging_scan/PRECISION_OVERLAP_SPEC_20261116.md']},
                 actual_overlap_estimated=False, overlap_bias_corrected=False,
                 independent_validation=False, new_MR_discovery=False)
    save(output/'RUN.json', audit)
    power_rows, covariance_rows = [], []
    counts, max_delta_diff = {}, 0.
    for policy, path in paths.items():
        written, estimated, seen, statuses = 0, 0, set(), {}
        out = output/f'pair_precision_{policy}.csv.gz'
        temporary = out.with_suffix('.tmp')
        with temporary.open('wb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as compressed, io.TextIOWrapper(compressed) as stream:
            for phenotype, data in groups(path, args.mode):
                if phenotype in seen or set(data.exposure_id) != set(assays) or len(data) != 2941:
                    raise ValueError('Incomplete or duplicate assay-by-trait group')
                seen.add(phenotype)
                data = data.set_index('exposure_id').loc[assays].reset_index()
                valid = data.status.eq('estimated')
                v = data[valid]
                n = len(v)
                bx, sx, by, sy = [v[c].to_numpy(float) for c in ['beta_x', 'se_x', 'beta_y_aligned', 'se_y']]
                s = sy/np.abs(bx)
                np.testing.assert_allclose(s, v.se, rtol=1e-12, atol=0)
                np.testing.assert_allclose(by/bx, v.beta, rtol=1e-12, atol=1e-15)
                np.testing.assert_allclose((bx/sx)**2, v.F, rtol=1e-12, atol=0)
                columns = data[['exposure_id', 'phenotype', 'family', 'status']].copy()
                columns['conditional_se'] = np.nan
                columns.loc[valid, 'conditional_se'] = s
                for label, alpha in ALPHAS.items():
                    record = dict(policy=policy, phenotype=phenotype, family=data.family.iloc[0],
                                  planned=2941, estimated=n, reference=label, alpha=alpha)
                    for power in [.8, .9]:
                        name = f'MDE{int(power*100)}_{label}'
                        values = detectable_effect(s, target_power=power, alpha=alpha)
                        columns[name] = np.nan
                        columns.loc[valid, name] = values
                        record.update(quantiles(values, f'MDE{int(power*100)}'))
                    power_rows.append(record)
                lo0, hi0 = bounded_fieller_interval(bx, sx, by, sy)
                delta0 = covariance_delta_se(bx, sx, by, sy)
                if n:
                    max_delta_diff = max(max_delta_diff, float(np.max(np.abs(delta0-v.se_delta))))
                np.testing.assert_allclose(delta0, v.se_delta, rtol=1e-12, atol=1e-15)
                lows, highs, deltas = [], [], []
                null_in = np.abs(by/sy) <= norm.isf(.025)
                for rho in RHOS:
                    lo, hi = bounded_fieller_interval(bx, sx, by, sy, rho=rho)
                    delta = covariance_delta_se(bx, sx, by, sy, rho=rho)
                    if not np.array_equal((lo <= 0) & (hi >= 0), null_in):
                        raise ValueError('Fieller null-inclusion identity failed')
                    record = dict(policy=policy, phenotype=phenotype, planned=2941, estimated=n,
                                  rho_error=float(rho), intervals_including_zero=int(null_in.sum()))
                    record.update(quantiles((hi-lo)/(hi0-lo0), 'Fieller_width_ratio_to_rho0'))
                    record.update(quantiles(delta/delta0, 'delta_SE_ratio_to_rho0'))
                    covariance_rows.append(record)
                    lows.append(lo); highs.append(hi); deltas.append(delta)
                derived = dict(Fieller95_rho0_lower=lo0, Fieller95_rho0_upper=hi0,
                               Fieller95_grid_lower_min=np.min(lows, axis=0),
                               Fieller95_grid_upper_max=np.max(highs, axis=0),
                               delta_SE_grid_min=np.min(deltas, axis=0), delta_SE_grid_max=np.max(deltas, axis=0))
                for name, values in derived.items():
                    columns[name] = np.nan
                    columns.loc[valid, name] = values
                if not columns.loc[~valid, columns.columns[4:]].isna().all().all():
                    raise ValueError('Missing slots were filled with estimates')
                columns.to_csv(stream, index=False, header=written == 0)
                written += len(data); estimated += n
                for status, count in data.status.value_counts().items():
                    statuses[status] = statuses.get(status, 0)+int(count)
        temporary.replace(out)
        assert written == 2941*ntraits and seen == set(metadata.phenotype)
        counts[policy] = dict(rows=written, estimated=estimated, status_counts=statuses, output_sha256=sha(out))
    pd.DataFrame(power_rows).to_csv(output/'conditional_MDE_by_phenotype.csv', index=False)
    pd.DataFrame(covariance_rows).to_csv(output/'covariance_sensitivity_by_phenotype.csv', index=False)
    audit.update(complete=True, finished_at_utc=datetime.now(timezone.utc).isoformat(), policies=counts,
                 delta_SE_rho0_max_difference=max_delta_diff,
                 all_nine_Fieller_null_inclusion_checks_passed=True,
                 missing_slots_preserved=True)
    save(output/'AUDIT.json', audit)
    save(output/'RUN.json', audit)
    save(output/'MANIFEST.json', {'artifacts_sha256': {str(p.relative_to(ROOT)): sha(p) for p in output.iterdir()
         if p.is_file() and p.name != 'MANIFEST.json'}})
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
