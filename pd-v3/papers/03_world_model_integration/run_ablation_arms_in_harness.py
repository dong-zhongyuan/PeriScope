#!/usr/bin/env python
"""Put the two ablation arms inside the frozen shared-representation benchmark harness.

Nothing under pd_product or pd_product_assets is written. The frozen harness modules are
imported unmodified; only the harness output root (run_benchmark.O) is rebound to this
output directory before current_benchmark is imported, so every file this run produces
lands here. The axis cache, seed representations and the published results are reached
through read-only symlinks.

Comparison logic
----------------
Published Table A uses readout `raw_shared` (paired_comparisons.csv): raw predictions, no
validation calibration. Direction rho is invariant to that calibration anyway - the
harness itself asserts the condition contrast is unchanged by it - so raw is both the
simplest and the matching choice.

For each arm the harness' frozen test-cell draws, axis cache, condition vectors and brain
observations are reused verbatim. Only the blood encodings are recomputed, with that arm's
own enc_blood / protein_head, because an arm has its own encoder; the equilibrium solve and
the brain decoder are then also the arm's own. That is exactly the path the harness uses
for `periscope`, so `full` doubles as a reproduction control and must return the published
periscope raw numbers to within float noise.
"""
import json, hashlib, os, sys, time, datetime
from pathlib import Path
import numpy as np

PRODUCT = Path('/public/home/mengxl/dzy/pd_product')
ASSETS = Path('/public/home/mengxl/dzy/pd_product_assets')
HARNESS = PRODUCT / '_work/benchmark_shared_20261003'
BENCH = ASSETS / 'results/benchmark_20261002'
OPT = ASSETS / 'results/optimization_20261003'
OUT = Path('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/ablation_in_harness_20261118')

SEEDS = [42, 43, 44]
# The published per-seed representations (benchmark_20261002/representations/seed_N) were
# built by re-encoding with the seed-N weights while KEEPING the seed-42 latent noise draw,
# so that only the weights differ between seeds. Verified numerically on seed 43:
# noise seed 42+100*bc reproduces the stored test_z to 7e-07, while 43+100*bc differs by 5.3.
# Any re-encoding here must use the same convention or it is not comparable to Table A.
PREPARE_SEED = 42
ARMS = {
    'full':           str(OPT / 'models/seed_{seed}/model.pt'),
    'coupling_off':   str(OPT / 'ablations/coupling_off/models/seed_{seed}/model.pt'),
    'random_pairing': str(OPT / 'ablations/random_pairing/models/seed_{seed}/model.pt'),
}
CEILING = {'astro': 0.5537, 'microglia_mhc2': 0.6116}

OUT.mkdir(parents=True, exist_ok=True)
# read-only views the frozen harness expects to find next to its output root
for name in ('cache', 'representations', 'prepared_eight.json'):
    link = OUT / name
    if not link.exists():
        link.symlink_to(BENCH / name)

os.environ.setdefault('BENCHMARK_DEVICE', 'cpu')
os.environ.setdefault('BENCHMARK_SEED', '42')
sys.path[:0] = [str(HARNESS), str(HARNESS / 'source')]
import run_benchmark as rb
rb.O = OUT                      # rebind BEFORE current_benchmark copies it
import current_benchmark as cb
assert cb.O == OUT, cb.O
import torch
import pandas as pd


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def write_overrides(arm):
    """Point the harness' own checkpoint resolver at this arm."""
    (OUT / 'checkpoint_overrides.json').write_text(json.dumps({
        'reason': f'ablation arm {arm} evaluated inside the shared-representation harness (20261118)',
        'checkpoints': {str(s): ARMS[arm].format(seed=s) for s in SEEDS},
    }, indent=2))


def arm_axis(arm, seed, axis, out_dir):
    """Raw periscope-path prediction for one arm/seed/axis, scored by the harness."""
    model = rb.load_model(seed, cb.DEVICE)
    for p in model.parameters():
        p.requires_grad_(False)
    k = cb.read(axis)
    # the arm has its own encoder: re-encode the SAME frozen test cells
    for bc in (0, 1):
        z, mu, u = cb.encode(model, k[f'test_x{bc}'], seed=PREPARE_SEED + 100 * bc)
        k[f'test_z{bc}'] = z
        k[f'test_u{bc}'] = u
    pred = {c: {} for c in (0, 1)}
    for c in (0, 1):
        for bc in (0, 1):
            pred[c][bc] = cb.full_predict(model, k, c, bc)
    keep = cb.RESULT_DIR
    cb.RESULT_DIR = out_dir
    cb.RESULT_DIR.mkdir(parents=True, exist_ok=True)
    cb.assess(arm, pred, k, seed, {'readout': 'raw', 'path': 'arm enc_blood + protein_head -> arm solve_equilibrium -> arm dec_brain',
                                   'test_cells': 'frozen harness draws, reused verbatim'})
    cb.RESULT_DIR = keep
    tag = f'{arm}__{str(k["blood_state"])}__{str(k["brain_state"])}'
    return json.loads((out_dir / f'{tag}.json').read_text())


def flatten(rec, arm, seed, axis):
    d = rec['distribution_average']
    row = dict(arm=arm, seed=seed, axis=axis,
               blood_state=rec['blood_state'], brain_state=rec['brain_state'],
               mean_mse=d['mean_mse'], energy_per_sqrt_gene=d['energy_per_sqrt_gene'],
               mmd=d['mmd'], variance_ratio=d['variance_ratio'])
    for task in ('condition_aware', 'fixed_control', 'fixed_PD'):
        row[f'{task}_rho'] = rec['direction'][task]['rho']
    row['condition_aware_delta_mse'] = rec['direction']['condition_aware']['delta_mse']
    row['condition_aware_predicted_rms'] = rec['direction']['condition_aware']['predicted_rms']
    row['condition_aware_observed_rms'] = rec['direction']['condition_aware']['observed_rms']
    row['ceiling'] = CEILING[rec['brain_state']]
    row['rho_over_ceiling'] = row['condition_aware_rho'] / row['ceiling']
    return row


def published_raw():
    """Published raw_shared per-cell values for every benchmark method."""
    p = pd.read_csv(BENCH / 'paired_comparisons.csv')
    p = p[p.readout.eq('raw_shared')]
    cols = ['mmd', 'mean_mse', 'energy_per_sqrt_gene', 'variance_ratio',
            'condition_aware_rho', 'fixed_control_rho', 'fixed_PD_rho']
    rows = []
    base = p.drop_duplicates(['seed', 'axis'])
    for _, r in base.iterrows():
        rows.append(dict(arm='periscope', seed=int(r.seed), axis=r.axis,
                         blood_state=r.blood_state_periscope, brain_state=r.brain_state_periscope,
                         **{c: float(r[f'{c}_periscope']) for c in cols}))
    for _, r in p.iterrows():
        rows.append(dict(arm=r.method_other, seed=int(r.seed), axis=r.axis,
                         blood_state=r.blood_state_other, brain_state=r.brain_state_other,
                         **{c: float(r[f'{c}_other']) for c in cols}))
    d = pd.DataFrame(rows).drop_duplicates(['arm', 'seed', 'axis'])
    d['ceiling'] = d.brain_state.map(CEILING)
    d['rho_over_ceiling'] = d.condition_aware_rho / d.ceiling
    return d


def per_seed(d):
    """Axis-average within a seed: the eight axes share five brain donors, so a per-cell
    paired test over axis x seed would treat correlated routes as independent."""
    g = d.groupby(['arm', 'seed'], as_index=False).agg(
        rho=('condition_aware_rho', 'mean'),
        norm=('rho_over_ceiling', 'mean'),
        fixed_control_rho=('fixed_control_rho', 'mean'),
        fixed_PD_rho=('fixed_PD_rho', 'mean'),
        mmd=('mmd', 'mean'), mean_mse=('mean_mse', 'mean'),
        n_axes=('axis', 'nunique'))
    return g


def main():
    t0 = time.time()
    rows = []
    for arm in ARMS:
        write_overrides(arm)
        for seed in SEEDS:
            ck = Path(ARMS[arm].format(seed=seed))
            assert rb.checkpoint_path(seed) == ck, (rb.checkpoint_path(seed), ck)
            out_dir = OUT / 'results' / arm / f'seed_{seed}'
            for axis in cb.AXES:
                rec = arm_axis(arm, seed, axis, out_dir)
                rows.append(flatten(rec, arm, seed, axis))
                print(f'DONE {arm} seed{seed} {axis} rho={rows[-1]["condition_aware_rho"]:+.4f} '
                      f'({time.time()-t0:.0f}s)', flush=True)
    mine = pd.DataFrame(rows)
    mine.to_csv(OUT / 'arm_per_cell.csv', index=False)

    pub = published_raw()
    pub.to_csv(OUT / 'published_raw_per_cell.csv', index=False)

    # --- reproduction control: arm "full" must equal published periscope ---
    a = mine[mine.arm.eq('full')].set_index(['seed', 'axis'])
    b = pub[pub.arm.eq('periscope')].set_index(['seed', 'axis'])
    chk = []
    for key in sorted(set(a.index) & set(b.index)):
        for col in ('condition_aware_rho', 'mmd', 'mean_mse', 'fixed_control_rho', 'fixed_PD_rho'):
            chk.append(dict(seed=key[0], axis=key[1], metric=col,
                            recomputed=float(a.loc[key, col]), published=float(b.loc[key, col]),
                            abs_diff=abs(float(a.loc[key, col]) - float(b.loc[key, col]))))
    chk = pd.DataFrame(chk)
    chk.to_csv(OUT / 'reproduction_check.csv', index=False)
    worst = chk.groupby('metric').abs_diff.max().to_dict()
    print('\nREPRODUCTION CONTROL (arm full vs published periscope, raw_shared)')
    for m, v in worst.items():
        print(f'  max |diff| {m:24s} {v:.3e}')
    rho_ok = worst['condition_aware_rho'] < 1e-6
    print('  VERDICT:', 'PASS' if rho_ok else 'FAIL - arm numbers are NOT comparable to Table A')

    # --- unified table ---
    combined = pd.concat([
        pub.assign(source='published_benchmark'),
        mine[mine.arm.ne('full')].assign(source='this_run'),
        mine[mine.arm.eq('full')].assign(arm='periscope_recomputed', source='this_run'),
    ], ignore_index=True)
    combined.to_csv(OUT / 'unified_per_cell.csv', index=False)

    ps = per_seed(combined)
    ps.to_csv(OUT / 'unified_per_seed.csv', index=False)

    summ = ps.groupby('arm', as_index=False).agg(
        rho_mean=('rho', 'mean'), rho_sd=('rho', 'std'),
        norm_mean=('norm', 'mean'), norm_sd=('norm', 'std'),
        fixed_control=('fixed_control_rho', 'mean'), fixed_PD=('fixed_PD_rho', 'mean'),
        mmd=('mmd', 'mean'), n_seeds=('seed', 'nunique'), n_axes=('n_axes', 'max'))
    summ = summ.sort_values('rho_mean', ascending=False)
    summ.to_csv(OUT / 'unified_summary.csv', index=False)

    # --- seed-matched paired contrast against periscope ---
    piv = ps.pivot(index='seed', columns='arm', values='rho')
    ref = 'periscope'
    pr = []
    for arm in piv.columns:
        if arm == ref:
            continue
        d = (piv[ref] - piv[arm]).dropna()
        pr.append(dict(comparator=arm, n_seeds=int(len(d)),
                       periscope_wins=int((d > 0).sum()),
                       margin_mean=float(d.mean()), margin_min=float(d.min()),
                       margin_max=float(d.max())))
    pr = pd.DataFrame(pr).sort_values('margin_mean')
    pr.to_csv(OUT / 'unified_paired_vs_periscope.csv', index=False)

    man = dict(created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               elapsed_s=round(time.time() - t0, 1), device=cb.DEVICE,
               readout='raw_shared (uncalibrated), matching published Table A',
               axes=cb.AXES, seeds=SEEDS,
               harness_files={p.name: sha(p) for p in sorted(HARNESS.glob('*.py'))},
               harness_source={p.name: sha(p) for p in sorted((HARNESS / 'source').glob('*.py'))},
               protocol_sha256=sha(HARNESS / 'protocol.txt'),
               checkpoints={arm: {str(s): sha(Path(ARMS[arm].format(seed=s))) for s in SEEDS}
                            for arm in ARMS},
               reproduction_control=dict(max_abs_diff=worst, rho_reproduced=bool(rho_ok)),
               ceiling=CEILING,
               note='frozen harness imported unmodified; only run_benchmark.O rebound to this '
                    'directory, so no file under pd_product or pd_product_assets is written')
    (OUT / 'MANIFEST.json').write_text(json.dumps(man, indent=2))

    print('\nUNIFIED SUMMARY (axis-averaged within seed, then across seeds)')
    print(summ.to_string(index=False))
    print('\nSEED-MATCHED PAIRED vs periscope')
    print(pr.to_string(index=False))
    print('\nPER SEED')
    print(piv.to_string())
    print(f'\nwrote {OUT}')


if __name__ == '__main__':
    main()
