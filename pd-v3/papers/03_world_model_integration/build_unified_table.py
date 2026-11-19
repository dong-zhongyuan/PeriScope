#!/usr/bin/env python
"""One table: the five shared-representation comparators, CellOT, PeriScope and the two
ablation arms, all on the benchmark harness' own frozen test draws and metric.

Sources
-------
baseline_metrics.csv  per-seed per-axis values for the seven published methods. Readout
                      `raw_shared` for every shared-representation method (matching the
                      published Table A); CellOT only exists as `original_cellot`, since it
                      keeps its own independent implementation - that is why it is absent
                      from paired_comparisons.csv, and it is labelled as such below.
arm_per_cell.csv      coupling_off / random_pairing, produced by
                      run_ablation_arms_in_harness.py. Its `full` arm reproduces the
                      published periscope raw numbers to 6.6e-07, which is what licenses
                      putting the arms in the same table.

Aggregation: the eight routes share five brain donors, so a per-cell paired test over
axis x seed would treat correlated routes as independent. Average the eight axes within a
seed first, then compare across the three seeds.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd

BENCH = Path('/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002')
OUT = Path('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/ablation_in_harness_20261118')
CEILING = {'astro': 0.5537, 'microglia_mhc2': 0.6116}
METRIC = 'condition_aware_rho'

base = pd.read_csv(BENCH / 'baseline_metrics.csv')
pub = pd.concat([
    base[base.readout.eq('raw_shared')],
    base[base.method.eq('cellot') & base.readout.eq('original_cellot')],
], ignore_index=True)
pub = pub.rename(columns={'method': 'arm'})
pub['source'] = 'published_benchmark'

arms = pd.read_csv(OUT / 'arm_per_cell.csv')
arms = arms[arms.arm.isin(['coupling_off', 'random_pairing'])].copy()
arms['readout'] = 'raw'
arms['source'] = 'this_run'

keep = ['arm', 'readout', 'source', 'seed', 'axis', 'blood_state', 'brain_state',
        METRIC, 'fixed_control_rho', 'fixed_PD_rho', 'mmd', 'mean_mse', 'variance_ratio']
cell = pd.concat([pub[keep], arms[keep]], ignore_index=True)
cell['ceiling'] = cell.brain_state.map(CEILING)
cell['rho_over_ceiling'] = cell[METRIC] / cell.ceiling
counts = cell.groupby(['arm', 'seed'])['axis'].nunique()
assert counts.eq(8).all(), counts
cell.to_csv(OUT / 'TABLE_UNIFIED_per_cell.csv', index=False)

seedwise = cell.groupby(['arm', 'seed'], as_index=False).agg(
    rho=(METRIC, 'mean'), norm=('rho_over_ceiling', 'mean'),
    fixed_control=('fixed_control_rho', 'mean'), fixed_PD=('fixed_PD_rho', 'mean'),
    mmd=('mmd', 'mean'), mean_mse=('mean_mse', 'mean'))
seedwise.to_csv(OUT / 'TABLE_UNIFIED_per_seed.csv', index=False)

summ = seedwise.groupby('arm', as_index=False).agg(
    rho_mean=('rho', 'mean'), rho_sd=('rho', 'std'),
    norm_mean=('norm', 'mean'),
    fixed_control=('fixed_control', 'mean'), fixed_PD=('fixed_PD', 'mean'),
    mmd=('mmd', 'mean'), mean_mse=('mean_mse', 'mean'), n_seeds=('seed', 'nunique'))
summ['pct_of_ceiling'] = 100 * summ.norm_mean
summ = summ.sort_values('rho_mean', ascending=False).reset_index(drop=True)
summ.to_csv(OUT / 'TABLE_UNIFIED_summary.csv', index=False)

piv = seedwise.pivot(index='seed', columns='arm', values='rho')
piv.to_csv(OUT / 'TABLE_UNIFIED_seed_matrix.csv')

ref = 'periscope'
rows = []
for arm in piv.columns:
    if arm == ref:
        continue
    d = (piv[ref] - piv[arm]).dropna()
    rows.append(dict(comparator=arm, n_seeds=int(len(d)),
                     periscope_wins=int((d > 0).sum()),
                     margin_mean=float(d.mean()), margin_min=float(d.min()),
                     margin_max=float(d.max())))
paired = pd.DataFrame(rows).sort_values('margin_mean').reset_index(drop=True)
paired.to_csv(OUT / 'TABLE_UNIFIED_paired_vs_periscope.csv', index=False)

# how a method ranks on MMD, to show the two axes disagree
mmd_rank = summ[['arm', 'rho_mean', 'mmd']].copy()
mmd_rank['rho_rank'] = mmd_rank.rho_mean.rank(ascending=False).astype(int)
mmd_rank['mmd_rank'] = mmd_rank.mmd.rank(ascending=True).astype(int)
mmd_rank.to_csv(OUT / 'TABLE_UNIFIED_rank_disagreement.csv', index=False)

lines = []
lines.append('# Unified Table A: PeriScope, five comparators, CellOT and both ablation arms')
lines.append('')
lines.append('Metric: condition_aware_rho -- Spearman across 4741 genes between the predicted')
lines.append('and the observed PD-minus-control brain direction. Eight blood-state x brain-state')
lines.append('routes, seeds 42/43/44, benchmark harness test draws, raw (uncalibrated) readout.')
lines.append('Ceiling: astro 0.5537, microglia_mhc2 0.6116 (between-donor reliability, measured).')
lines.append('')
lines.append('| method | rho | sd across seeds | % of ceiling | MMD | seed42 | seed43 | seed44 |')
lines.append('|---|---|---|---|---|---|---|---|')
for _, r in summ.iterrows():
    s = piv[r.arm]
    lines.append(f"| {r.arm} | {r.rho_mean:.4f} | {r.rho_sd:.4f} | {r.pct_of_ceiling:.1f}% | "
                 f"{r.mmd:.4f} | {s.get(42, float('nan')):+.4f} | {s.get(43, float('nan')):+.4f} | "
                 f"{s.get(44, float('nan')):+.4f} |")
lines.append('')
lines.append('## Seed-matched paired contrast against PeriScope')
lines.append('')
lines.append('| comparator | PeriScope wins | mean margin | min | max |')
lines.append('|---|---|---|---|---|')
for _, r in paired.iterrows():
    lines.append(f"| {r.comparator} | {r.periscope_wins}/{r.n_seeds} | {r.margin_mean:+.4f} | "
                 f"{r.margin_min:+.4f} | {r.margin_max:+.4f} |")
lines.append('')
lines.append('## Direction rank vs MMD rank')
lines.append('')
lines.append('| method | direction rank | MMD rank |')
lines.append('|---|---|---|')
for _, r in mmd_rank.sort_values('rho_rank').iterrows():
    lines.append(f"| {r.arm} | {r.rho_rank} | {r.mmd_rank} |")
lines.append('')
lines.append('Notes')
lines.append('')
lines.append('- CellOT is on readout `original_cellot`: it keeps its own independent implementation')
lines.append('  rather than the shared frozen representation, which is why it is missing from')
lines.append('  paired_comparisons.csv. Its per-seed values do exist in baseline_metrics.csv, so the')
lines.append('  seed-matched paired contrast above is computed here for the first time.')
lines.append('- The two arms were run by run_ablation_arms_in_harness.py. The `full` control in that')
lines.append('  run reproduces the published periscope raw numbers to max |diff| 6.6e-07 on this')
lines.append('  metric, which is what makes the arm rows comparable to the published rows.')
lines.append('- The arms keep the harness test draws and the seed-42 latent noise draw; only the')
lines.append('  weights (encoder, equilibrium solve, brain decoder) are the arm\'s own.')
(OUT / 'TABLE_UNIFIED.md').write_text('\n'.join(lines) + '\n')

print('\n'.join(lines))
print('\nwrote', OUT / 'TABLE_UNIFIED.md')
