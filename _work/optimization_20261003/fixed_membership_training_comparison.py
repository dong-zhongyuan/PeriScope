"""Isolate training changes while retaining previously defined gene membership.

Retrospective diagnostics only. These tables do not nominate targets or replace
the discovery/confirmation analysis of the optimized primary branch.
"""
from pathlib import Path
import json
import time
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

A = Path('/public/home/mengxl/dzy/pd_product_assets')
OLD = A / 'results/rescreen_20261002'
O = A / 'results/optimization_20261003'
W = Path(__file__).parent


def prepare(path):
    with np.load(path) as data:
        D, seeds = data['D'], data['seeds']
        bg = np.where(~data['is_control'])[0]
        genes = data['genes'].astype(str).tolist()
        proteins = data['proteins'].astype(str).tolist()
    mu = D[:, bg].mean(1, keepdims=True)
    sd = D[:, bg].std(1, keepdims=True)
    Z = (D - mu) / (sd + 1e-9)
    return dict(D=D, Z=Z, Z_mean=Z.mean(2), seeds=seeds, bg=bg,
                genes=genes, gi={g: i for i, g in enumerate(genes)},
                pi={p: i for i, p in enumerate(proteins)})


def metrics(z, protein, members, direction, gene_null=False):
    genes = z['genes']
    present = [g for g in members if g in z['gi']]
    if not present or protein not in z['pi']:
        return dict(status='missing_gene_or_antigen_coverage')
    ix = np.array([z['gi'][g] for g in present])
    pi = z['pi'][protein]
    conf = np.where(z['seeds'] >= 47)[0]
    assert len(conf) == 5
    D = z['D']
    sign = 1. if direction == 'up' else -1.
    bg = z['bg']
    score = sign * (z['Z'][:, :, ix].mean(2) - z['Z_mean'])
    panel = (score - score[:, bg].mean(1, keepdims=True)) / (score[:, bg].std(1, keepdims=True) + 1e-12)
    mean = panel[conf].mean(0)
    other = bg[bg != pi]
    values = sign * D[:, pi][:, ix].mean(1)
    row = dict(status='evaluated', n_defined_genes=len(members), n_present_genes=len(present),
               missing_genes=';'.join(g for g in members if g not in genes),
               confirmation_positive_fraction=float((values[conf] > 0).mean()),
               confirmation_signed_mean_slope=float(values[conf].mean()),
               protein_percentile=float((mean[other] < mean[pi]).mean()),
               p_protein_decoy=float((1 + (mean[other] >= mean[pi]).sum()) / (1 + len(other))))
    row.update({'signed_slope_seed_' + str(int(s)): float(v) for s, v in zip(z['seeds'], values)})
    if gene_null:
        rng = np.random.RandomState(7)
        null_index = np.stack([rng.choice(len(genes), len(ix), False) for _ in range(5000)])
        response = sign * D[conf, pi].mean(0)
        null = response[null_index].mean(1)
        observed = response[ix].mean()
        row['p_gene_decoy'] = float((1 + (null >= observed).sum()) / (1 + len(null)))
    return row


while not (O / 'final_quality_checks.json').exists():
    status = json.loads((O / 'status.json').read_text())
    if status.get('stage') == 'failed':
        raise RuntimeError(status)
    time.sleep(30)

legacy = json.loads((W / 'reference/legacy_programs.json').read_text())
previous = json.loads((OLD / 'programs.json').read_text())
families = [('legacy_manuscript', legacy), ('previous_full_rescreen', previous)]
rows = []
for axis_file in sorted((O / 'curves').glob('*__slopes.npz')):
    axis = axis_file.name.replace('__slopes.npz', '')
    source = {label: prepare(root / 'curves' / axis_file.name) for label, root in [('previous_training', OLD), ('corrected_training', O)]}
    for family, programs in families:
        for combo, members in programs.items():
            protein, remainder = combo.split('__', 1)
            program_axis, direction = (remainder, 'up') if family == 'legacy_manuscript' else remainder.rsplit('__', 1)
            if program_axis != axis:
                continue
            row = dict(family=family, combo=combo, protein=protein, axis=axis,
                       response_direction=direction, program_membership='held fixed', n_genes=len(members))
            for label, z in source.items():
                result = metrics(z, protein, members, direction, gene_null=family == 'legacy_manuscript')
                row.update({label + '_' + k: v for k, v in result.items()})
            rows.append(row)
    print('FIXED MEMBERSHIP', axis, len(rows), flush=True)
table = pd.DataFrame(rows)
for label in ['previous_training', 'corrected_training']:
    mask = table.family.eq('legacy_manuscript')
    table.loc[mask, label + '_q_gene_within_legacy_audit'] = false_discovery_control(table.loc[mask, label + '_p_gene_decoy'].to_numpy())
    table[label + '_direction_pass'] = table[label + '_confirmation_positive_fraction'] >= .8
    table[label + '_protein_percentile_pass'] = table[label + '_protein_percentile'] >= .95
table.to_csv(O / 'fixed_membership_training_comparison.csv', index=False)
summary = []
for family, d in table.groupby('family'):
    summary.append(dict(family=family, n_programs=len(d),
        previous_direction_pass=int(d.previous_training_direction_pass.sum()),
        corrected_direction_pass=int(d.corrected_training_direction_pass.sum()),
        previous_protein_percentile_pass=int(d.previous_training_protein_percentile_pass.sum()),
        corrected_protein_percentile_pass=int(d.corrected_training_protein_percentile_pass.sum()),
        rescued_direction=int((~d.previous_training_direction_pass & d.corrected_training_direction_pass).sum()),
        lost_direction=int((d.previous_training_direction_pass & ~d.corrected_training_direction_pass).sum())))
(O / 'fixed_membership_training_comparison.json').write_text(json.dumps(dict(
    status='complete', purpose='Retrospective isolation of training effects at fixed gene membership; not a nomination branch.',
    direction_gate_only='>=4/5 confirmation seeds; not the full monotonicity or biological gate.',
    definition='Same protein-panel calibration and signed slopes as main randKO; positive direction retained for legacy manuscript programs.',
    families=summary), indent=2))
print(json.dumps(summary), flush=True)
