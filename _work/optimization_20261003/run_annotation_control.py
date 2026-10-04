"""Isolate the external-cohort identifier repair on the prior frozen programs."""
from pathlib import Path
import json
import os
import subprocess
import sys
import pandas as pd
import numpy as np

W = Path(__file__).parent
A = Path('/public/home/mengxl/dzy/pd_product_assets')
O = A / 'results/optimization_20261003'
OLD = A / 'results/rescreen_20261002'
C = O / 'annotation_control'
S = W / 'annotation_control'
C.mkdir(exist_ok=True)
S.mkdir(exist_ok=True)
for name in ['programs.json', 'program_folds.json', 'program_gene_rankings.json']:
    p = C / name
    if not p.exists():
        p.symlink_to(OLD / name)
for name in ['downstream_input_checks.json']:
    p = C / name
    if not p.exists():
        p.symlink_to(O / name)
for name in ['independent_validation.py', 'evidence_cache.py']:
    text = (W / name).read_text().replace('results/optimization_20261003', 'results/optimization_20261003/annotation_control')
    (S / name).write_text(text)
for name in ['input_utils.py', 'reference']:
    p = S / name
    if not p.exists():
        p.symlink_to(W / name, target_is_directory=name == 'reference')
env = os.environ | dict(PYTHONPATH=str(W) + ':/public/home/mengxl/dzy/pd_product/src',
                        OPENBLAS_NUM_THREADS='2', OMP_NUM_THREADS='2')
with (C / 'independent_validation.log').open('w') as log:
    subprocess.run([sys.executable, '-u', str(S / 'independent_validation.py')], env=env,
                   stdout=log, stderr=subprocess.STDOUT, check=True)
e = pd.read_csv(OLD / 'all_program_evidence.csv').set_index('combo')
old_tier = e.tier.copy()
original = pd.read_csv(OLD / 'independent_validation_stats.csv').set_index('combo')
fixed = pd.read_csv(C / 'independent_validation_stats.csv').set_index('combo')
for column, external in [('independent_effect', 'effect_PD_minus_control'),
                         ('independent_p', 'p_exact_two_sided'), ('independent_q', 'q_BH_all_programs')]:
    e[column] = fixed[external].reindex(e.index)
direction = e.disease_direction.map({'up': 1., 'down': -1., 'unresolved': 0.})
e['independent_direction_pass'] = direction * e.independent_effect > 0
e['tier'] = 'D'
basic = e.technical_pass & e.reproducible_pass & (e.q_gene_decoy <= .05) & e.brain_pass
e.loc[basic, 'tier'] = 'C'
eligible = e.technical_pass & e.reproducible_pass & e.specificity_pass & e.brain_pass & e.independent_direction_pass
e.loc[eligible, 'tier'] = 'B'
e.loc[eligible, 'tier'] = 'A'
labels = [('gene_or_complex_identity', 'mapping_pass'), ('held_out_protein_prediction', 'protein_prediction_pass'),
          ('measured_dose_range', 'measured_dose_pass'), ('confirmation_seeds', 'reproducible_pass'),
          ('specificity', 'specificity_pass'), ('brain_disease_support', 'brain_pass'),
          ('independent_direction', 'independent_direction_pass')]
e['unmet_nomination_criteria'] = [';'.join(name for name, key in labels if not row[key])
                                 for row in e.to_dict('records')]
e.reset_index().to_csv(C / 'all_program_evidence.csv', index=False)
comparison = original.add_prefix('previous_').join(fixed.add_prefix('repaired_'), how='outer')
comparison.to_csv(C / 'external_annotation_effect_comparison.csv')
changed = original.effect_PD_minus_control.reindex(fixed.index) - fixed.effect_PD_minus_control
summary = dict(status='complete', source_programs=str(OLD / 'programs.json'),
               scope='Only external gene identifiers changed; previous trained models, programs, disease evidence and nomination rules retained.',
               previous_n_evaluable=len(original), repaired_n_evaluable=len(fixed),
               n_program_effect_changes=int((changed.abs() > 1e-10).sum()),
               previous_tiers=old_tier.value_counts().to_dict(), repaired_tiers=e.tier.value_counts().to_dict(),
               n_tier_changes=int((old_tier != e.tier).sum()))
(C / 'completed.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2), flush=True)
