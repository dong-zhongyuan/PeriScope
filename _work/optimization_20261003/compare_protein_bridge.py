"""Compare measured held-out protein prediction before and after DSM repair."""
from pathlib import Path
import hashlib
import json
import pandas as pd

A = Path('/public/home/mengxl/dzy/pd_product_assets/results')
O = A / 'optimization_20261003'
registry = json.loads((O / 'candidate_registry.json').read_text())
is_control = {r['target']: r['is_control'] for r in registry}
frames, sources = [], {}
for label, name in [('previous', 'rescreen_20261002'), ('corrected', 'optimization_20261003')]:
    root = A / name
    files = sorted((root / 'evaluation').glob('protein_accuracy_seed*.csv'))
    assert len(files) == 10
    a = pd.concat([pd.read_csv(p) for p in files])
    sources.update({str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    # Each seed table has one row per target, fine subtype and donor.
    # First average fine subtypes and models within donor, then weight donors equally.
    a = a.groupby(['target', 'blood_state', 'donor'], as_index=False)[['rho', 'mse']].mean()
    a['positive'] = a.rho > 0
    b = a.groupby(['target', 'blood_state']).agg(
        protein_rho=('rho', 'mean'), protein_mse=('mse', 'mean'),
        n_positive_donors=('positive', 'sum')).reset_index()
    b['prediction_pass'] = b.protein_rho.ge(.2) & b.n_positive_donors.ge(2)
    b['is_control'] = b.target.map(is_control)
    assert b.is_control.notna().all()
    b['branch'] = label
    frames.append(b)
d = pd.concat(frames, ignore_index=True)
d.to_csv(O / 'protein_bridge_training_comparison.csv', index=False)
(O / 'protein_bridge_comparison_definition.json').write_text(json.dumps(dict(
    cohorts='P6, P7, P8 held-out CITE-seq donors',
    summary='Equal fine-subtype/model weighting within donor and equal donor weighting, identical to nomination code.',
    control_rows='Retained with is_control flag; excluded from biological target counts.',
    source_sha256=sources), indent=2))
print(d.loc[~d.is_control].groupby('branch').prediction_pass.sum().to_dict())
