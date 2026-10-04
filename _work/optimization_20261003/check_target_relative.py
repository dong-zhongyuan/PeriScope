"""Validate vectorized leave-one-antigen-out contrasts against explicit decoys."""
from pathlib import Path
import json,numpy as np
A=Path('/public/home/mengxl/dzy/pd_product_assets')
p=next((A/'results/rescreen_20261002/curves').glob('*__slopes.npz'))
z=np.load(p);D=z['D'][:5].astype(float);bg=np.where(~z['is_control'])[0]
total=D[:,bg].sum(1,keepdims=True);squares=np.square(D[:,bg]).sum(1,keepdims=True);n=len(bg)-1
mu=(total-D)/n;var=np.maximum((squares-np.square(D))/n-mu**2,0);sd=np.sqrt(var)
relative=np.divide(D-mu,sd,out=np.zeros_like(mu),where=sd>1e-9)
rng=np.random.RandomState(701);errors=[]
for _ in range(100):
    si=rng.randint(5);pi=int(rng.choice(bg));gi=rng.randint(D.shape[-1]);pool=bg[bg!=pi]
    values=D[si,pool,gi];denom=values.std()
    brute=(D[si,pi,gi]-values.mean())/denom if denom>1e-9 else 0.
    errors.append(abs(brute-relative[si,pi,gi]))
assert max(errors)<1e-8 and np.isfinite(relative).all()
report=dict(status='passed',max_absolute_error=max(errors),n_random_coordinates=100,
    source=str(p),scope='Actual previous-run slopes used solely for arithmetic validation, not parameter or target selection',
    n_biological_antigens=len(bg),n_decoys_per_antigen=n,discovery_seed_indices=list(range(5)),
    controls_excluded_from_biological_pool=True)
(A/'results/optimization_20261003/target_relative_arithmetic_checks.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
