import json
from pathlib import Path
import numpy as np
from training_utils import BalancedSampler,mmd
A=Path('/public/home/mengxl/dzy/pd_product_assets');I=A/'interim/rescreen_20261002';O=A/'results/rescreen_20261002'
b=np.load(I/'brain_train.npz');X=b['X'];states=b['state_vocab'].astype(str).tolist();scale=json.loads((O/'pilot_beta5/training_inputs.json').read_text())['ot_scale'];rows=[];rng=np.random.RandomState(8801)
for state in ['astro','microglia_homeostatic','microglia_mhc2']:
 for cond in [0,1]:
  fit=np.where(~b['is_val']&(b['state']==states.index(state))&(b['cond2']==cond))[0];val=np.where(b['is_val']&(b['state']==states.index(state))&(b['cond2']==cond))[0]
  if len(fit)<25 or len(val)<25:continue
  tr=BalancedSampler(b['donor'],fit).draw(96,rng);va=BalancedSampler(b['donor'],val).draw(96,rng)
  rows.append(dict(state=state,condition=cond,label_marginal_mmd=mmd(X[tr],X[va],scale),label_mean_mmd=mmd(np.repeat(X[tr].mean(0)[None,:],96,axis=0),X[va],scale),mean_mse=float(np.mean((X[tr].mean(0)-X[va].mean(0))**2))))
(O/'validation_baselines.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows),flush=True)
