"""Accuracy check for OT iteration count on real repaired-input model readouts."""
from pathlib import Path
import numpy as np,torch,json
from ccwm import CCWM,CCWMConfig,sinkhorn_divergence
from training_utils import predict
A=Path('/public/home/mengxl/dzy/pd_product_assets');I=A/'interim/rescreen_20261002';O=A/'results/rescreen_20261002'
torch.set_num_threads(2)
ck=torch.load(O/'pilot_beta5/model.pt',map_location='cpu',weights_only=False)
m=CCWM(CCWMConfig(**ck['cfg']));m.load_state_dict(ck['state_dict']);m.eval().cuda()
b=np.load(I/'blood.npz');r=np.load(I/'brain_train.npz');Xb=b['X'];Xr=r['X'];states=r['state_vocab'].astype(str).tolist();rng=np.random.default_rng(76);records=[]
for state in ['astro','microglia_mhc2','microglia_homeostatic']:
 for cond in [0,1]:
  ii=np.where((r['state']==states.index(state))&(r['cond2']==cond)&~r['is_val'])[0]
  bi=np.where((b['state']==states.index('classical_mono'))&(b['cond2']==cond)&~(b['is_val']|b['is_test']))[0]
  if len(ii)<25:continue
  ri=rng.choice(ii,128,replace=len(ii)<128);bb=rng.choice(bi,128,replace=False)
  y=torch.from_numpy(Xr[ri]).cuda();c=torch.zeros(128,m.cfg.cond_dim,device='cuda');c[:,cond]=1;c[:,2+states.index(state)]=1;c[:,-1]=1
  x=predict(m,torch.from_numpy(Xb[bb]).cuda(),c).detach().requires_grad_(True)
  reference=sinkhorn_divergence(x,y,eps=m.cfg.sinkhorn_eps,iters=100,scale=m.cfg.ot_scale)
  rg=torch.autograd.grad(reference,x)[0]
  for it in [10,20,30,50]:
   v=sinkhorn_divergence(x,y,eps=m.cfg.sinkhorn_eps,iters=it,scale=m.cfg.ot_scale);g=torch.autograd.grad(v,x)[0]
   records.append(dict(state=state,condition=cond,iterations=it,cost_abs_error=abs(float(v-reference)),gradient_relative_error=float((g-rg).norm()/rg.norm().clamp_min(1e-12))))
(O/'sinkhorn_precision.json').write_text(json.dumps(records,indent=2));print(json.dumps(records),flush=True)
