"""Numerical regression checks for the substantive training repairs."""
import json
from pathlib import Path
import torch
from ccwm import CCWM,CCWMConfig,sinkhorn_divergence

torch.set_num_threads(2);torch.manual_seed(811)
cfg=CCWMConfig(n_genes=32,n_proteins=8,cond_dim=8,z_dim=8,hidden=32,depth=2,potential_beta=5,ot_scale=32)
m=CCWM(cfg)
x=torch.randn(24,32);y=torch.randn(24,32);c=torch.randn(24,8);u=torch.randn(24,8)
r=m.block_a(x,c,y,c);r['sinkhorn'].backward()
pg=float(m.u_brain.net[0].weight.grad[:,8:16].abs().max())
assert pg>0,pg
m.zero_grad();m.block_b(x,u,c)['dsm'].backward()
bg=float(m.u_blood.net[0].weight.grad[:,8:16].abs().max());assert bg>0
m.eval()
assert torch.allclose(m.protein_head(torch.randn(5,8)).sum(-1),torch.zeros(5),atol=1e-6)
with torch.no_grad():z,_=m.enc_blood(x)
_,z1,tr=m.solve_equilibrium(z,z.clone(),u,c,True)
_,z2=m.solve_equilibrium(z,z.clone(),u,c)
assert torch.equal(z1,z2)
_,single=m.solve_equilibrium(z[:1],z[:1].clone(),u[:1],c[:1])
assert torch.allclose(z1[:1],single,atol=1e-6), (z1[:1]-single).abs().max()
assert all(b<=a+1e-6 for a,b,t in tr),tr
same=float(sinkhorn_divergence(x,x,scale=32))
shift=float(sinkhorn_divergence(x,x+1,scale=32))
assert abs(same)<1e-6 and shift>same+.01,(same,shift)
truth=torch.tensor([[0.,1.,0.,2.]])
m.cfg.zero_weight=20
raw=m.zm_mse(torch.zeros_like(truth),truth,truth)
assert abs(float(raw)-(1+4)*20/4)<1e-6
out=dict(block_a_protein_gradient=pg,block_b_protein_gradient=bg,deterministic_solver=True,
         energy_nonincreasing=True,self_sinkhorn=same,shifted_sinkhorn=shift,raw_zero_mask=True)
p=Path('/public/home/mengxl/dzy/pd_product_assets/results/rescreen_20261002');p.mkdir(parents=True,exist_ok=True)
(p/'engineering_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
