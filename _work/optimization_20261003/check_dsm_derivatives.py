"""Check partial-gradient semantics and finite differences, without target outcomes."""
from pathlib import Path
import json,torch
from ccwm import CCWM,CCWMConfig

torch.set_num_threads(2);torch.manual_seed(319)
A=Path('/public/home/mengxl/dzy/pd_product_assets')
O=A/'results/optimization_20261003';O.mkdir(exist_ok=True,parents=True)
checkpoint=torch.load(A/'results/rescreen_20261002/models/seed_42/model.pt',map_location='cpu',weights_only=False)
m=CCWM(CCWMConfig(**checkpoint['cfg'])).double().eval();m.load_state_dict(checkpoint['state_dict'])
z=torch.randn(8,m.cfg.z_dim,dtype=torch.float64,requires_grad=True)
u=torch.randn(8,m.cfg.n_proteins,dtype=torch.float64);u-=u.mean(-1,keepdim=True)
c=torch.zeros(8,m.cfg.cond_dim,dtype=torch.float64);c[:,0]=1;c[:,-2]=1
b=z.detach().clone().requires_grad_(True)
e=m.potential(z,b,u,c);gb,gB=torch.autograd.grad(e.sum(),(z,b),create_graph=True)
gd=torch.autograd.grad(m.potential(z,z,u,c).sum(),z,create_graph=True)[0]
gp=torch.autograd.grad(m.potential(z,z.detach(),u,c).sum(),z,create_graph=True)[0]
v=torch.randn_like(z);v/=v.norm();eps=1e-5
with torch.no_grad():fd=float((m.potential(z+eps*v,b,u,c).sum()-m.potential(z-eps*v,b,u,c).sum())/(2*eps))
analytic=float((gp*v).sum());weights=list(m.u_brain.parameters())
old=torch.autograd.grad(gd.square().mean(),weights,allow_unused=True,retain_graph=True)
fixed=torch.autograd.grad(gp.square().mean(),weights,allow_unused=True)
norm=lambda gs:sum(float(g.square().sum()) for g in gs if g is not None)**.5
r=dict(status='passed',probe='engineering derivative identity on seed42 weights; random latent probe is not biological evidence',
    diagonal_gradient_is_sum_of_blood_and_brain_partials=bool(torch.allclose(gd,gb+gB,atol=1e-9,rtol=1e-9)),
    corrected_gradient_equals_blood_partial=bool(torch.allclose(gp,gb,atol=1e-9,rtol=1e-9)),
    partial_finite_difference_absolute_error=abs(fd-analytic),
    legacy_brain_potential_parameter_gradient_norm=norm(old),corrected_brain_potential_parameter_gradient_norm=norm(fixed),
    brain_partial_to_blood_partial_norm=float(gB.norm()/gb.norm()))
assert r['diagonal_gradient_is_sum_of_blood_and_brain_partials'] and r['corrected_gradient_equals_blood_partial']
assert r['partial_finite_difference_absolute_error']<1e-6 and norm(old)>0 and norm(fixed)==0
(O/'dsm_derivative_checks.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
