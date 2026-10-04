"""Confirm repaired density training and brain training remain correctly routed."""
from pathlib import Path
import json,torch,numpy as np
from ccwm import CCWM,CCWMConfig

torch.set_num_threads(2)
A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/optimization_20261003'
ck=torch.load(A/'results/rescreen_20261002/models/seed_42/model.pt',map_location='cpu',weights_only=False)
m=CCWM(CCWMConfig(**ck['cfg'])).eval();m.load_state_dict(ck['state_dict'])
cite=np.load(A/'interim/rescreen_20261002/citeseq.npz');bl=np.load(A/'interim/rescreen_20261002/blood.npz');br=np.load(A/'interim/rescreen_20261002/brain_train.npz')
ns=int(cite['n_states'])
def cond(z,ix,source):
    out=torch.zeros(len(ix),m.cfg.cond_dim)
    out[np.arange(len(ix)),z['cond2'][ix]]=1
    out[np.arange(len(ix)),2+z['state'][ix]]=1
    out[:,2+ns+source]=1
    return out
def norm(loss,module):
    gs=torch.autograd.grad(loss,list(module.parameters()),allow_unused=True,retain_graph=True)
    return sum(float(g.detach().square().sum()) for g in gs if g is not None)**.5
ix=np.random.RandomState(601).choice(np.where(~(cite['is_val']|cite['is_test']))[0],32,False)
x=torch.from_numpy(cite['X'][ix]);u=torch.from_numpy(cite['Y'][ix]);c=cond(cite,ix,0)
rows=[]
for mode in ['joint_legacy','partial_blood','blood_marginal']:
    m.cfg.dsm_mode=mode;torch.manual_seed(611)
    loss=m.block_b(x,u,c)['dsm']
    rows.append(dict(mode=mode,dsm=float(loss),blood_potential_grad=norm(loss,m.u_blood),brain_potential_grad=norm(loss,m.u_brain)))
assert rows[0]['brain_potential_grad']>0 and rows[1]['brain_potential_grad']==rows[2]['brain_potential_grad']==0
bix=np.where(~(bl['is_val']|bl['is_test'])&(bl['cond2']==1)&(bl['state']==list(bl['state_vocab']).index('classical_mono')))[0][:32]
rix=np.where(~br['is_val']&(br['cond2']==1)&(br['state']==list(br['state_vocab']).index('astro')))[0][:32]
torch.manual_seed(613)
loss=m.block_a(torch.from_numpy(bl['X'][bix]),cond(bl,bix,0),torch.from_numpy(br['X'][rix]),cond(br,rix,1))['sinkhorn']
brain_gradient=norm(loss,m.u_brain);assert brain_gradient>0
out=dict(status='passed',cite_density_gradients=rows,brain_distribution_loss_gradient=brain_gradient,
    scope='Engineering routing check using actual training batches. No candidate ranking or external validation used.')
(O/'real_batch_gradient_checks.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
