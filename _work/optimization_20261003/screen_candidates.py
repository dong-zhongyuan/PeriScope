"""All-antigen paired dose screening of the repaired CCWM checkpoints. Three donor-balanced PD blood draws, measured subtype-specific doses, CLR-preserving interventions."""
import os
import sys
import json
import glob
import time
import argparse

import numpy as np
from scipy.stats import combine_pvalues, wilcoxon

A = '/public/home/mengxl/dzy/pd_product_assets'
CACHE = A + '/results/optimization_20261003/curves'
QR = np.arange(5, dtype=float)
QUANTILES = [0.05, 0.25, 0.50, 0.75, 0.95]
DRAWS = [42, 43, 44]


def slope(C):
    x = QR - QR.mean()
    return np.tensordot(x, C, axes=(0, 2)) / (x * x).sum()


def bh(ps):
    ps = np.asarray(ps, float)
    m = len(ps)
    o = np.argsort(ps)
    q = np.empty(m)
    prev = 1.0
    for r_, i_ in zip(range(m, 0, -1), o[::-1]):
        prev = min(prev, ps[i_] * m / r_)
        q[i_] = prev
    return q


def emp_p(vr, vn):
    return float((1 + (vn >= vr).sum()) / (1 + len(vn)))


def build_cache(seed,output_root=None,limit_axes=None):
    from pathlib import Path
    import torch
    from ccwm import CCWM,CCWMConfig
    from training_utils import BalancedSampler
    torch.set_num_threads(2)
    I=Path(A)/'interim/rescreen_20261002';O=Path(output_root) if output_root else Path(A)/'results/optimization_20261003'
    out=O/'curves';out.mkdir(exist_ok=True)
    from model_registry import checkpoint, curve_provenance, valid_curve, stamp_curve
    provenance=curve_provenance(O,seed)
    ck=torch.load(checkpoint(O,seed),map_location='cpu',weights_only=False)
    m=CCWM(CCWMConfig(**ck['cfg']));m.load_state_dict(ck['state_dict']);m.eval().to('cuda')
    # Gradient flow is required for inference, higher-order graphs are not.
    original_grad=torch.autograd.grad
    def first_order(*a,**k):
        k['create_graph']=False
        return original_grad(*a,**k)
    torch.autograd.grad=first_order
    blood=np.load(I/'blood.npz');X=blood['X'];states=blood['state_vocab'].astype(str).tolist()
    registry=json.loads((O/'candidate_registry.json').read_text())
    dose=json.loads((O/'dose_definitions.json').read_text())
    training=json.loads((O/f'models/seed_{seed}/training_inputs.json').read_text())
    strata=training['strata'];pairs=sorted({(b,br) for c,b,br,*_ in strata if c==1})
    rng=np.random.RandomState(717);n=128;batch_targets=4
    if limit_axes is not None:pairs=pairs[:limit_axes]
    for bs,br in pairs:
        axis=states[bs]+'__x__'+states[br];dest=out/f'{axis}__seed{seed}.npz'
        if valid_curve(dest,provenance):continue
        if states[bs] not in dose['doses']:continue
        values=np.array(dose['doses'][states[bs]])
        ids=np.where(~(blood['is_test']|blood['is_val']) & (blood['state']==bs)&(blood['cond2']==1))[0]
        sampler=BalancedSampler(blood['donor'],ids)
        C=np.zeros((len(registry),5,X.shape[1]),np.float32)
        baseline=np.zeros(X.shape[1],np.float32); max_zero_sum=0.; diagnostics=[]
        for dr in DRAWS:
            sel=sampler.draw(n,np.random.RandomState(dr));xb=torch.from_numpy(X[sel]).cuda()
            c=torch.zeros(n,ck['cfg']['cond_dim'],device='cuda');c[:,1]=1;c[:,2+br]=1;c[:,-1]=1
            with torch.no_grad():
                mu,lv=m.enc_blood(xb);ub=m.protein_head(mu)
                if m.cfg.sample_latent_a:
                    generator=torch.Generator(device='cuda').manual_seed(dr)
                    mu=mu+torch.exp(.5*lv)*torch.randn(mu.shape,device='cuda',generator=generator)
            _,zb,tr=m.solve_equilibrium(mu,mu.clone(),ub,c,True)
            with torch.no_grad():baseline+=m.dec_brain(zb.detach()).mean(0).cpu().numpy()/len(DRAWS)
            diagnostics.append(dict(draw=dr,trajectory=tr,n_source_donors=len(np.unique(blood['donor'][sel]))))
            for start in range(0,len(registry),batch_targets):
                records=registry[start:start+batch_targets];controls=[]
                for r in records:
                    ii=r['channels'];other=[j for j in range(ub.shape[1]) if j not in ii]
                    for q in range(5):
                        uc=ub.clone();uc[:,ii]=torch.as_tensor(values[ii,q],device='cuda',dtype=ub.dtype)
                        delta=uc.sum(-1,keepdim=True)/len(other)
                        uc[:,other]-=delta
                        max_zero_sum=max(max_zero_sum,float(uc.sum(-1).abs().max()))
                        controls.append(uc)
                uu=torch.cat(controls);mm=mu.repeat(len(controls),1);cc=c.repeat(len(controls),1)
                _,zs=m.solve_equilibrium(mm,mm.clone(),uu,cc)
                with torch.no_grad():
                    # Decoder batches keep wide RNA outputs bounded in memory.
                    ys=[]
                    for k in range(len(controls)):
                        ys.append(m.dec_brain(zs[k*n:(k+1)*n].detach()).mean(0).cpu().numpy())
                    C[start:start+len(records)]+=np.stack(ys).reshape(len(records),5,-1)/len(DRAWS)
        assert np.isfinite(C).all() and max_zero_sum<.001
        # Mean clone CLR is the scalar dose coordinate; repeated clones move together.
        u=np.array([values[r['channels']].mean(0) for r in registry])
        np.savez_compressed(str(dest)+'.tmp.npz',C=C,baseline=baseline,proteins=np.array([r['target'] for r in registry]),
            genes=blood['genes'],u=u,seed=seed,quantiles=QUANTILES,is_control=np.array([r['is_control'] for r in registry]),
            diagnostics=json.dumps(diagnostics),clr_zero_sum_error=max_zero_sum)
        Path(str(dest)+'.tmp.npz').replace(dest)
        stamp_curve(dest,provenance)
        print(axis,'seed',seed,'all antigens',len(registry),'saved',flush=True)
    torch.autograd.grad=original_grad


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,required=True)
    ap.add_argument('--output-root');ap.add_argument('--limit-axes',type=int)
    args=ap.parse_args();build_cache(args.seed,args.output_root,args.limit_axes)


if __name__=='__main__':main()
