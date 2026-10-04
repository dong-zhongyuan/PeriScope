"""Donor-balanced sampling and validation; never reads locked test arrays."""
import numpy as np
import torch
from scipy.spatial.distance import cdist

class BalancedSampler:
    def __init__(self, labels, allowed=None):
        labels=np.asarray(labels)
        if allowed is None: allowed=np.arange(len(labels))
        self.groups=[allowed[labels[allowed]==s] for s in np.unique(labels[allowed])]
        self.groups=[g for g in self.groups if len(g)>=25]
        if not self.groups: raise ValueError('empty sampler')
    def draw(self,n,rng):
        group=rng.randint(len(self.groups),size=n)
        return np.array([self.groups[g][rng.randint(len(self.groups[g]))] for g in group])

def mmd(x,y,scale):
    xx=cdist(x,x,'sqeuclidean')/scale; yy=cdist(y,y,'sqeuclidean')/scale;xy=cdist(x,y,'sqeuclidean')/scale
    return float(np.mean([np.exp(-xx/t).mean()+np.exp(-yy/t).mean()-2*np.exp(-xy/t).mean() for t in [.25,1.,4.]]))

def predict(model,x,c):
    with torch.no_grad():
        mu,lv=model.enc_blood(x);u=model.protein_head(mu)
        source=mu
        if model.cfg.sample_latent_a:
            generator=torch.Generator(device=x.device).manual_seed(7701)
            source=mu+torch.exp(.5*lv)*torch.randn(mu.shape,device=x.device,generator=generator)
    with torch.enable_grad(): _,z=model.solve_equilibrium(source,source.clone(),u,c)
    with torch.no_grad(): return model.dec_brain(z.detach())

def validation(model,brain,blood,cite,X_brain,X_blood,X_cite,Y_cite,c_brain,brain_ids,blood_ids,prot_scale):
    model.eval();rng=np.random.RandomState(8801)
    # A fixed, subtype-balanced validation draw, chosen without examining results.
    labels=np.array([str(d)+'|'+str(t) for d,t in zip(cite['donor'],cite['subtype'])])
    ci=BalancedSampler(labels,np.where(cite['is_val'])[0]).draw(2048,rng)
    with torch.no_grad():
        mu,_=model.enc_blood(X_cite[ci]);ph=model.protein_head(mu)
        pmse=float((ph-Y_cite[ci]).square().mean())
    rows=[]
    for c in (0,1):
      for br in brain_ids:
        Br=np.where(brain['is_val'] & (brain['cond2']==c)&(brain['state']==br))[0]
        if len(Br)<25: continue
        ri=BalancedSampler(brain['donor'],Br).draw(96,rng)
        y=X_brain[ri].cpu().numpy()
        for bs in blood_ids:
          Bi=np.where(blood['is_val']&(blood['cond2']==c)&(blood['state']==bs))[0]
          if len(Bi)<25: continue
          bi=BalancedSampler(blood['donor'],Bi).draw(96,rng)
          p=predict(model,X_blood[bi],c_brain[ri]).cpu().numpy()
          rows.append({'cond':c,'blood':int(bs),'brain':int(br),'mmd':mmd(p,y,model.cfg.ot_scale),
                       'mean_mse':float(np.mean((p.mean(0)-y.mean(0))**2)),
                       'variance_ratio':float(p.var(0).mean()/max(y.var(0).mean(),1e-9))})
    if not rows: raise ValueError('no supported validation stratum')
    bscore=float(np.mean([r['mmd'] for r in rows]))
    result={'selection_score':pmse/prot_scale+bscore,'protein_mse':pmse,'brain_mmd':bscore,'strata':rows}
    model.train()
    return result
