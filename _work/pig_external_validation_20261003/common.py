from pathlib import Path
import hashlib,json,sys
import numpy as np
import torch

A=Path('/public/home/mengxl/dzy/pd_product_assets')
O=A/'results/pig_external_validation_20261003'
I=A/'interim/rescreen_20261002'
B=A/'results/benchmark_20261002'
W=Path('/public/home/mengxl/dzy/pd_product/_work/pig_external_validation_20261003')
sys.path.insert(0,str(W/'source'))
BS=['cDC','classical_mono','nonclassical_mono'];RS=['astro','microglia_mhc2']
AXES=[f'{b}__{r}' for b in BS for r in RS]

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for x in iter(lambda:f.read(4<<20),b''):h.update(x)
    return h.hexdigest()

def dump(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(value,indent=2,allow_nan=False));tmp.replace(path)

def normalize(v,mask):
    v=np.asarray(v,float);x=v[...,mask]
    return x/np.maximum(x.sum(axis=-1,keepdims=True),1e-12)*10000

def logfactor(target,reference,mask,alpha,cap):
    v=np.zeros(len(mask));v[mask]=np.clip(np.log((normalize(target,mask)+alpha)/(normalize(reference,mask)+alpha)),-cap,cap)
    return v

def adapt(linear,logfc,gain):
    old=np.asarray(linear,float);new=old*np.exp(gain*logfc)
    # Includes unmodelled genes as an unchanged remainder of the original CP10k library.
    total=np.maximum(10000+(new-old).sum(axis=-1,keepdims=True),1)
    return (new*(10000/total)).astype(np.float32)

def load_model(seed,device):
    from ccwm import CCWM,CCWMConfig
    f=A/f'results/optimization_20261003/models/seed_{seed}/model.pt'
    ck=torch.load(f,map_location=device,weights_only=False)
    m=CCWM(CCWMConfig(**ck['cfg'])).to(device);m.load_state_dict(ck['state_dict']);m.eval()
    for par in m.parameters():par.requires_grad_(False)
    return m

def condition(n,c,state,ns,device):
    v=torch.zeros((n,ns+4),device=device);v[:,c]=1;v[:,2+state]=1;v[:,-1]=1;return v
