"""Shared-representation benchmark. Frozen PeriScope encoding and decoding modules are reused.
CCWM prediction calls use the copied current training implementation. See protocol.txt.
"""
import argparse, hashlib, json, sys, time, copy, csv
from pathlib import Path
import numpy as np
import torch
from scipy.spatial.distance import cdist, pdist
from scipy.stats import spearmanr

ROOT=Path('/public/home/mengxl/dzy/pd_product')
A=ROOT.parent/'pd_product_assets'
W=ROOT/'_work/benchmark_shared_20261003'
O=A/'results/benchmark_20261002'
O.mkdir(parents=True,exist_ok=True)
sys.path[:0]=[str(W/'source')]
from ccwm import CCWM, CCWMConfig

def dump(path,obj):
 path=Path(path); tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False));tmp.replace(path)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<22),b''):h.update(b)
 return h.hexdigest()
def balanced(d,st,c,n,rng):
 ids=np.flatnonzero((d['state']==st)&(d['cond2']==c)); donors=np.unique(d['donor'][ids])
 if not len(donors):raise ValueError('empty stratum')
 chosen=np.resize(rng.permutation(donors),n)
 return np.array([rng.choice(ids[d['donor'][ids]==don]) for don in chosen])
def donor_profiles(d,st,c):
 ids=np.flatnonzero((d['state']==st)&(d['cond2']==c))
 return np.stack([d['X'][ids[d['donor'][ids]==don]].mean(0) for don in np.unique(d['donor'][ids])])

def mmd(x,y,bw):
 return float(np.exp(-cdist(x,x,'sqeuclidean')/(2*bw)).mean()+np.exp(-cdist(y,y,'sqeuclidean')/(2*bw)).mean()-2*np.exp(-cdist(x,y,'sqeuclidean')/(2*bw)).mean())
def measures(x,y,bw):
 g=x.shape[1]; xy=cdist(x,y);xx=cdist(x,x);yy=cdist(y,y)
 return dict(mean_mse=float(np.mean((x.mean(0)-y.mean(0))**2)),energy_per_sqrt_gene=float((2*xy.mean()-xx.mean()-yy.mean())/np.sqrt(g)),mmd=mmd(x,y,bw),variance_ratio=float(x.var(0).mean()/max(y.var(0).mean(),1e-12)))

def checkpoint_path(seed):
 overrides=O/'checkpoint_overrides.json'
 if overrides.exists():
  paths=json.loads(overrides.read_text()).get('checkpoints',{})
  if str(seed) in paths:return Path(paths[str(seed)])
 return A/f'results/optimization_20261003/models/seed_{seed}/model.pt'

def load_model(seed,device):
 ck=checkpoint_path(seed)
 d=torch.load(ck,map_location=device,weights_only=False)
 model=CCWM(CCWMConfig(**d['cfg'])).to(device);model.load_state_dict(d['state_dict']);model.eval();return model




if __name__=='__main__':
 from current_benchmark import main
 main()
