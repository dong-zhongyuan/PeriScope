"""Rerun only missing CellOT artifacts using the recovered original training function."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='1'
import argparse,copy,json,time,sys,subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
import torch
from sklearn.decomposition import PCA
sys.path.insert(0,'/public/home/mengxl/dzy/pd_product/_work/benchmark_shared_20261003')
import current_benchmark as b
sys.path.insert(0,str(Path(__file__).parent/'source'))
from cellot_losses import compute_loss_f,compute_loss_g
from cellot_icnns import ICNN
from run_benchmark import dump,mmd
DEVICE='cpu'
R=b.O/'cellot_recomputed_20261003'
MODEL_DIR=None
class Projection:
 def __init__(self,cache):self.mean_=cache['pca_mean'];self.components_=cache['pca_components']
 def transform(self,x):return (x-self.mean_)@self.components_.T
 def inverse_transform(self,x):return x@self.components_+self.mean_
def read(axis):
 k=dict(np.load(b.O/'cache'/f'{axis}.npz'))
 k.update(dict(np.load(R/'projection'/f'{axis}.npz')))
 return k
def assess(*args,**kwargs):return b.assess(*args,**kwargs)
def prepare():
 expected=[R/'projection'/f'{axis}.npz' for axis in b.AXES]
 if not any(p.exists() for p in expected):
  subprocess.run([sys.executable,str(Path(__file__).with_name('prepare_projection.py'))],check=True)
 assert all(p.exists() for p in expected), 'Incomplete PCA recovery; preserve and inspect before proceeding'
 assert (R/'projection_reconstruction.json').exists(), 'Verified PCA numerical configuration is required'

def cellot(seed,axis,max_steps=12000):
    k=read(axis);p=Projection(k);pred={};stats=[]
    for c in (0,1):
        ck=MODEL_DIR/f'cellot_{axis}_{c}.pt';log=MODEL_DIR/f'cellot_{axis}_{c}.json';torch.manual_seed(seed)
        init=lambda w:torch.nn.init.uniform_(w,0,.1)
        f=ICNN(64,[64]*4,kernel_init_fxn=init).to(DEVICE);g=ICNN(64,[64]*4,fnorm_penalty=1,kernel_init_fxn=init).to(DEVICE)
        x=p.transform(k[f'train_x{c}']);y=p.transform(k[f'train_y{c}']);scale=max(float(np.sqrt((np.mean(x*x)+np.mean(y*y))/2)),1e-6)
        if ck.exists() and log.exists() and json.loads(log.read_text()).get('complete'):
            state=torch.load(ck,map_location=DEVICE,weights_only=False);g.load_state_dict(state['state_dict']);stat=json.loads(log.read_text())
        else:
            sx=torch.as_tensor(x/scale,device=DEVICE);sy=torch.as_tensor(y/scale,device=DEVICE);vx=torch.as_tensor(p.transform(k[f'val_x{c}'])/scale,device=DEVICE)
            opts=[torch.optim.Adam(net.parameters(),lr=1e-4,betas=(.5,.9)) for net in [f,g]];rng=np.random.default_rng(seed);best=float('inf');bad=0;history=[];t0=time.time();beststate=None;reason='max_budget'
            for step in range(1,max_steps+1):
                target=sy[rng.integers(len(sy),size=256)]
                for _ in range(10):
                    source=sx[rng.integers(len(sx),size=256)].detach().requires_grad_(True);opts[1].zero_grad();loss=compute_loss_g(f,g,source).mean()+g.penalize_w();loss.backward();opts[1].step()
                source=sx[rng.integers(len(sx),size=256)].detach().requires_grad_(True);opts[0].zero_grad();lf=compute_loss_f(f,g,source,target).mean();lf.backward();opts[0].step();f.clamp_w()
                if step%1000==0 or step==max_steps:
                    q=g.transport(vx.detach().requires_grad_(True)).detach().cpu().numpy()*scale;val=p.inverse_transform(q);score=mmd(val,k[f'val_y{c}'],float(k['bw']))
                    if not np.isfinite(score):raise ValueError('nonfinite CellOT validation')
                    history.append({'step':step,'validation_mmd':score,'loss_f':float(lf.detach()),'loss_g':float(loss.detach()),'elapsed_s':time.time()-t0})
                    if score<best-1e-5:best=score;beststate=copy.deepcopy(g.state_dict());beststep=step;bad=0
                    else:bad+=1
                    dump(log,{'complete':False,'history':history,'best_step':beststep});print('CELLOT',axis,c,step,score,flush=True)
                    if step>=4000 and bad>=4:reason='validation_plateau';break
            g.load_state_dict(beststate);torch.save({'state_dict':beststate,'scale':scale},ck)
            stat={'complete':True,'history':history,'best_step':beststep,'stop_reason':reason,'budget':max_steps,'scale':scale};dump(log,stat)
        stats.append(stat);pred[c]={}
        for bc in (0,1):
            q=torch.as_tensor(p.transform(k[f'test_x{bc}'])/scale,device=DEVICE).requires_grad_(True)
            q=g.transport(q).detach().cpu().numpy()*scale;pred[c][bc]=p.inverse_transform(q).astype(np.float32)
    assess('cellot',pred,k,seed,stats)


def worker(seed,axis):
 env=os.environ.copy();env['BENCHMARK_SEED']=str(seed);env['BENCHMARK_DEVICE']='cpu'
 dest=R/'logs'/f'{seed}__{axis}.log';dest.parent.mkdir(parents=True,exist_ok=True)
 with dest.open('a') as log:
  subprocess.run([sys.executable,str(Path(__file__).resolve()),'--seed',str(seed),'--axis',axis],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 return (seed,axis)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int);ap.add_argument('--axis',choices=b.AXES);ap.add_argument('--workers',type=int,default=16);args=ap.parse_args()
 torch.set_num_threads(1)
 if args.seed is not None:
  MODEL_DIR=R/'models'/f'seed_{args.seed}';MODEL_DIR.mkdir(parents=True,exist_ok=True)
  b.RESULT_DIR=R/'results'/f'seed_{args.seed}';b.RESULT_DIR.mkdir(parents=True,exist_ok=True)
  cellot(args.seed,args.axis)
 else:
  prepare();done=[]
  with ThreadPoolExecutor(max_workers=args.workers) as pool:
   pairs=[(s,a) for s in [42,43,44] for a in b.AXES if s!=42 or a.startswith('pDC__')]
   jobs=[pool.submit(worker,s,a) for s,a in pairs]
   for job in as_completed(jobs):
    done.append(job.result());dump(R/'run_status.json',dict(complete=False,completed=done,total=len(pairs)));print('COMPLETE',done[-1],len(done),'/',len(pairs),flush=True)
  dump(R/'run_status.json',dict(complete=True,completed=done,total=len(pairs)))
