"""Short foreground-detached runs for reproduction and throughput checks.
usage: run_check.py <group> name:gpu:seed:epochs[:extra,args] ..."""
import os, subprocess, sys
from pathlib import Path
O=Path('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/random_pairing_20261119/training')
T=O/'code'; PY='/public/home/mengxl/dzy/envs/pd_bbm/bin/python'; DATA='/public/home/mengxl/dzy/pd_product_assets/interim/rescreen_20261002'
g=O/'checks'/sys.argv[1]; g.mkdir(parents=True,exist_ok=True)
pids={}
for spec in sys.argv[2:]:
    name,gpu,seed,ep,*rest=spec.split(':'); extra=rest[0].split(',') if rest and rest[0] else []
    out=g/name; out.mkdir(parents=True,exist_ok=False)
    cmd=['nice','-n','10',PY,'-u',str(T/'train_ccwm.py'),'--data-dir',DATA,'--out-dir',str(out),'--run-dir',str(out/'run'),'--potential-beta','5','--sample-latent-a','--dsm-mode','blood_marginal','--sinkhorn-iters','20','--seed',seed,'--epochs',ep,'--steps-per-epoch','200','--save-epochs','']+extra
    env=dict(os.environ,CUDA_VISIBLE_DEVICES=gpu,OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',PYTHONDONTWRITEBYTECODE='1')
    p=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=open(out/'console.log','w'),stderr=subprocess.STDOUT,env=env,cwd=T,start_new_session=True)
    pids[name]=p.pid
print(pids)
(g/'pids.txt').write_text('\n'.join(f'{k} {v}' for k,v in pids.items())+'\n')
