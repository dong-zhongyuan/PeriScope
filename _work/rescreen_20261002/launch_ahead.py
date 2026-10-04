"""Use available A6000 memory for two queued seeds; the main dispatcher adopts them."""
import subprocess,os,json,time
from pathlib import Path
import sys
W=Path(__file__).parent;A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/rescreen_20261002'
for seed,gpu in [(48,0),(49,1)]:
    out=O/f'models/seed_{seed}';out.mkdir(parents=True,exist_ok=True)
    if (out/'completed.json').exists() or (out/'adopted_training.json').exists():continue
    args=[sys.executable,'-u',str(W/'train_ccwm.py'),'--data-dir',str(A/'interim/rescreen_20261002'),'--out-dir',str(out),
          '--run-dir',str(out/'run'/time.strftime('%Y%m%d_%H%M%S')),'--run-id',f'rescreen-seed-{seed}',
          '--potential-beta','5','--sinkhorn-iters','20','--sample-latent-a','--seed',str(seed),'--epochs','200','--steps-per-epoch','200']
    env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
    with (out/'console.log').open('a') as f:p=subprocess.Popen(args,stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env,start_new_session=True)
    (out/'adopted_training.json').write_text(json.dumps(dict(pid=p.pid,reason='same frozen recipe, scheduled in available GPU memory')))
    print('LAUNCHED',seed,p.pid,flush=True)
