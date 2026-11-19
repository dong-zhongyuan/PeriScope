import json, os, subprocess, sys, time
from pathlib import Path
O=Path('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/random_pairing_20261119/training')
T=O/'code'; PY='/public/home/mengxl/dzy/envs/pd_bbm/bin/python'; DATA='/public/home/mengxl/dzy/pd_product_assets/interim/rescreen_20261002'
CFG=[('p3_t2',3,2),('p4_t1',4,1),('p6_t1',6,1)]
SPEC=[('full47',47,[]),('rp42',42,['--pairing','random']),('co44',44,['--coupling','off']),('ds42',42,['--pairing','disease_shuffled'])]
res=[]
for name,per,th in CFG:
    g=O/'checks'/f'bench_{name}'; g.mkdir(parents=True)
    ps=[]; k=0
    for i in range(per*2):
        gpu=i%2
        nm,seed,ex=SPEC[i] if i<len(SPEC) else (f'x{i}',100+i,[])
        out=g/nm; out.mkdir()
        cmd=['nice','-n','10',PY,'-u',str(T/'train_ccwm.py'),'--data-dir',DATA,'--out-dir',str(out),'--run-dir',str(out/'run'),'--potential-beta','5','--sample-latent-a','--dsm-mode','blood_marginal','--sinkhorn-iters','20','--seed',str(seed),'--epochs','5','--steps-per-epoch','200','--save-epochs','']+ex
        env=dict(os.environ,CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS=str(th),MKL_NUM_THREADS=str(th),OPENBLAS_NUM_THREADS=str(th),CCWM_NUM_THREADS=str(th),PYTHONDONTWRITEBYTECODE='1')
        ps.append((nm,subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=open(out/'console.log','w'),stderr=subprocess.STDOUT,env=env,cwd=T,start_new_session=True),out))
    (g/'pids.txt').write_text('\n'.join(f'{n} {p.pid}' for n,p,_ in ps)+'\n')
    t0=time.time(); [p.wait() for _,p,_ in ps]
    rates=[]
    for nm,_,out in ps:
        v=json.load(open(out/'validation.json')); rates.append(800/(v[-1]['seconds']-v[0]['seconds']))
    r=dict(config=name,per_gpu=per,threads_per_proc=th,n_procs=per*2,total_threads=per*2*th,wall_s=round(time.time()-t0,1),steps_s_per_proc=round(sum(rates)/len(rates),3),steps_s_total=round(sum(rates),2),min_rate=round(min(rates),3),max_rate=round(max(rates),3))
    res.append(r); (O/'checks'/'bench_summary.json').write_text(json.dumps(res,indent=1)); print(r,flush=True)
print('BENCH_DONE',flush=True)
