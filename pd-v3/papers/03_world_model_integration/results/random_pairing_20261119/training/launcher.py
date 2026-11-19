"""Priority queue launcher for the PeriScope training lane. Keeps QUEUE.csv current.
A failed run is left as failed; it is retried only by an explicit --retry arm/seed after diagnosis."""
import argparse, csv, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path
O=Path(__file__).resolve().parent; T=O/'code'; R=O/'runs'; Q=O/'QUEUE.csv'
PY='/public/home/mengxl/dzy/envs/pd_bbm/bin/python'; DATA='/public/home/mengxl/dzy/pd_product_assets/interim/rescreen_20261002'
BASE=['--potential-beta','5','--sample-latent-a','--dsm-mode','blood_marginal','--sinkhorn-iters','20','--epochs','200','--steps-per-epoch','200']
ARMS={'full':[], 'disease_shuffled':['--pairing','disease_shuffled'], 'random_pairing':['--pairing','random'], 'coupling_off':['--coupling','off']}
PLAN=[(1,'full',s) for s in range(52,62)]+[(2,'disease_shuffled',s) for s in (42,43,44)]+[(3,'random_pairing',s) for s in range(52,57)]+[(4,'coupling_off',s) for s in range(52,57)]
COLS=['order','priority','arm','seed','status','attempt','gpu','pid','start','end','seconds','best_epoch','epochs_run','steps','best_validation_score','checkpoint','model_sha256','note']
ap=argparse.ArgumentParser(); ap.add_argument('--per-gpu',type=int,required=True); ap.add_argument('--threads',type=int,required=True)
ap.add_argument('--min-free-mib',type=int,default=3000); ap.add_argument('--retry',action='append',default=[]); a=ap.parse_args()
now=lambda: time.strftime('%Y-%m-%d %H:%M:%S')
def save(rows):
    t=Q.with_suffix('.tmp')
    with t.open('w',newline='') as f:
        w=csv.DictWriter(f,COLS); w.writeheader(); w.writerows(rows)
    t.replace(Q)
def load():
    if Q.exists(): return list(csv.DictReader(Q.open()))
    return [dict({c:'' for c in COLS},order=i+1,priority=p,arm=m,seed=s,status='pending',attempt=0) for i,(p,m,s) in enumerate(PLAN)]
def out(r): return R/r['arm']/f"seed_{r['seed']}"
def free(g):
    x=subprocess.run(['nvidia-smi','-i',str(g),'--query-gpu=memory.free','--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=60)
    return int(x.stdout.strip())
def alive(r):
    try: return str(out(r)) in Path(f"/proc/{r['pid']}/cmdline").read_bytes().decode(errors='ignore').replace('\0',' ')
    except Exception: return False
def start(r,g):
    d=out(r); d.mkdir(parents=True,exist_ok=False)
    cmd=['nice','-n','10',PY,'-u',str(T/'train_ccwm.py'),'--data-dir',DATA,'--out-dir',str(d),'--run-dir',str(d/'run'),'--run-id',f"{r['arm']}-seed-{r['seed']}",'--seed',str(r['seed'])]+BASE+ARMS[r['arm']]
    t=str(a.threads); env=dict(os.environ,CUDA_VISIBLE_DEVICES=str(g),OMP_NUM_THREADS=t,MKL_NUM_THREADS=t,OPENBLAS_NUM_THREADS=t,CCWM_NUM_THREADS=t,PYTHONDONTWRITEBYTECODE='1')
    p=subprocess.Popen(cmd,stdin=subprocess.DEVNULL,stdout=open(d/'console.log','w'),stderr=subprocess.STDOUT,env=env,cwd=T,start_new_session=True)
    r.update(status='running',gpu=g,pid=p.pid,start=now(),end='',note='',attempt=int(r['attempt'] or 0)+1); return p
def finish(r,rc):
    d=out(r); r['end']=now(); r['pid']=r['pid']
    try:
        c=json.loads((d/'completed.json').read_text()); v=json.loads((d/'validation.json').read_text())
        e=[x for x in v if x.get('checkpoint_eligible')]; m=min(x['selection_score'] for x in e); b=next(x for x in e if x['epoch']==c['best_epoch'])
        assert rc in (0,None) and c['status']=='training_complete' and c['best_epoch']>=50 and all(x['epoch']>=50 for x in e) and b['selection_score']<=m+1e-5 and (d/'model.pt').exists()
        r.update(status='stopped_early' if (R/'STOP').exists() else 'done',seconds=round(c['seconds'],1),best_epoch=c['best_epoch'],epochs_run=c['epochs'],steps=c['steps'],best_validation_score=f"{c['best_validation_score']:.6f}",checkpoint=str(d/'model.pt'),model_sha256=hashlib.sha256((d/'model.pt').read_bytes()).hexdigest())
    except Exception as ex:
        tail=(d/'console.log').read_text(errors='ignore').strip().splitlines()[-1:] if (d/'console.log').exists() else ['']
        r.update(status='failed',note=f'rc={rc}; {type(ex).__name__}; '+(tail[0][:160] if tail else ''))
rows=load(); R.mkdir(exist_ok=True)
for k in a.retry:
    m,s=k.split('/'); r=next(x for x in rows if x['arm']==m and str(x['seed'])==s); assert r['status']=='failed' and int(r['attempt'])==1,(k,r['status'],r['attempt'])
    d=out(r); shutil.move(str(d),str(d)+'.failed_attempt1'); r.update(status='pending',note='retry after diagnosis',gpu='',pid='',start='',end='')
procs={}; save(rows)
(O/'launcher_settings.json').write_text(json.dumps(dict(per_gpu=a.per_gpu,threads=a.threads,min_free_mib=a.min_free_mib,started=now(),pid=os.getpid()),indent=1))
last=0
while True:
    for r in rows:
        if r['status']!='running': continue
        k=(r['arm'],r['seed']); p=procs.get(k)
        rc=p.poll() if p else (None if alive(r) else 1)
        if p is None and rc==1: rc=0 if (out(r)/'completed.json').exists() else 1
        if rc is not None: finish(r,rc); procs.pop(k,None); save(rows); print(now(),r['arm'],r['seed'],r['status'],r['best_epoch'],r['note'],flush=True)
    stop=(R/'STOP').exists(); pend=[r for r in rows if r['status']=='pending']
    if not stop and pend and time.time()-last>=20:
        used={g:sum(1 for r in rows if r['status']=='running' and str(r['gpu'])==str(g)) for g in (0,1)}
        for g in sorted(used,key=lambda x:(used[x],x)):
            if used[g]<a.per_gpu and free(g)>=a.min_free_mib:
                r=pend[0]; procs[(r['arm'],r['seed'])]=start(r,g); save(rows); last=time.time(); print(now(),'start',r['arm'],r['seed'],'gpu',g,'pid',r['pid'],flush=True); break
    if not any(r['status']=='running' for r in rows) and (stop or not any(r['status']=='pending' for r in rows)): break
    time.sleep(10)
print('QUEUE_FINISHED',{s:sum(1 for r in rows if r['status']==s) for s in sorted({r['status'] for r in rows})},flush=True)
