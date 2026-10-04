"""Repeat exactly seeds47/48/49 with post-warmup checkpoint eligibility."""
import concurrent.futures,fcntl,hashlib,json,os,subprocess,sys,threading,time
from pathlib import Path
W=Path(__file__).resolve().parent
A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/benchmark_20261002/seed_repair';M=A/'results/optimization_20261003/models'
O.mkdir(parents=True,exist_ok=True)
lock=(O/'launch.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
mu=threading.Lock();state={'status':'training','seeds':[47,48,49],'jobs':{}}
def save():
 with mu:
  p=O/'run_state.tmp';p.write_text(json.dumps(state,indent=2));p.replace(O/'run_state.json')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(seed,gpu):
 old=json.loads((M/f'seed_{seed}/training_inputs.json').read_text())['args'];args=old.copy()
 out=O/f'models/seed_{seed}';out.mkdir(parents=True,exist_ok=True)
 args.update(out_dir=str(out),run_dir=str(out/'run'),run_id=f'post-warmup-repair-seed-{seed}')
 cmd=[sys.executable,'-u',str(W/'train_ccwm.py')]
 for k,v in args.items():
  if isinstance(v,bool):
   if v:cmd+=['--'+k.replace('_','-')]
  elif v is not None:cmd+=['--'+k.replace('_','-'),str(v)]
 env=os.environ.copy();env.update(PYTHONPATH='/public/home/mengxl/dzy/pd_product/src',CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
 if (out/'completed.json').exists():
  state['jobs'][str(seed)]={'status':'complete','reused_completed_repair':True};save();return
 with (out/'console.log').open('w') as f:
  p=subprocess.Popen(cmd,stdout=f,stderr=subprocess.STDOUT,env=env,stdin=subprocess.DEVNULL)
  state['jobs'][str(seed)]={'status':'training','pid':p.pid,'gpu':gpu,'command':cmd,'started_unix':time.time()};save()
  code=p.wait()
 state['jobs'][str(seed)].update(status='complete' if code==0 else 'failed',returncode=code);save()
 if code:raise RuntimeError(f'seed{seed} failed')
 r=json.loads((out/'completed.json').read_text());assert r['best_epoch']>=50
 print('TRAINING_COMPLETE',seed,r,flush=True)

if __name__=='__main__':
 protocol={'seeds':[47,48,49],'change':'checkpoint eligibility starts after max(min_epochs, KL warmup); patience starts from eligible best; final training state saved','unchanged':'model architecture, data, losses, learning rate, batch sizes, original seed, maximum200 epochs x200steps','selection':'original composite internal validation score; no test selection','source_sha256':{p.name:sha(p) for p in [W/'train_ccwm.py',W/'ccwm.py',W/'training_utils.py']},'original_checkpoints':{str(s):sha(M/f'seed_{s}/model.pt') for s in range(42,52)},'output_root':str(O)}
 (O/'protocol.json').write_text(json.dumps(protocol,indent=2));save()
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
   jobs=[pool.submit(run,s,g) for s,g in [(47,0),(48,1),(49,0)]]
   for j in concurrent.futures.as_completed(jobs):j.result()
  state['status']='training_complete';save();print('ALL_THREE_REPAIRS_TRAINED',flush=True)
 except Exception as e:state['status']='failed';state['error']=str(e);save();raise
