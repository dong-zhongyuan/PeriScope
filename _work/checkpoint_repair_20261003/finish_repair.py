"""Finish current authorized training jobs, validate, and refresh the benchmark."""
import json,hashlib,os,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np,pandas as pd,torch
W=Path(__file__).resolve().parent
A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/benchmark_20261002';R=O/'seed_repair';B=W.parent/'benchmark_full_20261003';M=A/'results/optimization_20261003/models'
def dump(p,x):
 q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(x,indent=2));q.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(script,*args,seed=None):
 env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES='1',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',OMP_NUM_THREADS='2')
 if seed is not None:env['BENCHMARK_SEED']=str(seed)
 subprocess.run([sys.executable,'-u',str(B/script),*args],env=env,check=True)

def finish():
 protocol=json.loads((R/'protocol.json').read_text());records=[]
 # Original models stay available to the already completed downstream pipeline.
 for s in range(42,52):assert sha(M/f'seed_{s}/model.pt')==protocol['original_checkpoints'][str(s)]
 for s in [47,48,49]:
  p=R/f'models/seed_{s}';d=torch.load(p/'model.pt',map_location='cpu',weights_only=False)
  v=json.loads((p/'validation.json').read_text());eligible=[x for x in v if x['checkpoint_eligible']]
  assert eligible and all(x['epoch']>=50 for x in eligible)
  best=next(x for x in eligible if x['epoch']==d['best_epoch'])
  assert best['selection_score']<=min(x['selection_score'] for x in eligible)+1e-5
  assert d['cfg']['beta_kl']==.01 and d['best_epoch']>=50
  original=json.loads((M/f'seed_{s}/training_inputs.json').read_text())['args']
  actual=json.loads((p/'training_inputs.json').read_text())['args']
  assert {k:v for k,v in original.items() if k not in ['out_dir','run_dir','run_id']}=={k:v for k,v in actual.items() if k not in ['out_dir','run_dir','run_id']}
  records.append(dict(seed=s,best_epoch=d['best_epoch'],best_validation_score=d['best_validation_score'],model_sha256=sha(p/'model.pt'),steps=json.loads((p/'completed.json').read_text())['steps']))
 dump(R/'checkpoint_validation.json',{'status':'passed','records':records,'all_original_model_hashes_unchanged':True})
 previous=W/'previous_benchmark';previous.mkdir(exist_ok=True)
 # Store small before/after evidence and remove stale representation caches from active paths.
 if not (previous/'periscope_ten_seed_metrics.csv').exists():shutil.copy2(O/'periscope_ten_seed_metrics.csv',previous/'periscope_ten_seed_metrics.csv')
 overrides={'reason':'user-authorized post-warmup checkpoint repair for seeds47/48/49','checkpoints':{str(s):str(R/f'models/seed_{s}/model.pt') for s in [47,48,49]},'eligible_from_epoch':50}
 dump(O/'checkpoint_overrides.json',overrides)
 for s in [47,48,49]:
  oldcache=O/f'representations/seed_{s}';dst=previous/f'representations_seed_{s}'
  if oldcache.exists() and not dst.exists():shutil.move(str(oldcache),dst)
  out=O/f'results/seed_{s}';bak=previous/f'results_seed_{s}';bak.mkdir(exist_ok=True)
  for f in out.glob('periscope__*'):
   if not (bak/f.name).exists():shutil.move(str(f),bak/f.name)
  run('benchmark_worker.py','prepare','--seed',str(s),seed=s)
  run('benchmark_worker.py','stability','--seed',str(s),seed=s)
 # Full finalization recomputes metrics; the 132 primary comparison rows remain the same.
 run('finalize_benchmark.py')
 before=pd.read_csv(previous/'periscope_ten_seed_metrics.csv');after=pd.read_csv(O/'periscope_ten_seed_metrics.csv')
 paired=before.merge(after,on=['method','seed','blood_state','brain_state'],suffixes=('_before','_after'))
 paired.to_csv(R/'before_after_all_seeds.csv',index=False)
 fixed=paired[paired.seed.isin([47,48,49])];fixed.to_csv(R/'before_after_repaired_seeds.csv',index=False)
 for s in [42,43,44,45,46,50,51]:
  rr=paired[paired.seed==s]
  for m in ['mmd','condition_aware_rho','condition_aware_delta_mse']:assert np.array_equal(rr[m+'_before'].to_numpy(),rr[m+'_after'].to_numpy())
 lines=['47/48/49补训完成：预热后统一选模','', '架构、训练数据、损失权重、随机种子、学习率和最大训练预算均沿用原配置；只修正检查点资格与早停起点。']
 lines+=['seed | best epoch | validation score | steps']+[f"{r['seed']} | {r['best_epoch']} | {r['best_validation_score']:.6f} | {r['steps']}" for r in records]
 lines+=['','修复种子：疾病方向相关性，修复前→修复后']
 for _,r in fixed.iterrows():lines.append(f"{r['seed']} {r['blood_state']}→{r['brain_state']}: {r['condition_aware_rho_before']:.5f} → {r['condition_aware_rho_after']:.5f}")
 lines+=['','当前十种子方向稳定性：mean ± SD；min–max']
 for (b,r),df in after.groupby(['blood_state','brain_state']):
  x=df.condition_aware_rho;lines.append(f'{b}→{r}: {x.mean():.5f} ± {x.std():.5f}; {x.min():.5f}–{x.max():.5f}')
 lines+=['','其余七个种子的指标逐项保持一致。三种子主benchmark（42/43/44）以及两种方法的观测残差读出比较保持不变。','当前benchmark通过checkpoint_overrides.json使用补训模型；原模型保留供已经完成的靶点分析追溯。未新增或重跑训练消融。']
 (O/'补训结果说明.txt').write_text('\n'.join(lines)+'\n')
 readout=pd.read_csv(O/'readout_repair/metrics.csv').groupby('method').mmd.mean()
 header=f'补训更新：47/48/49已按预热后选模规则完成补训，十种子汇总已更新，见“补训结果说明.txt”。\n读出修复版：PeriScope平均MMD {readout["periscope"]:.5f}，CellOT同等修复后{readout["cellot"]:.5f}；见readout_repair/。下文为原始读出下的当前benchmark。\n\n'
 p=O/'评估与运行说明.txt';p.write_text(header+p.read_text())
 rep=R/'reproduce';rep.mkdir(exist_ok=True)
 for p in W.glob('*.py'):shutil.copy2(p,rep/p.name)
 status=json.loads((O/'run_status.json').read_text());status['checkpoint_repair']='complete: 47,48,49';dump(O/'run_status.json',status)
 dump(R/'completion.json',{'status':'complete','records':records,'finished_unix':time.time(),'validation':'passed; original seven seeds and primary benchmark unchanged'})
 print('REPAIR_EVALUATED_AND_DELIVERABLE',flush=True)

if __name__=='__main__':
 try:
  for _ in range(540):
   st=json.loads((R/'run_state.json').read_text())
   if st['status']=='training_complete':break
   if st['status']=='failed':raise RuntimeError(st)
   time.sleep(20)
  else:raise TimeoutError('Training exceeds three hours; inspect current run before continuing')
  finish()
 except Exception as e:
  dump(R/'finalization_error.json',{'error':str(e),'time':time.time()});raise
