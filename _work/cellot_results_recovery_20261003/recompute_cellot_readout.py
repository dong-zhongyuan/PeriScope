"""Frozen residual scale=1 chosen on seed42 validation; same repair for both methods."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd,torch
sys.path.insert(0,'/public/home/mengxl/dzy/pd_product/_work/benchmark_shared_20261003')
import current_benchmark as b
import recompute_cellot as recovery
b.Projection=recovery.Projection
W=Path(__file__).resolve().parent
D=recovery.R/'readout_repair'
D.mkdir(exist_ok=True)

def main():
 torch.set_num_threads(2);rows=[];valrows=[];checks=[]
 (D/'protocol.json').write_text(json.dumps({'scale':1.,'arithmetic':'float64 residual addition and contrasts to avoid float32 cancellation in large extrapolated outputs','selection':'smallest validation MMD on scale grid 0,0.5,1, evaluated before held-out scoring','pilot':'seed42; same scale grid and repair given to both methods','observation_model':'zero-mean empirical training reconstruction residuals pooled across conditions within each axis','same_noise_across_conditions_and_blood_inputs':True,'seed':20261003,'source_sha256':b.sha(Path(__file__))},indent=2))
 for seed in [42,43,44]:
  model=None;b.RESULT_DIR=D/'results'/f'seed_{seed}';b.RESULT_DIR.mkdir(parents=True,exist_ok=True)
  for axis in b.FOCUS_AXES:
   k=recovery.read(axis)
   p=b.Projection(k)
   for method in ['cellot']:
    bank=[]
    for c in (0,1):
     y=k[f'train_y{c}']
     if method=='periscope':rec=b.decode(model,b.encode(model,y,brain=True)[1])
     else:rec=p.inverse_transform(p.transform(y))
     r=y-rec;r-=r.mean(0);bank.append(r)
    bank=np.concatenate(bank);rng=np.random.default_rng(20261003);indices=rng.choice(len(bank),128,replace=False);noise=bank[indices]
    old=recovery.R/'results'/f'seed_{seed}'/f'{method}__{axis}'
    raw={n:v.astype(np.float64) for n,v in np.load(old.with_suffix('.npz')).items()};pred={c:{bc:raw[f'map{c}_blood{bc}']+noise for bc in (0,1)} for c in (0,1)}
    name=method+'_observation_residual'
    b.assess(name,pred,k,seed,{'readout':'centered training reconstruction residuals','scale':1.,'residual_bank_size':len(bank),'noise_indices':indices.tolist(),'raw_prediction_sha256':b.sha(old.with_suffix('.npz'))})
    f=b.RESULT_DIR/f'{name}__{axis}.json';r=json.loads(f.read_text());orig=json.loads(old.with_suffix('.json').read_text())
    for c1,b1,c0,b0 in [(1,1,0,0),(0,1,0,0),(1,1,1,0)]:assert np.allclose(pred[c1][b1].mean(0)-pred[c0][b0].mean(0),raw[f'map{c1}_blood{b1}'].mean(0)-raw[f'map{c0}_blood{b0}'].mean(0),atol=2e-6)
    row=dict(seed=seed,axis=axis,method=method,**r['distribution_average'])
    row.update({t+'_'+n:v for t,d in r['direction'].items() for n,v in d.items()})
    row.update({'raw_'+n:v for n,v in orig['distribution_average'].items()});rows.append(row)
    checks.append({'file':str(f.relative_to(D)),'sha256':b.sha(f),'prediction_sha256':b.sha(f.with_suffix('.npz')),'max_disease_contrast_change':float(np.max(np.abs((pred[1][1].mean(0)-pred[0][0].mean(0))-(raw['map1_blood1'].mean(0)-raw['map0_blood0'].mean(0)))))})
   print('DONE',seed,axis,flush=True)
 df=pd.DataFrame(rows);df.to_csv(D/'metrics.csv',index=False)
 df.groupby(['method','axis'])[['mmd','mean_mse','variance_ratio','condition_aware_rho','condition_aware_delta_mse']].agg(['mean','std']).to_csv(D/'axis_summary.csv')
 (D/'checks.json').write_text(json.dumps({'status':'passed','result_count':len(rows),'outputs':checks},indent=2))
 print('CELLOT_RESIDUAL_COMPLETE',len(rows),flush=True)
if __name__=='__main__':main()
