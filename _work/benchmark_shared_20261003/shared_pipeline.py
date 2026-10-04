"""All predictors share frozen encoders and brain decoder; validation-only readout tuning."""
import json,os
from pathlib import Path
import numpy as np,pandas as pd,torch
import current_benchmark as b
METHODS=['periscope','ot_ridge','conditional_mlp','conditional_icnn','shared_label_mean','shared_label_marginal']
D=b.O/'shared_readout';D.mkdir(exist_ok=True)
class Predictors:
 def __init__(self,model,k,axis,seed):
  self.model=model;self.k=k;self.axis=axis;self.seed=seed;self.ridges={};self.shared={}
  self.label_means={c:b.decode(model,k[f'train_yz{c}']).mean(0) for c in [0,1]}
  md=b.O/'models' if seed==42 else b.O/'models'/f'seed_{seed}'
  for c in [0,1]:self.ridges[c]=dict(np.load(md/f'ridge_{axis}_{c}.npz'))
  for m in ['conditional_mlp','conditional_icnn']:
   state=torch.load(md/f'{m}_eight_axis.pt',map_location=b.DEVICE,weights_only=False)
   net=b.ConditionalMap(m,model.cfg.z_dim,model.cfg.n_proteins+model.cfg.cond_dim).to(b.DEVICE);net.load_state_dict(state['state_dict']);net.eval();self.shared[m]=(net,state)
 def predict(self,m,z,u,c):
  n=len(z);model=self.model;k=self.k
  if m.startswith('shared_label'):
   if m.endswith('mean'):return np.repeat(self.label_means[c][None],n,0)
   yz=k[f'train_yz{c}'];return b.decode(model,yz[np.arange(n)%len(yz)])
  if m=='ot_ridge':
   r=self.ridges[c];f=(np.concatenate([z,u],1)-r['mean'])/r['sd'];return b.decode(model,f@r['coef']+r['target_mean'])
  cc=np.repeat(k[f'cond{c}'][None],n,0)
  if m=='periscope':
   zz=torch.as_tensor(z,device=b.DEVICE);uu=torch.as_tensor(u,device=b.DEVICE);ct=torch.as_tensor(cc,device=b.DEVICE)
   with torch.enable_grad():_,out=model.solve_equilibrium(zz,zz.clone(),uu,ct)
   return b.decode(model,out.detach().cpu().numpy())
  net,r=self.shared[m];zz=torch.as_tensor((z-r['zm'])/r['zs'],device=b.DEVICE);ctx=torch.as_tensor(np.concatenate([(u-r['um'])/r['us'],cc],1),device=b.DEVICE)
  with torch.enable_grad():out=net(zz,ctx).detach().cpu().numpy()*r['zs']+r['zm']
  return b.decode(model,out)

def noise_draw(bank,n):
 ix=np.random.default_rng(20261003).choice(len(bank),n,replace=False);noise=bank[ix].astype(float);return noise-noise.mean(0)

def main():
 torch.set_num_threads(2);seed=b.RUN_SEED;model=b.load_model(seed,b.DEVICE)
 for p in model.parameters():p.requires_grad_(False)
 rows=[];choices=[];selected=os.environ.get('BENCHMARK_METHODS',','.join(METHODS)).split(',')
 for axis in b.AXES:
  k=b.read(axis);predict=Predictors(model,k,axis,seed);bank=[]
  for c in [0,1]:
   y=k[f'train_y{c}'];r=y-b.decode(model,b.encode(model,y,brain=True)[1]);bank.append(r-r.mean(0))
  bank=np.concatenate(bank);record={}
  for m in selected:
   train=[predict.predict(m,k[f'train_z{c}'][:512],k[f'train_u{c}'][:512],c) for c in [0,1]]
   val=[predict.predict(m,k[f'val_z{c}'],k[f'val_u{c}'],c).astype(float) for c in [0,1]]
   offset=np.mean([k[f'train_y{c}'][:512].astype(float).mean(0)-train[c].astype(float).mean(0) for c in [0,1]],axis=0)
   # Condition-independent correction; disease contrasts are unchanged.
   alphas=[0.,.25,.5,1.];vs=[float(np.mean([np.mean((val[c].mean(0)+a*offset-k[f'val_y{c}'].mean(0))**2) for c in [0,1]])) for a in alphas]
   alpha=alphas[int(np.argmin(vs))];noisev=noise_draw(bank,len(val[0]));scales=[0.,.25,.5,1.]
   ms=[float(np.mean([b.mmd(val[c]+alpha*offset+t*noisev,k[f'val_y{c}'],float(k['bw'])) for c in [0,1]])) for t in scales];scale=scales[int(np.argmin(ms))]
   original=b.O/'results'/f'seed_{seed}'/f'{m}__{axis}.npz'
   if m.startswith('shared_label'):
    raw={c:{bc:predict.predict(m,k[f'test_z{bc}'],k[f'test_u{bc}'],c) for bc in [0,1]} for c in [0,1]}
    dest=b.RESULT_DIR;b.RESULT_DIR=b.O/'results'/f'seed_{seed}';b.assess(m,raw,k,seed,{'representation':'same frozen brain encoder/decoder; training latent reference, no blood input'});b.RESULT_DIR=dest
   else:
    x=np.load(original);raw={c:{bc:x[f'map{c}_blood{bc}'].astype(float) for bc in [0,1]} for c in [0,1]}
   noise=noise_draw(bank,128);pred={c:{bc:raw[c][bc].astype(float)+alpha*offset+scale*noise for bc in [0,1]} for c in [0,1]}
   before=raw[1][1].astype(float).mean(0)-raw[0][0].astype(float).mean(0);after=pred[1][1].mean(0)-pred[0][0].mean(0);err=float(np.max(abs(before-after)));assert err<1e-10
   dest=b.RESULT_DIR;b.RESULT_DIR=D/'results'/f'seed_{seed}';b.RESULT_DIR.mkdir(parents=True,exist_ok=True)
   rec=dict(offset_alpha=alpha,residual_scale=scale,validation_mean_mse_grid=vs,validation_mmd_grid=ms,alpha_grid=alphas,scale_grid=scales,selection='human validation only; identical grids for all methods',offset='pooled training-only, condition-independent',noise='same centered training brain reconstruction residual draw for all methods',contrast_error=err)
   b.assess(m,pred,k,seed,rec);b.RESULT_DIR=dest
   par=D/'parameters'/f'seed_{seed}';par.mkdir(parents=True,exist_ok=True);np.savez_compressed(par/f'{m}__{axis}.npz',offset=offset,alpha=alpha,residual_scale=scale,noise=noise)
   record[m]=rec;choices.append(dict(method=m,axis=axis,seed=seed,alpha=alpha,residual_scale=scale,validation_mean_mse_before=vs[0],validation_mean_mse_after=min(vs),validation_mmd_no_residual=ms[0],validation_mmd_after=min(ms),contrast_error=err))
  cp=D/f'choices_seed{seed}_{axis}.json';old=json.loads(cp.read_text()) if cp.exists() else {};old.update(record);cp.write_text(json.dumps(old,indent=2));print('OPTIMIZED',seed,axis,flush=True)
 cp=D/f'validation_choices_seed{seed}.csv';old=pd.read_csv(cp) if cp.exists() else pd.DataFrame();new=pd.DataFrame(choices)
 if len(old):new=pd.concat([old[~old.method.isin(selected)],new])
 new.to_csv(cp,index=False)
if __name__=='__main__':main()
