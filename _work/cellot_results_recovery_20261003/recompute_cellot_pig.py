"""Restore only CellOT pig predictions using the unchanged measured-blood adapter."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='2'
import sys,json
from pathlib import Path
import numpy as np,pandas as pd,torch
sys.path.insert(0,'/public/home/mengxl/dzy/pd_product/_work/pig_external_validation_20261003')
import common as c
import evaluate_external as e
from predict_external import to_linear_mean
import recompute_cellot as r
from cellot_icnns import ICNN
D=r.R/'pig';D.mkdir(exist_ok=True)

def main():
 torch.set_num_threads(2);device='cuda:1'
 ref=dict(np.load(c.O/'human_reference_cells.npz'));adapted=dict(np.load(c.O/'adapted_pig_blood_reference_inputs.npz'))
 genes=ref['genes'];mask=ref['mask'];names=adapted['sample_ids'];assert np.array_equal(genes,adapted['genes'])
 model=c.load_model(42,device);center=model.brain_center.cpu().numpy();del model
 preds=np.empty((2,len(c.AXES),len(names),len(genes)),np.float32);diag=[]
 for ai,axis in enumerate(c.AXES):
  k=r.read(axis);p=r.Projection(k);bs,_=axis.split('__')
  for condition in [0,1]:
   ck=torch.load(r.R/'models/seed_42'/f'cellot_{axis}_{condition}.pt',map_location=device,weights_only=False)
   net=ICNN(64,[64]*4,fnorm_penalty=1,kernel_init_fxn=lambda w:torch.nn.init.uniform_(w,0,.1)).to(device);net.load_state_dict(ck['state_dict']);net.eval()
   for param in net.parameters():param.requires_grad_(False)
   for j,name in enumerate(names):
    x=adapted[bs][j];xp=p.transform(x)
    source=torch.as_tensor(xp/ck['scale'],device=device,dtype=torch.float32).requires_grad_(True)
    with torch.enable_grad():cp=net.transport(source).detach().cpu().numpy()*ck['scale']
    cp=p.inverse_transform(cp);preds[condition,ai,j]=to_linear_mean(cp,center)
    diag.append(dict(condition=condition,axis=axis,sample_id=name,method='cellot',fraction_negative_implied_log_expression=float(np.mean(cp+center<0)),max_implied_log_expression=float(np.max(cp+center)),fraction_above_log10001=float(np.mean(cp+center>np.log(10001))),predicted_mean_linear=float(preds[condition,ai,j].mean())))
 np.savez_compressed(D/'cellot_predictions.npz',linear_profiles=preds,axes=np.array(c.AXES),sample_ids=names,genes=genes)
 pd.DataFrame(diag).to_csv(D/'prediction_diagnostics.csv',index=False)
 # All predictions saved before opening measured brain.
 meta=pd.read_csv(c.O/'inputs/matched_animals.csv').set_index('sample_id').loc[names].reset_index();group=meta.group.eq('LPS').to_numpy();ctrl=~group
 brain=pd.read_csv(c.O/'inputs/brain_one2one_sum_provided_transcript_TPM.csv',index_col=0).loc[genes,names].to_numpy().T
 obs=e.response(brain,mask,ctrl);mapped=genes[mask];audit=pd.read_csv(c.O/'inputs/existing_TPM_feature_audit.csv');hc=set(audit[(audit.tissue=='Brain')&audit.high_confidence].model_gene)
 markers=pd.read_csv(c.O/'human_training_marker_sets.csv');subsets={'all_one2one':np.ones(mask.sum(),bool),'high_confidence':np.isin(mapped,list(hc))}
 for state in c.RS:subsets[state+'_markers']=np.isin(mapped,markers[markers.receiver.eq(state)].gene)
 rows=[];axisrows=[];primary=None
 for condition in [0,1]:
  pp=np.stack([e.response(preds[condition,ai],mask,ctrl) for ai in range(len(c.AXES))]);aggregate=pp.mean(0)
  if condition==0:
   primary=aggregate;pd.DataFrame(aggregate.T,index=mapped,columns=names).to_csv(D/'predicted_response_cellot.csv',index_label='gene')
   pd.DataFrame({'gene':mapped,'cellot':e.effect(aggregate,group)}).to_csv(D/'gene_level_group_effects.csv',index=False)
  for subset,keep in subsets.items():
   row=dict(method='cellot',condition=condition,subset=subset,**e.metrics(aggregate[:,keep],obs[:,keep],group))
   if condition==0 and subset=='all_one2one':
    row.update(e.exact_group_p(aggregate,obs,group));pr,null=e.paired_test(aggregate,obs,group);row.update(pr);np.savez_compressed(D/'within_group_pairing_nulls.npz',cellot=null)
   rows.append(row)
  for ai,axis in enumerate(c.AXES):axisrows.append(dict(method='cellot',condition=condition,axis=axis,**e.metrics(pp[ai],obs,group)))
 pd.DataFrame(rows).to_csv(D/'response_metrics.csv',index=False);pd.DataFrame(axisrows).to_csv(D/'all_six_axis_metrics.csv',index=False)
 rng=np.random.default_rng(20261003);cc=np.flatnonzero(ctrl);ll=np.flatnonzero(group);bg=np.r_[np.zeros(6,bool),np.ones(4,bool)];boot=[]
 for rep in range(1000):
  ix=np.r_[rng.choice(cc,6,replace=True),rng.choice(ll,4,replace=True)];m=e.metrics(primary[ix],obs[ix],bg)
  boot.append(dict(replicate=rep,method='cellot',spearman=m['spearman'],improvement_over_zero=m['improvement_over_zero']))
 boot=pd.DataFrame(boot);boot.to_csv(D/'paired_animal_bootstrap.csv',index=False)
 ci=dict(method='cellot',spearman_ci_low=float(boot.spearman.quantile(.025)),spearman_ci_high=float(boot.spearman.quantile(.975)),improvement_ci_low=float(boot.improvement_over_zero.quantile(.025)),improvement_ci_high=float(boot.improvement_over_zero.quantile(.975)))
 pd.DataFrame([ci]).to_csv(D/'primary_bootstrap_intervals.csv',index=False)
 primaryrow=next(x for x in rows if x['condition']==0 and x['subset']=='all_one2one');primaryrow.update(ci);pd.DataFrame([primaryrow]).to_csv(D/'primary_summary.csv',index=False)
 loo=[]
 for j,name in enumerate(names):
  keep=np.arange(len(names))!=j;loo.append(dict(removed_animal=name,removed_group=meta.group.iloc[j],removed_duration_min=int(meta.duration_min.iloc[j]),method='cellot',**e.metrics(primary[keep],obs[keep],group[keep])))
 pd.DataFrame(loo).to_csv(D/'leave_one_animal_out.csv',index=False)
 print('CELLOT_PIG_COMPLETE',primaryrow,flush=True)
if __name__=='__main__':main()
