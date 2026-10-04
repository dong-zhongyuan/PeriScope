"""Compare recomputed CellOT with recovered original numeric outputs."""
import json
import numpy as np,pandas as pd
import recompute_cellot as r
b=r.b;R=r.R;old=b.O/'cellot_restored'
METRICS=['mean_mse','mmd','energy_per_sqrt_gene','variance_ratio','condition_aware_rho','condition_aware_delta_mse','condition_aware_predicted_rms','fixed_control_rho','fixed_control_delta_mse','fixed_control_predicted_rms','fixed_PD_rho','fixed_PD_delta_mse','fixed_PD_predicted_rms']

def main():
 rows=[];hashes={}
 for seed in [42,43,44]:
  for axis in b.AXES:
   path=R/'results'/f'seed_{seed}'/f'cellot__{axis}.json';z=json.loads(path.read_text())
   row=dict(method='cellot',seed=seed,axis=axis,blood_state=axis.split('__')[0],brain_state=axis.split('__')[1],**z['distribution_average'])
   for task,d in z['direction'].items():row.update({task+'_'+key:v for key,v in d.items()})
   rows.append(row)
   pred=dict(np.load(path.with_suffix('.npz')));assert len(pred)==4
   assert all(v.shape==(128,4741) and np.isfinite(v).all() for v in pred.values())
   hashes[str(path.relative_to(R))]=b.sha(path);hashes[str(path.with_suffix('.npz').relative_to(R))]=b.sha(path.with_suffix('.npz'))
 df=pd.DataFrame(rows);assert len(df)==24;df.to_csv(R/'cellot_metrics.csv',index=False)
 means=df.groupby('seed')[METRICS].mean();summary=means.agg(['mean','std']);summary.to_csv(R/'cellot_summary.csv');checks=[]
 def check(label,actual,expected,atol=2e-6,rtol=1e-5):
  checks.append(dict(item=label,actual=float(actual),expected=float(expected),absolute_error=float(abs(actual-expected)),passed=bool(np.isclose(actual,expected,atol=atol,rtol=rtol))))
 original=pd.read_csv(old/'comparison_summary.csv',header=[0,1],index_col=0).loc['cellot']
 for metric,stat in original.index:check('macro/'+metric+'/'+stat,summary.loc[stat,metric],original[(metric,stat)])
 original=pd.read_csv(old/'direction_by_axis.csv').set_index('axis');axisstats=df.groupby('axis').condition_aware_rho.agg(['mean','std'])
 for axis,z in original.iterrows():
  for stat in ['mean','std']:check('axis_direction/'+axis+'/'+stat,axisstats.loc[axis,stat],z['condition_aware_rho_'+stat],atol=5.1e-6,rtol=0)
 original=pd.read_csv(old/'raw_axis_summary_partial.csv',header=[0,1],index_col=[0,1,2]);new=df.groupby(['method','blood_state','brain_state'])[METRICS].agg(['mean','std','min','max'])
 for ix,z in original.iterrows():
  for metric,stat in z.index:check('raw_axis/'+'/'.join(ix)+'/'+metric+'/'+stat,new.loc[ix,(metric,stat)],z[(metric,stat)])
 residual=pd.read_csv(R/'readout_repair/metrics.csv');assert len(residual)==24
 original=pd.read_csv(old/'residual_axis_summary_partial.csv',header=[0,1],index_col=[0,1]);new=residual.groupby(['method','axis'])[['mmd','mean_mse','variance_ratio','condition_aware_rho','condition_aware_delta_mse']].agg(['mean','std'])
 for ix,z in original.iterrows():
  for metric,stat in z.index:check('residual_axis/'+'/'.join(ix)+'/'+metric+'/'+stat,new.loc[ix,(metric,stat)],z[(metric,stat)])
 original=pd.read_csv(old/'pig_primary_summary.csv').iloc[0];new=pd.read_csv(R/'pig/primary_summary.csv').iloc[0]
 for key in original.index:
  if key in ['method','condition','subset','BH_q_four_models'] or key not in new:continue
  check('pig_primary/'+key,new[key],original[key])
 original=pd.read_csv(old/'pig_response_metrics.csv').set_index(['condition','subset']);new=pd.read_csv(R/'pig/response_metrics.csv').set_index(['condition','subset'])
 for ix,z in original.iterrows():
  for key in z.index:
   if key in ['method','BH_q_four_models'] or key not in new or pd.isna(z[key]):continue
   check('pig_response/'+str(ix)+'/'+key,new.loc[ix,key],z[key])
 for p in (R/'models').rglob('*.pt'):hashes[str(p.relative_to(R))]=b.sha(p)
 assert len(list((R/'models').rglob('*.pt')))==48
 pd.DataFrame(checks).to_csv(R/'comparison_with_original.csv',index=False)
 result=dict(status='passed' if all(x['passed'] for x in checks) else 'differences_found',checks=len(checks),passed=sum(x['passed'] for x in checks),failed=sum(not x['passed'] for x in checks),raw_results=24,models=48,residual_results=24,pig_axes=6,files=hashes)
 (R/'recovery_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='files'}),flush=True)
 if result['failed']:print(pd.DataFrame(checks).query('not passed').head(30).to_string(index=False),flush=True)
if __name__=='__main__':main()
