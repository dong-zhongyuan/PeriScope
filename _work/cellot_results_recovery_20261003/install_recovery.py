"""Install verified CellOT-only artifacts; preserve every existing non-CellOT row."""
import csv,io,json,shutil
from pathlib import Path
import pandas as pd
import recompute_cellot as r
B=r.b.O;R=r.R;P=B.parent/'pig_external_validation_20261003';installed=[]

def copy_missing(src,dest):
 dest.parent.mkdir(parents=True,exist_ok=True)
 if dest.exists():
  assert r.b.sha(src)==r.b.sha(dest),f'Refusing to overwrite different existing file: {dest}'
 else:shutil.copy2(src,dest)
 installed.append(str(dest.relative_to(B.parent)))

def append_table(path,rows,index_columns,metric_columns=None):
 old=path.read_bytes();lines=old.decode().splitlines();header=list(csv.reader(lines))
 assert old.endswith(b'\n');records=[]
 if metric_columns is None:
  names=header[0]
  for row in rows:records.append([str(row.get(k,'')) if pd.notna(row.get(k,'')) else '' for k in names])
 else:
  for keys,stats in rows:
   records.append(list(keys)+[str(stats[(header[0][i],header[1][i])]) for i in range(index_columns,len(header[0]))])
 # Idempotent: the exact CellOT identity has one record per readout/axis/seed.
 existing={tuple(x[:index_columns]):x for x in header[(1 if metric_columns is None else 3):]}
 new=[]
 for row in records:
  key=tuple(row[:index_columns])
  if key in existing:
   assert existing[key]==row,f'Existing CellOT row differs: {path} {key}'
  else:new.append(row)
 if not new:return
 backup=R/'before_install'/path.name;backup.parent.mkdir(exist_ok=True)
 if not backup.exists():backup.write_bytes(old)
 stream=io.StringIO();csv.writer(stream,lineterminator='\n').writerows(new)
 path.write_bytes(old+stream.getvalue().encode());assert path.read_bytes()[:len(old)]==old
 installed.append(str(path.relative_to(B.parent)))

def main():
 check=json.loads((R/'recovery_validation.json').read_text());assert check['status']=='passed', 'Original-result comparisons must pass before installation'
 # Preserve an independent copy of the recovered weights and predictions first.
 backup=Path('/public/home/mengxl/dzy/pd_recovery_backup_20261003_cellot/recomputed_cellot')
 if not backup.exists():shutil.copytree(R,backup)
 for seed in [42,43,44]:
  modeldir=B/'models' if seed==42 else B/'models'/f'seed_{seed}'
  for src in (R/'models'/f'seed_{seed}').glob('cellot_*'):copy_missing(src,modeldir/src.name)
  for src in (R/'results'/f'seed_{seed}').glob('cellot__*'):copy_missing(src,B/'results'/f'seed_{seed}'/src.name)
  for src in (R/'readout_repair/results'/f'seed_{seed}').glob('cellot_observation_residual__*'):copy_missing(src,B/'readout_repair/results'/f'seed_{seed}'/src.name)
 for src in (R/'projection').glob('*.npz'):copy_missing(src,B/'cellot_projection'/src.name)
 raw=pd.read_csv(R/'cellot_metrics.csv');raw['readout']='original_cellot'
 residual=pd.read_csv(R/'readout_repair/metrics.csv');residual['readout']='original_cellot_observation_residual'
 for df in [raw,residual]:
  df['blood_state']=df.axis.str.split('__').str[0];df['brain_state']=df.axis.str.split('__').str[1]
 df=pd.concat([raw,residual],ignore_index=True)
 append_table(B/'baseline_metrics.csv',df.to_dict('records'),6)
 metrics=['mmd','mean_mse','variance_ratio','condition_aware_rho','condition_aware_delta_mse','fixed_control_rho','fixed_PD_rho']
 seed=df.groupby(['readout','method','seed'])[metrics].mean().reset_index();append_table(B/'seed_summary.csv',seed.to_dict('records'),3)
 axis=df.groupby(['readout','method','axis'])[metrics].agg(['mean','std']);append_table(B/'axis_summary.csv',[(ix,s.to_dict()) for ix,s in axis.iterrows()],3,metrics)
 # The previously recovered exact original raw summary remains as it is.
 summary=seed[seed.readout.eq('original_cellot_observation_residual')].groupby(['readout','method'])[metrics].agg(['mean','std'])
 append_table(B/'comparison_summary.csv',[(ix,s.to_dict()) for ix,s in summary.iterrows()],2,metrics)
 for name in ['cellot_metrics.csv','cellot_summary.csv','comparison_with_original.csv','recovery_validation.json','recovery_environment.json','execution_modes.json','projection_reconstruction.json']:
  copy_missing(R/name,B/'cellot_restored'/name)
 copy_missing(R/'readout_repair/metrics.csv',B/'cellot_restored/residual_metrics.csv')
 for src in (R/'pig').glob('*'):
  if src.is_file():copy_missing(src,P/'cellot_restored'/src.name)
 copy_missing(R/'pig/predicted_response_cellot.csv',P/'predicted_response_cellot.csv')
 training=[]
 for seed_value in [42,43,44]:
  for axis_value in r.b.AXES:
   for condition_value in [0,1]:
    log=json.loads((R/'models'/f'seed_{seed_value}'/f'cellot_{axis_value}_{condition_value}.json').read_text())
    chosen=next(h for h in log['history'] if h['step']==log['best_step'])
    training.append(dict(method='cellot',seed=seed_value,axis=axis_value,condition=condition_value,best_step=log['best_step'],steps_run=log['history'][-1]['step'],validation_mmd=chosen['validation_mmd'],stop_reason=log['stop_reason'],elapsed_s=log['history'][-1]['elapsed_s']))
 append_table(B/'training_summary.csv',training,4)
 # Restore CellOT detail rows while retaining existing other-method bytes.
 for name,nkeys in [('all_six_axis_metrics.csv',3),('leave_one_animal_out.csv',4),('paired_animal_bootstrap.csv',2),('primary_bootstrap_intervals.csv',1),('prediction_diagnostics.csv',4)]:
  append_table(P/name,pd.read_csv(R/'pig'/name).to_dict('records'),nkeys)
 bootstrap=pd.read_csv(P/'paired_animal_bootstrap.csv')
 paired=bootstrap[bootstrap.method.eq('periscope')][['replicate','spearman']].merge(bootstrap[bootstrap.method.eq('cellot')][['replicate','spearman']],on='replicate',suffixes=('_periscope','_cellot'),validate='one_to_one')
 assert len(paired)==1000
 delta=(paired.spearman_periscope-paired.spearman_cellot).dropna()
 append_table(P/'periscope_comparator_bootstrap.csv',[dict(comparator='cellot',delta_spearman_ci_low=float(delta.quantile(.025)),delta_spearman_ci_high=float(delta.quantile(.975)),bootstrap_fraction_periscope_higher=float((delta>0).mean()))],1)
 evaluation=[]
 for (readout,method,axis),group in df.groupby(['readout','method','axis']):
  for item,names in [(1,['mmd','energy_per_sqrt_gene','variance_ratio']),(2,['mean_mse']),(3,['condition_aware_rho','condition_aware_delta_mse']),(4,['fixed_control_rho','fixed_PD_rho'])]:
   for name in names:evaluation.append(dict(evaluation_item=item,biological_axis=axis,cohort='human_PD',method=method,readout=readout,metric=name,estimate=group[name].mean(),seed_sd=group[name].std(),n_seeds=3,interpretation='blood contribution diagnostic' if item==4 else 'prediction',source='baseline_metrics.csv'))
 oldpig=pd.read_csv(P/'primary_summary.csv').query("method == 'cellot'").iloc[0]
 for name in ['spearman','response_RMSE','improvement_over_zero','exact_group_p_positive','within_group_pairing_p']:
  evaluation.append(dict(evaluation_item=5,biological_axis='bulk_external_cohort',cohort='pig_LPS',method='cellot',readout='original_cellot',metric=name,estimate=oldpig[name],n_seeds=1,interpretation='external bulk response',source='../pig_external_validation_20261003/primary_summary.csv'))
 append_table(B/'benchmark_evaluation_summary.csv',evaluation,6)
 # Refresh hashes only for files that were changed or added by this installation.
 manifest=B/'manifest.json';obj=json.loads(manifest.read_text())
 for name in ['baseline_metrics.csv','seed_summary.csv','axis_summary.csv','comparison_summary.csv','benchmark_evaluation_summary.csv','training_summary.csv']:
  if name in obj.get('core_files',{}):obj['core_files'][name]=r.b.sha(B/name)
 obj['cellot_recovery_validation_sha256']=r.b.sha(B/'cellot_restored/recovery_validation.json');manifest.write_text(json.dumps(obj,indent=2));installed.append(str(manifest.relative_to(B.parent)))
 manifest=P/'output_manifest.json';obj=json.loads(manifest.read_text())
 for name in installed:
  if name.startswith(P.name+'/'):
   rel=name[len(P.name)+1:];obj['files'][rel]=r.b.sha(P/rel)
 manifest.write_text(json.dumps(obj,indent=2));installed.append(str(manifest.relative_to(B.parent)))
 (R/'installed_files_sha256.json').write_text(json.dumps({x:r.b.sha(B.parent/x) for x in sorted(set(installed))},indent=2))
 (R/'installed_files.txt').write_text('\n'.join(sorted(set(installed)))+'\n')
 (R/'installation_status.json').write_text(json.dumps(dict(complete=True,installed_files=len(set(installed)),only_cellot_results_added=True,existing_other_rows_preserved=True),indent=2))
 print('CELLOT_INSTALLED',len(set(installed)),flush=True)
if __name__=='__main__':main()
