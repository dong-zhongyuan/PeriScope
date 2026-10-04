"""Validate and publish only the shared-representation comparison."""
from pathlib import Path
import json,hashlib,csv
import numpy as np,pandas as pd
import current_benchmark as b
from shared_pipeline import METHODS,D

def main():
 # Preserve the restored original CellOT records when refreshing shared methods.
 cellot_rows={}
 for name in ['baseline_metrics.csv','seed_summary.csv','comparison_summary.csv','axis_summary.csv']:
  path=b.O/name
  cellot_rows[name]=[line for line in path.read_bytes().splitlines(keepends=True) if len(next(csv.reader([line.decode()])))>1 and next(csv.reader([line.decode()]))[1]=='cellot'] if path.exists() else []
 rows=[];checks=[]
 oldfile=b.O/"validation_checks.json";old=json.loads(oldfile.read_text()) if oldfile.exists() else {};prior={r["path"]:r for r in old.get("outputs",[]) if "path" in r}
 for readout,root in [('raw_shared',b.O),('validation_calibrated_shared',D)]:
  for seed in [42,43,44]:
   for axis in b.AXES:
    k=dict(np.load(b.O/'cache'/f'{axis}.npz'))
    for m in METHODS:
     p=root/'results'/f'seed_{seed}'/f'{m}__{axis}.json';r=json.loads(p.read_text());z=dict(np.load(p.with_suffix('.npz')))
     assert all(a.shape==(128,len(k['genes'])) and np.isfinite(a).all() for a in z.values())
     donor=pd.DataFrame(r['distribution_by_donor']);balanced=donor.groupby('condition')[['mmd','mean_mse','energy_per_sqrt_gene','variance_ratio']].mean().mean()
     previous=prior.get(str(p.relative_to(b.O)),{})
     unchanged=previous.get('sha256')==b.sha(p) and previous.get('prediction_sha256')==b.sha(p.with_suffix('.npz'))
     if not unchanged:
      for c in [0,1]:
       for j,don in enumerate(k[f'donors{c}']):
        got=b.measures(z[f'map{c}_blood{c}'],k[f'obs{c}'][j],float(k['bw']));orig=donor[donor.condition.eq(c)&donor.donor.eq(int(don))].iloc[0]
        for key,val in got.items():assert np.isclose(val,orig[key],rtol=1e-5,atol=1e-8)
     row=dict(readout=readout,method=m,seed=seed,axis=axis,blood_state=axis.split('__')[0],brain_state=axis.split('__')[1],**balanced.to_dict())
     for task,v in r['direction'].items():row.update({task+'_'+n:x for n,x in v.items()})
     rows.append(row);checks.append(dict(path=str(p.relative_to(b.O)),sha256=b.sha(p),prediction_sha256=b.sha(p.with_suffix('.npz'))))
 df=pd.DataFrame(rows);df.to_csv(b.O/'baseline_metrics.csv',index=False);metrics=['mmd','mean_mse','variance_ratio','condition_aware_rho','condition_aware_delta_mse','fixed_control_rho','fixed_PD_rho']
 perseed=df.groupby(['readout','method','seed'])[metrics].mean();perseed.to_csv(b.O/'seed_summary.csv');perseed.groupby(['readout','method']).agg(['mean','std']).to_csv(b.O/'comparison_summary.csv');df.groupby(['readout','method','axis'])[metrics].agg(['mean','std']).to_csv(b.O/'axis_summary.csv')
 for name,lines in cellot_rows.items():
  if lines:
   with (b.O/name).open('ab') as output:output.write(b''.join(lines))
 choices=pd.concat([pd.read_csv(D/f'validation_choices_seed{s}.csv') for s in [42,43,44]]);assert len(choices)==144
 assert (choices.validation_mean_mse_after<=choices.validation_mean_mse_before+1e-10).all();assert (choices.validation_mmd_after<=choices.validation_mmd_no_residual+1e-10).all();assert choices.contrast_error.max()<1e-10
 choices.to_csv(D/'validation_choices.csv',index=False)
 pairs=[]
 for readout,g in df.groupby('readout'):
  per=g[g.method.eq('periscope')]
  for m in METHODS:
   if m=='periscope':continue
   q=per.merge(g[g.method.eq(m)],on=['seed','axis'],suffixes=('_periscope','_other'));q['comparator']=m;q['readout']=readout;pairs.append(q)
 pd.concat(pairs).to_csv(b.O/'paired_comparisons.csv',index=False)
 (b.O/'validation_checks.json').write_text(json.dumps(dict(status='passed',shared_methods=METHODS,primary_rows=len(df),unique_validated_outputs=len(checks),original_prediction_rows=144,calibrated_prediction_rows=144,validation_only_selection=True,maximum_contrast_change=float(choices.contrast_error.max()),outputs=checks),indent=2))
 report=['共享表示Benchmark（2026-10-03）','', '正式方法为PeriScope、共享Ridge、条件MLP、条件ICNN；另有共享脑编码器/解码器的无血液均值与分布对照。全部基因映射经过相同冻结表示与脑解码器。', '8条轴、3个种子；每种方法保留原始读出和统一验证集校准读出。校准以训练集计算条件无关偏差和重构残差，所有方法使用同一参数网格和验证集；测试集不参与选择。此为观察初轮结果后的方法修订。', '原始均值读出与分布采样分别报告；校准采用零均值配对噪声，疾病响应差异保持不变。表达指标按PD/control各1/2、组内供者等权汇总；种子SD是训练波动。', '共享表示比较控制了编码器和解码器，不等于所有方法的训练预算和目标完全相同；PeriScope仍是原完整多任务模型。条件ICNN是适配的条件势函数基线。', '', 'method | readout | MMD | mean MSE | disease rho']
 for (rd,m),g in df.groupby(['readout','method']):
  z=g[metrics].mean();report.append(f'{m} | {rd} | {z.mmd:.6f} | {z.mean_mse:.6f} | {z.condition_aware_rho:.6f}')
 report+=['','固定条件指标为血液贡献/疾病方向对齐诊断，不是实测靶点干预精度。外部猪队列独立报告，不与人类八轴合成总分。MR是靶点筛选门槛，不是benchmark指标。']
 (b.O/'评估与运行说明.txt').write_text('\n'.join(report)+'\n')
 (b.O/'benchmark_scope_index.json').write_text(json.dumps(dict(formal_axes=b.AXES,n_formal_benchmark_biological_axes=8,methods=METHODS,comparison_seeds=[42,43,44],readouts=['raw_shared','validation_calibrated_shared'],n_primary_method_axis_seed_results=144,n_readout_results=288,shared_representation=True,evaluation_items=['expression distribution','mean expression','disease response','blood contribution diagnostic','independent pig cohort'],no_pooled_five_axis_score=True,source_directory='benchmark_shared_20261003'),indent=2))
 (b.O/'run_status.json').write_text(json.dumps(dict(complete=True,methods=METHODS,shared_representation=True,validated_rows=288),indent=2))
 print(df.groupby(['readout','method'])[metrics[:5]].mean().to_string(),flush=True)
if __name__=='__main__':main()
