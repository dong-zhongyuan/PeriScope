"""Verify every direction-comparison metric with separate scalar formulas."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import pearsonr,spearmanr
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'papers/03_world_model_integration/results/perturbation_readout_20261117'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def check(v,w):
 if np.isnan(v) or np.isnan(w):assert np.isnan(v) and np.isnan(w);return 0.
 d=abs(float(v)-float(w));assert d<1e-11,(v,w,d);return d
def main():
 assert json.loads((O/'COMPARISON_RUN.json').read_text())['complete']
 src=json.loads((O/'COMPARISON_SOURCES.json').read_text())
 with np.load(O/'frozen_model_decrease_responses.npz',allow_pickle=False) as z:
  genes=z['gene'].astype(str);keys=list(zip(z['target'].astype(str),z['axis'].astype(str),z['seed'].tolist()));pred=z['low_minus_high'].copy()
 registry=json.loads((ROOT/'papers/03_world_model_integration/results/cross_domain_inputs_20261117/original_all_217_model_targets.json').read_text())
 lookup={(t['genes'][0] if len(t['genes'])==1 else ''):i for i,t in enumerate(registry)}
 idx={key:i for i,key in enumerate(keys)}
 for item in src['model_curves']:
  p=Path(item['path']);assert sha(p)==item['sha256'];axis,seed=p.stem.rsplit('__seed',1)
  with np.load(p,allow_pickle=False) as z:
   assert np.array_equal(z['genes'].astype(str),genes)
   for target in ['CD86','MERTK','IFNGR1']:
    j=lookup[target];v=z['C'][j].astype(np.float64)
    assert np.array_equal(v[0]-v[4],pred[idx[(target,axis,int(seed))]])
 rows=pd.read_csv(O/'all_direction_comparisons.csv.gz',float_precision='round_trip');assert len(rows)==26880
 maxdiff=0.;verified=0
 for item in src['readouts']:
  p=Path(item['path']);assert sha(p)==item['sha256'];sample=p.name.split('_')[0]
  with np.load(p,allow_pickle=False) as z:
   sourcegenes=z['gene'].astype(str);counts=pd.Series(sourcegenes).value_counts();unique={g:i for i,g in enumerate(sourcegenes) if counts[g]==1}
   common=[(i,unique[g],g) for i,g in enumerate(genes) if g in unique];gm={v:i for i,v in enumerate(z['group_id'].astype(str))}
   for key,block in rows[rows['sample']==sample].groupby(['assignment','kind','label','target_gene','gene_scope'],sort=False):
    assignment,kind,label,target,scope=key;gi=gm[assignment+'|'+kind+'|'+label];ci=gm[assignment+'|NTC|NTC']
    use=[(i,j) for i,j,g in common if scope=='all_common_genes' or g!=target];mi=[i for i,j in use];ei=[j for i,j in use]
    x=z['mean_log1p_10k'][gi,ei]-z['mean_log1p_10k'][ci,ei]
    for row in block.itertuples():
     assert row.n_cells==z['n_cells'][gi] and row.n_control_cells==z['n_cells'][ci] and row.n_common_genes==len(use)
     y=pred[idx[(target,row.axis,int(row.seed))]][mi]
     if row.n_cells==0 or row.n_control_cells==0:
      assert np.isnan(row.spearman) and np.isnan(row.pearson);verified+=1;continue
     pr=pearsonr(x,y).statistic if np.ptp(x)>0 and np.ptp(y)>0 else np.nan
     sr=spearmanr(x,y).statistic if np.ptp(x)>0 and np.ptp(y)>0 else np.nan
     denom=np.sqrt(np.sum(x*x)*np.sum(y*y));cos=np.sum(x*y)/denom if denom else np.nan
     b=(x!=0)&(y!=0);agree=float(((x[b]>0)==(y[b]>0)).sum()/b.sum()) if b.any() else np.nan
     for a,v in [(pr,row.pearson),(sr,row.spearman),(cos,row.cosine),(agree,row.sign_agreement)]:maxdiff=max(maxdiff,check(a,v))
     assert (x!=0).sum()==row.observed_nonzero_genes and (y!=0).sum()==row.model_nonzero_genes and b.sum()==row.joint_nonzero_genes
     verified+=1
  print('VERIFIED_COMPARISON',sample,verified,flush=True)
 keys=['assignment','kind','label','target_gene','model_target','axis','gene_scope']
 sums=pd.read_csv(O/'direction_comparison_summary.csv',float_precision='round_trip').set_index(keys)
 for key,group in rows.groupby(keys):
  s=sums.loc[key];assert s.comparisons==len(group) and s.defined_comparisons==group.spearman.notna().sum() and s.wells==group['sample'].nunique() and s.model_seeds==group.seed.nunique()
  for col,fun in [('spearman_median',np.nanmedian),('spearman_min',np.nanmin),('spearman_max',np.nanmax),('pearson_median',np.nanmedian),('sign_agreement_median',np.nanmedian)]:
   source=col.rsplit('_',1)[0];maxdiff=max(maxdiff,check(fun(group[source].to_numpy()),s[col]))
 assert verified==len(rows)
 report=dict(complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),comparison_rows=verified,summary_rows=len(sums),all_80_source_curves_exact=True,maximum_metric_difference=maxdiff,independent_functions='scipy pearsonr/spearmanr versus production centered dot/rank formulas; separate source membership reconstruction',model_validation_complete=False,script_sha256=sha(Path(__file__)),comparison_sha256=sha(O/'all_direction_comparisons.csv.gz'))
 (O/'INDEPENDENT_COMPARISON_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
