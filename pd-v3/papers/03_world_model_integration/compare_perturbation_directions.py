"""Compare frozen model directions with public microglial knockdown readouts; no equivalence claim."""
from pathlib import Path
import json,hashlib,os
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import rankdata
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'papers/03_world_model_integration/results/perturbation_readout_20261117'
OLD=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(n,x):(OUT/n).write_text(json.dumps(x,indent=2)+'\n')
def metric(x,y):
 xc=x-x.mean();yc=y-y.mean();den=np.linalg.norm(xc)*np.linalg.norm(yc);rawden=np.linalg.norm(x)*np.linalg.norm(y)
 both=(x!=0)&(y!=0);rx=rankdata(x);ry=rankdata(y);rx-=rx.mean();ry-=ry.mean();rden=np.linalg.norm(rx)*np.linalg.norm(ry)
 return dict(pearson=float(xc@yc/den) if den else np.nan,spearman=float(rx@ry/rden) if rden else np.nan,cosine=float(x@y/rawden) if rawden else np.nan,observed_nonzero_genes=int((x!=0).sum()),model_nonzero_genes=int((y!=0).sum()),joint_nonzero_genes=int(both.sum()),sign_agreement=float(np.mean(np.sign(x[both])==np.sign(y[both]))) if both.any() else np.nan,observed_positive_fraction=float(np.mean(x[x!=0]>0)) if np.any(x!=0) else np.nan,model_positive_fraction=float(np.mean(y[y!=0]>0)) if np.any(y!=0) else np.nan)
def main():
 assert json.loads((OUT/'RUN.json').read_text())['complete']
 if (OUT/'COMPARISON_RUN.json').exists():raise RuntimeError('Inspect previous comparison')
 run=dict(status='running',complete=False,pid=os.getpid(),started_at_utc=datetime.now(timezone.utc).isoformat(),scope='Descriptive direction comparison across explicitly unmatched interventions. Not a validation-pass decision.')
 dump('COMPARISON_RUN.json',run)
 reg=json.loads((OLD/'candidate_registry.json').read_text());cov=pd.read_csv(OUT/'all_217_model_target_coverage.csv',keep_default_na=False)
 overlap=cov[cov.exact_TSS_gene_overlap!=''];assert set(overlap.exact_TSS_gene_overlap)=={'CD86','IFNGR1','MERTK'}
 links={r.exact_TSS_gene_overlap:r.model_target for r in overlap.itertuples()}
 targets=[i for i,t in enumerate(reg) if t['target'] in set(links.values())];assert len(targets)==3
 identities=json.loads((ROOT/'papers/03_world_model_integration/results/cross_domain_inputs_20261117/SOURCE_IDENTITIES.json').read_text())
 source={Path(x['path']).name:x['sha256'] for x in identities['curves']}
 predictions={};predrows=[];curvefiles=[]
 for p in sorted((OLD/'curves').glob('*__seed*.npz')):
  axis,seed=p.stem.rsplit('__seed',1)
  if '__x__microglia_' not in axis:continue
  assert sha(p)==source[p.name]
  with np.load(p,allow_pickle=False) as z:
   modelgenes=z['genes'].astype(str);assert len(set(modelgenes))==len(modelgenes)
   delta=z['C'][targets,0].astype(np.float64)-z['C'][targets,-1].astype(np.float64)
   for j,ti in enumerate(targets):
    gene=reg[ti]['genes'][0];predictions[(gene,axis,int(seed))]=delta[j]
    predrows.append(dict(gene=gene,model_target=reg[ti]['target'],axis=axis,seed=int(seed),dose_low=float(z['u'][ti,0]),dose_high=float(z['u'][ti,-1]),dose_constant=bool(z['u'][ti,0]==z['u'][ti,-1])))
  curvefiles.append(dict(path=str(p),sha256=source[p.name]))
 assert len(predictions)==240 and len(curvefiles)==80
 pd.DataFrame(predrows).to_csv(OUT/'model_comparison_coverage.csv',index=False)
 ordered=sorted(predictions)
 np.savez_compressed(OUT/'frozen_model_decrease_responses.npz',gene=modelgenes,target=np.array([k[0] for k in ordered]),axis=np.array([k[1] for k in ordered]),seed=np.array([k[2] for k in ordered]),low_minus_high=np.stack([predictions[k] for k in ordered]))
 metadata=pd.read_csv(OUT/'group_coverage.csv',keep_default_na=False)
 records=[];targeteffects=[];sampleidentities=[]
 for path in sorted((OUT/'samples').glob('*_readouts.npz')):
  sample=path.name.split('_')[0];sampleidentities.append(dict(path=str(path),sha256=sha(path)))
  with np.load(path,allow_pickle=False) as z:
   genes=z['gene'].astype(str);ids=z['gene_id'].astype(str);groups=z['group_id'].astype(str);n=z['n_cells'];values=z['mean_log1p_10k'];raw=z['raw_pseudobulk'];det=z['detected_cells']
   # Exact unique gene symbols only; ambiguous symbols remain in source coverage.
   unique=pd.Series(genes).value_counts();gmap={g:i for i,g in enumerate(genes) if unique[g]==1}
   common=[g for g in modelgenes if g in gmap];mi=np.array([np.flatnonzero(modelgenes==g)[0] for g in common]);ei=np.array([gmap[g] for g in common])
   assert len(common)>0
   if not (OUT/'comparison_gene_coverage.csv').exists():
    pd.DataFrame(dict(model_gene=modelgenes,in_source=[g in set(genes) for g in modelgenes],unique_source_symbol=[g in gmap for g in modelgenes])).to_csv(OUT/'comparison_gene_coverage.csv',index=False)
   for i,gid in enumerate(groups):
    scheme,kind,label=gid.split('|');target=label if kind=='gene' else label.rsplit('_tss_g',1)[0] if kind=='guide' and '_tss_g' in label else ''
    if target not in links:continue
    control=np.flatnonzero(groups==scheme+'|NTC|NTC');assert len(control)==1;ci=int(control[0])
    valid=n[i]>0 and n[ci]>0
    if target in gmap:
     t=gmap[target]
     targeteffects.append(dict(sample=sample,assignment=scheme,kind=kind,label=label,target_gene=target,model_target=links[target],n_cells=int(n[i]),n_control_cells=int(n[ci]),target_mean_log1p_10k=float(values[i,t]),control_mean_log1p_10k=float(values[ci,t]),target_gene_delta=float(values[i,t]-values[ci,t]) if valid else np.nan,target_detected_cells=int(det[i,t]),control_detected_cells=int(det[ci,t]),target_raw_count=int(raw[i,t]),control_raw_count=int(raw[ci,t])))
    observed=values[i,ei]-values[ci,ei]
    for scope in ['all_common_genes','exclude_perturbed_gene']:
     keep=np.ones(len(common),dtype=bool) if scope=='all_common_genes' else np.array([g!=target for g in common])
     x=observed[keep]
     for (gene,axis,seed),pred in predictions.items():
      if gene!=target:continue
      y=pred[mi][keep]
      result=metric(x,y) if valid and np.isfinite(x).all() else {k:np.nan for k in ['pearson','spearman','cosine','observed_nonzero_genes','model_nonzero_genes','joint_nonzero_genes','sign_agreement','observed_positive_fraction','model_positive_fraction']}
      records.append(dict(sample=sample,assignment=scheme,kind=kind,label=label,target_gene=target,model_target=links[target],axis=axis,seed=seed,gene_scope=scope,n_cells=int(n[i]),n_control_cells=int(n[ci]),n_common_genes=int(keep.sum()),comparison_defined=bool(valid),intervention_cell_matched=False,recipient_subtype_matched=False,**result))
  print('COMPARED',sample,len(records),flush=True)
 frame=pd.DataFrame(records);frame.to_csv(OUT/'all_direction_comparisons.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 pd.DataFrame(targeteffects).to_csv(OUT/'target_gene_readouts.csv',index=False)
 keys=['assignment','kind','label','target_gene','model_target','axis','gene_scope']
 summary=frame.groupby(keys).agg(comparisons=('spearman','size'),defined_comparisons=('spearman','count'),wells=('sample','nunique'),model_seeds=('seed','nunique'),spearman_median=('spearman','median'),spearman_min=('spearman','min'),spearman_max=('spearman','max'),pearson_median=('pearson','median'),sign_agreement_median=('sign_agreement','median')).reset_index()
 summary.to_csv(OUT/'direction_comparison_summary.csv',index=False)
 # Median/min/max describe technical wells and trained seeds; neither is a biological confidence interval.
 dump('COMPARISON_SOURCES.json',dict(model_curves=curvefiles,readouts=sampleidentities,analysis_specification_sha256=sha(OUT/'ANALYSIS_SPECIFICATION.json'),script_sha256=sha(Path(__file__))))
 run.update(status='complete',complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),source_samples=7,matched_genes=3,model_seeds=10,axes=8,comparison_rows=len(frame),summary_rows=len(summary),model_validation_complete=False)
 dump('COMPARISON_RUN.json',run)
if __name__=='__main__':main()
