"""Connect full assay coverage to all original model antigens and preserve ten-seed dose summaries."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2]
OLD=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
OUT=ROOT/'papers/03_world_model_integration/results/cross_domain_inputs_20261117'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(n,j):
 f=OUT/n;t=f.with_suffix(f.suffix+'.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(f)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'RUN.json').exists():raise RuntimeError('Existing run; inspect before restarting')
 run=dict(status='running',pid=os.getpid(),started_at_utc=datetime.now(timezone.utc).isoformat(),complete=False);dump('RUN.json',run)
 mreg=json.loads((OLD/'model_registry.json').read_text());registry=json.loads((OLD/'candidate_registry.json').read_text());dose=json.loads((OLD/'dose_definitions.json').read_text())
 assay_path=ROOT/'data/ukb_ppp/extracted/assay_coverage.csv';assays=pd.read_csv(assay_path,keep_default_na=False);assert len(registry)==217 and len(assays)==2941 and assays.exposure_id.is_unique
 sources={str(p):sha(p) for p in [OLD/'model_registry.json',OLD/'candidate_registry.json',OLD/'dose_definitions.json',assay_path,Path(__file__)]}
 # Exact author-supplied gene symbols only; no assumption of plasma/ADT assay equivalence.
 links=[];matched=set()
 for index,entry in enumerate(registry):
  genes=list(entry['genes']);hits=assays[assays.gene.isin(genes)] if not entry['is_control'] else assays.iloc[:0]
  base=dict(model_target_index=index,model_target=entry['target'],model_genes=';'.join(genes),model_mapping=entry['mapping'],is_control=entry['is_control'],model_channels=';'.join(map(str,entry['channels'])),n_assay_gene_matches=len(hits),identity_evidence='exact_author_gene_symbol_only',units_equivalent=False)
  if len(hits):
   for a in hits.to_dict('records'):
    links.append(dict(**base,exposure_id=a['exposure_id'],UKB_gene=a['gene'],UKB_uniprot=a['uniprot'],UKB_assay=a['assay'],has_cis_primary=a['has_cis_primary']));matched.add(a['exposure_id'])
  else:links.append(dict(**base,exposure_id='',UKB_gene='',UKB_uniprot='',UKB_assay='',has_cis_primary=''))
 mapping=pd.DataFrame(links);mapping.to_csv(OUT/'all_model_target_assay_links.csv',index=False)
 assaycoverage=assays.copy();assaycoverage['model_exact_gene_label_match']=assaycoverage.exposure_id.isin(matched);assaycoverage.to_csv(OUT/'all_2941_assay_model_coverage.csv',index=False)
 pd.DataFrame(registry).to_json(OUT/'original_all_217_model_targets.json',orient='records',indent=2)
 # The two astro targets remain separate and unresolved by the original astro model state.
 pd.DataFrame([
  {'paper2_state':'astro_fibrous','model_brain_state':'astro','relationship':'model_state_not_subtype_resolved','output_compatible':False},
  {'paper2_state':'astro_protoplasmic','model_brain_state':'astro','relationship':'model_state_not_subtype_resolved','output_compatible':False},
  {'paper2_state':'microglia_homeostatic','model_brain_state':'microglia_homeostatic','relationship':'same_label_different_measurement','output_compatible':False},
  {'paper2_state':'microglia_mhc2','model_brain_state':'microglia_mhc2','relationship':'same_label_different_measurement','output_compatible':False}
 ]).to_csv(OUT/'original_four_state_crosswalk.csv',index=False)
 curves=sorted((OLD/'curves').glob('*__seed*.npz'));assert len(curves)==120
 rows=[];identities=[];genes_reference=None;all_axes=set()
 for number,path in enumerate(curves,1):
  axis,seedtxt=path.stem.rsplit('__seed',1);seed=int(seedtxt);assert seed in mreg['seeds'];all_axes.add(axis)
  provenance=json.loads(path.with_suffix('.provenance.json').read_text());h=sha(path)
  assert provenance['output_sha256']==h and provenance['checkpoint_sha256']==mreg['models'][str(seed)]['sha256']
  assert provenance['candidate_registry_sha256']==sources[str(OLD/'candidate_registry.json')] and provenance['dose_definitions_sha256']==sources[str(OLD/'dose_definitions.json')]
  with np.load(path,allow_pickle=False) as z:
   C=z['C'].astype(np.float64);baseline=z['baseline'].astype(np.float64);u=z['u'].astype(np.float64);genes=z['genes'].astype(str)
   assert z['proteins'].astype(str).tolist()==[x['target'] for x in registry]
   assert np.array_equal(z['is_control'],np.array([x['is_control'] for x in registry]))
   assert int(z['seed'])==seed and np.array_equal(z['quantiles'],dose['quantiles'])
   assert C.shape==(217,5,len(genes)) and u.shape==(217,5) and np.isfinite(C).all() and np.isfinite(u).all() and np.isfinite(baseline).all()
   if genes_reference is None:genes_reference=genes
   else:assert np.array_equal(genes_reference,genes)
   centered=u-u.mean(1,keepdims=True);den=(centered**2).sum(1);delta=C[:,-1]-C[:,0]
   slopes=np.full((217,len(genes)),np.nan);good=den>0
   slopes[good]=np.einsum('tq,tqg->tg',centered[good],C[good])/den[good,None]
   for i,entry in enumerate(registry):
    rows.append(dict(seed=seed,axis=axis,model_target_index=i,model_target=entry['target'],is_control=entry['is_control'],genes=len(genes),dose_min=float(u[i,0]),dose_max=float(u[i,-1]),dose_coordinate_variance=float(den[i]),dose_constant=bool(not good[i]),delta_95_vs_05_RNA_RMS=float(np.sqrt(np.mean(delta[i]**2))),delta_95_vs_05_RNA_mean=float(delta[i].mean()),slope_per_CLR_RMS=float(np.sqrt(np.mean(slopes[i]**2))) if good[i] else np.nan,slope_per_CLR_mean=float(slopes[i].mean()) if good[i] else np.nan,baseline_RNA_RMS=float(np.sqrt(np.mean(baseline**2)))))
  identities.append(dict(path=str(path),sha256=h,seed=seed,axis=axis,provenance_sha256=sha(path.with_suffix('.provenance.json'))))
  run.update(verified_curve_files=number);dump('RUN.json',run)
  if number%12==0:print('VERIFIED_CURVES',number,flush=True)
 assert len(all_axes)==12 and len(rows)==26040
 frame=pd.DataFrame(rows);frame.to_csv(OUT/'all_targets_dose_response_by_seed.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 summary=frame.groupby(['axis','model_target_index','model_target','is_control'],sort=False).agg(n_seeds=('seed','nunique'),mean_response_RMS=('delta_95_vs_05_RNA_RMS','mean'),min_response_RMS=('delta_95_vs_05_RNA_RMS','min'),max_response_RMS=('delta_95_vs_05_RNA_RMS','max'),mean_slope_per_CLR_RMS=('slope_per_CLR_RMS','mean'),n_constant_dose=('dose_constant','sum')).reset_index();assert (summary.n_seeds==10).all();summary.to_csv(OUT/'all_targets_dose_response_across_seeds.csv',index=False)
 pd.Series(genes_reference).to_csv(OUT/'model_gene_universe.csv',header=['gene'],index=False)
 for f,h in sources.items():assert sha(Path(f))==h,f
 dump('SOURCE_IDENTITIES.json',dict(files=sources,curves=identities))
 run.update(status='complete',complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),all_model_targets=217,all_UKB_assays=2941,UKB_assays_with_exact_model_gene_label=len(matched),model_targets_with_exact_UKB_gene_label=int(mapping.loc[mapping.exposure_id!='','model_target_index'].nunique()),curve_files=120,seed_axis_target_rows=len(frame),axes=sorted(all_axes),scientific_scope='All original model targets retained; descriptive model-space gene-expression response. Exact gene labels do not equate plasma pQTL and cellular ADT assays. No new MR branch, candidate decisions, state-proportion or MRI predictions.')
 dump('RUN.json',run);dump('MANIFEST.json',{'artifacts_sha256':{str(p.relative_to(ROOT)):sha(p) for p in OUT.iterdir() if p.is_file() and p.name not in ['MANIFEST.json','controller.log']}})
if __name__=='__main__':main()
