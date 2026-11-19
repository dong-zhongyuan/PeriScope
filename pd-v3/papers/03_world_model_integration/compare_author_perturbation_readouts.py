"""Use all reported author CSF1R RNA effects, without significance filtering, for direction comparison."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,os
import numpy as np,pandas as pd
from scipy.stats import rankdata,pearsonr,spearmanr
ROOT=Path(__file__).resolve().parents[2];SRC=ROOT/'data/perturbation_sources_20261117/GSE178317';OLD=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003');OUT=ROOT/'papers/03_world_model_integration/results/perturbation_author_readouts_20261117'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 OUT.mkdir(exist_ok=True)
 assert not (OUT/'RUN.json').exists()
 specification=dict(study='GSE178317',source='Author Supplementary Table 9, all reported genes for every exact registry overlap; no added P/FDR/effect-size selection.',frozen_predictions='All original 8 microglial axes and 10 seeds, q05 minus q95 RNA response, no refit. Full shared gene set and sensitivity excluding perturbed gene.',scope='Descriptive signed response agreement; microglial CRISPRi differs from blood-cell ADT model control; source has no matched model substate or donor-heldout estimate. This does not complete causal blood-to-brain validation.',target_selection='Exact intersection of Supplementary Table 5 target genes and original all-217 antigen genes, before reading effect values.')
 (OUT/'ANALYSIS_SPECIFICATION.json').write_text(json.dumps(specification,indent=2)+'\n')
 coverage=pd.read_csv(SRC/'all_217_target_coverage.csv',keep_default_na=False);match=coverage[coverage.exact_author_guide_gene_overlap!=''];assert len(match)==1 and match.iloc[0].model_target=='CD115' and match.iloc[0].exact_author_guide_gene_overlap=='CSF1R'
 registry=json.loads((OLD/'candidate_registry.json').read_text());targetindex=next(i for i,t in enumerate(registry) if t['target']=='CD115')
 source_table=SRC/'all_author_reported_RNA_effects.csv.gz';table=pd.read_csv(source_table,float_precision='round_trip');obs=table[table.TargetGene=='CSF1R'].copy();assert len(obs)==8610 and obs.Gene.is_unique
 obs.to_csv(OUT/'all_8610_author_CSF1R_RNA_effects.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 identity=json.loads((ROOT/'papers/03_world_model_integration/results/cross_domain_inputs_20261117/SOURCE_IDENTITIES.json').read_text());expected={Path(x['path']).name:x['sha256'] for x in identity['curves']}
 records=[];vectors=[];sources=[];largest=0.
 for p in sorted((OLD/'curves').glob('*__seed*.npz')):
  axis,seed=p.stem.rsplit('__seed',1)
  if '__x__microglia_' not in axis:continue
  assert sha(p)==expected[p.name];sources.append(dict(path=str(p),sha256=expected[p.name]))
  with np.load(p,allow_pickle=False) as z:
   genes=z['genes'].astype(str);delta=z['C'][targetindex,0].astype(float)-z['C'][targetindex,-1].astype(float)
  shared=obs.set_index('Gene').reindex(genes);mask=shared.log2FC.notna().to_numpy()
  if not vectors:
   pd.DataFrame(dict(model_gene=genes,reported_for_CSF1R=mask)).to_csv(OUT/'all_model_gene_coverage.csv',index=False)
   shared.loc[mask].to_csv(OUT/'shared_gene_observed_effects.csv',index_label='gene')
  vectors.append(delta)
  for scope in ['all_shared_genes','exclude_perturbed_gene']:
   use=mask if scope=='all_shared_genes' else mask&(genes!='CSF1R');x=shared.log2FC.to_numpy()[use];y=delta[use]
   xc=x-x.mean();yc=y-y.mean();pr=float(xc@yc/np.sqrt((xc@xc)*(yc@yc)))
   rx=rankdata(x);ry=rankdata(y);rx-=rx.mean();ry-=ry.mean();sr=float(rx@ry/np.sqrt((rx@rx)*(ry@ry)))
   check_pr=pearsonr(x,y).statistic;check_sr=spearmanr(x,y).statistic;largest=max(largest,abs(pr-check_pr),abs(sr-check_sr));assert largest<1e-12
   nz=(x!=0)&(y!=0)
   records.append(dict(study='GSE178317',target_gene='CSF1R',model_target='CD115',axis=axis,seed=int(seed),gene_scope=scope,n_shared_genes=int(use.sum()),pearson=pr,spearman=sr,sign_agreement=float((np.sign(x[nz])==np.sign(y[nz])).mean()),joint_nonzero_genes=int(nz.sum()),observed_positive_fraction=float((x[x!=0]>0).mean()),intervention_cell_matched=False,recipient_subtype_matched=False))
 assert len(records)==160 and len(sources)==80
 frame=pd.DataFrame(records);frame.to_csv(OUT/'all_direction_comparisons.csv',index=False)
 summary=frame.groupby(['axis','gene_scope']).agg(seeds=('seed','nunique'),spearman_median=('spearman','median'),spearman_min=('spearman','min'),spearman_max=('spearman','max'),pearson_median=('pearson','median'),sign_agreement_median=('sign_agreement','median')).reset_index();summary.to_csv(OUT/'direction_comparison_summary.csv',index=False)
 np.savez_compressed(OUT/'frozen_model_decrease_responses.npz',gene=genes,low_minus_high=np.stack(vectors),axis=np.array([Path(x['path']).stem.rsplit('__seed',1)[0] for x in sources]),seed=np.array([int(Path(x['path']).stem.rsplit('__seed',1)[1]) for x in sources]))
 sources_meta=dict(workbook_sha256=sha(SRC/'41593_2022_1131_MOESM4_ESM.xlsx'),author_extraction_sha256=sha(SRC/'AUTHOR_TABLE_EXTRACTION.json'),reported_RNA_table_sha256=sha(source_table),curves=sources,script_sha256=sha(Path(__file__)))
 (OUT/'SOURCE_IDENTITIES.json').write_text(json.dumps(sources_meta,indent=2)+'\n')
 report=dict(complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),source_RNA_rows=8610,unique_source_genes=8610,all_217_model_targets_retained_in_source_coverage=True,comparison_rows=160,all_80_curve_hashes_verified=True,independent_scalar_formula_max_difference=largest,model_validation_complete=False)
 (OUT/'RUN.json').write_text(json.dumps(report,indent=2)+'\n')
 (OUT/'MANIFEST.json').write_text(json.dumps({'files_sha256':{p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name not in ['MANIFEST.json','comparison.log']}},indent=2)+'\n')
 print(json.dumps(report),flush=True);print(summary.to_string(index=False),flush=True)
if __name__=='__main__':main()
