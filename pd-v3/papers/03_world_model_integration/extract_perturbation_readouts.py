"""Extract public perturbation responses; preserve unmatched intervention domains explicitly."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,gzip,hashlib,re
import numpy as np,pandas as pd
from scipy import sparse,io
ROOT=Path(__file__).resolve().parents[2]
SRC=ROOT/'data/perturbation_sources_20261117/GSE311359'
OUT=ROOT/'papers/03_world_model_integration/results/perturbation_readout_20261117'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def now():return datetime.now(timezone.utc).isoformat()
def dump(n,v):
 p=OUT/n;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
def load_sample(sample):
 fs=list(SRC.glob(sample+'_*_features.tsv.gz'));assert len(fs)==1
 prefix=fs[0].name.removesuffix('features.tsv.gz')
 features=pd.read_csv(fs[0],sep='\t',header=None,names=['id','gene','type'],keep_default_na=False)
 barcodes=pd.read_csv(SRC/(prefix+'barcodes.tsv.gz'),sep='\t',header=None)[0].astype(str).to_numpy()
 with gzip.open(SRC/(prefix+'matrix.mtx.gz'),'rb') as f:
  matrix=io.mmread(f).tocsr()
 assert matrix.shape==(len(features),len(barcodes)) and len(set(barcodes))==len(barcodes)
 assert np.isfinite(matrix.data).all() and (matrix.data>=0).all() and np.equal(matrix.data,np.floor(matrix.data)).all()
 matrix.sum_duplicates();matrix.eliminate_zeros();matrix.sort_indices()
 return features,barcodes,matrix
def main():
 if (OUT/'RUN.json').exists():raise RuntimeError('Inspect existing extraction before restart')
 acquisition=json.loads((SRC/'COUNT_ACQUISITION.json').read_text());assert acquisition['complete'] and len(acquisition['files'])==21
 spec_hash=sha(OUT/'ANALYSIS_SPECIFICATION.json')
 run=dict(status='running',complete=False,pid=os.getpid(),started_at_utc=now(),samples_completed=0,analysis_specification_sha256=spec_hash)
 dump('RUN.json',run)
 expected={Path(x['file']).name:x for x in acquisition['files']}
 sampleids=sorted({x['sample'] for x in acquisition['files']});summary=[];coverage=[];identities=[]
 (OUT/'samples').mkdir(exist_ok=True)
 reference=None
 for sample in sampleids:
  for path in SRC.glob(sample+'_*'):
   if path.name.endswith(('.tsv.gz','.mtx.gz')):assert sha(path)==expected[path.name]['sha256']
  feat,bar,m=load_sample(sample)
  if reference is None:reference=feat.copy();feat.to_csv(OUT/'all_source_features.csv',index=False)
  else:assert feat.equals(reference)
  rnaidx=np.flatnonzero(feat.type.to_numpy()=='Gene Expression');gidx=np.flatnonzero(feat.type.to_numpy()=='CRISPR Guide Capture')
  assert len(gidx)==381 and len(rnaidx)+len(gidx)==len(feat)
  rna=m[rnaidx].T.tocsr().astype(np.float64);guides=m[gidx].T.tocsr()
  gnames=feat.iloc[gidx].gene.to_numpy().astype(str);rnames=feat.iloc[rnaidx].gene.to_numpy().astype(str)
  total=np.asarray(rna.sum(axis=1)).ravel();detected=np.diff(rna.indptr)
  mito=np.asarray(rna[:,np.char.startswith(rnames,'MT-')].sum(axis=1)).ravel()
  mitofrac=np.divide(mito,total,out=np.full_like(total,np.nan),where=total>0)
  complexity=np.full_like(total,np.nan);ok=(total>1)&(detected>0);complexity[ok]=np.log10(detected[ok])/np.log10(total[ok])
  qc=(detected>=250)&(total>=500)&(mitofrac<=.2)&(complexity>=.8)
  strict=(detected>250)&(total>500)&(mitofrac<.2)&(complexity>.8)
  positive=guides.copy();positive.data=(positive.data>=3).astype(np.int8);positive.eliminate_zeros()
  positive=positive.astype(np.int32);npos=np.asarray(positive.sum(axis=1)).ravel()
  isnt=np.char.startswith(gnames,'non-targeting_');assert isnt.sum()==17
  nnt=np.asarray(positive[:,isnt].sum(axis=1)).ravel();ntc=(npos>0)&(npos==nnt)
  singleguide=np.full(len(bar),'',dtype=object)
  one=np.flatnonzero(npos==1)
  singleguide[one]=gnames[positive[one].indices]
  pd.DataFrame(dict(sample=sample,barcode=bar,RNA_UMI=total.astype(np.int64),RNA_genes=detected,mitochondrial_fraction=mitofrac,log10_genes_per_log10_UMI=complexity,manuscript_QC=qc,GEO_strict_QC=strict,positive_guides=npos,positive_NTC_guides=nnt,NTC_only=ntc,single_positive_guide=singleguide)).to_csv(OUT/'samples'/f'{sample}_cell_coverage.csv.gz',index=False,compression={'method':'gzip','mtime':0})
  tss={i:re.fullmatch(r'(.+)_tss_g([0-9]+)',g) for i,g in enumerate(gnames)};tss={i:m.group(1) for i,m in tss.items() if m}
  assert len(tss)==31 and len(set(tss.values()))==10
  # Membership uses guides only, never measured RNA response or target knockdown.
  basegroups=[('NTC','NTC','',ntc)]
  for i in range(len(gnames)):
   if i in tss or isnt[i]:
    basegroups.append(('guide',gnames[i],tss.get(i,''),np.asarray(positive[:,i].toarray()).ravel()>0 if i in tss else (ntc & (np.asarray(positive[:,i].toarray()).ravel()>0))))
  for gene in sorted(set(tss.values())):
   cols=[i for i,g in tss.items() if g==gene]
   basegroups.append(('gene',gene,gene,np.asarray(positive[:,cols].sum(axis=1)).ravel()>0))
  groupids=[];groupmeta=[];members=[]
  for scheme in ['author_inclusive','single_positive_guide']:
   for kind,label,gene,mask in basegroups:
    ids=np.flatnonzero(qc & mask & ((npos==1) if scheme=='single_positive_guide' else True))
    gid=scheme+'|'+kind+'|'+label;groupids.append(gid);members.append(ids)
    groupmeta.append(dict(sample=sample,group_id=gid,assignment=scheme,kind=kind,label=label,target_gene=gene,n_cells=len(ids),co_guide_cells=int((npos[ids]>1).sum()),GEO_boundary_cells=int((~strict[ids]).sum())))
  row=np.concatenate([np.full(len(ids),i,dtype=np.int32) for i,ids in enumerate(members)])
  col=np.concatenate(members)
  A=sparse.coo_matrix((np.ones(len(row)),(row,col)),shape=(len(members),len(bar))).tocsr()
  raw=np.asarray((A@rna).toarray());n=np.array([len(x) for x in members],dtype=np.int64)
  norm=rna.copy();norm.data*=np.repeat(np.divide(10000,total,out=np.zeros_like(total),where=total>0),np.diff(norm.indptr));norm.data=np.log1p(norm.data)
  sums=np.asarray((A@norm).toarray())
  mean=np.divide(sums,n[:,None],out=np.full_like(sums,np.nan),where=n[:,None]>0)
  detect=np.asarray((A@(rna>0)).toarray()).astype(np.int64)
  assert np.isfinite(raw).all() and (raw>=0).all() and np.equal(raw,np.floor(raw)).all()
  np.savez_compressed(OUT/'samples'/f'{sample}_readouts.npz',group_id=np.array(groupids),n_cells=n,raw_pseudobulk=raw.astype(np.int64),mean_log1p_10k=mean,detected_cells=detect,gene=feat.iloc[rnaidx].gene.to_numpy(dtype=str),gene_id=feat.iloc[rnaidx].id.to_numpy(dtype=str))
  summary.extend(groupmeta)
  coverage.append(dict(sample=sample,input_barcodes=len(bar),manuscript_QC=int(qc.sum()),GEO_strict_QC=int(strict.sum()),QC_any_guide=int((qc&(npos>0)).sum()),QC_no_guide=int((qc&(npos==0)).sum()),QC_single_guide=int((qc&(npos==1)).sum()),QC_multi_guide=int((qc&(npos>1)).sum()),QC_NTC_only=int((qc&ntc).sum()),RNA_features=len(rnaidx),guide_features=len(gidx),input_nonzero=m.nnz))
  identities.append(dict(sample=sample,outputs_sha256={x.name:sha(x) for x in (OUT/'samples').glob(sample+'_*')}))
  run['samples_completed']+=1;dump('RUN.json',run);print('COMPLETE_SAMPLE',sample,coverage[-1],flush=True)
 pd.DataFrame(summary).to_csv(OUT/'group_coverage.csv',index=False);pd.DataFrame(coverage).to_csv(OUT/'sample_coverage.csv',index=False)
 guide_table=[]
 for i,name in enumerate(gnames):
  match=re.fullmatch(r'(.+)_tss_g([0-9]+)',name)
  guide_table.append(dict(guide=name,kind='TSS' if match else 'NTC' if name.startswith('non-targeting_') else 'other_locus_guide',target_gene=match.group(1) if match else '',gene_identity_evidence='author_feature_TSS_label' if match else 'no_direct_gene_assignment'))
 pd.DataFrame(guide_table).to_csv(OUT/'all_381_guide_identity.csv',index=False)
 registry=json.loads((ROOT/'papers/03_world_model_integration/results/cross_domain_inputs_20261117/original_all_217_model_targets.json').read_text())
 pd.DataFrame([dict(model_target=t['target'],model_genes=';'.join(t['genes']),exact_TSS_gene_overlap=';'.join(sorted(set(t['genes'])&set(tss.values()))),is_control=t['is_control'],intervention_cell_matches=False,recipient_subtype_resolved=False) for t in registry]).to_csv(OUT/'all_217_model_target_coverage.csv',index=False)
 assert sha(OUT/'ANALYSIS_SPECIFICATION.json')==spec_hash
 run.update(status='complete',complete=True,finished_at_utc=now(),all_source_samples=len(sampleids),all_source_barcodes=sum(x['input_barcodes'] for x in coverage),source_QC_guide_positive_cells=sum(x['QC_any_guide'] for x in coverage),source_author_reported_QC_guide_positive_cells=65906,source_author_reported_NTC_cells=1021,observed_NTC_cells=sum(x['QC_NTC_only'] for x in coverage),model_validation_complete=False)
 dump('SOURCE_IDENTITIES.json',dict(acquisition_sha256=sha(SRC/'COUNT_ACQUISITION.json'),source_files=expected,sample_outputs=identities,script_sha256=sha(Path(__file__))))
 dump('RUN.json',run)
 dump('MANIFEST.json',{'artifacts_sha256':{str(p.relative_to(ROOT)):sha(p) for p in OUT.glob('*') if p.is_file() and p.name not in ['MANIFEST.json','extraction.log']}})
if __name__=='__main__':main()
