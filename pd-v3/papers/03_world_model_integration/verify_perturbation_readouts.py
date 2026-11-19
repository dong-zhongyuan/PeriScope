"""Independent source-matrix reconstruction of every exported perturbation aggregate."""
from pathlib import Path
from datetime import datetime,timezone
import json,gzip,hashlib,os
import numpy as np,pandas as pd
from scipy.io import mmread
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'data/perturbation_sources_20261117/GSE311359';O=ROOT/'papers/03_world_model_integration/results/perturbation_readout_20261117'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 assert json.loads((O/'RUN.json').read_text())['complete']
 receipts=json.loads((D/'COUNT_ACQUISITION.json').read_text())['files'];expected={Path(r['file']).name:r for r in receipts}
 groups=pd.read_csv(O/'group_coverage.csv',keep_default_na=False);sample_records=[];globaldiff=0.;totalgroups=0
 for p in sorted((O/'samples').glob('*_readouts.npz')):
  sample=p.name.split('_')[0];ff=list(D.glob(sample+'_*_features.tsv.gz'));assert len(ff)==1
  prefix=ff[0].name.removesuffix('features.tsv.gz')
  for suffix in ['features.tsv.gz','barcodes.tsv.gz','matrix.mtx.gz']:assert sha(D/(prefix+suffix))==expected[prefix+suffix]['sha256']
  f=pd.read_csv(ff[0],header=None,sep='\t',keep_default_na=False)
  bc=pd.read_csv(D/(prefix+'barcodes.tsv.gz'),header=None)[0].astype(str).to_numpy()
  with gzip.open(D/(prefix+'matrix.mtx.gz'),'rb') as stream:matrix=mmread(stream).tocsc()
  genes=f.loc[f[2]=='Gene Expression',1].astype(str).to_numpy();geneids=f.loc[f[2]=='Gene Expression',0].astype(str).to_numpy()
  rna=matrix[(f[2]=='Gene Expression').to_numpy()].astype(np.float64).tocsc()
  g=matrix[(f[2]=='CRISPR Guide Capture').to_numpy()].toarray()
  names=f.loc[f[2]=='CRISPR Guide Capture',1].astype(str).to_numpy();pos=g>=3;ng=pos.sum(axis=0);nt=np.array([name.startswith('non-targeting_') for name in names])
  ntc=(pos[nt].sum(axis=0)>0)&(~pos[~nt].any(axis=0))
  sums=np.asarray(rna.sum(axis=0)).ravel();nfeatures=np.asarray((rna>0).sum(axis=0)).ravel();mt=np.array([s.startswith('MT-') for s in genes])
  mfrac=np.asarray(rna[mt].sum(axis=0)).ravel()/sums
  with np.errstate(divide='ignore',invalid='ignore'):complexity=np.log10(nfeatures)/np.log10(sums)
  mask=(nfeatures>=250)&(sums>=500)&(mfrac<=.2)&(complexity>=.8)
  strict=(nfeatures>250)&(sums>500)&(mfrac<.2)&(complexity>.8)
  cells=pd.read_csv(O/'samples'/f'{sample}_cell_coverage.csv.gz',float_precision='round_trip',keep_default_na=False)
  assert np.array_equal(cells.barcode.to_numpy(),bc)
  for column,v in [('RNA_UMI',sums),('RNA_genes',nfeatures),('positive_guides',ng),('positive_NTC_guides',pos[nt].sum(0)),('manuscript_QC',mask),('GEO_strict_QC',strict),('NTC_only',ntc)]:assert np.array_equal(cells[column].to_numpy(),v),(sample,column)
  single=np.array([names[np.flatnonzero(pos[:,i])[0]] if ng[i]==1 else '' for i in range(len(bc))])
  assert np.array_equal(cells.single_positive_guide.to_numpy(),single)
  # Independently scale columns of gene-by-cell CSC and use direct selected-cell reductions.
  norm=rna.copy();norm.data=np.log1p(norm.data*np.repeat(10000/sums,np.diff(norm.indptr)))
  with np.load(p,allow_pickle=False) as out:
   assert np.array_equal(out['gene'],genes) and np.array_equal(out['gene_id'],geneids)
   for i,gid in enumerate(out['group_id'].astype(str)):
    scheme,kind,label=gid.split('|')
    if kind=='NTC':membership=ntc
    elif kind=='gene':
     selected=np.array([n.startswith(label+'_tss_g') for n in names]);assert selected.any();membership=pos[selected].any(axis=0)
    else:
     selected=np.flatnonzero(names==label);assert len(selected)==1;membership=pos[selected[0]].copy()
     if label.startswith('non-targeting_'):membership &= ntc
    select=mask&membership
    if scheme=='single_positive_guide':select&=(ng==1)
    else:assert scheme=='author_inclusive'
    count=int(select.sum());assert out['n_cells'][i]==count
    metadata=groups[(groups['sample']==sample)&(groups.group_id==gid)];assert len(metadata)==1
    assert int(metadata.n_cells.iloc[0])==count
    assert int(metadata.co_guide_cells.iloc[0])==int((select&(ng>1)).sum())
    assert int(metadata.GEO_boundary_cells.iloc[0])==int((select&~strict).sum())
    counts=np.asarray(rna[:,select].sum(axis=1)).ravel().astype(np.int64)
    detect=np.asarray((rna[:,select]>0).sum(axis=1)).ravel().astype(np.int64)
    assert np.array_equal(counts,out['raw_pseudobulk'][i]) and np.array_equal(detect,out['detected_cells'][i])
    if count:
     means=np.asarray(norm[:,select].sum(axis=1)).ravel()/count;diff=float(np.max(np.abs(means-out['mean_log1p_10k'][i])));assert diff<1e-11,(sample,gid,diff);globaldiff=max(globaldiff,diff)
    else:assert np.isnan(out['mean_log1p_10k'][i]).all()
    totalgroups+=1
  sample_records.append(dict(sample=sample,barcodes=len(bc),groups=len(groups[groups['sample']==sample]),author_QC_positive=int((mask&(ng>0)).sum()),NTC_cells=int((mask&ntc).sum()),readout_sha256=sha(p)))
  print('VERIFIED_SAMPLE',sample,'maxdiff',globaldiff,flush=True)
 report=dict(complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat(),samples=sample_records,groups_reconstructed=totalgroups,every_raw_pseudobulk_count_exact=True,every_detected_cell_count_exact=True,every_cell_assignment_exact=True,maximum_mean_log1p_difference=globaldiff,source_reported_QC_positive=65906,observed_QC_positive=sum(x['author_QC_positive'] for x in sample_records),source_reported_NTC=1021,observed_NTC=sum(x['NTC_cells'] for x in sample_records),model_validation_complete=False,verifier_sha256=sha(Path(__file__)))
 (O/'INDEPENDENT_READOUT_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
