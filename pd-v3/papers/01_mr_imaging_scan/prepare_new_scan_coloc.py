"""Exact-identity regional joins for the current new-MR public acquisition batch."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];DATA=ROOT/'data/new_scan_regional_20261117';OUT=ROOT/'papers/01_mr_imaging_scan/results/new_scan_regions_20261117'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,j):p.write_text(json.dumps(j,indent=2)+'\n')
def key(row):
 a,b=row.effect_allele,row.other_allele
 if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b:return None
 return f'{int(row.chromosome)}:{int(row.base_pair_location)}:'+ '_'.join(sorted([a,b]))
def write(frame,p):
 frame.to_csv(p,index=False,compression={'method':'gzip','mtime':0})
 return dict(rows=len(frame),sha256=sha(p),path=str(p.relative_to(ROOT)))
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'INPUT_AUDIT.json').exists():raise RuntimeError('Existing completed input audit; inspect before rerunning')
 for name in ['RUN.json','IMAGING_RUN.json']:
  assert json.loads((DATA/name).read_text())['complete'],name
 records=json.loads((DATA/'FIRST_DOWNLOAD_BATCH.json').read_text())['records'];audits=[]
 for rec in records:
  accession=rec['accession'];folder=DATA/accession;d=OUT/accession;d.mkdir(exist_ok=True)
  acquired=json.loads((folder/'ACQUISITION.json').read_text());imaging=json.loads((folder/'IMAGING_REGION.json').read_text());assert acquired['complete'] and imaging['complete']
  xp=folder/'regional_statistics.csv.gz';yp=folder/'imaging_region_GRCh37.csv.gz';xh=sha(xp);yh=sha(yp)
  assert acquired['download']['validation']['outputs_sha256'][str(xp.relative_to(ROOT))]==xh
  assert imaging['output_sha256']==yh
  x=pd.read_csv(xp,float_precision='round_trip');y=pd.read_csv(yp,float_precision='round_trip')
  sources={str(f.relative_to(ROOT)):sha(f) for f in [xp,yp,folder/'ACQUISITION.json',folder/'IMAGING_REGION.json']}
  for frame in [x,y]:
   assert frame.chromosome.nunique()==1 and frame.base_pair_location.gt(0).all()
   assert np.isfinite(frame[['beta','standard_error','effect_allele_frequency']].to_numpy()).all()
   assert frame.standard_error.gt(0).all() and frame.effect_allele_frequency.between(0,1).all()
   frame['variant_key']=frame.apply(key,axis=1)
   assert not frame.loc[frame.variant_key.notna(),'variant_key'].duplicated().any(),'Ambiguous regional SNP identity'
  cols=['variant_key','chromosome','base_pair_location','variant_id','effect_allele','other_allele','effect_allele_frequency','beta','standard_error','p_value']
  a=x.loc[x.variant_key.notna(),cols+['n_metadata']].rename(columns={'n_metadata':'n'})
  b=y.loc[y.variant_key.notna(),cols+['n_y']].rename(columns={'n_y':'n'})
  a=a.rename(columns={c:c+'_x' for c in a.columns if c!='variant_key'});b=b.rename(columns={c:c+'_y' for c in b.columns if c!='variant_key'})
  union=a.merge(b,on='variant_key',how='outer',validate='one_to_one',indicator=True);counts=union['_merge'].value_counts().to_dict()
  shared=union[union['_merge'].eq('both')].copy().drop(columns='_merge')
  assert len(shared)>0 and shared.chromosome_x.eq(shared.chromosome_y).all() and shared.base_pair_location_x.eq(shared.base_pair_location_y).all()
  shared['pos37']=shared.base_pair_location_x.astype(int);shared['accession']=accession;shared['gene']=rec['gene'];shared['phenotype']=rec['best_phenotype']
  shared['palindromic']=shared.variant_key.str.endswith(('A_T','C_G'))
  shared['direction_status']=np.where(shared.palindromic,'not_resolved_for_signed_use','exact_allele_alignment_required_for_signed_use')
  shared=shared.sort_values(['pos37','variant_key'])
  region=acquired['download']['validation']['region'];assert shared.n_x.eq(rec['sample_size']).all() and shared.n_y.gt(0).all()
  audit=dict(**rec,region=region,protein_source_rows=len(x),imaging_source_rows=len(y),unsupported_protein_rows=int(x.variant_key.isna().sum()),unsupported_imaging_rows=int(y.variant_key.isna().sum()),match_counts=counts,shared_snps=len(shared),input_sha256=sources,outputs={})
  audit['outputs']['union']=write(union,d/'variant_identity_audit.csv.gz');audit['outputs']['shared']=write(shared,d/'shared_snps.csv.gz')
  audit['outputs']['unsupported_protein']=write(x[x.variant_key.isna()],d/'unsupported_protein_rows.csv.gz');audit['outputs']['unsupported_imaging']=write(y[y.variant_key.isna()],d/'unsupported_imaging_rows.csv.gz')
  save(d/'INPUT_AUDIT.json',audit);audits.append(audit);print('COLOC_INPUT',rec['gene'],len(shared),counts,flush=True)
 pd.DataFrame([dict(accession=a['accession'],gene=a['gene'],phenotype=a['best_phenotype'],center=a['region']['center']) for a in audits]).to_csv(OUT/'batch_regions.csv',index=False)
 save(OUT/'INPUT_AUDIT.json',dict(complete=True,created_at_utc=datetime.now(timezone.utc).isoformat(),datasets=audits,code_sha256=sha(Path(__file__)),batch_sha256=sha(DATA/'FIRST_DOWNLOAD_BATCH.json'),matching='GRCh37 chromosome, coordinate, exact unordered alleles; biallelic SNVs only',p_or_maf_selection=None,protein_SD='Estimated explicitly from standard errors, allele frequencies and reported assay N using pinned coloc sdY.est; not assumed 1; N is not variant-specific',imaging_SD='Original Oxford deconfounded SD units: 1',batch_scope='First acquisition batch in full new-MR descriptive order, not a candidate or support cutoff',cross_platform_scope='AGES aptamer protein and UKB Olink assay share author gene labels; assay equivalence not established',model='Single-causal-variant ABF within each region; cannot resolve multiple causal signals',palindromic_use='Included in squared-Z ABF only; no signed-direction conclusion',windows_bp=[500000,1000000],prior_sensitivity='Existing p12 grid 1e-6 through 3e-5 plus 1e-5; p1=p2=1e-4'))
if __name__=='__main__':main()
