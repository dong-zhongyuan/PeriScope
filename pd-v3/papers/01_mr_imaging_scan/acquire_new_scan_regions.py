"""First public regional-validation download batch ordered by the complete replacement MR scan."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import sys,hashlib,json,os,csv,math
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(Path(__file__).parent))
from s07_acquire_big40 import acquire
OUT=ROOT/'data/new_scan_regional_20261117'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(p,j):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def normchr(v):
 v=str(v).removeprefix('chr').lstrip('0');return 'X' if v in ['23','X'] else v
def key(chrom,pos,a,b):return (normchr(chrom),int(pos),tuple(sorted([a,b])))
def main():
 if (OUT/'RUN.json').exists():raise RuntimeError('Inspect existing run before restarting')
 plan=OUT/'FIRST_DOWNLOAD_BATCH.json';records=json.loads(plan.read_text())['records'];instruments=pd.read_csv(ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv')
 peaks_path=ROOT/'data/big40_peak_catalog/all_10889_source_peaks_with_frequencies.csv';peaks=pd.read_csv(peaks_path,dtype={'idp':str,'chromosome':str})
 inputs={str(p.relative_to(ROOT)):sha(p) for p in [plan,peaks_path,ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv',Path(__file__),ROOT/'papers/01_mr_imaging_scan/s07_acquire_big40.py']}
 run=dict(status='running',pid=os.getpid(),complete=False,started_at_utc=datetime.now(timezone.utc).isoformat(),input_sha256=inputs,completed=[],failures={});dump(OUT/'RUN.json',run)
 def one(rec):
  acc=rec['accession'];folder=OUT/acc;lead=instruments[instruments.exposure_id.eq(rec['exposure_id'])];assert len(lead)==1;lead=lead.iloc[0]
  chrom=normchr(lead.chrom);center=int(lead.pos37);low=max(1,center-1000000);high=center+1000000
  pp=peaks[peaks.phenotype==rec['best_phenotype']].copy();peakkeys={key(x.chromosome,x.position_grch37,x.effect_allele,x.other_allele) for x in pp.itertuples()};pp.to_csv(folder/'all_source_imaging_peaks.csv',index=False)
  def validate(path):
   h=hashlib.md5()
   with path.open('rb') as f:
    for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
   assert h.hexdigest()==rec['published_md5'],'Published original MD5 differs'
   regional=[];reverse=[];total=0
   with path.open() as f:
    reader=csv.DictReader(f,delimiter='\t');assert reader.fieldnames=='variant_id p_value chromosome base_pair_location effect_allele other_allele effect_allele_frequency beta standard_error'.split()
    for total,row in enumerate(reader,1):
     if None in row or any(v is None for v in row.values()):raise ValueError('Malformed source row')
     c=normchr(row['chromosome']);pos=int(row['base_pair_location']);a=row['effect_allele'];b=row['other_allele']
     in_region=c==chrom and low<=pos<=high;reverse_match=key(c,pos,a,b) in peakkeys
     if not in_region and not reverse_match:continue
     values=[float(row[k]) for k in ['p_value','effect_allele_frequency','beta','standard_error']];assert all(math.isfinite(v) for v in values) and 0<=values[0]<=1 and 0<=values[1]<=1 and values[3]>0
     saved=dict(**row,source_row=total,gene=rec['gene'],assay=rec['assay'],accession=acc,build='GRCh37',n_metadata=rec['sample_size'],n_scope='reported_assay_sample_size_not_per_variant_N')
     if in_region:regional.append(saved)
     if reverse_match:reverse.append(saved)
   assert total>0 and regional
   for name,rows in [('regional_statistics.csv.gz',regional),('imaging_peak_outcomes.csv.gz',reverse)]:pd.DataFrame(rows,columns=list(regional[0])).to_csv(folder/name,index=False,compression={'method':'gzip','mtime':0})
   return dict(published_md5_verified=True,md5=rec['published_md5'],source_rows=total,source_eof_verified=True,region={'chrom':chrom,'start':low,'end':high,'center':center},regional_rows=len(regional),source_imaging_peaks=len(pp),imaging_peak_outcome_rows=len(reverse),full_summary_statistics_retained=True,protein_SD_not_assumed=True,outputs_sha256={str((folder/n).relative_to(ROOT)):sha(folder/n) for n in ['regional_statistics.csv.gz','imaging_peak_outcomes.csv.gz','all_source_imaging_peaks.csv']})
  result=acquire(rec['url'],OUT/'download'/rec['file_name'],validator=validate)
  receipt=dict(**rec,complete=True,download=result,metadata_sha256=sha(folder/(acc+'_buildGRCh37.tsv-meta.yaml')))
  dump(folder/'ACQUISITION.json',receipt);print('ACQUIRED',acc,rec['gene'],flush=True);return acc
 with ThreadPoolExecutor(max_workers=3) as pool:
  futures={pool.submit(one,rec):rec for rec in records}
  for future in as_completed(futures):
   try:run['completed'].append(future.result())
   except Exception as e:run['failures'][futures[future]['accession']]=repr(e)
   dump(OUT/'RUN.json',run)
 for f,h in inputs.items():assert sha(ROOT/f)==h,f
 run.update(status='complete' if not run['failures'] else 'complete_with_download_errors',complete=not bool(run['failures']),finished_at_utc=datetime.now(timezone.utc).isoformat());dump(OUT/'RUN.json',run)
if __name__=='__main__':main()
