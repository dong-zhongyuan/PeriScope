"""Extract complete original-GRCh37 imaging regions for the new-scan public batch."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib,json,gzip,os,sys
import pandas as pd,numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from pdv3.native_imaging import stream_oxford
OUT=ROOT/'data/new_scan_regional_20261117'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(p,j):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def norm(v):return str(v).lstrip('0').replace('X','23')
def one(task):
 rec,meta,frequency=task;folder=OUT/rec['accession'];points=[(x['chromosome'],x['base_pair_location']) for x in frequency]
 url=f'https://open.oxcin.ox.ac.uk/ukbiobank/big40/release2/stats/{meta["idp"]}.txt.gz'
 rows,receipt=stream_oxford(url,points,meta);frame=pd.DataFrame(rows);freq=pd.DataFrame(frequency)
 keys=['chromosome','base_pair_location','effect_allele','other_allele','rsid']
 assert not frame.duplicated(keys).any() and not freq.duplicated(keys).any()
 frame=frame.merge(freq[keys+['effect_allele_frequency','imputation_info']],on=keys,how='outer',validate='one_to_one',indicator=True)
 assert frame['_merge'].eq('both').all(),'Native statistics and full variant catalog regional rows differ'
 frame=frame.drop(columns='_merge');assert len(frame)==len(frequency)
 output=folder/'imaging_region_GRCh37.csv.gz';frame.to_csv(output,index=False,compression={'method':'gzip','mtime':0})
 receipt.update(complete=True,phenotype_record=meta,full_regional_variant_count=len(frame),frequency_source='Unfiltered official full variant table, exact position/rsID/allele join',allele_filters=None,output=str(output.relative_to(ROOT)),output_sha256=sha(output))
 dump(folder/'IMAGING_REGION.json',receipt);print('IMAGING_REGION',rec['gene'],len(frame),flush=True);return rec['accession']
def main():
 status=OUT/'IMAGING_RUN.json'
 if status.exists():raise RuntimeError('Inspect existing run before restarting')
 records=json.loads((OUT/'FIRST_DOWNLOAD_BATCH.json').read_text())['records'];inst=pd.read_csv(ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv');meta=pd.read_csv(ROOT/'data/big40_atlas/phenotypes_3935.csv',dtype={'idp':str})
 varpath=ROOT/'data/big40_variants/download/variants.txt.gz';expected='152517ada2fb5c572b26959d184f099369e521eeac681c01b3a3ca94b1e7299a';assert sha(varpath)==expected
 run=dict(status='preparing_regions',pid=os.getpid(),complete=False,started_at_utc=datetime.now(timezone.utc).isoformat(),completed=[],failures={},input_sha256={str(varpath.relative_to(ROOT)):expected,str(Path(__file__).relative_to(ROOT)):sha(Path(__file__)),str((ROOT/'pdv3/native_imaging.py').relative_to(ROOT)):sha(ROOT/'pdv3/native_imaging.py')});dump(status,run)
 regions={};frequency={}
 for rec in records:
  lead=inst[inst.exposure_id==rec['exposure_id']];assert len(lead)==1;lead=lead.iloc[0];regions[rec['accession']]=(norm(int(lead.chrom)),max(1,int(lead.pos37)-1000000),int(lead.pos37)+1000000);frequency[rec['accession']]=[]
 with gzip.open(varpath,'rt') as f:
  assert f.readline().split()==['chr','rsid','pos','a1','a2','af','info']
  for total,line in enumerate(f,1):
   v=line.split();assert len(v)==7;c=norm(v[0]);position=int(v[2])
   for accession,(chrom,low,high) in regions.items():
    if c==chrom and low<=position<=high:
     frequency[accession].append(dict(chromosome=c,base_pair_location=position,rsid=v[1],other_allele=v[3],effect_allele=v[4],effect_allele_frequency=float(v[5]),imputation_info=float(v[6])))
 assert total==17103079
 tasks=[]
 for rec in records:
  selected=meta[meta.phenotype==rec['best_phenotype']];assert len(selected)==1;m=selected.iloc[0].to_dict();f=frequency[rec['accession']];assert f
  pd.DataFrame(f).to_csv(OUT/rec['accession']/'complete_regional_variant_frequencies.csv.gz',index=False,compression={'method':'gzip','mtime':0});tasks.append((rec,m,f))
 run.update(status='streaming_three_original_imaging_sources',regional_variant_counts={k:len(v) for k,v in frequency.items()});dump(status,run)
 with ProcessPoolExecutor(max_workers=3) as pool:
  futures={pool.submit(one,t):t[0]['accession'] for t in tasks}
  for future in as_completed(futures):
   try:run['completed'].append(future.result())
   except Exception as error:run['failures'][futures[future]]=repr(error)
   dump(status,run)
 for f,h in run['input_sha256'].items():assert sha(ROOT/f)==h
 run.update(status='complete' if not run['failures'] else 'complete_with_source_errors',complete=not bool(run['failures']),finished_at_utc=datetime.now(timezone.utc).isoformat());dump(status,run)
if __name__=='__main__':main()
