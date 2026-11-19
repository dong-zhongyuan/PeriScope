"""Public EUR genotypes for all reported reverse-direction peaks in the first new-scan batch."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import os,json,hashlib,re
import requests,numpy as np,pandas as pd,pysam
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/'data/new_scan_regional_20261117';OUT=BASE/'reverse_LD'
URL='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/';PANEL='integrated_call_samples_v3.20130502.ALL.panel'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,j):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'RUN.json').exists():raise RuntimeError('Existing run; inspect before restarting')
 run=dict(status='running',pid=os.getpid(),complete=False,started_at_utc=datetime.now(timezone.utc).isoformat(),completed_chromosomes=[],failures={});dump(OUT/'RUN.json',run)
 records=json.loads((BASE/'FIRST_DOWNLOAD_BATCH.json').read_text())['records'];allpeaks=pd.concat([pd.read_csv(BASE/x['accession']/'all_source_imaging_peaks.csv').assign(accession=x['accession']) for x in records],ignore_index=True)
 allpeaks['variant_key']=allpeaks.apply(lambda x:str(x.chromosome)+':'+str(x.position_grch37)+':'+ '_'.join(sorted([x.effect_allele,x.other_allele])),axis=1)
 allpeaks.to_csv(OUT/'all_requested_peaks.csv',index=False);unique=allpeaks.drop_duplicates('variant_key')
 panel_path=OUT/PANEL
 if not panel_path.exists():
  response=requests.get(URL+PANEL,timeout=45);response.raise_for_status();panel_path.write_bytes(response.content)
 panel=pd.read_csv(panel_path,sep='\t',usecols=['sample','pop','super_pop']);eur=panel[panel.super_pop=='EUR'];samples=eur['sample'].tolist();assert len(samples)==len(set(samples))==503
 response=requests.get(URL,timeout=45);response.raise_for_status();(OUT/'source_listing.html').write_text(response.text);names=re.findall(r'href="([^"]+)"',response.text)
 def chromosome(chrom):
  d=OUT/('chr'+str(chrom));d.mkdir(exist_ok=True);entries=unique[unique.chromosome==chrom]
  choices=[x for x in names if x.startswith(f'ALL.chr{chrom}.') and x.endswith('.genotypes.vcf.gz')];assert len(choices)==1;name=choices[0];url=URL+name;index=d/(name+'.tbi')
  if not index.exists():
   response=requests.get(url+'.tbi',timeout=60);response.raise_for_status();index.write_bytes(response.content)
  response=requests.head(url,timeout=45);response.raise_for_status();identity={k:response.headers.get(k) for k in ['ETag','Last-Modified','Content-Length']}
  rows=[];dosages=[]
  with pysam.VariantFile(url,index_filename=str(index)) as source:
   source.subset_samples(samples);actual_samples=list(source.header.samples);assert set(actual_samples)==set(samples);order=[actual_samples.index(x) for x in samples]
   for entry in entries.itertuples():
    candidates=[]
    for variant in source.fetch(str(chrom),int(entry.position_grch37)-1,int(entry.position_grch37)):
     if variant.pos!=entry.position_grch37 or not variant.alts or len(variant.alts)!=1:continue
     if sorted([variant.ref,variant.alts[0]])!=sorted([entry.effect_allele,entry.other_allele]):continue
     if len(variant.ref)!=1 or len(variant.alts[0])!=1:continue
     candidates.append(variant)
    assert len(candidates)<=1,'Ambiguous reference allele identity'
    if not candidates:
     rows.append(dict(variant_key=entry.variant_key,chromosome=chrom,position_grch37=entry.position_grch37,status='not_matched_biallelic_SNV_in_reference',ref='',alt='',effect_allele=entry.effect_allele,dosage_row=-1));continue
    v=candidates[0];g=[]
    for sample in samples:
     gt=v.samples[sample]['GT'];g.append(-1 if len(gt)!=2 or any(a is None for a in gt) else sum(gt))
    a=np.array(g,dtype=np.int8);valid=a>=0
    if entry.effect_allele==v.ref:a[valid]=2-a[valid]
    row=dict(variant_key=entry.variant_key,chromosome=chrom,position_grch37=entry.position_grch37,status='matched',ref=v.ref,alt=v.alts[0],effect_allele=entry.effect_allele,dosage_row=len(dosages),missing_calls=int((a<0).sum()),effect_allele_frequency=float(a[a>=0].mean()/2) if valid.any() else np.nan)
    rows.append(row);dosages.append(a)
  response=requests.head(url,timeout=45);response.raise_for_status();assert identity=={k:response.headers.get(k) for k in identity}
  table=pd.DataFrame(rows);table.to_csv(d/'coverage.csv',index=False)
  g=np.stack(dosages) if dosages else np.empty((0,503),dtype=np.int8);np.savez_compressed(d/'effect_allele_dosages.npz',dosages=g,samples=np.array(samples),variant_keys=table.loc[table.status=='matched','variant_key'].to_numpy(str))
  dump(d/'MANIFEST.json',dict(complete=True,url=url,http_identity=identity,index_sha256=sha(index),panel_sha256=sha(panel_path),source_reference='1000 Genomes Phase3 EUR503 GRCh37',cohort_matched=False,whole_remote_file_checksum_verified=False,requested=len(entries),matched=len(dosages),outputs_sha256={n:sha(d/n) for n in ['coverage.csv','effect_allele_dosages.npz']}));print('LD_CHR',chrom,len(dosages),'of',len(entries),flush=True);return int(chrom)
 with ThreadPoolExecutor(max_workers=6) as pool:
  futures={pool.submit(chromosome,c):c for c in sorted(unique.chromosome.unique())}
  for future in as_completed(futures):
   try:run['completed_chromosomes'].append(future.result())
   except Exception as error:run['failures'][str(futures[future])]=repr(error)
   dump(OUT/'RUN.json',run)
 run.update(status='complete' if not run['failures'] else 'complete_with_chromosome_errors',complete=not bool(run['failures']),finished_at_utc=datetime.now(timezone.utc).isoformat(),all_source_peak_rows=len(allpeaks),unique_requested_variants=len(unique),samples=503,pysam_version=pysam.__version__,code_sha256=sha(Path(__file__)));dump(OUT/'RUN.json',run)
if __name__=='__main__':main()
