"""Resolve every exact-gene AGES source for all 2,941 registered UKB assays."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import csv,hashlib,json,re,time,os,threading
import requests
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'data/all_target_followup_20261118'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):
 p.parent.mkdir(parents=True,exist_ok=True)
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n');t.replace(p)
_tls=threading.local()
def get(url,path):
 receipt=path.with_suffix(path.suffix+'.receipt.json')
 if path.exists() and receipt.exists():
  j=json.loads(receipt.read_text());assert sha(path)==j['sha256'] and j['url']==url
  return path.read_text()
 if not hasattr(_tls,'session'):_tls.session=requests.Session()
 for attempt in range(4):
  try:
   r=_tls.session.get(url,timeout=(10,30));r.raise_for_status()
   assert len(r.content)<200000
   path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(r.content)
   dump(receipt,dict(url=url,status=r.status_code,headers=dict(r.headers),sha256=sha(path),fetched_at=datetime.now(timezone.utc).isoformat()))
   return r.text
  except Exception:
   if attempt==3:raise
   time.sleep(2**attempt)
def resolve(row):
 acc=row['accession'];d=OUT/'sources'/acc
 if (d/'SOURCE.json').exists():
  old=json.loads((d/'SOURCE.json').read_text())
  if old.get('complete'):return old
 n=int(acc[4:]);lo=((n-1)//1000)*1000+1
 base=f'http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{lo}-GCST{lo+999}/{acc}/'
 name=acc+'_buildGRCh37.tsv'
 raw=get(base+name+'-meta.yaml',d/(name+'-meta.yaml'))
 def field(pattern,text=raw):
  m=re.search(pattern,text,re.M);assert m,pattern
  return m.group(1).strip()
 assert field(r'^gwas_id:\s*(\S+)$')==acc
 assert field(r'^genome_assembly:\s*(\S+)$')=='GRCh37'
 assert field(r'^data_file_name:\s*(\S+)$')==name
 trait=field(r'^trait_description:\s*\n\s*-\s*(.+)$')
 assert trait==row['reported_trait'],(acc,trait,row['reported_trait'])
 ns=re.findall(r'^\s+sample_size:\s*(\d+)\s*$',raw,re.M);assert len(ns)==1
 md5=field(r'^data_file_md5sum:\s*(\w+)$');assert len(md5)==32
 directory=get(base+'harmonised/',d/'harmonised_index.html')
 names=re.findall(r'href="([^"]+-Build37\.f\.tsv\.gz)"',directory)
 assert len(names)<=1 and all(acc in n for n in names),(acc,names)
 fname=names[0] if names else ''
 sums=get(base+'harmonised/md5sum.txt',d/'harmonised_md5sum.txt')
 checks={v[1]:v[0] for line in sums.splitlines() if len(v:=line.split())==2}
 assert all(re.fullmatch('[0-9a-f]{32}',v) for v in checks.values())
 fmd5=checks.get(fname,'')
 files=[d/(name+'-meta.yaml'),d/'harmonised_index.html',d/'harmonised_md5sum.txt']
 if fname and f'href="{fname}-meta.yaml"' in directory:
  fmeta=get(base+'harmonised/'+fname+'-meta.yaml',d/(fname+'-meta.yaml'))
  side_md5=field(r'^data_file_md5sum:\s*(\w+)$',fmeta)
  if fmd5:assert side_md5==fmd5
  else:fmd5=side_md5
  fbuild=field(r'^genome_assembly:\s*(\S+)$',fmeta);files.append(d/(fname+'-meta.yaml'))
 else:fbuild='formatted_sidecar_not_published'
 result=dict(complete=True,accession=acc,ages_assay=row['ages_assay'],gene=row['ages_target'],
 original_url=base+name,original_name=name,original_published_md5=md5,
 original_build='GRCh37',sample_size=int(ns[0]),sample_size_scope='reported_assay_not_per_variant',
 formatted_url=base+'harmonised/'+fname if fname else '',formatted_name=fname,formatted_published_md5=fmd5,
 formatted_metadata_build=fbuild,
 formatted_checksum_status='published_MD5_available' if fmd5 else 'no_published_compressed_MD5;require_reconstructed_original_MD5',
 indexed_harmonized_files=[k for k in checks if k.endswith('.tbi')],
 formatted_use_status='not_approved_pending_full_original_content_comparison' if fname else 'formatted_file_not_listed',
 source_scope='exact_author_gene_label_not_assay_equivalence',
 input_sha256={str(f.relative_to(ROOT)):sha(f) for f in files})
 dump(d/'SOURCE.json',result);return result
def main():
 OUT.mkdir(exist_ok=True)
 mapping=ROOT/'data/ages_catalog/ukb_all_assays_to_ages.csv'
 rows=list(csv.DictReader(mapping.open()));unique={}
 for r in rows:
  if r['accession']:
   if r['accession'] in unique:assert unique[r['accession']]['ages_target']==r['ages_target']
   unique[r['accession']]=r
 run=dict(status='running',complete=False,pid=os.getpid(),planned_ukb_assays=2941,
 planned_unique_ages_sources=len(unique),matched_ukb_assays=len({r['exposure_id'] for r in rows if r['accession']}),
 completed=[],failures={},workers=8,started_at=datetime.now(timezone.utc).isoformat(),
 input_sha256={str(mapping.relative_to(ROOT)):sha(mapping)},code_sha256=sha(Path(__file__)))
 dump(OUT/'METADATA_RUN.json',run)
 with ThreadPoolExecutor(max_workers=8) as pool:
  jobs={pool.submit(resolve,r):a for a,r in sorted(unique.items())}
  for f in as_completed(jobs):
   a=jobs[f]
   try:f.result();run['completed'].append(a)
   except Exception as e:run['failures'][a]=repr(e)
   run['updated_at']=datetime.now(timezone.utc).isoformat();dump(OUT/'METADATA_RUN.json',run)
   if len(run['completed'])%50==0:print('SOURCES',len(run['completed']),len(unique),'errors',len(run['failures']),flush=True)
 assert sha(mapping)==run['input_sha256'][str(mapping.relative_to(ROOT))]
 run.update(complete=not run['failures'],status='complete' if not run['failures'] else 'finished_with_source_errors',
 finished_at=datetime.now(timezone.utc).isoformat());dump(OUT/'METADATA_RUN.json',run)
if __name__=='__main__':main()
