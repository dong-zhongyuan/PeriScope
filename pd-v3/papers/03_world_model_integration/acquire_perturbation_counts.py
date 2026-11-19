"""Acquire the seven public GSE311359 count matrices directly on A6000."""
from pathlib import Path
import gzip,json,hashlib,os,time,concurrent.futures
from datetime import datetime,timezone
import requests
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'data/perturbation_sources_20261117/GSE311359'
def now():return datetime.now(timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,j):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def main():
 runpath=OUT/'COUNT_ACQUISITION.json'
 if runpath.exists():raise RuntimeError('Inspect existing acquisition before restarting')
 run=dict(status='running',complete=False,pid=os.getpid(),started_at_utc=now(),files=[])
 write(runpath,run)
 soft=OUT/'GSE311359_family.soft.gz';run['source_SOFT_sha256']=sha(soft)
 jobs=[]
 with gzip.open(soft,'rt') as f:
  sample=None
  for line in f:
   if line.startswith('^SAMPLE = '):sample=line.strip().split(' = ')[1]
   if line.startswith('!Sample_supplementary_file_'):
    url=line.strip().split(' = ',1)[1].replace('ftp://','https://');name=url.rsplit('/',1)[1]
    jobs.append((sample,url,name))
 assert len(jobs)==21 and len(set(j[0] for j in jobs))==7
 def fetch(job):
  sample,url,name=job;p=OUT/name;receipt=p.with_name(name+'.acquisition.json')
  started=time.monotonic();reused=p.exists()
  if not reused:
   part=p.with_suffix(p.suffix+'.part')
   if part.exists():raise RuntimeError('Partial already exists: '+str(part))
   with requests.get(url,stream=True,timeout=(15,60)) as r:
    r.raise_for_status();headers={k:v for k,v in r.headers.items() if k.lower() in ['content-length','last-modified','etag','content-type']}
    with part.open('wb') as f:
     for b in r.iter_content(1024*1024):f.write(b)
    if 'Content-Length' in r.headers:assert part.stat().st_size==int(r.headers['Content-Length'])
   part.replace(p)
  else:headers={}
  # Read to EOF so the gzip trailer and CRC are checked independently of transfer length.
  plain_bytes=0
  with gzip.open(p,'rb') as f:
   for b in iter(lambda:f.read(4*1024*1024),b''):plain_bytes+=len(b)
  m=dict(sample=sample,url=url,file=str(p),sha256=sha(p),bytes=p.stat().st_size,uncompressed_bytes=plain_bytes,gzip_crc_to_eof=True,reused=reused,headers=headers,finished_at_utc=now(),seconds=time.monotonic()-started)
  write(receipt,m);return m
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=7) as pool:
   futs=[pool.submit(fetch,j) for j in jobs]
   for f in concurrent.futures.as_completed(futs):
    m=f.result();run['files'].append(m);write(runpath,run);print('COMPLETE',m['sample'],Path(m['file']).name,m['bytes'],flush=True)
  assert sha(soft)==run['source_SOFT_sha256']
  run.update(status='complete',complete=True,finished_at_utc=now(),total_bytes=sum(x['bytes'] for x in run['files']))
 except Exception as e:run.update(status='failed',error=repr(e),finished_at_utc=now());raise
 finally:write(runpath,run)
if __name__=='__main__':main()
