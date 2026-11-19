"""Resume NCBI matrices with validated parallel HTTP ranges, preserving serial prefixes."""
from pathlib import Path
import requests,json,hashlib,gzip,concurrent.futures,os,time
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[2];D=ROOT/'data/perturbation_sources_20261117/GSE311359'
def now():return datetime.now(timezone.utc).isoformat()
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def write(p,j):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def main():
 old=json.loads((D/'COUNT_ACQUISITION_serial.json').read_text());assert len(old['files'])==14
 run=dict(status='running',complete=False,pid=os.getpid(),started_at_utc=now(),source_SOFT_sha256=old['source_SOFT_sha256'],files=old['files'],workers=8,chunk_bytes=2*1024*1024,serial_prefixes=[],chunks_complete=0)
 assert sha(D/'GSE311359_family.soft.gz')==run['source_SOFT_sha256']
 write(D/'COUNT_ACQUISITION.json',run)
 chunks=D/'matrix_range_chunks';chunks.mkdir(exist_ok=True)
 urls=[]
 with gzip.open(D/'GSE311359_family.soft.gz','rt') as f:
  sample=None
  for l in f:
   if l.startswith('^SAMPLE = '):sample=l.strip().split(' = ')[1]
   if l.startswith('!Sample_supplementary_file_') and l.strip().endswith('matrix.mtx.gz'):
    urls.append((sample,l.strip().split(' = ',1)[1].replace('ftp://','https://')))
 plans=[]
 for sample,url in urls:
  p=D/url.rsplit('/',1)[1];prefix=p.with_suffix(p.suffix+'.serial-prefix')
  assert prefix.exists() and not p.exists()
  with requests.head(url,timeout=(15,30)) as r:
   r.raise_for_status();size=int(r.headers['Content-Length']);lm=r.headers['Last-Modified']
  assert 0<=prefix.stat().st_size<size
  plan=dict(sample=sample,url=url,file=str(p),prefix=str(prefix),prefix_bytes=prefix.stat().st_size,prefix_sha256=sha(prefix),total_bytes=size,last_modified=lm)
  plan['chunks']=[dict(start=start,end=min(start+run['chunk_bytes'],size)-1) for start in range(plan['prefix_bytes'],size,run['chunk_bytes'])]
  plans.append(plan)
 run['serial_prefixes']=plans;write(D/'COUNT_ACQUISITION.json',run)
 def fetch(plan,b):
  start,end=b['start'],b['end'];p=chunks/(Path(plan['file']).name+'.'+str(start)+'.bin')
  prior=json.loads((D/'PRESERVED_RANGE_CHUNKS.json').read_text())
  if p.name in prior:
   assert p.stat().st_size==prior[p.name]['bytes']==end-start+1 and sha(p)==prior[p.name]['sha256']
   return dict(sample=plan['sample'],start=start,end=end,file=str(p),sha256=sha(p),bytes=p.stat().st_size,reused_prior_validated_range=True)
  for attempt in range(8):
   try:
    h={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'}
    with requests.get(plan['url'],headers=h,timeout=(15,60),stream=True) as r:
     assert r.status_code==206,(r.status_code,plan['url'])
     assert r.headers['Content-Range']==f"bytes {start}-{end}/{plan['total_bytes']}"
     assert r.headers['Last-Modified']==plan['last_modified']
     with p.with_suffix('.part').open('wb') as f:
      for block in r.iter_content(256*1024):f.write(block)
    assert p.with_suffix('.part').stat().st_size==end-start+1
    p.with_suffix('.part').replace(p)
    return dict(sample=plan['sample'],start=start,end=end,file=str(p),sha256=sha(p),bytes=p.stat().st_size)
   except Exception as e:
    print('RANGE_RETRY',plan['sample'],start,attempt+1,repr(e),flush=True)
    if attempt==7:raise
    time.sleep(min(60,2**(attempt+2)))
 receipts=[];failures=[]
 try:
  with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
   # Round-robin across files avoids monopolizing the first sample.
   jobs=[]
   for index in range(max(len(p['chunks']) for p in plans)):
    for p in plans:
     if index<len(p['chunks']):jobs.append(pool.submit(fetch,p,p['chunks'][index]))
   for future in concurrent.futures.as_completed(jobs):
    try:receipts.append(future.result())
    except Exception as e:failures.append(repr(e))
    run['chunks_complete']=len(receipts);run['range_failures']=failures;write(D/'COUNT_ACQUISITION.json',run)
    if len(receipts)%16==0:print('CHUNKS',len(receipts),'/',len(jobs),flush=True)
  write(D/'MATRIX_RANGE_RECEIPTS.json',receipts)
  assert not failures,failures
  for plan in plans:
   p=Path(plan['file']);prefix=Path(plan['prefix']);assert sha(prefix)==plan['prefix_sha256']
   parts=[prefix]+[chunks/(p.name+'.'+str(b['start'])+'.bin') for b in plan['chunks']]
   with p.with_suffix(p.suffix+'.assembling').open('wb') as out:
    for part in parts:
     with part.open('rb') as f:
      for block in iter(lambda:f.read(4*1024*1024),b''):out.write(block)
   assembled=p.with_suffix(p.suffix+'.assembling');assert assembled.stat().st_size==plan['total_bytes']
   plainbytes=0
   with gzip.open(assembled,'rb') as f:
    for block in iter(lambda:f.read(4*1024*1024),b''):plainbytes+=len(block)
   with requests.head(plan['url'],timeout=(15,30)) as r:
    r.raise_for_status();assert r.headers['Last-Modified']==plan['last_modified'] and int(r.headers['Content-Length'])==plan['total_bytes']
   assembled.replace(p)
   receipt=dict(sample=plan['sample'],url=plan['url'],file=str(p),sha256=sha(p),bytes=p.stat().st_size,uncompressed_bytes=plainbytes,gzip_crc_to_eof=True,serial_prefix_bytes=plan['prefix_bytes'],validated_ranges=len(plan['chunks']),last_modified=plan['last_modified'],finished_at_utc=now())
   write(p.with_name(p.name+'.acquisition.json'),receipt);run['files'].append(receipt);write(D/'COUNT_ACQUISITION.json',run);print('COMPLETE_MATRIX',plan['sample'],flush=True)
  assert len(run['files'])==21
  run.update(status='complete',complete=True,finished_at_utc=now(),total_bytes=sum(x['bytes'] for x in run['files']),preserved_prefix_bytes=sum(x['prefix_bytes'] for x in plans))
 except Exception as e:run.update(status='failed',error=repr(e),finished_at_utc=now());raise
 finally:write(D/'COUNT_ACQUISITION.json',run)
if __name__=='__main__':main()
