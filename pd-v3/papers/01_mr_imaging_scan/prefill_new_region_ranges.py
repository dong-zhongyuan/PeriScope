"""Parallel immutable HTTP ranges; preserve the original contiguous prefixes and validate full published MD5."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib,json,time,os,subprocess,requests,shutil
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'data/new_scan_regional_20261117'
def digest(p,kind='sha256'):
 h=hashlib.new(kind)
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(p,j):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def main():
 plan=json.loads((OUT/'PARALLEL_RANGE_PLAN.json').read_text());chunks=[];prefixes={};started=time.monotonic()
 for rec in plan['records']:
  part=OUT/'download'/(rec['file_name']+'.part');assert part.stat().st_size==rec['prefix_size'] and digest(part)==rec['prefix_sha256'];prefixes[rec['accession']]=part
  d=OUT/'range_chunks'/rec['accession'];d.mkdir(parents=True,exist_ok=True)
  for start in range(rec['prefix_size'],rec['size'],plan['chunk_bytes']):chunks.append((rec,start,min(rec['size']-1,start+plan['chunk_bytes']-1),d))
 # Interleave files so concurrency is spread over all three current inputs.
 chunks.sort(key=lambda x:(x[1]-x[0]['prefix_size'],x[0]['accession']))
 run=dict(status='running',pid=os.getpid(),complete=False,workers=plan['workers'],planned_chunks=len(chunks),completed_chunks=0,downloaded_bytes=0,failures=[],started_at_utc=datetime.now(timezone.utc).isoformat());dump(OUT/'PARALLEL_RANGE_RUN.json',run)
 def one(task):
  rec,start,end,d=task;path=d/f'{start:012d}-{end:012d}.bin';size=end-start+1;last=None
  if path.exists() and path.with_suffix('.json').exists():
   receipt=json.loads(path.with_suffix('.json').read_text())
   if path.stat().st_size==size and digest(path)==receipt['sha256']:return size,0.0
   raise ValueError('Existing range failed identity verification')
  for attempt in range(5):
   try:
    t=time.monotonic();tmp=path.with_suffix('.part')
    headers={'Range':f'bytes={start}-{end}','If-Range':rec['http']['etag'] or rec['http']['last_modified'],'Accept-Encoding':'identity'}
    with requests.get(rec['url'],headers=headers,stream=True,timeout=(20,60)) as response:
     response.raise_for_status();assert response.status_code==206 and response.headers.get('Content-Range')==f'bytes {start}-{end}/{rec["size"]}'
     for key,header in [('etag','ETag'),('last_modified','Last-Modified')]:
      if rec['http'][key]:assert response.headers.get(header)==rec['http'][key]
     count=0;h=hashlib.sha256()
     with tmp.open('wb') as f:
      for b in response.iter_content(1024*1024):
       count+=len(b);assert count<=size;f.write(b);h.update(b)
     assert count==size;tmp.replace(path)
    elapsed=time.monotonic()-t;dump(path.with_suffix('.json'),dict(start=start,end=end,size=count,sha256=h.hexdigest(),url=rec['url'],seconds=elapsed,source_identity=rec['http']));return size,elapsed
   except Exception as e:last=repr(e);time.sleep(min(2**attempt,10))
  raise RuntimeError((rec['accession'],start,last))
 with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
  futures={pool.submit(one,x):x for x in chunks}
  for future in as_completed(futures):
   try:
    size,elapsed=future.result();run['downloaded_bytes']+=size;run['completed_chunks']+=1
   except Exception as e:run['failures'].append(repr(e))
   run['elapsed_seconds']=time.monotonic()-started;run['aggregate_MiB_per_second']=run['downloaded_bytes']/1048576/run['elapsed_seconds'];dump(OUT/'PARALLEL_RANGE_RUN.json',run)
   if run['completed_chunks']%12==0:print('RANGES',run['completed_chunks'],run['aggregate_MiB_per_second'],'MiB/s',flush=True)
 if run['failures']:raise RuntimeError('Range failures preserved; do not expose incomplete assembled input')
 assembled=[]
 for rec in plan['records']:
  part=prefixes[rec['accession']];assert part.stat().st_size==rec['prefix_size'] and digest(part)==rec['prefix_sha256']
  temp=part.with_suffix('.assembled')
  with temp.open('wb') as f:
   with part.open('rb') as src:shutil.copyfileobj(src,f,8*1024*1024)
   relevant=sorted([x for x in chunks if x[0]['accession']==rec['accession']],key=lambda x:x[1]);offset=rec['prefix_size']
   for _,start,end,d in relevant:
    assert start==offset;path=d/f'{start:012d}-{end:012d}.bin';receipt=json.loads(path.with_suffix('.json').read_text());assert digest(path)==receipt['sha256']
    with path.open('rb') as src:shutil.copyfileobj(src,f,8*1024*1024)
    offset=end+1
   f.flush();os.fsync(f.fileno())
  assert offset==rec['size']==temp.stat().st_size and digest(temp,'md5')==rec['published_md5'];temp.replace(part);assembled.append(rec['accession'])
 run.update(status='complete',complete=True,assembled=assembled,published_MD5_all_verified=True,finished_at_utc=datetime.now(timezone.utc).isoformat());dump(OUT/'PARALLEL_RANGE_RUN.json',run)
 # The unchanged scientific validator sees a fully prefilled contiguous part and validates every source row.
 with (OUT/'controller.log').open('a') as log:
  job=subprocess.Popen(['/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python','-u',str(ROOT/'papers/01_mr_imaging_scan/acquire_new_scan_regions.py')],cwd=str(ROOT),stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 dump(OUT/'JOB.json',dict(pid=job.pid,prefilled_by='PARALLEL_RANGE_RUN.json'));print('VALIDATOR_STARTED',job.pid,flush=True)
if __name__=='__main__':main()
