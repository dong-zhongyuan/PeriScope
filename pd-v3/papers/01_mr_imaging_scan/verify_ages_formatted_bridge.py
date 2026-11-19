"""Validate a smaller published representation against an already MD5-verified original."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import requests,hashlib,json,gzip,time,os,itertools,math
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'data/all_target_followup_20261118'
def sha(p,kind='sha256'):
 h=hashlib.new(kind)
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,v):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(v,indent=2)+'\n');t.replace(p)
def main():
 d=OUT/'formatted_bridge';d.mkdir(exist_ok=True)
 name='35078996-GCST90090368-EFO_0007937-Build37.f.tsv.gz'
 url='http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90090001-GCST90091000/GCST90090368/harmonised/'+name
 size=177389960;md5='e365929edd45154e894a864acbce3127'
 run=dict(status='downloading',complete=False,pid=os.getpid(),started_at=datetime.now(timezone.utc).isoformat(),url=url,downloaded_bytes=0)
 dump(d/'RUN.json',run);t0=time.monotonic()
 probe=json.loads((OUT/'probe/probe.json').read_text());identity=probe['headers']
 first=OUT/'probe'/(name+'.first1MiB');assert first.stat().st_size==1048576 and sha(first)==probe['sha256']
 parts=[(0,1048575,first)]
 for start in range(1048576,size,4<<20):parts.append((start,min(start+(4<<20)-1,size-1),d/f'{start}.bin'))
 def download(task):
  start,end,path=task
  if start==0:return end-start+1
  receipt=path.with_suffix('.json')
  if path.exists() and receipt.exists():
   info=json.loads(receipt.read_text());assert path.stat().st_size==end-start+1 and sha(path)==info['sha256'];return end-start+1
  for attempt in range(5):
   try:
    with requests.get(url,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},stream=True,timeout=(10,45)) as r:
     assert r.status_code==206 and r.headers['Content-Range']==f'bytes {start}-{end}/{size}'
     for k in ['ETag','Last-Modified']:assert r.headers[k]==identity[k]
     tmp=path.with_suffix('.part');count=0
     with tmp.open('wb') as f:
      for b in r.iter_content(1<<20):count+=len(b);assert count<=end-start+1;f.write(b)
     assert count==end-start+1;tmp.replace(path)
    dump(receipt,dict(start=start,end=end,url=url,sha256=sha(path),bytes=count));return count
   except Exception:
    if attempt==4:raise
    time.sleep(min(2**attempt,10))
 with ThreadPoolExecutor(max_workers=8) as pool:
  for future in as_completed([pool.submit(download,x) for x in parts]):
   run['downloaded_bytes']+=future.result();run['elapsed_seconds']=time.monotonic()-t0;dump(d/'RUN.json',run)
 target=d/name
 with target.with_suffix('.assembled').open('wb') as out:
  for _,_,part in parts:
   with part.open('rb') as f:
    for b in iter(lambda:f.read(8<<20),b''):out.write(b)
 target.with_suffix('.assembled').replace(target)
 assert target.stat().st_size==size and sha(target,'md5')==md5
 run.update(status='comparing_every_source_row',download_seconds=time.monotonic()-t0,compressed_sha256=sha(target));dump(d/'RUN.json',run)
 original=ROOT/'data/new_scan_regional_20261117/download/GCST90090368_buildGRCh37.tsv'
 assert sha(original,'md5')=='9fd25b6cecb70b234f9207b2f47057d6'
 projected=hashlib.md5();ident=hashlib.sha256();mismatch=[];different_rows=0;numeric_different=0;extras=set()
 with original.open('rb') as a,gzip.open(target,'rb') as b:
  ah=a.readline().rstrip(b'\r\n').split(b'\t');bh=b.readline().rstrip(b'\r\n').split(b'\t');assert bh[:9]==ah
  projected.update(b'\t'.join(bh[:9])+b'\n')
  for count,(ra,rb) in enumerate(itertools.zip_longest(a,b),1):
   assert ra is not None and rb is not None
   av=ra.rstrip(b'\r\n').split(b'\t');bv=rb.rstrip(b'\r\n').split(b'\t')
   assert len(av)==9 and len(bv)==12
   ident.update(b'\t'.join([bv[i] for i in [0,2,3,4,5]])+b'\n')
   projected.update(b'\t'.join(bv[:9])+b'\n')
   extras.add(tuple(x.decode() for x in bv[9:]))
   if av!=bv[:9]:
    different_rows+=1
    equal=all(av[i]==bv[i] if i in [0,2,3,4,5] else float(av[i])==float(bv[i]) for i in range(9))
    if not equal:numeric_different+=1
    if len(mismatch)<10:mismatch.append(dict(row=count,original=ra.decode().strip(),formatted=rb.decode().strip(),all_values_equal=equal))
 assert count==7506463
 result=dict(complete=True,rows=count,original_published_md5_verified=True,formatted_published_md5_verified=True,
 original_sha256=sha(original),formatted_sha256=sha(target),first_nine_projected_md5=projected.hexdigest(),
 first_nine_text_different_rows=different_rows,first_nine_value_different_rows=numeric_different,
 difference_examples=mismatch,appended_field_values=sorted(extras),original_variant_identity_sha256=ident.hexdigest(),
 original_build='GRCh37',formatted_metadata_claim='GRCh38 conflicting with source content',
 permitted_use='Only after complete source identity and published MD5 checks; cannot infer build from filename or conflicting sidecar',
 elapsed_seconds=time.monotonic()-t0)
 dump(d/'FULL_CONTENT_COMPARISON.json',result)
 run.update(complete=True,status='complete',finished_at=datetime.now(timezone.utc).isoformat());dump(d/'RUN.json',run)
 print(json.dumps(result,indent=2),flush=True)
if __name__=='__main__':main()
