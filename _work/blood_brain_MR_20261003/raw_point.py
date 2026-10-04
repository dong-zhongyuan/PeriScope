"""Recover exact records from sorted public plaintext tables via HTTP byte ranges."""
from range_tabix import get
import urllib.request,hashlib,json
from pathlib import Path
class RawPoint:
 def __init__(self,url,cache):
  self.url=url;self.cache=Path(cache);self.cache.mkdir(parents=True,exist_ok=True)
  r=urllib.request.urlopen(urllib.request.Request(url,headers={'Range':'bytes=0-4095'}),timeout=40);assert r.status==206
  b=r.read();self.size=int(r.headers['Content-Range'].split('/')[-1]);self.header=b.decode().splitlines()[0].split('\t')
 def fetch(self,chrom,pos,rsid):
  target=(int(chrom),int(pos));left=0;right=self.size-1
  for _ in range(40):
   start=max(0,(left+right)//2-32768);end=min(self.size-1,start+65535)
   p=self.cache/(hashlib.sha256(f'{self.url}:{start}:{end}'.encode()).hexdigest()+'.tsv.part')
   if p.exists():b=p.read_bytes()
   else:b=get(self.url,start,end);p.write_bytes(b)
   rows=[dict(zip(self.header,l.split('\t'))) for l in b.decode().splitlines()[1:-1]]
   def key(r):return (({'X':23,'Y':24,'MT':25}.get(r['chromosome'].replace('chr','')) or int(r['chromosome'].replace('chr',''))),int(r['base_pair_location']))
   keys=[key(r) for r in rows];assert keys==sorted(keys),'unsorted source chunk'
   found=[r for r in rows if key(r)==target and r['variant_id']==rsid]
   if found:return found
   if target<keys[0]:right=start-1
   elif target>keys[-1]:left=end+1
   else:return []
   if left>right:return []
  raise RuntimeError('point lookup did not converge')
