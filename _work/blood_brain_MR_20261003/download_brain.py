from pathlib import Path
import urllib.request,hashlib,json,time,concurrent.futures,gzip
W=Path(__file__).resolve().parent
O=Path('/public/home/mengxl/dzy/pd_product_assets/raw/blood_brain_MR_20261003');O.mkdir(parents=True,exist_ok=True)
plan=json.loads((W/'brain_panel_lock.json').read_text())
base='https://open.oxcin.ox.ac.uk/ukbiobank/big40/release2/'
tasks=[('BIG40_variants.txt.gz',base+'variants.txt.gz')]+[(x['id']+'.txt.gz',base+'stats33k/'+x['id']+'.txt.gz') for x in plan['outcomes']]
def fetch(item):
 name,url=item;p=O/name;info=O/(name+'.json')
 if p.exists() and info.exists():return json.loads(info.read_text())
 for attempt in range(3):
  try:
   tmp=p.with_suffix('.partial');h=hashlib.sha256();n=0
   with urllib.request.urlopen(url,timeout=90) as r,tmp.open('wb') as f:
    total=int(r.headers.get('Content-Length',0))
    while True:
     b=r.read(2**20)
     if not b:break
     f.write(b);h.update(b);n+=len(b)
   assert n==total or total==0
   tmp.rename(p)
   with gzip.open(p,'rt') as f:header=f.readline().strip();first=f.readline().strip()
   data=dict(name=name,url=url,bytes=n,sha256=h.hexdigest(),header=header,first_row=first,completed=time.time())
   info.write_text(json.dumps(data,indent=2));print('DOWNLOADED',name,n,flush=True);return data
  except Exception as e:
   print('RETRY',name,attempt,str(e),flush=True);time.sleep(2)
 raise RuntimeError(name)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(fetch,tasks))
(O/'download_manifest.json').write_text(json.dumps(dict(complete=True,files=results,protocol=plan),indent=2))
print('ALL_BRAIN_OUTCOMES_DOWNLOADED',flush=True)
