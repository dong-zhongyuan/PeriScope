from pathlib import Path
import urllib.request,json,csv,time,concurrent.futures,hashlib
W=Path(__file__).resolve().parent;O=W.parents[2]/'delivery/blood_brain_MR_20261003';raw=W/'AGES_raw_leads';raw.mkdir(exist_ok=True)
rows=list(csv.DictReader((O/'AGES_published_cis_instruments.csv').open()))
def retrieve(url,start,end):
 for i in range(3):
  try:
   r=urllib.request.urlopen(urllib.request.Request(url,headers={'Range':f'bytes={start}-{end}'}),timeout=25)
   assert r.status==206
   b=r.read();return b,int(r.headers['Content-Range'].split('/')[-1])
  except Exception:
   if i==2:raise
   time.sleep(1)
def run(x):
 acc=x['Study Accession'];out=raw/(acc+'.json')
 if out.exists():return json.loads(out.read_text())
 n=int(acc[4:]);lo=(n-1)//1000*1000+1
 url=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{lo:08d}-GCST{lo+999:08d}/{acc}/{acc}_buildGRCh37.tsv'
 b,size=retrieve(url,0,4095);header=b.decode().splitlines()[0].split('\t');left=0;right=size-1;target=(int(x['Chr']),int(x['Pos (GRCh37)']))
 trace=[]
 for it in range(35):
  start=max(0,(left+right)//2-16384);end=min(size-1,start+32767)
  block,_=retrieve(url,start,end);lines=block.decode().splitlines()[1:-1]
  parsed=[dict(zip(header,l.split('\t'))) for l in lines]
  def key(r):return (int(r['chromosome']),int(r['base_pair_location']))
  keys=[key(r) for r in parsed];assert keys==sorted(keys),(acc,'unsorted block')
  trace.append(dict(start=start,end=end,first=keys[0],last=keys[-1]))
  matched=[r for r in parsed if key(r)==target and r['variant_id']==x['rsID']]
  if matched:
   assert len(matched)==1;v=matched[0];assert v['effect_allele']==x['EA'];assert abs(float(v['beta'])-float(x['beta']))<1e-5
   result=dict(accession=acc,gene=x['Protein (Entrez symbol)'],source_url=url,source_build='GRCh37',row=v,search_trace=trace,original_table_instrument=x,source_range_sha256=hashlib.sha256(block).hexdigest())
   out.write_text(json.dumps(result,indent=2));print('FOUND',acc,x['Protein (Entrez symbol)'],flush=True);return result
  if target<keys[0]:right=start-1
  elif target>keys[-1]:left=end+1
  else:raise ValueError((acc,'target absent within bracket'))
  assert left<=right
 raise ValueError((acc,'search did not converge'))
results=[];errors=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 jobs={pool.submit(run,x):x for x in rows}
 for j in concurrent.futures.as_completed(jobs):
  try:results.append(j.result())
  except Exception as e:errors.append(dict(accession=jobs[j]['Study Accession'],error=str(e)));print('ERROR',errors[-1],flush=True)
(O/'AGES_full_summary_lead_rows.json').write_text(json.dumps(dict(rows=results,errors=errors),indent=2))
print('COMPLETE',len(results),len(errors),flush=True)
