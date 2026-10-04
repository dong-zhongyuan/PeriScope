import urllib.request,concurrent.futures,gzip,json,csv,hashlib,os,time
from range_tabix import RangeTabix,get
from pathlib import Path
W=Path(__file__).resolve().parent;R=Path('/public/home/mengxl/dzy/pd_product_assets/raw/blood_brain_MR_20261003/indexed_checked') if Path('/public/home/mengxl').exists() else W/'raw/indexed_checked';R.mkdir(parents=True,exist_ok=True)
O=R/'stage_summary';O.mkdir(parents=True,exist_ok=True);os.chdir(R)
plan=json.loads((W/'indexed_plan.json').read_text());ex=list(csv.DictReader((W/'Sun_exposure_instruments.csv').open()))+list(csv.DictReader((W/'AGES_exposure_instruments.csv').open()))

def work(x):
 acc=x['accession'];p=R/(acc+'.sentinels.json')
 if p.exists():
  saved=json.loads(p.read_text())
  if len(saved.get('audit',[]))==len(ex) and not saved.get('errors') and saved.get('engine')=='checked_http_range':return saved
 url=x['url'];meta=get(x['metadata_url']);(R/(acc+'.yaml')).write_bytes(meta);assert 'GRCh38' in meta.decode(),meta
 t=RangeTabix(url,R/'blocks');header=t.header;out=[];audit=[]
 for e in ex:
  pos=int(e['pos']);lines=t.fetch(str(e['chrom']),pos)
  match=[]
  for line in lines:
   r=line
   if e['rsid'] in [r.get('rsid'),r.get('variant_id'),r.get('hm_rsid')]:match.append(r)
  for r in match:out.append(dict(exposure_assay=e['assay'],exposure_gene=e['gene'],candidate_antigens=e['candidate_antigens'],outcome_id=x['id'],outcome_trait=x['trait'],outcome_accession=acc,**r))
  audit.append(dict(assay=e['assay'],rsid=e['rsid'],position_rows=len(lines),rsid_matched_rows=len(match)))
 data=dict(outcome=x,header=header,rows=out,audit=audit,source_url=url,unfiltered_genome_wide_source=True,engine='checked_http_range')
 p.write_text(json.dumps(data,indent=2));print('EXTRACTED',acc,len(out),'of',len(ex),flush=True);return data
results=[];errors=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 jobs={pool.submit(work,x):x for x in plan['outcomes'] if x['accession']!='GCST90002436'}
 for j in concurrent.futures.as_completed(jobs):
  try:results.append(j.result())
  except Exception as e:errors.append(dict(accession=jobs[j]['accession'],error=str(e)));print('ERROR',errors[-1],flush=True)
rows=[r for d in results for r in d['rows']]
if rows:
 with (O/'blood_protein_brain_sentinel_statistics.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary=dict(complete=len(results)==16 and not errors,n_outcomes_completed=len(results),n_requested_assays=len(ex),n_matched_rows=len(rows),expected_pairs=len(ex)*16,errors=errors,outcomes=[dict(accession=x['outcome']['accession'],n_matched=len(x['rows']),n_requested=len(x['audit']),raw_path=str(R/(x['outcome']['accession']+'.sentinels.json'))) for x in results],plan=plan,raw_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in R.glob('*.sentinels.json')})
(O/'brain_data_acquisition.json').write_text(json.dumps(summary,indent=2));print('INDEXED_BRAIN_COMPLETE',summary['complete'],len(rows),flush=True)
