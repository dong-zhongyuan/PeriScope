from pathlib import Path
import csv,json,concurrent.futures,os
from raw_point import RawPoint
W=Path(__file__).resolve().parent;R=Path('/public/home/mengxl/dzy/pd_product_assets/raw/blood_brain_MR_20261003/indexed_checked') if Path('/public/home/mengxl').exists() else W/'raw/indexed_checked';R.mkdir(parents=True,exist_ok=True)
x=json.loads((W/'indexed_plan.json').read_text())['outcomes'][0];assert x['accession']=='GCST90002436'
u='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90002001-GCST90003000/GCST90002436/GCST90002436_buildGRCh37.tsv'
t=RawPoint(u,R/'raw_blocks');ex=list(csv.DictReader((W/'Sun_exposure_instruments.csv').open()))+list(csv.DictReader((W/'AGES_exposure_instruments.csv').open()))
cache=R/'raw_single';cache.mkdir(exist_ok=True)
def f(e):
 p=cache/(e['rsid']+'.json')
 if p.exists():rr=json.loads(p.read_text())
 else:rr=t.fetch(e['chrom'],e['pos37'],e['rsid']);p.write_text(json.dumps(rr))
 out=[]
 for r in rr:
  out.append(dict(exposure_assay=e['assay'],exposure_gene=e['gene'],candidate_antigens=e['candidate_antigens'],outcome_id=x['id'],outcome_trait=x['trait'],outcome_accession=x['accession'],outcome_source_build='GRCh37',exposure_pos37=e['pos37'],**r))
 print('RAW_POINT',e['gene'],len(out),flush=True);return out,dict(assay=e['assay'],rsid=e['rsid'],rsid_matched_rows=len(out))
rows=[];audit=[];errors=[]
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
 jobs={pool.submit(f,e):e for e in ex}
 for j in concurrent.futures.as_completed(jobs):
  try:r,a=j.result();rows+=r;audit.append(a)
  except Exception as e:errors.append(dict(assay=jobs[j]['assay'],error=str(e)));print('ERROR',errors[-1],flush=True)
(R/'GCST90002436.sentinels.json').write_text(json.dumps(dict(outcome=x,header=t.header,rows=rows,audit=audit,source_url=u,unfiltered_genome_wide_source=True,source_build='GRCh37',errors=errors),indent=2));print('RAW_COMPLETE',len(rows),errors,flush=True)
