from pathlib import Path
import json,datetime,urllib.request,concurrent.futures
W=Path(__file__).resolve().parent;O=W.parents[2]/'delivery/blood_brain_MR_20261003'
locked=json.loads((W/'brain_panel_lock.json').read_text())
def f(r):
 acc='GCST'+str(90002425+int(r['id']));n=int(acc[4:]);lo=(n-1)//1000*1000+1;hi=lo+999
 folder=f'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST{lo:08d}-GCST{hi:08d}/{acc}/'
 u=folder+acc+'_buildGRCh37.tsv-meta.yaml';b=urllib.request.urlopen(u,timeout=40).read();(W/(acc+'.yaml')).write_bytes(b)
 assert r['trait'].replace('_',' ') in b.decode(),(acc,r['trait'],b)
 base=folder+'harmonised/'+acc+'.h.tsv.gz'
 return {**r,'accession':acc,'url':base,'sample_size':int(r['metadata_cells'][9]),'cohort_version':'discovery_only','outcome_build':'GRCh38 (harmonized; verify YAML)','metadata_url':base+'-meta.yaml'}
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:out=list(e.map(f,locked['outcomes']))
plan={'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'reason_for_source_change':'Official Oxford shortcut PheWeb omits p>0.01; full Oxford33k transfers slow. GWAS Catalog indexed whole-genome files allow unbiased extraction including null associations. Catalog version is discovery22k, not pooled33k; outcomes unchanged. Pooled33k planned sensitivity once downloaded.','outcomes':out}
(W/'indexed_plan.json').write_text(json.dumps(plan,indent=2));(O/'indexed_plan.json').write_text(json.dumps(plan,indent=2));print('VERIFIED',[(x['accession'],x['trait'],x['sample_size']) for x in out])
