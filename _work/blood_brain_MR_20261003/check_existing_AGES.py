import json,gzip,csv,concurrent.futures
from pathlib import Path
p=Path('/public/home/mengxl/dzy/pd_product_assets/raw/mr');co=json.loads((p/'gene_coords.json').read_text())['grch38']
def f(g):
 name=next(p.glob('pqtl_'+g+'_*.h.tsv.gz'));c=co[g];best=None;n=0
 with gzip.open(name,'rt') as h:
  for r in csv.DictReader(h,delimiter='\t'):
   if r['hm_chrom']!=c['chr']:continue
   try:pos=int(r['hm_pos']);pv=float(r['p_value'])
   except:continue
   if c['start']-1000000<=pos<=c['end']+1000000:
    n+=1
    if best is None or pv<best['p']:best=dict(gene=g,p=pv,row=r,file=str(name))
 return best,n
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as e:out=list(e.map(f,['CD22','CD4']))
Path('/public/home/mengxl/dzy/pd_product_assets/raw/blood_brain_MR_20261003/AGES_CD22_CD4_local_check.json').write_text(json.dumps(out,indent=2))
for r,n in out:print(r['gene'],r['p'],r['row']['hm_rsid'],r['row']['beta'],r['row']['standard_error'],n)
