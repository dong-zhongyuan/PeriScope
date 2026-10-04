from pathlib import Path
import json,urllib.request,concurrent.futures,csv
W=Path(__file__).resolve().parent;O=W.parents[2]/'delivery/blood_brain_MR_20261003'
d=json.loads((O/'AGES_full_summary_lead_rows.json').read_text())['rows']
def f(x):
 r=x['row'];chrom=r['chromosome'];pos=r['base_pair_location'];u=f'https://rest.ensembl.org/map/human/GRCh37/{chrom}:{pos}..{pos}:1/GRCh38?content-type=application/json'
 p=W/'AGES_raw_leads'/(x['accession']+'.mapping.json')
 if p.exists():m=json.loads(p.read_text())
 else:b=urllib.request.urlopen(u,timeout=30).read();p.write_bytes(b);m=json.loads(b)
 assert len(m['mappings'])==1
 h=m['mappings'][0]['mapped'];assert h['start']==h['end'] and h['strand']==1 and str(h['seq_region_name'])==str(chrom)
 return dict(assay=x['accession'],gene=x['gene'],candidate_antigens=x['original_table_instrument']['candidate_antigens'],chrom=chrom,pos=h['start'],pos37=pos,rsid=r['variant_id'],effect_allele=r['effect_allele'],other_allele=r['other_allele'],beta_exposure=r['beta'],se_exposure=r['standard_error'],eaf_exposure=r['effect_allele_frequency'],p_exposure=r['p_value'],F_exposure=(float(r['beta'])/float(r['standard_error']))**2,source='Gudjonsson2022_AGES_serum',source_url=x['source_url'])
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:out=list(e.map(f,d))
with (O/'AGES_exposure_instruments.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
print('MAPPED',len(out))
