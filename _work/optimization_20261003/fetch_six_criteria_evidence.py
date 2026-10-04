"""Cache public target identity and clinical mechanism evidence for all antigens."""
from pathlib import Path
import argparse, concurrent.futures, hashlib, json, time
import urllib.request, urllib.parse, urllib.error

def main(result_dir):
    root = Path(result_dir)
    out = root / 'six_criteria_20261003'
    cache = out / 'reference_cache'
    cache.mkdir(parents=True, exist_ok=True)
    log = []
    def get(url):
        p = cache / (hashlib.sha256(url.encode()).hexdigest()+'.json')
        if p.exists():
            data = json.loads(p.read_text())
        else:
            for attempt in range(3):
                try:
                    req = urllib.request.Request(url, headers={'User-Agent':'PeriScope-evidence-review/1.0'})
                    with urllib.request.urlopen(req, timeout=55) as r:
                        data = json.load(r)
                    p.write_text(json.dumps(data))
                    break
                except Exception:
                    if attempt == 2: raise
                    time.sleep(2)
        log.append(dict(url=url, file=str(p.relative_to(out)), sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
        return data
    def pages(endpoint,params,key):
        url='https://www.ebi.ac.uk/chembl/api/data/'+endpoint+'.json?'+urllib.parse.urlencode(dict(params,limit=1000))
        result=[]
        while url:
            d=get(url);result.extend(d[key]);n=d['page_meta']['next']
            url='https://www.ebi.ac.uk'+n if n and n.startswith('/') else n
        return result
    registry=json.loads((root/'candidate_registry.json').read_text())
    genes=sorted({g for a in registry if not a['is_control'] for g in a['genes']})
    def uniprot(chunk):
        q='('+' OR '.join('gene_exact:'+g for g in chunk)+') AND organism_id:9606 AND reviewed:true'
        u='https://rest.uniprot.org/uniprotkb/search?'+urllib.parse.urlencode(dict(query=q,format='json',size=500,
          fields='accession,gene_primary,cc_subcellular_location,ft_topo_dom,protein_name,cc_tissue_specificity'))
        d=get(u)
        assert len(d['results'])<500
        return d['results']
    chunks=lambda x,n:[x[i:i+n] for i in range(0,len(x),n)]
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        uni=sum(pool.map(uniprot,chunks(genes,45)),[])
    uni=list({x['primaryAccession']:x for x in uni}.values())
    (out/'uniprot_targets.json').write_text(json.dumps(uni))
    print('UNIPROT',len(uni),flush=True)
    acc=[x['primaryAccession'] for x in uni]
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        targets=sum(pool.map(lambda x:pages('target',{'target_components__accession__in':','.join(x)},'targets'),chunks(acc,35)),[])
    targets=list({x['target_chembl_id']:x for x in targets if x.get('tax_id')==9606}.values())
    (out/'chembl_targets.json').write_text(json.dumps(targets))
    print('TARGETS',len(targets),flush=True)
    tids=[x['target_chembl_id'] for x in targets]
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        mechs=sum(pool.map(lambda x:pages('mechanism',{'target_chembl_id__in':','.join(x)},'mechanisms'),chunks(tids,35)),[])
    mechs=list({x['mec_id']:x for x in mechs}.values())
    (out/'chembl_mechanisms.json').write_text(json.dumps(mechs))
    print('MECHANISMS',len(mechs),flush=True)
    mids=sorted({x['molecule_chembl_id'] for x in mechs})
    with concurrent.futures.ThreadPoolExecutor(4) as pool:
        mols=sum(pool.map(lambda x:pages('molecule',{'molecule_chembl_id__in':','.join(x)},'molecules'),chunks(mids,35)),[])
    (out/'chembl_molecules.json').write_text(json.dumps(mols))
    (out/'reference_manifest.json').write_text(json.dumps(dict(retrieved_date='2026-10-03',sources=sorted(log,key=lambda x:x['url'])),indent=2))
    print('MOLECULES',len(mols),'REQUESTS',len(log),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--result-dir',required=True)
    main(ap.parse_args().result_dir)
