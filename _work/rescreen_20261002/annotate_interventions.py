"""Refresh public drug-gene evidence for every mapped panel gene."""
import json,time,urllib.request,hashlib
from pathlib import Path
import pandas as pd

O=Path('/public/home/mengxl/dzy/pd_product_assets/results/rescreen_20261002');URL='https://dgidb.org/api/graphql'

def main():
    registry=json.loads((O/'candidate_registry.json').read_text());genes=sorted({g for r in registry if not r['is_control'] for g in r['genes']})
    raw=O/'drug_database_queries';raw.mkdir(exist_ok=True);rows=[];statuses=[];used=[]
    for start in range(0,len(genes),20):
        batch=genes[start:start+20];query='{ genes(names: '+json.dumps(batch)+') { nodes { name conceptId interactions { drug { name conceptId approved } interactionScore interactionTypes { type directionality } publications { pmid } sources { sourceDbName } } } } }'
        query_hash=hashlib.sha256(query.encode()).hexdigest()
        dest=raw/f'query_{query_hash}.json';used.append(dest)
        if dest.exists():response=json.loads(dest.read_text())
        else:
            try:
                req=urllib.request.Request(URL,data=json.dumps({'query':query}).encode(),headers={'Content-Type':'application/json','User-Agent':'PeriScope-research-reanalysis'})
                response=json.load(urllib.request.urlopen(req,timeout=45))
                if response.get('errors'):raise RuntimeError(response['errors'])
                dest.write_text(json.dumps(response,indent=2));time.sleep(.3)
            except Exception as e:
                statuses.extend(dict(gene=g,status='query_failed',error=str(e)) for g in batch);continue
        nodes=response.get('data',{}).get('genes',{}).get('nodes',[]);found={n['name'] for n in nodes}
        for n in nodes:
            interactions=n.get('interactions',[]);statuses.append(dict(gene=n['name'],status='returned',n_interactions=len(interactions)))
            for r in interactions:
                rows.append(dict(gene=n['name'],drug=r['drug']['name'],drug_id=r['drug']['conceptId'],drug_approved=r['drug'].get('approved'),
                    interaction_score=r.get('interactionScore'),interaction_types=';'.join(t['type'] for t in r.get('interactionTypes',[])),
                    directionality=';'.join(str(t.get('directionality')) for t in r.get('interactionTypes',[])),
                    pmids=';'.join(str(p['pmid']) for p in r.get('publications',[])),sources=';'.join(t['sourceDbName'] for t in r.get('sources',[]))))
        statuses.extend(dict(gene=g,status='not_returned') for g in batch if g not in found)
    pd.DataFrame(rows).to_csv(O/'drug_gene_interactions.csv',index=False);pd.DataFrame(statuses).to_csv(O/'drug_gene_query_status.csv',index=False)
    (O/'drug_annotation_source.json').write_text(json.dumps(dict(url=URL,documentation='https://dgidb.org/api',queried=time.strftime('%Y-%m-%d %H:%M:%S'),n_genes=len(genes),n_interactions=len(rows),
        approval_semantics='drug regulatory status recorded by DGIdb, not approval for Parkinson disease',
        registry_sha256=hashlib.sha256((O/'candidate_registry.json').read_bytes()).hexdigest(),
        query_files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in used if p.exists()}),indent=2))
    print('DRUG ANNOTATION',len(genes),len(rows),flush=True)

if __name__=='__main__':main()
