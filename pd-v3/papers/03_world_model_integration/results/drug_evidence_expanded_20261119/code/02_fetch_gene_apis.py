#!/usr/bin/env python
"""Fetch Open Targets and DGIdb records per gene; raw responses kept. Usage: 02_fetch_gene_apis.py genes.txt"""
import json, os, sys, time, gzip, urllib.request, urllib.error
from concurrent.futures import ThreadPoolExecutor
O='/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119'
H={'User-Agent':'research-client/1.0','Content-Type':'application/json','Accept':'application/json'}
OT='https://api.platform.opentargets.org/api/v4/graphql'; DG='https://dgidb.org/api/graphql'
os.makedirs(O+'/raw/opentargets/genes',exist_ok=True); os.makedirs(O+'/raw/dgidb/genes',exist_ok=True)
def post(u,q,v=None,tries=6):
    for i in range(tries):
        try:
            b=urllib.request.urlopen(urllib.request.Request(u,json.dumps({'query':q,'variables':v or {}}).encode(),H),timeout=120).read()
            return json.loads(b)
        except Exception as e:
            err=repr(e)
            if isinstance(e,urllib.error.HTTPError):
                err+=' '+e.read()[:300].decode('utf8','replace')
            time.sleep(4*(i+1))
    return {'errors':[{'message':err}]}
Q_MAP='query($t:[String!]!){mapIds(queryTerms:$t,entityNames:["target"]){mappings{term hits{id name entity score}}}}'
Q_T='''query($e:String!){target(ensemblId:$e){id approvedSymbol approvedName biotype
 proteinIds{id source}
 drugAndClinicalCandidates{count rows{id maxClinicalStage
  drug{id name drugType maximumClinicalStage
   crossReferences{source ids}
   parentMolecule{id name}
   mechanismsOfAction{rows{mechanismOfAction actionType targetName references{source ids urls} targets{id approvedSymbol}}}
   drugWarnings{warningType toxicityClass country year description chemblIds references{id source url}}
  }}}}}'''
Q_D='''query($n:[String!]!){genes(names:$n){nodes{name conceptId longName geneCategories{name}
 interactions{interactionScore interactionTypes{type directionality definition}
  interactionAttributes{name value} publications{pmid citation} sources{sourceDbName fullName}
  drug{name conceptId approved immunotherapy antiNeoplastic}}}}}'''
genes=[g.strip() for g in open(sys.argv[1]) if g.strip()]
mp=O+'/raw/opentargets/symbol_to_ensembl.json'; M=json.load(open(mp)) if os.path.exists(mp) else {}
todo=[g for g in genes if g not in M]
for i in range(0,len(todo),50):
    t=todo[i:i+50]; r=post(OT,Q_MAP,{'t':t})
    for m in (r.get('data') or {}).get('mapIds',{}).get('mappings',[]):
        hs=[h for h in m['hits'] if h['entity']=='target']
        ex=[h for h in hs if h['name'].upper()==m['term'].upper()]
        M[m['term']]=(ex or [{}])[0].get('id')
json.dump(M,open(mp,'w'),indent=0)
print('mapped',sum(1 for g in genes if M.get(g)),'of',len(genes),'unmapped',[g for g in genes if not M.get(g)],flush=True)
def ot(g):
    p=f'{O}/raw/opentargets/genes/{g}.json.gz'
    if os.path.exists(p) or not M.get(g): return 0
    r=post(OT,Q_T,{'e':M[g]})
    if r.get('errors') and not r.get('data'): print('OTERR',g,str(r['errors'])[:200],flush=True); return 0
    with gzip.open(p+'.tmp','wt') as f: json.dump({'symbol':g,'ensembl':M[g],'retrieved':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'response':r},f)
    os.replace(p+'.tmp',p); return 1
def dg(g):
    p=f'{O}/raw/dgidb/genes/{g}.json.gz'
    if os.path.exists(p): return 0
    r=post(DG,Q_D,{'n':[g]})
    if r.get('errors') and not r.get('data'): print('DGERR',g,str(r['errors'])[:200],flush=True); return 0
    with gzip.open(p+'.tmp','wt') as f: json.dump({'symbol':g,'retrieved':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'response':r},f)
    os.replace(p+'.tmp',p); return 1
with ThreadPoolExecutor(4) as ex: a=sum(ex.map(ot,genes))
print('OT fetched',a,flush=True)
with ThreadPoolExecutor(3) as ex: b=sum(ex.map(dg,genes))
print('DGIdb fetched',b,flush=True)
meta=post(OT,'{meta{name apiVersion{x y z} dataVersion{year month iteration}}}')
json.dump({'retrieved':time.strftime('%Y-%m-%dT%H:%M:%S%z'),'meta':meta},open(O+'/raw/opentargets/meta.json','w'))
print('DONE',flush=True)
