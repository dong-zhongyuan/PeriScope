#!/usr/bin/env python
"""Fresh ChEMBL REST pull: targets containing antigen, ligand or partner proteins, their mechanisms and molecules."""
import sys; sys.path.insert(0,'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119/code')
from common import *
import gzip, glob, os, time, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor
R=O+'/raw/chembl/rest_20261119/'; os.makedirs(R,exist_ok=True); B='https://www.ebi.ac.uk/chembl/api/data/'; fails=[]
def get(ep,q):
    u=B+ep+'.json?'+urllib.parse.urlencode(q)
    for i in range(7):
        try: return json.loads(urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'research-client/1.0','Accept':'application/json'}),timeout=120).read())
        except Exception as e: err=repr(e)[:100]; time.sleep(4*(i+1))
    fails.append((ep,str(q)[:100],err)); print('FAIL',ep,str(q)[:80],err,flush=True); return None
def pages(ep,key,q):
    out=[]; off=0
    while True:
        r=get(ep,{**q,'limit':500,'offset':off})
        if r is None: break
        out+=r.get(key,[])
        if not r['page_meta'].get('next'): break
        off+=500
    return out
L=pd.read_csv(O+'/curation/ligand_pairs.csv'); P=pd.read_csv(O+'/curation/complex_pairs.csv')
genes=sorted({g for a in antigens() for g in a['genes']}|set(L.ligand_gene)|set(P.partner_gene))
acc={}
for g in genes:
    f=O+f'/raw/opentargets/genes/{g}.json.gz'
    if not os.path.exists(f): continue
    t=(json.load(gzip.open(f,'rt'))['response'].get('data') or {}).get('target') or {}
    a=[x['id'] for x in t.get('proteinIds') or [] if x.get('source')=='uniprot_swissprot']
    if a: acc[g]=a[0]
print('genes',len(genes),'accessions',len(acc),flush=True)
T={}
with ThreadPoolExecutor(4) as ex:
    for res in ex.map(lambda a:pages('target','targets',{'target_components__accession':a}),sorted(set(acc.values()))):
        for t in res:
            if t.get('organism')=='Homo sapiens': T[t['target_chembl_id']]=t
print('targets',len(T),flush=True); json.dump(list(T.values()),gzip.open(R+'targets.json.gz','wt'))
ids=sorted(T); ch=[ids[i:i+40] for i in range(0,len(ids),40)]; M=[]
with ThreadPoolExecutor(4) as ex:
    for res in ex.map(lambda c:pages('mechanism','mechanisms',{'target_chembl_id__in':','.join(c)}),ch): M+=res
print('mechanisms',len(M),flush=True); json.dump(M,gzip.open(R+'mechanisms.json.gz','wt'))
mids=sorted({m['molecule_chembl_id'] for m in M}|{m.get('parent_molecule_chembl_id') for m in M if m.get('parent_molecule_chembl_id')}); ch=[mids[i:i+40] for i in range(0,len(mids),40)]; X=[]
with ThreadPoolExecutor(4) as ex:
    for res in ex.map(lambda c:pages('molecule','molecules',{'molecule_chembl_id__in':','.join(c)}),ch): X+=res
print('molecules',len(X),flush=True); json.dump(X,gzip.open(R+'molecules.json.gz','wt'))
st=get('status',{}); json.dump(dict(retrieved=time.strftime('%Y-%m-%dT%H:%M:%S%z'),status=st,n_genes=len(genes),n_accessions=len(acc),n_targets=len(T),n_mechanisms=len(M),n_molecules=len(X),failures=fails),open(R+'meta.json','w'),indent=1)
print('DONE failures',len(fails),flush=True)
