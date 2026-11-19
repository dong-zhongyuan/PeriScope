#!/usr/bin/env python
"""openFDA drug labels (mechanism of action section) for approved biologics in scope; raw fields saved."""
import sys; sys.path.insert(0,'/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119/code')
from common import *
import urllib.request, urllib.parse, urllib.error, time, os, gzip
I=O+'/interim/'; C=O+'/curation/'; R=O+'/raw/openfda/'; os.makedirs(R,exist_ok=True)
EV=pd.read_csv(I+'evidence_long.csv.gz',low_memory=False).fillna({'target_genes':''}); DM=pd.read_csv(I+'drug_master.csv.gz',low_memory=False)
L=pd.read_csv(C+'ligand_pairs.csv'); P=pd.read_csv(C+'complex_pairs.csv')
ag={g for a in antigens() for g in a['genes']}; scope=ag|set(L.ligand_gene)|set(P.partner_gene)
e=EV[[bool(scope&set(s.split(';'))) for s in EV.target_genes]]
DM=DM[DM.drug_key.isin(e.drug_key)]
ap=(DM.ot_stage>=4)|(DM.chembl_phase>=4)|DM.gtop_approved.fillna(False)
bio=DM.drug_type.astype(str).str.contains('Antibody|Protein|Enzyme|Peptide|Oligonucleotide|Unknown|Cell|Gene',case=False)
D=DM[ap&bio]; print('drugs to query',len(D),flush=True)
out=R+'labels.jsonl'; done=set()
if os.path.exists(out):
    for l in open(out): done.add(json.loads(l)['drug_key'])
def q(n):
    n=n.replace('"','')
    s=f'openfda.generic_name:"{n}" OR openfda.substance_name:"{n}" OR openfda.brand_name:"{n}"'
    u='https://api.fda.gov/drug/label.json?'+urllib.parse.urlencode({'search':s,'limit':5})
    for k in range(4):
        try: return json.loads(urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'research-client/1.0'}),timeout=60).read()).get('results',[])
        except urllib.error.HTTPError as ex:
            if ex.code==404: return []
            time.sleep(5*(k+1))
        except Exception: time.sleep(5*(k+1))
    return None
from concurrent.futures import ThreadPoolExecutor
import threading
lock=threading.Lock(); f=open(out,'a')
def work(r):
    if r.drug_key in done: return
    hits=q(str(r.name).lower())
    if hits is None: print('fail',r.name,flush=True); return
    kept=[]
    for h in hits:
        o=h.get('openfda',{}); sid=o.get('spl_set_id') or h.get('set_id') or ['']
        sid=sid[0] if isinstance(sid,list) else sid
        kept.append(dict(set_id=sid,effective=h.get('effective_time'),generic=o.get('generic_name'),brand=o.get('brand_name'),moa=(h.get('mechanism_of_action') or [''])[0][:3000]))
    with lock:
        f.write(json.dumps(dict(drug_key=r.drug_key,name=r.name,retrieved=time.strftime('%Y-%m-%d'),hits=kept))+'\n'); f.flush()
with ThreadPoolExecutor(5) as ex: list(ex.map(work,list(D.itertuples())))
print('DONE',flush=True)
