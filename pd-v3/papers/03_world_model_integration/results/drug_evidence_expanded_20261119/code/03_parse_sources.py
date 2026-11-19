#!/usr/bin/env python
"""Parse every raw drug source into one long evidence table plus a drug master table."""
import pandas as pd, numpy as np, json, glob, gzip, re, os, time, urllib.request, collections
O='/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119'
R=O+'/raw/'; I=O+'/interim/'; os.makedirs(I,exist_ok=True)
A='/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/'; S6=A+'six_criteria_20261003/'
H={'User-Agent':'research-client/1.0','Content-Type':'application/json'}
OTU='https://api.platform.opentargets.org/api/v4/graphql'
def nk(s): return re.sub('[^A-Z0-9]','',str(s).upper())
def lst(x):
    if x is None: return []
    if isinstance(x,float) and np.isnan(x): return []
    return list(x)
def pq(d): return pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(R+f'opentargets/bulk/{d}/*.parquet'))],ignore_index=True)
moa=pq('drug_mechanism_of_action'); mol=pq('drug_molecule'); ct=pq('clinical_target'); dw=pq('drug_warning')
# ensembl -> symbol
ens=sorted({e for t in moa.targets for e in lst(t)}|set(ct.targetId))
mp=R+'opentargets/ensembl_to_symbol.json'; E=json.load(open(mp)) if os.path.exists(mp) else {}
todo=[e for e in ens if e not in E]
for i in range(0,len(todo),200):
    q='query($e:[String!]!){targets(ensemblIds:$e){id approvedSymbol}}'
    for k in range(5):
        try:
            r=json.loads(urllib.request.urlopen(urllib.request.Request(OTU,json.dumps({'query':q,'variables':{'e':todo[i:i+200]}}).encode(),H),timeout=120).read()); break
        except Exception as ex: time.sleep(5*(k+1)); r=None
    for t in ((r or {}).get('data') or {}).get('targets',[]): E[t['id']]=t['approvedSymbol']
json.dump(E,open(mp,'w')); print('symbols',len(E),'missing',len([e for e in ens if e not in E]),flush=True)
# ---- drug master from OT
STAGE={'APPROVAL':4,'PHASE_4':4,'WITHDRAWAL':4,'PREAPPROVAL':3.5,'PHASE_3':3,'PHASE_2_3':2.5,'PHASE_2':2,'PHASE_1_2':1.5,'PHASE_1':1,'EARLY_PHASE_1':.5,'IND':0,'PRECLINICAL':0,'UNKNOWN':np.nan}
parent={r.id:(r.parentId if isinstance(r.parentId,str) and r.parentId else r.id) for r in mol.itertuples()}
def root(c):
    s=set()
    while c in parent and parent[c]!=c and c not in s: s.add(c); c=parent[c]
    return c
D={}
def rec(k):
    return D.setdefault(k,dict(drug_key=k,chembl_id=k if k.startswith('CHEMBL') else '',name='',drug_type='',ot_stage=np.nan,ot_stage_label='',chembl_phase=np.nan,gtop_approved=False,dgidb_approved=False,drugcentral=False,withdrawn=False,withdrawn_detail=set(),names=set(),description=''))
for r in mol.itertuples():
    k=root(r.id); d=rec(k)
    if r.id==k or not d['name']: d['name']=r.name; d['drug_type']=r.drugType; d['description']=str(r.description or '')[:300]
    s=STAGE.get(r.maximumClinicalStage,np.nan)
    if not np.isnan(s) and not (d['ot_stage']>=s): d['ot_stage']=s; d['ot_stage_label']=r.maximumClinicalStage
    d['names'].add(r.name)
    for f in (r.synonyms,r.tradeNames):
        for y in lst(f):
            l=y.get('label') if isinstance(y,dict) else y
            if l: d['names'].add(l)
# ---- ChEMBL REST: snapshot of 2026-10-03 plus the fresh pull of 2026-11-19 when present
acc2gene={}
for f in glob.glob(R+'opentargets/genes/*.json.gz'):
    o=json.load(gzip.open(f,'rt')); t=(o['response'].get('data') or {}).get('target') or {}
    for x in t.get('proteinIds') or []:
        if x.get('source')=='uniprot_swissprot': acc2gene[x['id']]=t.get('approvedSymbol') or o['symbol']
CH=[('ChEMBL_REST_20261003',json.load(open(S6+'chembl_mechanisms.json')),json.load(open(S6+'chembl_molecules.json')),json.load(open(S6+'chembl_targets.json')))]
if os.path.exists(R+'chembl/rest_20261119/meta.json'):
    rd=lambda n:json.load(gzip.open(R+f'chembl/rest_20261119/{n}.json.gz','rt'))
    CH.append(('ChEMBL_REST_20261119',rd('mechanisms'),rd('molecules'),rd('targets')))
for _,cm,cmo,ctg in CH:
    for m in cmo:
        h=m.get('molecule_hierarchy') or {}; k=root(h.get('parent_chembl_id') or m['molecule_chembl_id']); d=rec(k)
        if not d['name']: d['name']=m.get('pref_name') or k
        if not d['drug_type']: d['drug_type']=m.get('molecule_type') or ''
        p=m.get('max_phase'); p=float(p) if p not in (None,'') else np.nan
        if not np.isnan(p) and not (d['chembl_phase']>=p): d['chembl_phase']=p
        if m.get('withdrawn_flag'): d['withdrawn']=True; d['withdrawn_detail'].add('ChEMBL withdrawn_flag')
        if m.get('pref_name'): d['names'].add(m['pref_name'])
        for y in m.get('molecule_synonyms') or []:
            if y.get('molecule_synonym'): d['names'].add(y['molecule_synonym'])
# withdrawn from OT warnings and target-level stage
for r in dw.itertuples():
    if str(r.warningType).lower().startswith('withdraw'):
        for c in lst(r.chemblIds):
            d=rec(root(c)); d['withdrawn']=True; d['withdrawn_detail'].add(f'OT warning {r.country} {"" if pd.isna(r.year) else int(r.year)}'.strip())
ev=[]
def add(**k): ev.append(k)
# ---- GtoPdb ligand identities
gl=pd.read_csv(R+'gtopdb/ligands.csv',skiprows=1,low_memory=False,dtype=str).fillna(''); gm=pd.read_csv(R+'gtopdb/ligand_id_mapping.csv',skiprows=1,low_memory=False,dtype=str).fillna('')
strip=lambda s:re.sub('<[^>]+>','',s)
name2={}
def build_names():
    name2.clear()
    for k,d in D.items():
        for n in d['names']:
            x=nk(n)
            if len(x)>=4: name2.setdefault(x,set()).add(k)
build_names()
def resolve(chembl=None,names=()):
    if chembl and str(chembl).startswith('CHEMBL'):
        c=root(chembl); rec(c); return c
    for n in names:
        h=name2.get(nk(n))
        if not h: continue
        if len(h)==1: return next(iter(h))
        c=[k for k in h if k.startswith('CHEMBL')]
        e=[k for k in c if nk(D[k]['name'])==nk(n)]
        if e: return max(e,key=lambda k:(0 if np.isnan(D[k]['ot_stage']) else D[k]['ot_stage']))
    return None
GK={}
glc=gl.set_index('Ligand ID')
for lid,r in glc.iterrows():
    n=strip(r['Name']); inn=strip(r['INN']); syn=[strip(x) for x in r['Synonyms'].split('|') if x]
    k=resolve(r['ChEMBL ID'],[inn,n]+syn) or f'GTOP:{lid}'
    d=rec(k); GK[lid]=k
    if not d['name']: d['name']=inn or n
    if not d['drug_type']: d['drug_type']=r['Type']
    d['names'].update(x for x in [n,inn] if x)
    if r['Approved'].lower() in('yes','true','t'): d['gtop_approved']=True
    if r['Withdrawn'].lower() in('yes','true','t'): d['withdrawn']=True; d['withdrawn_detail'].add('GtoPdb withdrawn')
build_names()
# ---- Open Targets mechanisms (ChEMBL-derived, release 26.09)
def refs(x):
    out=[]
    for y in lst(x):
        for u in lst(y.get('urls')):
            out.append(u)
    return list(dict.fromkeys(out))
for r in moa.itertuples():
    gs=sorted({E.get(e,e) for e in lst(r.targets)})
    for c in lst(r.chemblIds):
        k=root(c); rec(k); u=refs(r.references)
        add(source='OpenTargets_26.09',drug_key=k,src_drug=c,target_genes=';'.join(gs),target_type=str(r.targetType).lower(),target_name=r.targetName or '',action=str(r.actionType or '').upper(),mechanism=r.mechanismOfAction or '',primary=True,pact=np.nan,act_type='',refs=';'.join(u[:6]),link=f'https://platform.opentargets.org/drug/{c}')
# ---- ChEMBL REST mechanisms; genes come from the component accession when it is known
def comps(t):
    o=set()
    for c in t['target_components']:
        g=acc2gene.get(c.get('accession'))
        if g: o.add(g); continue
        sy=[s['component_synonym'] for s in c['target_component_synonyms'] if s['syn_type']=='GENE_SYMBOL']
        if sy: o.add(sy[0])
    return sorted(o)
for src,cm,cmo,ctg in CH:
    tg={x['target_chembl_id']:x for x in ctg}
    for m in cm:
        t=tg.get(m['target_chembl_id'])
        if not t: continue
        k=root(m['molecule_chembl_id']); rec(k)
        add(source=src,drug_key=k,src_drug=m['molecule_chembl_id'],target_genes=';'.join(comps(t)),target_type=str(t['target_type']).lower(),target_name=t.get('pref_name') or '',action=str(m.get('action_type') or '').upper(),mechanism=m.get('mechanism_of_action') or '',primary=m.get('direct_interaction')==1 and m.get('molecular_mechanism')==1,pact=np.nan,act_type='',refs=';'.join(x['ref_url'] for x in m.get('mechanism_refs',[]) if x.get('ref_url'))[:400],link='https://www.ebi.ac.uk/chembl/explore/compound/'+m['molecule_chembl_id'],extra=str(m.get('mechanism_comment') or '')[:200])
# ---- GtoPdb interactions
g=pd.read_csv(R+'gtopdb/interactions.csv',skiprows=1,low_memory=False,dtype=str).fillna('')
gg=g[(g['Target Species'].str.lower().isin(['human','']))]
for r in gg.to_dict('records'):
    tgt=r['Target Gene Symbol'] or r['Target Ligand Gene Symbol']
    if not tgt: continue
    if r['Endogenous'].lower() in('true','t','yes'): continue
    if r['Target Species'].lower()!='human' and not r['Target Ligand Gene Symbol']: continue
    k=GK.get(r['Ligand ID']) or f"GTOP:{r['Ligand ID']}"; rec(k)
    genes=sorted(set(tgt.replace('|',';').split(';')))
    if not D[k]['name']: D[k]['name']=strip(r['Ligand'])
    aff=np.nan
    for c in ('Affinity Median','Affinity High','Affinity Low'):
        try:
            aff=float(r[c]); break
        except: pass
    tl=bool(r['Target Ligand Gene Symbol']) and not r['Target Gene Symbol']
    add(source='GtoPdb_2026.3.1',drug_key=k,src_drug='GTOP:'+r['Ligand ID'],target_genes=';'.join(genes),target_type='ligand_target' if tl else ('protein complex' if len(genes)>1 else 'single protein'),target_name=strip(r['Target'] or r['Target Ligand']),action=('GTOP|'+r['Type']+'|'+r['Action']).upper(),mechanism=strip(r['Action comment']),primary=r['Primary Target'].lower() in('true','t','yes'),pact=aff,act_type=r['Affinity Units'],refs=';'.join('https://pubmed.ncbi.nlm.nih.gov/'+p for p in re.split(r'[|; ]+',r['PubMed ID']) if p)[:300],link='https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId='+r['Ligand ID'])
# ---- DrugCentral
dc=pd.read_csv(R+'drugcentral/drug.target.interaction.tsv.gz',sep='\t',low_memory=False)
dc=dc[dc.ORGANISM.eq('Homo sapiens')&dc.GENE.notna()]
gd=dict(zip(gm['Drug Central ID'],gm['Ligand id']))
for r in dc.to_dict('records'):
    lid=gd.get(str(r['STRUCT_ID'])); k=GK.get(lid) if lid else None
    k=k or resolve(None,[r['DRUG_NAME']]) or f"DC:{r['STRUCT_ID']}"; d=rec(k); d['drugcentral']=True
    if not d['name']: d['name']=r['DRUG_NAME']
    d['names'].add(r['DRUG_NAME'])
    genes=sorted(set(str(r['GENE']).split('|')))
    v=r['ACT_VALUE']; add(source='DrugCentral_2026-09',drug_key=k,src_drug=f"DC:{r['STRUCT_ID']}",target_genes=';'.join(genes),target_type='protein complex' if len(genes)>1 else 'single protein',target_name=r['TARGET_NAME'] or '',action=str(r['ACTION_TYPE'] if pd.notna(r['ACTION_TYPE']) else '').upper(),mechanism='' ,primary=r['MOA']==1,pact=float(v) if pd.notna(v) else np.nan,act_type=str(r['ACT_TYPE'] or ''),refs=str(r['MOA_SOURCE_URL'] if pd.notna(r['MOA_SOURCE_URL']) else (r['ACT_SOURCE_URL'] if pd.notna(r['ACT_SOURCE_URL']) else ''))[:300],link=f"https://drugcentral.org/drugcard/{r['STRUCT_ID']}")
build_names()
# ---- DGIdb
for f in sorted(glob.glob(R+'dgidb/genes/*.json.gz')):
    o=json.load(gzip.open(f,'rt'))
    for n in ((o['response'].get('data') or {}).get('genes') or {}).get('nodes',[]):
        for x in n['interactions']:
            dr=x['drug']; cid=str(dr.get('conceptId') or ''); c=cid.split(':',1)[1].upper() if cid.lower().startswith('chembl:') else None
            k=resolve(c,[dr['name']]) or 'DGI:'+nk(dr['name']); d=rec(k)
            if not d['name']: d['name']=dr['name']
            d['names'].add(dr['name'])
            if dr.get('approved'): d['dgidb_approved']=True
            ty=[(t.get('type') or '').upper() for t in x['interactionTypes']]; di=sorted({(t.get('directionality') or '').upper() for t in x['interactionTypes']}-{''})
            so=sorted({s['sourceDbName'] for s in x['sources']}); pm=[str(p['pmid']) for p in x['publications'] if p.get('pmid')]
            add(source='DGIdb_live',drug_key=k,src_drug=cid or dr['name'],target_genes=n['name'],target_type='single protein',target_name=n.get('longName') or '',action='DGI|'+';'.join(sorted(set(ty)))+'|'+';'.join(di),mechanism='sources='+','.join(so),primary=False,pact=np.nan,act_type='',refs=';'.join('https://pubmed.ncbi.nlm.nih.gov/'+p for p in pm[:6]),link='https://dgidb.org/drugs/'+str(cid or dr['name']).replace(' ','%20'))
EV=pd.DataFrame(ev); EV['extra']=EV.get('extra','').fillna('') if 'extra' in EV else ''
# target-level OT stage
ctg_stage=collections.defaultdict(lambda:np.nan); ct_w=set()
for r in ct.itertuples():
    k=root(r.drugId); g=E.get(r.targetId,r.targetId); s=STAGE.get(r.maxClinicalStage,np.nan)
    if not np.isnan(s) and not (ctg_stage[(k,g)]>=s): ctg_stage[(k,g)]=s
    if r.maxClinicalStage=='WITHDRAWAL': ct_w.add(k)
pd.DataFrame([dict(drug_key=k[0],gene=k[1],ot_target_stage=v) for k,v in ctg_stage.items()]).to_csv(I+'ot_target_stage.csv.gz',index=False)
for k in ct_w:
    d=rec(k); d['withdrawn_detail'].add('OT clinical_target WITHDRAWAL')
DM=pd.DataFrame([{**{a:b for a,b in d.items() if a not in('names','withdrawn_detail')},'withdrawn_detail':';'.join(sorted(d['withdrawn_detail'])),'names':'|'.join(sorted(d['names']))[:600]} for d in D.values()])
DM.to_csv(I+'drug_master.csv.gz',index=False); EV.to_csv(I+'evidence_long.csv.gz',index=False)
print(DM.shape,EV.shape); print(EV.source.value_counts().to_dict()); print(EV[EV.source.str.startswith('Open')].target_type.value_counts().to_dict()); print(EV[EV.source.str.startswith('ChEMBL')].target_type.value_counts().to_dict())
