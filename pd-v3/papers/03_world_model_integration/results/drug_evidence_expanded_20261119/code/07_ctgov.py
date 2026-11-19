#!/usr/bin/env python
"""ClinicalTrials.gov v2 intervention lookups for manually curated priority drugs."""
import urllib.request, urllib.parse, json, os, time, sys
O='/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119'
os.makedirs(O+'/raw/clinicaltrials',exist_ok=True)
names=[l.strip() for l in open(sys.argv[1]) if l.strip()]; out=O+'/raw/clinicaltrials/lookup.json'
R=json.load(open(out)) if os.path.exists(out) else {}
F='NCTId,BriefTitle,Phase,OverallStatus,InterventionName,StartDate'
for n in names:
    if n in R: continue
    u='https://clinicaltrials.gov/api/v2/studies?'+urllib.parse.urlencode({'query.intr':n,'fields':F,'pageSize':50,'format':'json'})
    for k in range(4):
        try:
            d=json.loads(urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'research-client/1.0','Accept':'application/json'}),timeout=60).read()); break
        except Exception as e:
            d=None; time.sleep(5*(k+1))
    if d is None: print('fail',n,flush=True); continue
    S=[]
    for s in d.get('studies',[]):
        p=s['protocolSection']; ids=p['identificationModule']; st=p['statusModule']; ds=p.get('designModule',{}); iv=p.get('armsInterventionsModule',{}).get('interventions',[])
        S.append(dict(nct=ids['nctId'],title=ids.get('briefTitle','')[:90],phase=ds.get('phases',[]),status=st.get('overallStatus'),start=(st.get('startDateStruct') or {}).get('date'),interventions=[i.get('name') for i in iv][:6]))
    R[n]=dict(retrieved=time.strftime('%Y-%m-%d'),n=len(S),studies=S); json.dump(R,open(out,'w')); time.sleep(.4)
    ph=sorted({x for s in S for x in s['phase']}); print(n,'|',len(S),ph,[s['nct'] for s in S[:4]],flush=True)
print('DONE')
