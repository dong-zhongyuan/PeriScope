"""Full ADT registry, clone groups and measured training-subtype doses."""
import re,json
from pathlib import Path
import numpy as np
import pandas as pd

A=Path('/public/home/mengxl/dzy/pd_product_assets')
I=A/'interim/rescreen_20261002';O=A/'results/rescreen_20261002'
QUANTILES=[.05,.25,.5,.75,.95]

def weighted_quantiles(x,donor):
    w=np.zeros(len(x))
    for d in np.unique(donor): w[donor==d]=1/np.sum(donor==d)/len(np.unique(donor))
    o=np.argsort(x)
    return np.interp(QUANTILES,np.cumsum(w[o])-w[o]/2,x[o])

def atomic_text(path,text):
    temporary=path.with_name(path.name+'.pending');temporary.write_text(text);temporary.replace(path)

def main():
    z=np.load(I/'citeseq.npz');proteins=z['proteins'].astype(str);states=z['state_vocab'].astype(str)
    h=pd.read_csv(Path(__file__).parent/'reference/hgnc_complete_set.txt',sep='\t',dtype=str).fillna('')
    lookup={}
    for _,r in h.iterrows():
        for a in [r.symbol]+r.alias_symbol.split('|')+r.prev_symbol.split('|'):
            if a: lookup.setdefault(a.upper().replace('-',''),set()).add(r.symbol)
    groups={}
    for i,p in enumerate(proteins):
        # Only remove suffixes when the matched pair of antibody clone channels exists.
        root=re.sub(r'-[12]$','',p) if re.match(r'^CD[0-9]+[A-Za-z]*-[12]$',p) else p
        if not (root+'-1' in proteins and root+'-2' in proteins): root=p
        groups.setdefault(root,[]).append(i)
    curated={r['target']:r for r in json.loads((Path(__file__).parent/'reference/antigen_identity_curation.json').read_text())['records']}
    records=[]
    for name,ii in groups.items():
        control=bool(re.match(r'Ra[gt]-IgG',name))
        exact=h[h.symbol.str.upper()==name.upper().replace('-','')].symbol.tolist()
        match=set(exact) if exact else lookup.get(name.upper().replace('-',''),set())
        if name=='GP130':match={'IL6ST'}
        established={'Galectin-9':['LGALS9'],'Podoplanin':['PDPN'],'IgD':['IGHD'],'IgM':['IGHM'],'CD199':['CCR9'],'CD307c/FcRL3':['FCRL3']}
        complexes={'HLA-DR':['HLA-DRA','HLA-DRB1','HLA-DRB3','HLA-DRB4','HLA-DRB5'],'CD3':['CD3D','CD3E','CD3G','CD247'],'CD11a/CD18':['ITGAL','ITGB2']}
        if name in established:match=set(established[name])
        if name in complexes:match=set(complexes[name])
        if name in ['CD15','CD57']:match=set()
        if name in ['CD45RA','CD45RB','CD45RO']: match=lookup.get('CD45',set())
        if name in curated:match=set(curated[name]['genes'])
        records.append(dict(target=name,channels=ii,channel_names=proteins[ii].tolist(),is_control=control,
                            genes=sorted(match),mapping='defined_complex' if name in complexes else 'unique' if len(match)==1 else 'multiple' if match else 'unresolved',
                            clone_group_size=len(ii)))
        if name in curated:
            records[-1].update(mapping=curated[name]['mapping'],antibody_catalog=curated[name]['catalog'],antibody_clone=curated[name]['clone'],identity_note=curated[name]['note'])
    O.mkdir(parents=True,exist_ok=True);atomic_text(O/'candidate_registry.json',json.dumps(records,indent=2))
    train=~(z['is_test']|z['is_val']);Y=z['Y'];doses={};rows=[]
    for st in ['cDC','classical_mono','nonclassical_mono','pDC']:
        ix=np.where(train & (z['state']==list(states).index(st)))[0]
        if len(ix)<25: continue
        q=np.stack([weighted_quantiles(Y[ix,j],z['donor'][ix]) for j in range(len(proteins))])
        doses[st]=q.tolist()
        for r in records:
            ii=r['channels'];iqr=float(np.mean(q[ii,3]-q[ii,1]))
            rows.append(dict(target=r['target'],blood_state=st,n_cells=len(ix),n_donors=len(np.unique(z['donor'][ix])),
                dose_iqr=iqr,dose_range=float(np.mean(q[ii,-1]-q[ii,0])),is_control=r['is_control'],
                eligibility='measured_range' if iqr>1e-6 else 'no_measured_range',genes=';'.join(r['genes']),mapping=r['mapping']))
    atomic_text(O/'dose_definitions.json',json.dumps(dict(quantiles=QUANTILES,proteins=proteins.tolist(),doses=doses,
        training_donors=z['donor_names'][np.unique(z['donor'][train])].tolist(),weighting='equal donors',
        intervention='target CLR clamp, equal compensation of other channels to preserve zero sum'),indent=2))
    atomic_text(O/'all_candidate_input_support.csv',pd.DataFrame(rows).to_csv(index=False))
    print('REGISTRY',len(records),'antigens; biological',sum(not r['is_control'] for r in records),flush=True)

if __name__=='__main__': main()
