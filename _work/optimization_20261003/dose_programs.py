"""Signed programs from regenerated five-point dose slopes. Full candidate family; training gene universe; complete enrichment overlaps; discovery-only definitions."""
import os
import sys
import json

import numpy as np

A = '/public/home/mengxl/dzy/pd_product_assets'
CACHE = A + '/results/optimization_20261003/curves'
sys.path.insert(0, '/public/home/mengxl/dzy/pd_product/src')
PMAX = int(os.environ.get('RECIPE_PURE_MAX', '51'))
SEEDS = [int(x) for x in os.environ.get('CCWM_SEEDS', '42,43,44,45,46,47,48,49,50,51').split(',')]
SEEDS = [s for s in SEEDS if s <= PMAX]
QUANTILES = [0.05, 0.25, 0.50, 0.75, 0.95]

import gseapy as gp

# The original top-three loop below is generalized to the complete candidate family.
from pathlib import Path
import hashlib,urllib.request
from scipy.stats import hypergeom
from input_utils import canonical_symbols

O=Path(A)/'results/optimization_20261003';I=Path(A)/'interim/rescreen_20261002'
REF=Path(__file__).parent/'reference'
LIBRARIES=['GO_Biological_Process_2023','KEGG_2021_Human']

def bh(ps):
    ps=np.asarray(ps,float)
    if not len(ps):return ps.copy()
    order=np.argsort(ps);q=np.empty(len(ps))
    q[order]=np.minimum(1.,np.minimum.accumulate((ps[order]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1])
    return q

def atomic_json(path,value):
    temporary=path.with_name(path.name+'.pending');temporary.write_text(json.dumps(value,indent=2));temporary.replace(path)

def load_libraries(genes):
    background=set(genes);out=[];provenance=[]
    for library in LIBRARIES:
        path=REF/(library+'.gmt');url='https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName='+library
        if not path.exists(): urllib.request.urlretrieve(url,path)
        for line in path.read_text().splitlines():
            x=line.split('\t');members=set(canonical_symbols(x[2:]))&background
            if 5<=len(members)<=500:out.append((library,x[0],members))
        provenance.append(dict(library=library,url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    (O/'enrichment_sources.json').write_text(json.dumps(provenance,indent=2))
    return out


def main():
    registry=json.loads((O/'candidate_registry.json').read_text());programs={};folds={};all_rows=[];gene_rankings={}
    paths=sorted(Path(CACHE).glob('*__seed42.npz'))
    if not paths:raise RuntimeError('No regenerated curves')
    genes=np.load(paths[0])['genes'].astype(str).tolist();gidx={g:i for i,g in enumerate(genes)}
    libraries=load_libraries(genes)
    sizes=np.array([len(t[2]) for t in libraries]);term_ix=[np.array([gidx[g] for g in t[2]]) for t in libraries]
    from scipy import sparse
    term_matrix=sparse.csr_matrix((np.ones(sum(len(ix) for ix in term_ix),np.int32),(np.repeat(np.arange(len(term_ix)),[len(ix) for ix in term_ix]),np.concatenate(term_ix))),shape=(len(term_ix),len(genes)))
    for f in paths:
        axis=f.name.rsplit('__seed',1)[0]
        arr=[np.load(Path(CACHE)/f'{axis}__seed{s}.npz') for s in SEEDS]
        C=np.stack([z['C'] for z in arr]);u=arr[0]['u'];ps=arr[0]['proteins'].astype(str).tolist();control=arr[0]['is_control']
        xx=u-u.mean(1,keepdims=True)
        D=np.einsum('spqg,pq->spg',C,xx)/np.maximum((xx**2).sum(1)[None,:,None],1e-12)
        np.savez_compressed(O/'curves'/f'{axis}__slopes.npz',D=D,genes=genes,proteins=ps,seeds=SEEDS,is_control=control,u=u)
        fold_sets={}
        for fold,si in [('A',[0,1,2]),('B',[3,4])]:
            ens=D[si].mean(0)
            for pi,pname in enumerate(ps):
                if control[pi]:continue
                for direction,sgn in [('up',1),('down',-1)]:
                    rd=sgn*ens[pi];gi=np.argsort(-rd)[:200];gi=gi[rd[gi]>0]
                    if len(gi)<5:continue
                    mask=np.zeros(len(genes),bool);mask[gi]=True
                    overlap=term_matrix@mask.astype(np.int32)
                    p=hypergeom.sf(overlap-1,len(genes),sizes,len(gi));q=bh(p)
                    selected=[j for j in np.argsort(q) if q[j]<.05 and overlap[j]>=3][:6]
                    combo=pname+'__'+axis+'__'+direction
                    members=sorted(set.union(*(libraries[j][2]&{genes[g] for g in gi} for j in selected))) if selected else []
                    if len(members)>=5:fold_sets.setdefault(combo,{})[fold]=members
                    for j in selected:
                        all_rows.append(dict(combo=combo,fold=fold,library=libraries[j][0],term=libraries[j][1],
                            p=float(p[j]),q=float(q[j]),overlap=int(overlap[j]),background=len(genes),query=len(gi),
                            genes=';'.join(sorted(libraries[j][2]&{genes[g] for g in gi}))))
        # Discovery union uses only seeds 42..46; confirmation seeds 47..51 remain unused here.
        for combo,fs in fold_sets.items():
            programs[combo]=sorted(set.union(*(set(v) for v in fs.values())));folds[combo]=fs
            pname=combo.split('__',1)[0];direction=combo.rsplit('__',1)[1];sign=1 if direction=='up' else -1
            response=D[[i for i,seed in enumerate(SEEDS) if seed<=46],ps.index(pname)].mean(0)
            gene_rankings[combo]=sorted(programs[combo],key=lambda g:-sign*response[gidx[g]])
        print('programs',axis,len(fold_sets),flush=True)
    atomic_json(O/'programs.json',programs)
    atomic_json(O/'program_gene_rankings.json',gene_rankings)
    atomic_json(O/'program_folds.json',folds)
    fingerprint=O/'discovery_definition_fingerprint.json'
    if not fingerprint.exists():atomic_json(fingerprint,{n:hashlib.sha256((O/n).read_bytes()).hexdigest() for n in ['programs.json','program_gene_rankings.json','program_folds.json']})
    import pandas as pd
    pd.DataFrame(all_rows,columns=['combo','fold','library','term','p','q','overlap','background','query','genes']).to_csv(O/'enrichment_terms.csv',index=False)
    (O/'program_definition.json').write_text(json.dumps(dict(
        discovery_seeds=[42,43,44,45,46],confirmation_seeds=[47,48,49,50,51],
        score='five measured-dose-point OLS slope; mean across discovery models',
        enrichment='two signed top-200 queries per candidate; hypergeometric ORA on model gene universe; BH across GO BP+KEGG terms',
        membership='all overlap genes from at most six q<0.05 terms per discovery fold; union of folds',
        n_programs=len(programs),candidate_restriction='all biological antigen groups'),indent=2))

if __name__=='__main__':
    from evidence_cache import acquire
    lock=acquire(__file__)
    main()
