"""Read counts without loading duplicate h5ad layers; fixed donor splits."""
from pathlib import Path
import json
import h5py
import anndata as ad
from anndata.experimental import read_elem, sparse_dataset
import numpy as np
import pandas as pd
from scipy import sparse

A = Path('/public/home/mengxl/dzy/pd_product_assets')
ROOT = A / 'interim/rescreen_20261002'
REPORT = json.loads((A / 'interim/ccwm_v23/build_report.json').read_text())
BRAIN_VAL = ['SN-2569', 'SN-4956']
BLOOD_VAL = ['b223:BL-E2', 'b223:BL-N2', 'bmb:CTRL2', 'bmb:PD2']
BLOOD_TEST = ['bmb:CTRL1', 'bmb:PD1']
COND = {'Disease':1, 'PD':1, 'PDD':1, 'Early PD':1, 'Late PD':1,
        'Ctrl':0, 'Control':0, 'Normal':0, 'CTRL':0}

def metadata(path):
    with h5py.File(path) as f:
        return read_elem(f['obs']), read_elem(f['var'])

_HGNC_RENAME = None

def canonical_symbols(values):
    global _HGNC_RENAME
    if _HGNC_RENAME is None:
        h=pd.read_csv(Path(__file__).parent/'reference/hgnc_complete_set.txt',sep='\t',dtype=str).fillna('')
        approved=set(h.symbol);previous={}
        for _,r in h.iterrows():
            for old in r.prev_symbol.split('|'):
                if old and old not in approved:previous.setdefault(old,set()).add(r.symbol)
        _HGNC_RENAME={old:next(iter(new)) for old,new in previous.items() if len(new)==1}
    return np.array([_HGNC_RENAME.get(str(v),str(v)) for v in values],dtype=str)

def symbols(var):
    return canonical_symbols(var['gene_symbol'].astype(str).to_numpy() if 'gene_symbol' in var else var.index.astype(str).to_numpy())

def log_count_blocks(path, rows=None, selected=None, batch=4096):
    obs, var = metadata(path)
    sy = symbols(var)
    unique = sorted(set(sy)) if selected is None else list(selected)
    pos = {g:i for i,g in enumerate(unique)}
    ri = np.array([i for i,g in enumerate(sy) if g in pos])
    ci = np.array([pos[sy[i]] for i in ri])
    # Sum duplicate gene symbols at count level before normalization/log1p.
    collapse = sparse.csr_matrix((np.ones(len(ri),np.float32),(ri,ci)),shape=(len(sy),len(unique)))
    rows = np.arange(len(obs)) if rows is None else np.asarray(rows)
    with h5py.File(path) as f:
        raw = sparse_dataset(f['layers/counts'])
        for start in range(0,len(rows),batch):
            rr = rows[start:start+batch]
            x = raw[rr].astype(np.float32)
            lib = np.asarray(x.sum(1)).ravel().clip(1)
            x = (x @ collapse).multiply((1e4/lib)[:,None]).tocsr()
            x.data = np.log1p(x.data)
            yield rr, x, unique

def load_expression(path, selected):
    obs,var = metadata(path)
    mat = np.empty((len(obs),len(selected)),np.float32)
    for rr,x,_ in log_count_blocks(path,selected=selected): mat[rr]=x.toarray()
    return ad.AnnData(X=mat,obs=obs,var=pd.DataFrame({'gene_symbol':selected},index=selected),layers={'lognorm':mat})

def donor_mean(x, donors, mask):
    ds=np.unique(donors[mask])
    return np.stack([x[mask & (donors==d)].mean(0) for d in ds]).mean(0).astype(np.float32)

def prepare_features():
    ROOT.mkdir(parents=True,exist_ok=True)
    split=json.loads(Path(REPORT['inputs']['split']['path']).read_text())['donors']
    sets={}; detail={}
    for key,cohort,n in [('brain','c178',2000),('blood223','b223',1500),('bloodmb','bmb',1500)]:
        path=REPORT['inputs'][key]['path'];obs,var=metadata(path)
        dn=obs.donor.astype(str).to_numpy()
        if key=='brain':
            keep=np.isin(dn,[d for d,v in split.items() if v=='train' and d not in BRAIN_VAL])
            keep &= obs.cell_type.astype(str).isin(['Astrocyte','Microglia']).to_numpy()
        else:
            dn=np.array([cohort+':'+d for d in dn])
            keep=obs.condition.astype(str).isin(COND).to_numpy() & ~np.isin(dn,BLOOD_VAL+BLOOD_TEST)
            keep &= obs.cell_type.astype(str).isin(['Monocyte_CD14','Monocyte_CD16','DC','pDC']).to_numpy()
        rng=np.random.default_rng(9102); rows=[]
        # Equal maximum cells per training donor; no validation/test expression used.
        for d in np.unique(dn[keep]):
            ii=np.where(keep&(dn==d))[0];rows.extend(rng.choice(ii,min(2000,len(ii)),replace=False))
        rows=np.sort(rows); s1=s2=None;cnt=0
        for rr,x,gs in log_count_blocks(path,rows=rows):
            a=np.asarray(x.sum(0)).ravel();b=np.asarray(x.power(2).sum(0)).ravel()
            s1=a if s1 is None else s1+a;s2=b if s2 is None else s2+b;cnt+=len(rr)
        if not cnt: raise ValueError(('empty feature selection',key))
        mu=s1/cnt;var=s2/cnt-mu**2;disp=var/np.maximum(mu,1e-8)
        # Standardize dispersion within expression bins to avoid pure abundance selection.
        bins=pd.qcut(mu,20,labels=False,duplicates='drop'); score=np.zeros(len(mu))
        for b in np.unique(bins):
            m=bins==b;v=disp[m];score[m]=(v-np.mean(v))/(np.std(v)+1e-8)
        score[mu<=.01]=-np.inf
        sets[key]=[gs[i] for i in np.argsort(-score)[:n]]
        detail[key]={'n_train_cells':cnt,'donors':sorted(set(dn[rows])),'n_hvg':n}
        print('HVG',key,detail[key],flush=True)
    (ROOT/'tissue_features.json').write_text(json.dumps({'genes':sorted(set.union(*(set(v) for v in sets.values()))),'by_source':sets,'selection':detail},indent=2))

if __name__=='__main__': prepare_features()
