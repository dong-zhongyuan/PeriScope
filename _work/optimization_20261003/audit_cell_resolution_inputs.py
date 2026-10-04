"""Read actual cellular units and donor coverage; do not modify model inputs."""
from pathlib import Path
import hashlib, json
import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

A=Path('/public/home/mengxl/dzy/pd_product_assets')
O=A/'results/optimization_20261003/cell_resolution_audit_20261003'
O.mkdir(exist_ok=True)
I=A/'interim/rescreen_20261002'
F=A/'interim/v0.1/purification'
files={'c178':'gse178265/v0.1/GSE178265_sn_annotated.h5ad',
       'c157':'gse157783/v0.1/gse157783_qc.h5ad',
       'b223':'gse223138/v0.1/GSE223138_pbmc_annotated.h5ad',
       'bmb':'moquin_beaudry2025/v0.1/MB2025_pbmc_annotated.h5ad',
       'ma':'gse253975/v0.1/GSE253975_geomx.h5ad'}
split=json.loads((A/'interim/v0.1/splits/GSE178265_sn_donor_split_v1.json').read_text())['donors']
counts=[];clusters=[];sources={};metadata={}
for cohort,rel in files.items():
    path=A/'processed'/rel
    with h5py.File(path) as f:
        obs=read_elem(f['obs'])
        metadata[cohort]={'path':str(path),'n_observations':len(obs),'obs_columns':list(obs.columns),
                          'obsm_keys':list(f.get('obsm',{})),'uns_keys':list(f.get('uns',{}))}
    if cohort=='ma':
        obs.groupby(['donor','condition'],observed=True).size().rename('n_spots').reset_index().to_csv(O/'spatial_donor_units.csv',index=False)
        continue
    families=['blood_myeloid'] if cohort.startswith('b') else ['brain_astro','brain_microglia']
    for fam in families:
        path=F/f'{cohort}_{fam}.csv';m=pd.read_csv(path,index_col=0)
        sources[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        assert m.index.is_unique
        q=m.join(obs[['donor','condition']+(['tissue'] if 'tissue' in obs else [])],how='left',validate='one_to_one')
        assert q.donor.notna().all()
        q['split']=q.donor.astype(str).map(split).fillna('not_model_brain_split') if cohort=='c178' else 'not_model_brain_split'
        for (cl,st,kept),g in q.groupby(['cl','state_pure','kept'],observed=True):
            clusters.append(dict(cohort=cohort,family=fam,cluster=cl,state=st,kept=bool(kept),n_cells=len(g),n_donors=g.donor.nunique()))
        q=q[q.kept]
        if cohort=='c178':q=q[q.tissue.astype(str).eq('substantia_nigra')]
        for (st,d,cond,sp),g in q.groupby(['state_pure','donor','condition','split'],observed=True):
            counts.append(dict(cohort=cohort,state=st,donor=d,condition=cond,split=sp,n_cells=len(g),meets_25_cells=len(g)>=25))
pd.DataFrame(counts).to_csv(O/'cell_state_donor_coverage.csv',index=False)
pd.DataFrame(clusters).to_csv(O/'purification_cluster_assignments.csv',index=False)
z=np.load(I/'citeseq.npz')
v=z['state_vocab'].astype(str)
cite=pd.DataFrame({'donor':z['donor_names'][z['donor']],'subtype':z['subtype'],
    'model_state':v[z['state']], 'split':np.where(z['is_test'],'test',np.where(z['is_val'],'validation','train'))})
cite.groupby(['donor','subtype','model_state','split']).size().rename('n_cells').reset_index().to_csv(O/'cite_subtype_to_model_state.csv',index=False)
for name in ['brain_train','brain_locked_test','brain_calibration','blood']:
    z=np.load(I/(name+'.npz'));v=z['state_vocab'].astype(str)
    d=pd.DataFrame({'donor':z['donor_names'][z['donor']],'state':v[z['state']],'condition':z['cond2']})
    if 'is_val' in z:d['is_val']=z['is_val']
    if 'is_test' in z:d['is_test']=z['is_test']
    d.groupby(list(d.columns)).size().rename('n_cells').reset_index().to_csv(O/(name+'_units.csv'),index=False)
(O/'input_metadata.json').write_text(json.dumps(metadata,indent=2))
(O/'purification_map_hashes.json').write_text(json.dumps(sources,indent=2))
print(json.dumps({'status':'complete','output':str(O),'cohorts':{k:v['n_observations'] for k,v in metadata.items()}},indent=2))
