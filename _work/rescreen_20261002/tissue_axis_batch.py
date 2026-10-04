"""Retain the original sequencing-batch sensitivity with donor/state lane proportions."""
from pathlib import Path
import json,time,h5py
import numpy as np,pandas as pd
from anndata.io import read_elem
from concurrent.futures import ProcessPoolExecutor
from tissue_axis import O,A,D,compare

def main():
    while not (O/'tissue_axis_associations.json').exists():time.sleep(10)
    with h5py.File(A/'processed/moquin_beaudry2025/v0.1/MB2025_pbmc_annotated.h5ad') as f:obs=read_elem(f['obs'])
    pure=pd.read_csv(A/'interim/v0.1/purification/bmb_blood_myeloid.csv',index_col=0);pure=pure[pure.kept]
    typed=pd.Series('other',index=obs.index);ix=typed.index.intersection(pure.index);typed.loc[ix]=pure.loc[ix,'state_pure'].astype(str)
    lanes=obs.index.str.split('_').str[0];lane_names=sorted(set(lanes));blocks={};batch_rows=[]
    for file in D.glob('*.npz'):
        cohort,state=file.stem.split('__',1);z=np.load(file);C=z['covariates'];cov=C if C.shape[1] else None
        if cohort=='bmb':
            extra=[]
            for donor in z['donors']:
                mask=(obs.donor.astype(str)==donor)&(typed==state);fractions=np.array([(lanes[mask]==lane).mean() for lane in lane_names])
                assert np.isclose(fractions.sum(),1)
                extra.append(fractions[:-1])
                for lane,fraction in zip(lane_names,fractions):batch_rows.append(dict(donor=donor,blood_state=state,lane=lane,cell_fraction=float(fraction)))
            extra=np.array(extra);extra=extra[:,extra.std(0)>1e-12];extra=(extra-extra.mean(0))/(extra.std(0)+1e-12)
            cov=np.column_stack([C,extra])
        blocks.setdefault(cohort,{})[state]=dict(means=z['means'],conds=z['conditions'],donors=z['donors'].astype(str).tolist(),cov=cov)
    pd.DataFrame(batch_rows).to_csv(O/'blood_lane_fractions_by_donor.csv',index=False)
    tasks=[]
    for blood,a in blocks['bmb'].items():
        for brain in sorted(set(blocks.get('c178',{}))|set(blocks.get('c157',{}))):
            bs=[blocks[c][brain] for c in ['c178','c157'] if brain in blocks.get(c,{})]
            tasks.append((blood+'__x__brain_meta_'+brain,a,bs,True))
        for c in ['ma','pdd']:
            for name,b in blocks.get(c,{}).items():tasks.append((blood+'__x__'+c+'_'+name,a,[b],True))
    with ProcessPoolExecutor(4) as p:new=list(p.map(compare,tasks))
    source=O/'tissue_axis_associations.json';data=json.loads(source.read_text());old=[r for r in data['results'] if r.get('covariate_model')!='age_sex_batch']
    for row in old:row['covariate_model']='age_sex' if row['covariate_adjusted'] else 'unadjusted'
    for row in new:row['covariate_model']='age_sex_batch'
    data['results']=old+new;data['batch_sensitivity']='MB2025 age, sex and five independent multiplex-lane fractions within donor/state; Kamath age, sex and PMI; adjusted bootstrap requires an identifiable design'
    data['primary_covariate_model']='age_sex_batch; age_sex and unadjusted retained as sensitivities, without choosing the model by correlation strength'
    source.write_text(json.dumps(data,indent=2));pd.DataFrame(data['results']).to_csv(O/'tissue_axis_associations.csv',index=False)
    print('TISSUE BATCH SENSITIVITY COMPLETE',len(new),flush=True)
if __name__=='__main__':main()
