"""Donor-level age/sex and lane-composition sensitivity for the measured MB cohort."""
from pathlib import Path
import json,time,h5py
import numpy as np,pandas as pd
from anndata.io import read_elem
from scipy.stats import t
from randko import bh

A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/optimization_20261003'

def fit(y,condition,cov):
    keep=cov.std(0)>1e-12
    cov=cov[:,keep];cov=(cov-cov.mean(0))/(cov.std(0)+1e-12)
    X=np.column_stack([np.ones(len(y)),condition,cov]);rank=np.linalg.matrix_rank(X)
    basic=dict(n_donors=len(y),n_PD=int(condition.sum()),n_control=int((1-condition).sum()),n_parameters=X.shape[1])
    if rank<X.shape[1] or len(y)<=rank+2 or min(basic['n_PD'],basic['n_control'])<2:
        return dict(**basic,status='nonidentified_or_insufficient_donors')
    bread=np.linalg.pinv(X.T@X);beta=bread@X.T@y;res=y-X@beta;hat=np.sum((X@bread)*X,1)
    if np.any(hat>=1-1e-8):return dict(**basic,status='unit_leverage_no_HC3_estimate')
    covariance=bread@(X.T@((res/(1-hat))[:,None]**2*X))@bread
    se=float(np.sqrt(max(covariance[1,1],0)))
    return dict(**basic,status='ok',effect_PD_minus_control=float(beta[1]),se_HC3=se,
                p=float(2*t.sf(abs(beta[1])/max(se,1e-12),len(y)-rank)))

def main():
    while not (O/'blood_disease_donor_values.csv').exists():
        if json.loads((O/'status.json').read_text()).get('stage')=='failed':raise RuntimeError('Main run failed')
        time.sleep(15)
    with h5py.File(A/'processed/moquin_beaudry2025/v0.1/MB2025_pbmc_annotated.h5ad') as f:obs=read_elem(f['obs'])
    # Use all retained fine-state cells to form lane proportions, including sparse pDC donors.
    pure=pd.read_csv(A/'interim/v0.1/purification/bmb_blood_myeloid.csv',index_col=0);pure=pure[pure.kept]
    typed=pd.Series('other',index=obs.index);ix=typed.index.intersection(pure.index);typed.loc[ix]=pure.loc[ix,'state_pure'].astype(str)
    lanes=obs.index.str.split('_').str[0];lane_names=sorted(set(lanes));covariates={}
    for (donor,state),indices in obs.assign(state_pure=typed).groupby(['donor','state_pure'],observed=True).groups.items():
        if state=='other':continue
        rows=obs.loc[indices];row=rows.iloc[0]
        assert rows.Age.nunique()==1 and rows.Sex.nunique()==1
        age=float(row.Age);sex=1. if row.Sex=='M' else 0. if row.Sex=='F' else np.nan
        mask=obs.index.isin(indices);fraction=np.array([(lanes[mask]==lane).mean() for lane in lane_names])
        assert np.isfinite([age,sex]).all() and np.isclose(fraction.sum(),1)
        covariates['bmb:'+str(donor),str(state)]=np.r_[age,sex,fraction[:-1]]
    donor=pd.read_csv(O/'blood_disease_donor_values.csv');donor=donor[donor.cohort==1];results=[]
    for (assay,target,state),df in donor.groupby(['assay','target','blood_state']):
        C=np.stack([covariates[d,state] for d in df.donor])
        for model,width in [('age_sex',2),('age_sex_batch',C.shape[1])]:
            row=fit(df.value.to_numpy(),df.condition.to_numpy(),C[:,:width])
            results.append(dict(assay=assay,target=target,blood_state=state,covariate_model=model,**row))
    df=pd.DataFrame(results);df['q_BH']=np.nan
    for _,ix in df[df.status=='ok'].groupby(['assay','covariate_model']).groups.items():df.loc[ix,'q_BH']=bh(df.loc[ix,'p'])
    df.to_csv(O/'blood_disease_covariate_sensitivity.csv',index=False)
    (O/'blood_covariate_definition.json').write_text(json.dumps(dict(cohort='Moquin-Beaudry 2025',
        models=['age and sex','age, sex and five independent multiplex-lane cell fractions within donor/state'],
        estimator='OLS disease coefficient, HC3 covariance, donor-level t reference; rank and leverage checked',
        multiple_testing='BH across all candidates and states separately within each assay and covariate model',
        status_counts=df.status.value_counts().to_dict(),scope='prespecified sensitivity; no change to frozen nomination criteria; GSE223138 lacks source age/sex metadata'),indent=2))
    print('BLOOD COVARIATE SENSITIVITY COMPLETE',len(df),flush=True)

if __name__=='__main__':main()
