"""Full-candidate donor-level blood disease association, cohort-adjusted HC3."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from randko import bh

A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/rescreen_20261002';I=A/'interim/rescreen_20261002'

def test(y,pdlabel,cohort):
    columns=[np.ones(len(y)),np.asarray(pdlabel)]
    if len(set(cohort))>1:columns.append(np.asarray(cohort))
    X=np.stack(columns,axis=1);rank=np.linalg.matrix_rank(X)
    if rank<X.shape[1] or min(sum(pdlabel==0),sum(pdlabel==1))<2 or len(y)<=rank+1:return {}
    bread=np.linalg.pinv(X.T@X);beta=bread@X.T@y;res=y-X@beta;h=np.sum((X@bread)*X,axis=1)
    meat=X.T@((res/np.maximum(1-h,1e-8))[:,None]**2*X);cov=bread@meat@bread
    se=np.sqrt(max(cov[1,1],0));stat=beta[1]/max(se,1e-12)
    return dict(effect_PD_minus_control_adjusted=float(beta[1]),se_HC3=float(se),p=float(2*student_t.sf(abs(stat),len(y)-rank)),n_donors=len(y),n_PD=int(sum(pdlabel==1)),n_control=int(sum(pdlabel==0)))

def main():
    records=json.loads((O/'candidate_registry.json').read_text());z=np.load(I/'blood.npz');X=z['X'];genes=z['genes'].astype(str).tolist();gidx={g:i for i,g in enumerate(genes)}
    states=z['state_vocab'].astype(str).tolist();donor_rows=[];tested=[];gene_rows=[]
    candidate_genes=sorted({g for r in records for g in r['genes'] if g in gidx});candidate_ix=[gidx[g] for g in candidate_genes]
    for state in ['cDC','classical_mono','nonclassical_mono','pDC']:
        for d in np.unique(z['donor'][z['state']==states.index(state)]):
            ix=np.where((z['state']==states.index(state))&(z['donor']==d))[0]
            if len(ix)<25:continue
            mean=X[ix].mean(0)+z['gene_mu']
            detected=(X[ix][:,candidate_ix]+z['gene_mu'][candidate_ix]>1e-6).mean(0)
            for j,g in enumerate(candidate_genes):gene_rows.append(dict(gene=g,blood_state=state,donor=z['donor_names'][d],condition=int(z['cond2'][ix[0]]),cohort=int(z['cohort'][ix[0]]),n_cells=len(ix),mean_log1p_CP10k=float(mean[gidx[g]]),fraction_nonzero=float(detected[j])))
            for r in records:
                ii=[gidx[g] for g in r['genes'] if g in gidx]
                if not ii or r['is_control']:continue
                donor_rows.append(dict(target=r['target'],blood_state=state,donor=z['donor_names'][d],condition=int(z['cond2'][ix[0]]),cohort=int(z['cohort'][ix[0]]),value=float(mean[ii].mean()),assay='observed_RNA',n_cells=len(ix),n_genes_present=len(ii),genes=';'.join(r['genes'])))
    pd.DataFrame(gene_rows).to_csv(O/'target_gene_expression_by_donor.csv',index=False)
    predicted=pd.concat([pd.read_csv(p) for p in sorted((O/'evaluation').glob('blood_donor_predictions_seed*.csv'))],ignore_index=True)
    predicted=predicted.groupby(['target','blood_state','donor','condition','cohort'],as_index=False).agg(value=('predicted_ADT','mean'),n_cells=('n','first'))
    predicted['assay']='inferred_ADT';valid={r['target'] for r in records if not r['is_control']};predicted=predicted[predicted.target.isin(valid)]
    donor=pd.concat([pd.DataFrame(donor_rows),predicted],ignore_index=True)
    pending=O/'blood_disease_donor_values.csv.pending';donor.to_csv(pending,index=False);pending.replace(O/'blood_disease_donor_values.csv')
    for (assay,target,state),df in donor.groupby(['assay','target','blood_state']):
        for cohort in ['combined',0,1]:
            d=df if cohort=='combined' else df[df.cohort==cohort]
            r=test(d.value.to_numpy(),d.condition.to_numpy(),d.cohort.to_numpy())
            if r:tested.append(dict(assay=assay,target=target,blood_state=state,cohort=cohort,**r))
    df=pd.DataFrame(tested);df['q_BH']=np.nan
    for keys,ix in df.groupby(['assay','cohort']).groups.items():df.loc[ix,'q_BH']=bh(df.loc[ix,'p'])
    df.to_csv(O/'blood_disease_associations.csv',index=False)
    (O/'disease_association_definition.json').write_text(json.dumps(dict(unit='donor',observed='donor-average log1p(CP10k) RNA',
        inferred='mean of ten-model ADT predictions, aggregated by donor and blood subtype',
        model='OLS disease coefficient adjusted for blood cohort, HC3 covariance; cohort-specific estimates also reported',
        family='all candidate x blood-state associations within each assay and cohort analysis',
        n_tests=len(df),scope='within-project disease association; disease-cohort surface ADT was not directly measured'),indent=2))

if __name__=='__main__':main()
