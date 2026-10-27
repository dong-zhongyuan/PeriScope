"""Compartment-matched sensitivity: Kamath and GSE157783 fine brain states.

Spatial mixed-tissue and single-cell directions need not be identical because
cell composition can differ. This analysis reports matched-state evidence and
does not change the frozen nomination criteria or remove spatial results.
"""
from pathlib import Path
import json,time,hashlib,h5py
import numpy as np,pandas as pd
from anndata.io import read_elem,sparse_dataset
from scipy import sparse,stats
from input_utils import canonical_symbols
from external_gene_symbols import external_symbols
from brain_competitive import _competitive_test
from scipy.stats import false_discovery_control

A=Path('/public/home/mengxl/dzy/pd_product_assets');P=A/'results/optimization_20261003';S=P/'cell_state_sensitivity';S.mkdir(exist_ok=True)
source=A/'processed/gse157783/v0.1/gse157783_qc.h5ad'
with h5py.File(source,'r') as f:
    obs=read_elem(f['obs']);var=read_elem(f['var'])
    symbols=external_symbols(var)
    genes=sorted(set(symbols));gi={g:i for i,g in enumerate(genes)}
    collapse=sparse.csr_matrix((np.ones(len(symbols)),(np.arange(len(symbols)),[gi[g] for g in symbols])),shape=(len(symbols),len(genes)))
    state=pd.Series('other',index=obs.index)
    for family in ['brain_astro','brain_microglia']:
        mapping=pd.read_csv(A/f'interim/v0.1/purification/c157_{family}.csv',index_col=0);mapping=mapping[mapping.kept]
        common=state.index.intersection(mapping.index);state.loc[common]=mapping.loc[common,'state_pure'].astype(str)
    selected=np.where(state!='other')[0];raw=sparse_dataset(f['X'])[selected].astype(float)
raw=(raw@collapse).multiply((1e6/np.asarray(raw.sum(1)).ravel().clip(1))[:,None]).tocsr();raw.data=np.log1p(raw.data)
meta=obs.iloc[selected].copy();meta['state']=state.iloc[selected].to_numpy();profiles={};support=[]
for name in ['astro','microglia_mhc2','microglia_homeostatic']:
    cells=meta.state.to_numpy()==name;donors=[];means=[];conditions=[]
    for d in sorted(set(meta.loc[cells,'donor'].astype(str))):
        mask=cells&(meta.donor.astype(str).to_numpy()==d)
        if mask.sum()<25:continue
        donors.append(d);means.append(np.asarray(raw[mask].mean(0)).ravel());conditions.append(meta.loc[mask,'condition'].iloc[0]=='PD')
        support.append(dict(cohort='GSE157783',state=name,donor=d,condition='PD' if conditions[-1] else 'Control',n_cells=int(mask.sum())))
    if len(means) and min(sum(conditions),len(conditions)-sum(conditions))>=2:
        profiles[name]=(np.stack(means),np.array(conditions),donors)

O=P/'figure_inputs_20261004';O.mkdir(exist_ok=True)
pd.DataFrame(support).to_csv(O/'matched_donor_support.csv',index=False)
programs=json.loads((P/'joint_analysis/programs.json').read_text());e=pd.read_csv(P/'mr_parallel_gate_20261003/all_programs_parallel_gate.csv');e=e[e.priority_pass | e.C7.eq('pass')];rows=[];checks=[]
from brain_competitive import _welch_t
for row in e.itertuples():
    if row.brain_state not in profiles:continue
    M,is_pd,donors=profiles[row.brain_state];ok=np.isfinite(_welch_t(M,is_pd));cols=np.array([gi[g] for g in programs[row.combo] if g in gi and ok[gi[g]]])
    if len(cols)<4:continue
    a=M[:,cols];a=(a-a.mean(0))/(a.std(0)+1e-12);score=a.mean(1)
    eff=float(score[is_pd].mean()-score[~is_pd].mean());expected=row.cell_state_sensitivity_gse157783_effect
    assert np.isclose(eff,expected,rtol=1e-8,atol=1e-10),(row.combo,eff,expected)
    for i,d in enumerate(donors):rows.append(dict(combo=row.combo,brain_state=row.brain_state,donor=d,condition='PD' if is_pd[i] else 'Control',score=score[i],n_genes=len(cols)))
    checks.append(dict(combo=row.combo,calculated_effect=eff,recorded_effect=expected))
pd.DataFrame(rows).to_csv(O/'matched_independent_donor_scores.csv',index=False);pd.DataFrame(checks).to_csv(O/'matched_score_checks.csv',index=False)
print('MATCHED DONORS',len(rows),'checks',len(checks),flush=True)
