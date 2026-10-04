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
pd.DataFrame(support).to_csv(S/'donor_support.csv',index=False)
protocol=dict(scope='Compartment-matched sensitivity only; no nomination gates changed',
    cohorts=['locked Kamath donors','independent GSE157783 donors'],unit='donor mean of cell log1p(CPM), within existing fine state labels',
    minimum_cells_per_donor=25,minimum_donors_per_arm=2,
    statistic='Same mean-gene-Welch-t competitive test and 5000 expression-matched random gene sets as the main analysis',
    multiplicity='BH across all programs and both directions; missing-state tests enter as P=1 and remain labeled unavailable',
    meta='Sample-size-weighted Stouffer, requiring both cell-state cohorts to be evaluable',
    unchanged='All Ma spatial and original brain-gate results retained. This comparison does not relabel A/B candidates.')
(S/'protocol.json').write_text(json.dumps(protocol,indent=2))
for label,O in [('raw_discovery',P),('target_specific_discovery',P/'target_specific')]:
    while not (O/'discovery_evidence_completed.json').exists():
        if json.loads((P/'status.json').read_text()).get('stage')=='failed':raise RuntimeError('Main run failed')
        time.sleep(20)
    programs=json.loads((O/'programs.json').read_text());brain=json.loads((O/'brain_spatial_competitive.json').read_text());rows=[]
    for state in ['astro','microglia_mhc2','microglia_homeostatic']:
        for combo,geneset in programs.items():
            if combo.split('__x__')[1].split('__')[0]!=state:continue
            kam=brain['kamath_sets']['gsea::'+combo]
            if state in profiles:
                M,is_pd,donors=profiles[state]
                c157=_competitive_test(M,is_pd,[gi[g] for g in geneset if g in gi])
                c157.update(n_donors=len(donors),n_PD=int(is_pd.sum()),n_Control=int((~is_pd).sum()))
            else:c157=dict(status='insufficient_donors')
            available=kam.get('status')=='ok' and c157.get('status')=='ok'
            record=dict(combo=combo,brain_state=state,kamath_status=kam.get('status'),gse157783_status=c157.get('status'),both_available=available)
            for prefix,result in [('kamath',kam),('gse157783',c157)]:
                record[prefix+'_effect']=result.get('donor_score_effect_PD_minus_control',np.nan)
                record[prefix+'_relative_effect']=result.get('stat',np.nan)-result.get('null_mean',np.nan)
                for direction in ['up','down']:record[prefix+'_p_'+direction]=result.get('p_'+direction,1.)
            record['actual_directions_agree']=bool(record['kamath_effect']*record['gse157783_effect']>0)
            record['relative_directions_agree']=bool(record['kamath_relative_effect']*record['gse157783_relative_effect']>0)
            for direction in ['up','down']:
                if available:
                    w=np.sqrt([kam['n_donors'],c157['n_donors']]);pp=np.clip([kam['p_'+direction],c157['p_'+direction]],1e-300,1-1e-16)
                    record['meta_p_'+direction]=float(stats.norm.sf(np.dot(stats.norm.isf(pp),w)/np.linalg.norm(w)))
                else:record['meta_p_'+direction]=1.
            rows.append(record)
    table=pd.DataFrame(rows)
    for source_name in ['kamath','gse157783','meta']:
        q=false_discovery_control(table[[source_name+'_p_up',source_name+'_p_down']].to_numpy().ravel()).reshape(-1,2)
        table[source_name+'_q_up']=q[:,0];table[source_name+'_q_down']=q[:,1]
    direction=np.where(table.meta_q_up<=table.meta_q_down,1.,-1.)
    table['matched_state_support']=table.both_available&table.actual_directions_agree&(direction*table.kamath_relative_effect>0)&(direction*table.gse157783_relative_effect>0)&((table[['meta_q_up','meta_q_down']].min(axis=1))<=.05)
    table.to_csv(S/(label+'.csv'),index=False)
    print('CELL-STATE SENSITIVITY',label,len(table),int(table.matched_state_support.sum()),flush=True)
(S/'completed.json').write_text(json.dumps(dict(status='complete',nomination_changes=False),indent=2))
