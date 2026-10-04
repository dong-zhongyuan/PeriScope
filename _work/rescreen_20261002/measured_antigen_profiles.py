"""Observed antigen and mapped-gene profiles, retaining donor and fine subtype."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
O=Path('/public/home/mengxl/dzy/pd_product_assets/results/rescreen_20261002');I=O.parent.parent/'interim/rescreen_20261002'
z=np.load(I/'citeseq.npz');X=z['X'];Y=z['Y'];idx={g:i for i,g in enumerate(z['genes'].astype(str))};states=z['state_vocab'].astype(str).tolist();registry=json.loads((O/'candidate_registry.json').read_text());rows=[]
for st in ['cDC','classical_mono','nonclassical_mono','pDC']:
    for d in np.unique(z['donor'][z['state']==states.index(st)]):
        for subtype in np.unique(z['subtype'][(z['state']==states.index(st))&(z['donor']==d)]):
            ix=np.where((z['donor']==d)&(z['subtype']==subtype))[0]
            if len(ix)<25:continue
            raw=X[ix]+z['gene_mu'];detected=(raw>1e-6).mean(0);mean_rna=raw.mean(0)
            split='test' if z['is_test'][ix[0]] else 'validation' if z['is_val'][ix[0]] else 'train'
            for r in registry:
                y=Y[ix][:,r['channels']].mean(1);gi=[idx[g] for g in r['genes'] if g in idx]
                rows.append(dict(target=r['target'],blood_state=st,citeseq_subtype=str(subtype),donor=z['donor_names'][d],split=split,n_cells=len(ix),
                    mean_measured_ADT_CLR=float(y.mean()),median_measured_ADT_CLR=float(np.median(y)),sd_measured_ADT_CLR=float(y.std()),
                    q05_measured_ADT_CLR=float(np.quantile(y,.05)),q95_measured_ADT_CLR=float(np.quantile(y,.95)),
                    genes=';'.join(r['genes']),n_model_genes_present=len(gi),RNA_fraction_nonzero_mean_components=float(detected[gi].mean()) if gi else np.nan,
                    RNA_mean_log1p_CP10k_components=float(mean_rna[gi].mean()) if gi else np.nan,is_control=r['is_control']))
pd.DataFrame(rows).to_csv(O/'measured_antigen_profiles.csv',index=False)
(O/'measured_antigen_profile_definition.json').write_text(json.dumps(dict(source=str(I/'citeseq.npz'),scope='observed RNA and ADT; donors and fine subtypes retained',
    ADT='zero-sum CLR across the measured panel; across-antibody CLR magnitude is not an absolute protein-concentration calibration',
    RNA='mean component-gene values for complexes; missing genes retained as missing, not set to zero',n_profiles=len(rows)),indent=2))
print('MEASURED ANTIGEN PROFILES',len(rows),flush=True)
