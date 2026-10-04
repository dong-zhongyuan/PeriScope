"""Post-screen identity diagnostics; does not alter frozen nominations or thresholds."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata

O=Path('/public/home/mengxl/dzy/pd_product_assets/results/rescreen_20261002')
I=O.parent.parent/'interim/rescreen_20261002'

def rho(a,b):
    return float(spearmanr(a,b).statistic) if np.std(a)>1e-10 and np.std(b)>1e-10 else np.nan

def partial(a,b,c):
    x=np.column_stack([np.ones(len(c)),rankdata(c)])
    ar=rankdata(a);br=rankdata(b)
    ar=ar-x@np.linalg.lstsq(x,ar,rcond=None)[0]
    br=br-x@np.linalg.lstsq(x,br,rcond=None)[0]
    return float(np.corrcoef(ar,br)[0,1]) if min(np.std(ar),np.std(br))>1e-10 else np.nan

def main():
    registry=json.loads((O/'candidate_registry.json').read_text());by={r['target']:r for r in registry}
    z=np.load(I/'citeseq.npz');X=z['X'];Y=z['Y'];genes=z['genes'].astype(str).tolist()
    states=z['state_vocab'].astype(str).tolist();state=z['state'];donor=z['donor'];subtype=z['subtype']
    markers=[g for g in ['PF4','PPBP','GP9','GP1BB','ITGA2B'] if g in genes]
    peers=[t for t in ['CD41','CD61','CD42a'] if t in by]
    raw_gp=X[:,genes.index('GP1BA')]+z['gene_mu'][genes.index('GP1BA')]
    score=(X[:,[genes.index(g) for g in markers]]+z['gene_mu'][[genes.index(g) for g in markers]]).mean(1)
    y=Y[:,by['CD42b']['channels']].mean(1)
    controls=np.median(np.column_stack([Y[:,r['channels']].mean(1) for r in registry if r['is_control']]),axis=1)
    rows=[]
    for st in ['classical_mono','nonclassical_mono']:
        for d in np.unique(donor[state==states.index(st)]):
            for sub in np.unique(subtype[(state==states.index(st))&(donor==d)]):
                ix=np.where((state==states.index(st))&(donor==d)&(subtype==sub))[0]
                if len(ix)<25:continue
                base=dict(target='CD42b',blood_state=st,citeseq_subtype=str(sub),donor=str(z['donor_names'][d]),n_cells=len(ix),
                    split='test' if z['is_test'][ix[0]] else 'validation' if z['is_val'][ix[0]] else 'train',
                    GP1BA_RNA_detection=float((raw_gp[ix]>1e-6).mean()),GP1BA_RNA_mean=float(raw_gp[ix].mean()))
                for feature,v in [('GP1BA_RNA',raw_gp),('platelet_RNA_score_excluding_GP1BA',score),('isotype_median',controls)]+[(t+'_ADT',Y[:,by[t]['channels']].mean(1)) for t in peers]:
                    rows.append(base|dict(feature=feature,rho=rho(y[ix],v[ix]),rho_adjusted_for_isotype=np.nan if feature=='isotype_median' else partial(y[ix],v[ix],controls[ix])))
    pd.DataFrame(rows).to_csv(O/'priority_identity_qc.csv',index=False)
    definition=dict(status='complete',timing='post-screen interpretation only; no change to frozen candidate ranking',
        candidate='CD42b / GP1BA',RNA_context_genes=markers,ADT_context_antigens=peers,
        method='within donor and fine subtype Spearman correlations; partial rank correlation after residualizing both ranks on median isotype rank; descriptive, no p values',
        interpretation='CD42b on monocyte-labelled events may reflect platelet association; these correlations do not distinguish attachment, protein transfer, or ambient signal',
        sources=[dict(url='https://www.ncbi.nlm.nih.gov/gene/2811',supports='GP1BA encodes the platelet VWF-receptor alpha subunit'),
                 dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC11430373/',supports='Primary CITE-seq study used imaging flow cytometry to identify CD42b-positive monocyte-platelet aggregates; this is contextual evidence, not validation in the current PD cohort')])
    (O/'priority_identity_qc_definition.json').write_text(json.dumps(definition,indent=2))
    print('PRIORITY IDENTITY QC',len(rows),flush=True)

if __name__=='__main__':main()
