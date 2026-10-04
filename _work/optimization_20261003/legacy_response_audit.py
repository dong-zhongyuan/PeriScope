"""Keep legacy gene membership fixed to separate model changes from set changes.

These previously selected hypotheses are retrospective diagnostics, not a new
nomination family. The legacy positive-response direction is retained.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

A=Path('/public/home/mengxl/dzy/pd_product_assets')
O=A/'results/optimization_20261003'
legacy=json.loads((Path(__file__).parent/'reference/legacy_programs.json').read_text())
rows=[];perseed=[];draws={}
for file in sorted((A/'results/rescreen_20261002/curves').glob('*__slopes.npz')):
    axis=file.name.replace('__slopes.npz','');z=np.load(file)
    genes=z['genes'].astype(str).tolist();ps=z['proteins'].astype(str).tolist()
    conf=np.where(z['seeds']>=47)[0];D=z['D'];bg=np.where(~z['is_control'])[0]
    mu=D[:,bg].mean(1,keepdims=True);sd=D[:,bg].std(1,keepdims=True)
    Z=(D-mu)/(sd+1e-9)
    for combo,members in legacy.items():
        protein,caxis=combo.split('__',1)
        if caxis!=axis:continue
        present=[g for g in members if g in genes]
        if not present:continue
        ix=np.array([genes.index(g) for g in present]);pi=ps.index(protein)
        values=D[:,pi][:,ix].mean(1);response=D[conf,pi].mean(0)
        key=(len(genes),len(ix))
        if key not in draws:
            rng=np.random.RandomState(7)
            draws[key]=np.stack([rng.choice(len(genes),len(ix),False) for _ in range(5000)])
        null=response[draws[key]].mean(1);obs=response[ix].mean()
        score=Z[:,:,ix].mean(2)-Z.mean(2)
        panel=(score-score[:,bg].mean(1,keepdims=True))/(score[:,bg].std(1,keepdims=True)+1e-12)
        mean=panel[conf].mean(0);other=bg[bg!=pi]
        rows.append(dict(combo=combo,n_legacy_genes=len(members),n_present=len(ix),
            missing_genes=';'.join(g for g in members if g not in genes),
            confirmation_positive_fraction=float((values[conf]>0).mean()),
            confirmation_mean_slope=float(values[conf].mean()),
            p_gene_decoy=float((1+(null>=obs).sum())/(1+len(null))),
            protein_percentile=float((mean[other]<mean[pi]).mean()),
            hypothesis_direction='legacy positive-response direction',status='retrospective diagnostic only'))
        perseed.extend(dict(combo=combo,seed=int(s),mean_slope=float(v)) for s,v in zip(z['seeds'],values))
df=pd.DataFrame(rows);df['q_within_legacy_diagnostic_family']=false_discovery_control(df.p_gene_decoy)
df.to_csv(O/'legacy_membership_current_model_audit.csv',index=False)
pd.DataFrame(perseed).to_csv(O/'legacy_membership_current_seed_scores.csv',index=False)
print(df[df.combo.str.startswith(('CD22__cDC__x__astro','CD35__nonclassical_mono__x__astro'))].to_string(index=False))
