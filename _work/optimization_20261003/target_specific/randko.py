"""Full-family randKO on regenerated measured-dose responses. Discovery seeds 42-46 define programs; confirmation seeds 47-51 test gene selectivity and protein-panel calibration."""
import os
import sys
import json
import glob
import time
import argparse

import numpy as np
from scipy.stats import combine_pvalues, wilcoxon

A = '/public/home/mengxl/dzy/pd_product_assets'
CACHE = A + '/results/optimization_20261003/target_specific/curves'
QR = np.arange(5, dtype=float)
QUANTILES = [0.05, 0.25, 0.50, 0.75, 0.95]
DRAWS = [42, 43, 44]


def slope(C):
    x = QR - QR.mean()
    return np.tensordot(x, C, axes=(0, 2)) / (x * x).sum()


def bh(ps):
    ps = np.asarray(ps, float)
    m = len(ps)
    o = np.argsort(ps)
    q = np.empty(m)
    prev = 1.0
    for r_, i_ in zip(range(m, 0, -1), o[::-1]):
        prev = min(prev, ps[i_] * m / r_)
        q[i_] = prev
    return q


def emp_p(vr, vn):
    return float((1 + (vn >= vr).sum()) / (1 + len(vn)))


def main():
    from pathlib import Path
    import pandas as pd
    O=Path(A)/'results/optimization_20261003/target_specific';I=Path(A)/'interim/rescreen_20261002'
    programs=json.loads((O/'programs.json').read_text());rows=[];seed_rows=[];loso=[];gene_rows=[];null_indices={}
    for file in sorted((O/'curves').glob('*__slopes.npz')):
        axis=file.name.replace('__slopes.npz','');z=np.load(file);D=z['D'];genes=z['genes'].astype(str).tolist();ps=z['proteins'].astype(str).tolist()
        background=np.where(~z['is_control'])[0];gidx={g:i for i,g in enumerate(genes)}
        # Retain the prior protein-panel calibration, covering every biological antigen.
        mu=D[:,background].mean(1,keepdims=True);sd=D[:,background].std(1,keepdims=True)
        Z=(D-mu)/(sd+1e-9)
        conf=[i for i,s in enumerate(z['seeds']) if s>=47];disc=[i for i,s in enumerate(z['seeds']) if s<=46]
        ens=D[conf].mean(0)
        curves=np.stack([np.load(O/'curves'/f'{axis}__seed{s}.npz')['C'] for s in z['seeds']])
        for combo,gs in programs.items():
            pname,rest=combo.split('__',1);paxis,direction=rest.rsplit('__',1)
            if paxis!=axis:continue
            ri=ps.index(pname);si=np.array([gidx[g] for g in gs]);sgn=1 if direction=='up' else -1
            pool=np.array([i for i in background if i!=ri])
            # The all-gene mean is the exact expectation of the mean of uniformly sampled fixed-size gene sets; subtract it without finite Monte Carlo noise.
            vset=sgn*Z[:,:,si].mean(2);Iscore=vset-sgn*Z.mean(2)
            panel=(Iscore-Iscore[:,background].mean(1,keepdims=True))/(Iscore[:,background].std(1,keepdims=True)+1e-12)
            mz=panel[conf].mean(0);p1=emp_p(mz[ri],mz[pool])
            isotype=np.where(z['is_control'])[0]
            control_max=float(mz[isotype].max()) if len(isotype) else np.nan
            dR=sgn*ens[ri];denom=np.mean(np.abs(dR))+1e-12
            observed=dR[si].mean()/denom
            # Original size-matched random gene sets; never reuse old q-values.
            key=(len(genes),len(si))
            if key not in null_indices:
                rng=np.random.RandomState(7)
                null_indices[key]=np.stack([rng.choice(len(genes),len(si),replace=False) for _ in range(5000)]).astype(np.int32)
            null=dR[null_indices[key]].mean(1)/denom
            p2=emp_p(observed,null)
            values=sgn*D[:,:,si].mean(2)[:,ri]
            for ii,seed in enumerate(z['seeds']):seed_rows.append(dict(combo=combo,seed=int(seed),signed_program_slope=float(values[ii]),protein_selectivity=float(panel[ii,ri]),split='discovery' if seed<=46 else 'confirmation'))
            for omitted in conf:
                others=[i for i in conf if i!=omitted];v=panel[others].mean(0)
                loso.append(dict(combo=combo,omitted_seed=int(z['seeds'][omitted]),p_protein_decoy=emp_p(v[ri],v[pool]),signed_program_slope=float(values[others].mean())))
            # Top genes selected only in discovery, tested against exhaustive protein decoys in confirmation.
            top=np.argsort(-np.abs(D[disc,ri].mean(0)))[:10]
            for gi in top:gene_rows.append(dict(combo=combo,gene=genes[gi],slope=float(ens[ri,gi]),p_protein_decoy=emp_p(abs(ens[ri,gi]),np.abs(ens[pool,gi]))))
            cv=curves[:,ri][:,:,si].mean(2)
            # Curves retain the original dose axis.
            monotonic=(sgn*np.diff(cv,axis=1)>=-1e-7).mean(1)
            rows.append(dict(combo=combo,protein=pname,axis=axis,response_direction=direction,n_genes=len(si),
                p_gene_decoy=p2,p_protein_decoy=p1,protein_selectivity=float(mz[ri]),isotype_max_selectivity=control_max,above_isotype_controls=bool(mz[ri]>control_max),
                protein_percentile=float((mz[pool]<mz[ri]).mean()),signed_program_slope=float(values[conf].mean()),
                confirmation_positive_fraction=float((values[conf]>0).mean()),
                confirmation_monotonic_fraction=float(monotonic[conf].mean()),all_seed_positive_fraction=float((values>0).mean()),
                gene_selectivity=float(observed),null_gene_mean=float(null.mean()),null_gene_sd=float(null.std()),
                n_confirmation_seeds=len(conf),n_protein_decoys=len(pool)))
        print('randKO',axis,flush=True)
    df=pd.DataFrame(rows)
    if not len(df):df=pd.DataFrame(columns=['combo','protein','axis','response_direction','p_gene_decoy','p_protein_decoy','q_gene_decoy','q_protein_decoy'])
    if len(df):
        df['q_gene_decoy']=bh(df.p_gene_decoy);df['q_protein_decoy']=bh(df.p_protein_decoy)
    df.to_csv(O/'randko_all_programs.csv',index=False)
    pd.DataFrame(seed_rows).to_csv(O/'program_seed_scores.csv',index=False)
    pd.DataFrame(loso).to_csv(O/'leave_one_seed_out.csv',index=False)
    gdf=pd.DataFrame(gene_rows)
    if len(gdf):gdf['q_BH_all_reported_genes']=bh(gdf.p_protein_decoy)
    gdf.to_csv(O/'gene_level_randko.csv',index=False)
    (O/'randko_definition.json').write_text(json.dumps(dict(discovery_seeds=list(range(42,47)),confirmation_seeds=list(range(47,52)),
        gene_decoys=5000,protein_decoys='all other biological antigens',family='all new programs across all targets, axes and both response directions',
        gene_test='signed confirmation response on frozen discovery program vs size-matched random gene sets',
        protein_test='per-gene protein panel z; program minus genome-wide mean; per-seed panel standardization; confirmation mean vs protein pool',
        n_tested=len(df)),indent=2))

if __name__=='__main__':main()
