"""CP-derived randKO: recompute gene-null FDR after omitting each current confirmation seed."""
import os
import sys
import json
import glob
import time
import argparse

import numpy as np
from scipy.stats import combine_pvalues, wilcoxon

A = '/public/home/mengxl/dzy/pd_product_assets'
CACHE = A + '/results/optimization_20261003/curves'
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
    from model_registry import ROOT as O, sha256, atomic_json
    out = O / 'seed_robustness'
    out.mkdir(exist_ok=True)
    registry_hash = sha256(O/'model_registry.json')
    verified = json.loads((O/'current_pipeline_verification.json').read_text())
    assert verified['model_registry_sha256'] == registry_hash
    programs = json.loads((O/'joint_analysis/programs.json').read_text())
    evidence = pd.read_csv(O/'joint_analysis/all_program_evidence.csv').set_index('combo')
    priority = pd.read_csv(O/'mr_parallel_gate_20261003/priority_all_programs.csv')
    assert set(programs) == set(evidence.index)
    manifest = json.loads((O/'server_binary_manifest.json').read_text())
    binary_hashes = {r['path']:r['sha256'] for r in manifest['artifacts']}
    hashes = {name:sha256(O/name) for name in ['model_registry.json','joint_analysis/programs.json',
        'joint_analysis/all_program_evidence.csv','mr_parallel_gate_20261003/priority_all_programs.csv']}
    rows = []
    null_indices = {}
    full_checks = []
    atomic_json(out/'status.json', dict(status='running', model_registry_sha256=registry_hash))
    for file in sorted((O/'curves').glob('*__slopes.npz')):
        hashes[str(file.relative_to(O))] = sha256(file)
        # Check against the refreshed binary inventory before evaluating the arrays.
        assert hashes[str(file.relative_to(O))] == binary_hashes[str(file)]
        axis = file.name.replace('__slopes.npz','')
        with np.load(file) as z:
            D=z['D'];genes=z['genes'].astype(str).tolist();ps=z['proteins'].astype(str).tolist();seeds=z['seeds']
        assert list(seeds)==list(range(42,52)) and np.isfinite(D).all()
        gidx={g:i for i,g in enumerate(genes)}
        conf=[i for i,s in enumerate(seeds) if s>=47]
        ensembles={0:D[conf].mean(0)}
        for omitted in conf:
            ensembles[int(seeds[omitted])] = D[[i for i in conf if i!=omitted]].mean(0)
        for combo,gs in programs.items():
            branch,source_combo=combo.split('::',1)
            pname,rest=source_combo.split('__',1);paxis,direction=rest.rsplit('__',1)
            if paxis!=axis:continue
            ri=ps.index(pname);si=np.array([gidx[g] for g in gs]);sgn=1 if direction=='up' else -1
            key=(len(genes),len(si))
            if key not in null_indices:
                rng=np.random.RandomState(7)
                null_indices[key]=np.stack([rng.choice(len(genes),len(si),replace=False) for _ in range(5000)]).astype(np.int32)
            for omitted,ens in ensembles.items():
                # Identical signed statistic and Monte Carlo null to randko.py.
                dR=sgn*ens[ri];denom=np.mean(np.abs(dR))+1e-12
                observed=dR[si].mean()/denom
                null=dR[null_indices[key]].mean(1)/denom
                p2=emp_p(observed,null)
                rec=dict(combo=combo,program_definition=branch,protein=pname,axis=axis,
                    response_direction=direction,n_genes=len(si),omitted_seed=omitted,
                    n_confirmation_seeds=5 if omitted==0 else 4,selectivity=float(observed),
                    signed_program_slope=float(dR[si].mean()),p_gene_decoy=p2)
                rows.append(rec)
                if omitted==0:
                    np.testing.assert_allclose(p2,evidence.loc[combo,'p_gene_decoy'],atol=1e-15,rtol=0)
                    np.testing.assert_allclose(observed,evidence.loc[combo,'gene_selectivity'],atol=1e-12,rtol=1e-12)
                    full_checks.append(combo)
        print('SEED_ROBUSTNESS',axis,len(rows),flush=True)
        atomic_json(out/'status.json',dict(status='running',last_axis=axis,rows=len(rows),model_registry_sha256=registry_hash))
    df=pd.DataFrame(rows)
    assert len(full_checks)==len(programs) and len(df)==6*len(programs)
    assert not df.duplicated(['combo','omitted_seed']).any()
    for omitted,ix in df.groupby('omitted_seed').groups.items():
        assert set(df.loc[ix,'combo'])==set(programs)
        df.loc[ix,'q_within_omission']=bh(df.loc[ix,'p_gene_decoy'].to_numpy())
    full=df[df.omitted_seed.eq(0)].set_index('combo')
    np.testing.assert_allclose(full.q_within_omission,evidence.q_gene_decoy.reindex(full.index),atol=1e-12,rtol=1e-12)
    loo=df[df.omitted_seed.ne(0)].copy()
    loo['FDR_pass']=loo.q_within_omission.le(.05)
    summary=loo.groupby('combo').agg(retained_runs=('FDR_pass','sum'),n_omissions=('omitted_seed','size'),
        min_selectivity=('selectivity','min'),max_q=('q_within_omission','max')).reset_index()
    summary['retention_percent']=100*summary.retained_runs/summary.n_omissions
    summary['full_q']=summary.combo.map(full.q_within_omission)
    summary['full_FDR_pass']=summary.full_q.le(.05)
    summary['priority_program']=summary.combo.isin(priority.combo)
    summary['protein']=summary.combo.map(evidence.protein)
    summary['axis']=summary.combo.map(evidence.axis)
    summary['program_definition']=summary.combo.map(evidence.program_definition)
    summary=summary.sort_values(['protein','axis','combo'])
    full.reset_index().to_csv(out/'all_five_confirmation.csv',index=False)
    loo.to_csv(out/'leave_one_seed_out.csv',index=False)
    summary.to_csv(out/'leave_one_seed_summary.csv',index=False)
    # All final programs selected by the unchanged MR/perturbation criteria,
    # irrespective of the leave-one-seed-out result.
    selected=summary[summary.priority_program]
    selected.to_csv(out/'priority_program_retention.csv',index=False)
    plot=loo[loo.combo.isin(priority.combo)].copy()
    plot['p']=plot.p_gene_decoy
    plot.to_csv(out/'priority_plot_input.csv',index=False)
    statistics=[]
    for label,t in [('all_programs',summary),('full_FDR_programs',summary[summary.full_FDR_pass]),('priority_programs',selected)]:
        for k in range(6):
            statistics.append(dict(population=label,retained_runs=k,n_programs=int(t.retained_runs.eq(k).sum()),total_programs=len(t)))
    pd.DataFrame(statistics).to_csv(out/'retention_distribution.csv',index=False)
    omission_summary=loo.groupby('omitted_seed').agg(n_programs=('combo','size'),FDR_pass=('FDR_pass','sum')).reset_index()
    omission_summary.to_csv(out/'omission_summary.csv',index=False)
    for name,h in hashes.items():assert sha256(O/name)==h,name
    definition=dict(status='complete',model_registry_sha256=registry_hash,n_programs=len(programs),n_gene_decoys=5000,
        discovery_seeds=list(range(42,47)),confirmation_seeds=list(range(47,52)),omitted_seeds=list(range(47,52)),
        random_seed=7,statistic='Signed mean confirmation slope within frozen program / mean absolute whole-gene slope',
        null='5000 uniformly sampled gene sets without replacement, matched to program size; identical draws for each gene universe and set size in every omission',
        multiplicity='Separate BH over all 6849 raw+relative programs in each omission; never BH on the displayed subset',
        threshold=.05,role='Model seed sensitivity evaluation only; no new nomination gate',
        plot_selection='All programs already passing the current MR+perturbation and C1/C2/C4/C5/C6 criteria, without selection on this robustness result',
        original_plot='panel_grid/source/round2_panels.R: original 4E, delivered 4f, planned new 5f',
        all_five_p_and_q_reproduce_current_tables=True,
        input_sha256=hashes,source_sha256=sha256(__file__))
    atomic_json(out/'definition.json',definition)
    report=dict(status='passed',n_programs=len(programs),n_omission_tests=len(loo),n_full_tests=len(full),
        full_FDR_programs=int(summary.full_FDR_pass.sum()),all_five_omissions_retained=int(summary.retained_runs.eq(5).sum()),
        full_FDR_and_all_five_retained=int((summary.full_FDR_pass & summary.retained_runs.eq(5)).sum()),
        priority_programs=len(selected),priority_all_five_retained=int(selected.retained_runs.eq(5).sum()),
        priority_retention_histogram={str(k):int(v) for k,v in selected.retained_runs.value_counts().sort_index().items()},
        model_registry_sha256=registry_hash,all_five_matches_current_statistics=True)
    atomic_json(out/'verification.json',report)
    atomic_json(out/'manifest.json',{p.name:sha256(p) for p in out.iterdir() if p.is_file() and p.name not in ['manifest.json','status.json']})
    atomic_json(out/'status.json',dict(status='complete',model_registry_sha256=registry_hash))
    active=json.loads((O/'current_results.json').read_text())
    for role,name in [('seed_omissions','seed_robustness/leave_one_seed_out.csv'),('seed_robustness_summary','seed_robustness/leave_one_seed_summary.csv'),('seed_robustness_plot','seed_robustness/priority_plot_input.csv')]:
        active['files'][role]=name;active['sha256'][role]=sha256(O/name)
    atomic_json(O/'current_results.json',active)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
