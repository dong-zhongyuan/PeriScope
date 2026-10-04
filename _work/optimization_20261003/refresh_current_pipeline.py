"""Propagate the authoritative model set into all current candidate evidence."""
from pathlib import Path
import os, sys, json, time, subprocess, traceback, fcntl, shutil
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control
from model_registry import ROOT as O, sha256, atomic_json, curve_provenance, valid_curve

W=Path(__file__).resolve().parent
STATUS=O/'pipeline_refresh_status.json'
ENV=os.environ|{'OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2',
               'PYTHONPATH':'/public/home/mengxl/dzy/pd_product/src'}

def run(script,args=(),env=None):
    p=Path(script);p=p if p.is_absolute() else W/p
    atomic_json(STATUS,{'stage':'running','step':str(p),'time':time.time()})
    with (O/('refresh_'+p.stem+'.log')).open('w') as f:
        subprocess.run([sys.executable,'-u',str(p),*map(str,args)],env=ENV|(env or {}),
                       stdout=f,stderr=subprocess.STDOUT,check=True)
    print('REFRESHED',p.name,flush=True)

def rebuild_slopes():
    seeds=list(range(42,52))
    for source in sorted((O/'curves').glob('*__seed42.npz')):
        axis=source.stem.rsplit('__seed',1)[0]
        arrays=[]
        for seed in seeds:
            p=O/'curves'/f'{axis}__seed{seed}.npz'
            assert valid_curve(p,curve_provenance(O,seed)),p
            with np.load(p) as z:arrays.append({k:z[k] for k in ['C','u','genes','proteins','is_control']})
        u=arrays[0]['u'];C=np.stack([z['C'] for z in arrays]);xx=u-u.mean(1,keepdims=True)
        D=np.einsum('spqg,pq->spg',C,xx)/np.maximum((xx**2).sum(1)[None,:,None],1e-12)
        assert np.isfinite(D).all()
        p=O/'curves'/f'{axis}__slopes.npz'
        old=np.load(p)
        # Discovery and the two unchanged confirmation models must remain identical.
        for seed in [42,43,44,45,46,50,51]:
            assert np.array_equal(D[seeds.index(seed)],old['D'][list(old['seeds']).index(seed)]), (axis,seed)
        tmp=p.with_name(p.stem+'.pending.npz')
        np.savez_compressed(tmp,D=D,u=u,seeds=seeds,genes=arrays[0]['genes'],
                            proteins=arrays[0]['proteins'],is_control=arrays[0]['is_control'])
        tmp.replace(p)
        target=O/'target_specific/curves'/p.name
        if target.is_symlink():
            assert target.resolve()==p.resolve()
        else:shutil.copy2(p,target)
        for seed in seeds:
            target=O/'target_specific/curves'/f'{axis}__seed{seed}.npz'
            source=O/'curves'/target.name
            if target.is_symlink():assert target.resolve()==source.resolve()
            elif sha256(target)!=sha256(source):shutil.copy2(source,target)
        print('SLOPES',axis,flush=True)

def update_columns(table,source,columns):
    assert table.index.is_unique and source.index.is_unique
    for c in columns:
        if c in source:table[c]=source[c].reindex(table.index)
    return table

def refresh_evidence():
    ac=pd.concat([pd.read_csv(p) for p in sorted((O/'evaluation').glob('protein_accuracy_seed*.csv'))])
    assert set(ac.seed)==set(range(42,52))
    ac=ac.groupby(['target','blood_state','donor'],as_index=False)[['rho','mse']].mean()
    ac['positive']=ac.rho.gt(0)
    bridge=ac.groupby(['target','blood_state']).agg(protein_rho=('rho','mean'),protein_mse=('mse','mean'),
                    positive_protein_donors=('positive','sum'),n_protein_donors=('donor','size'))
    assoc=pd.read_csv(O/'blood_disease_associations.csv')
    assoc=assoc[assoc.cohort.astype(str).eq('combined')].set_index(['target','blood_state','assay'])
    joined=[]
    for label,root in [('raw',O),('relative',O/'target_specific')]:
        r=pd.read_csv(root/'randko_all_programs.csv').set_index('combo')
        effect=pd.read_csv(root/'program_effect_sizes.csv').set_index('combo')
        e=pd.read_csv(root/'all_program_evidence.csv').set_index('combo')
        assert set(e.index)==set(r.index)
        e=update_columns(e,r,list(r));e=update_columns(e,effect,list(effect))
        for idx,row in e.iterrows():
            key=(row.protein,row.blood_state)
            if key in bridge.index:
                for col,val in bridge.loc[key].items():e.at[idx,col]=val
            for assay in ['observed_RNA','inferred_ADT']:
                ak=(*key,assay)
                if ak in assoc.index:
                    e.at[idx,'blood_'+assay+'_effect']=assoc.loc[ak,'effect_PD_minus_control_adjusted']
                    e.at[idx,'blood_'+assay+'_q']=assoc.loc[ak,'q_BH']
        # Diagnostic flags remain descriptive and never enter active C1--C6/MR eligibility.
        e['protein_prediction_pass']=e.protein_rho.ge(.2)&e.positive_protein_donors.ge(2)
        e['reproducible_pass']=e.confirmation_positive_fraction.ge(.8)&e.confirmation_monotonic_fraction.ge(.75)
        e['specificity_pass']=e.q_gene_decoy.le(.05)
        e['technical_pass']=e.mapping_pass&e.measured_dose_pass
        e.to_csv(root/'all_program_evidence.csv',index_label='combo')
        x=e.copy();x.index=label+'::'+x.index
        joined.append(x)
    dynamic=pd.concat(joined)
    joint_path=O/'joint_analysis/all_program_evidence.csv'
    e=pd.read_csv(joint_path).set_index('combo')
    assert set(e.index)==set(dynamic.index) and len(e)==6849
    raw_columns=list(pd.read_csv(O/'randko_all_programs.csv',nrows=0).columns)
    effect_columns=list(pd.read_csv(O/'program_effect_sizes.csv',nrows=0).columns)
    cols=[c for c in raw_columns+effect_columns if c not in ['combo','protein','axis','response_direction','n_genes']]
    cols+=list(bridge)+['blood_'+a+'_'+v for a in ['observed_RNA','inferred_ADT'] for v in ['effect','q']]
    cols+=['protein_prediction_pass','reproducible_pass','specificity_pass','technical_pass']
    update_columns(e,dynamic,cols)
    for p,q in [('p_gene_decoy','q_gene_decoy'),('p_protein_decoy','q_protein_decoy')]:
        e[q+'_within_branch']=dynamic[q].reindex(e.index)
        e[q]=false_discovery_control(e[p].to_numpy())
    e['specificity_pass']=e.q_gene_decoy.le(.05)
    e.to_csv(joint_path,index_label='combo')
    for name in ['program_seed_scores.csv','leave_one_seed_out.csv','gene_level_randko.csv']:
        parts=[]
        for label,root in [('raw',O),('relative',O/'target_specific')]:
            part=pd.read_csv(root/name);part['combo']=label+'::'+part.combo;parts.append(part)
        d=pd.concat(parts,ignore_index=True)
        if name=='gene_level_randko.csv':d['q_BH_all_reported_genes']=false_discovery_control(d.p_protein_decoy.to_numpy())
        d.to_csv(O/'joint_analysis'/name,index=False)
    # Do not recalculate expression/spatial evidence: fixed memberships are unchanged.
    atomic_json(O/'joint_analysis/model_refresh.json',{'status':'updated','programs':len(e),
        'model_registry_sha256':sha256(O/'model_registry.json'),'table_sha256':sha256(joint_path),
        'gene_decoy_FDR_programs':int(e.q_gene_decoy.le(.05).sum()),
        'protein_panel_metrics_role':'evaluation_only'})

def rejoin_spatial():
    source=O/'six_criteria_20261003/all_6849_programs_six_criteria.csv'
    six=pd.read_csv(source).set_index('combo');S=O/'seven_criteria_20261003'
    old=pd.read_csv(S/'all_6849_programs_seven_criteria.csv').set_index('combo')
    assert set(six.index)==set(old.index)
    e=old.copy()
    for c in six:e[c]=six[c]
    e['seven_criteria_pass']=e.six_criteria_pass&e.C7.eq('pass')
    e['seven_criteria_approved_subset']=e.six_criteria_approved_subset&e.C7.eq('pass')
    e['criteria_pass_prefix_length_seven']=e[[f'C{i}' for i in range(1,8)]].eq('pass').astype(int).cumprod(axis=1).sum(axis=1)
    if 'unmet_seven_criteria' in e:e['unmet_seven_criteria']=[';'.join(f'{c}:{r[c]}' for c in [f'C{i}' for i in range(1,8)] if r[c]!='pass') for r in e.to_dict('records')]
    e=e.reset_index().sort_values(['seven_criteria_pass','q_gene_decoy','matched_brain_q','protein','combo'],ascending=[False,True,True,True,True])
    selections={'all_6849_programs_seven_criteria.csv':e,'all_213_antigen_dispositions.csv':e.drop_duplicates('protein'),
        'A7_all_programs.csv':e[e.seven_criteria_pass],'A7_antigens.csv':e[e.seven_criteria_pass].drop_duplicates('protein'),
        'A7_approved_drug_subset.csv':e[e.seven_criteria_approved_subset].drop_duplicates('protein'),
        'six_to_seven_changes.csv':e[e.six_criteria_pass],'CD22_CD35_seven_criteria.csv':e[e.protein.isin(['CD22','CD35'])],
        'five_biological_plus_spatial_without_drug_requirement.csv':e[e.five_criteria_pass&e.C7.eq('pass')],
        'C7_unresolved_programs.csv':e[e.C7.eq('unresolved')],
        'A6_pending_spatial_resolution.csv':e[e.six_criteria_pass&e.C7.eq('unresolved')],
        'A6_tested_without_spatial_support.csv':e[e.six_criteria_pass&e.C7.eq('not_supported')],
        'spatial_supported_previously_excluded_programs.csv':e[e.C7.eq('pass')&~e.six_criteria_pass]}
    for name,d in selections.items():d.to_csv(S/name,index=False)
    mixed=pd.read_csv(S/'mixed_tissue_vs_matched_state.csv');e[list(mixed)].to_csv(S/'mixed_tissue_vs_matched_state.csv',index=False)
    count=lambda d:{'programs':len(d),'antigens':int(d.protein.nunique())}
    flow=[dict(step='all',**count(e))];keep=pd.Series(True,index=e.index)
    for i in range(1,8):
        before=e[keep];keep &= e[f'C{i}'].eq('pass')
        flow.append(dict(step=f'C{i}',**count(e[keep]),removed_programs=len(before)-int(keep.sum()),removed_antigens=int(before.protein.nunique()-e[keep].protein.nunique())))
    pd.DataFrame(flow).to_csv(S/'seven_criteria_funnel.csv',index=False)
    summary=json.loads((S/'completed.json').read_text());a=e[e.seven_criteria_pass]
    summary.update(source_six_criteria_sha256=sha256(source),A6=count(e[e.six_criteria_pass]),A7=count(a),
        approved_subset=count(e[e.seven_criteria_approved_subset]),flow=flow,A7_antigens=sorted(a.protein.unique()),
        A6_C7_status=e[e.six_criteria_pass].C7.value_counts().to_dict(),CD22=count(a[a.protein.eq('CD22')]),
        CD35=count(a[a.protein.eq('CD35')]),model_registry_sha256=sha256(O/'model_registry.json'))
    for p in [S/'completed.json',O/'spatial_subtype_author_20261003/completed.json']:atomic_json(p,summary)
    summary['A7_previously_opposed_in_mixed_tissue']=count(a[a.old_mixed_spatial_direction_opposed_matched_brain])
    atomic_json(S/'completed.json',summary)
    # Compact summary is regenerated from the same current rows.
    brief=a.drop_duplicates('protein')[['protein','genes','blood_state','brain_state','q_gene_decoy','matched_brain_q','spatial_q','effect','drug_examples','approved_examples']]
    brief.to_csv(S/'七项通过候选汇总.csv',index=False,encoding='utf-8-sig')

def main():
    lock=(O/'pipeline_refresh.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    protected=[O/'joint_analysis/programs.json',O/'programs.json',O/'target_specific/programs.json',
        O.parent/'blood_brain_MR_20261003/all_brain_MR_estimates.csv',
        O/'cell_state_sensitivity/raw_discovery.csv',O/'cell_state_sensitivity/target_specific_discovery.csv',
        O/'spatial_subtype_author_20261003/all_leaf_program_spatial_tests.csv']
    fixed={str(p):sha256(p) for p in protected}
    try:
        while True:
            status=json.loads((O/'model_refresh_status.json').read_text())
            if status['stage']=='failed':raise RuntimeError(status)
            if status['stage']=='model_outputs_complete':break
            time.sleep(20)
        if '--verify-only' not in sys.argv:
            atomic_json(STATUS,{'stage':'running','step':'slopes'})
            rebuild_slopes()
            for p in ['randko.py','target_specific/randko.py','program_effect_sizes.py','target_specific/program_effect_sizes.py',
                      'disease_association.py','blood_covariate_sensitivity.py']:
                run(p)
            refresh_evidence()
            run('screen_six_criteria.py',['--result-dir',O])
            rejoin_spatial()
            run(W.parent/'spatial_subtype_20261003/verify_spatial_delivery.py',env={'SPATIAL_OUTPUT_DIR':str(O/'spatial_subtype_author_20261003')})
            run('apply_mr_parallel_gate.py',['--results-root',O.parent])
            run('refresh_pig_seeds.py')
        run('verify_current_pipeline.py')
        run('randko_seed_robustness.py')
        for p,h in fixed.items():assert sha256(p)==h,p
        registry=json.loads((O/'model_registry.json').read_text())
        for rel,h in registry['discovery_files'].items():assert sha256(O/rel)==h,rel
        pointer=json.loads((O/'active_nomination.json').read_text());pointer.update(status='complete',model_registry='model_registry.json',
            model_registry_sha256=sha256(O/'model_registry.json'),pipeline_entrypoint=str(W/'refresh_current_pipeline.py'))
        atomic_json(O/'active_nomination.json',pointer)
        atomic_json(STATUS,{'stage':'complete','model_registry_sha256':sha256(O/'model_registry.json'),
            'unchanged_evidence_hashes':fixed,'completed':time.time(),
            'summary':json.loads((O/'mr_parallel_gate_20261003/summary.json').read_text())})
        print('PIPELINE_REFRESH_COMPLETE',flush=True)
    except Exception:
        atomic_json(STATUS,{'stage':'failed','traceback':traceback.format_exc()});raise

if __name__=='__main__':main()
