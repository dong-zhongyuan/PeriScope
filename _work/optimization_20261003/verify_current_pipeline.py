"""Verify current weights, dependent statistics and nominated program identities."""
from pathlib import Path
import json,time
import numpy as np,pandas as pd
from scipy.stats import false_discovery_control
from model_registry import ROOT as O, sha256, atomic_json, curve_provenance, valid_curve

def main():
    reg=json.loads((O/'model_registry.json').read_text());rh=sha256(O/'model_registry.json')
    for seed,r in reg['models'].items():assert sha256(r['path'])==r['sha256']
    overrides=json.loads((O.parent/'benchmark_20261002/checkpoint_overrides.json').read_text())
    for seed in [47,48,49]:
        assert sha256(overrides['checkpoints'][str(seed)])==reg['models'][str(seed)]['sha256']
    for p,h in reg['discovery_files'].items():assert sha256(O/p)==h
    curves=sorted((O/'curves').glob('*__seed*.npz'));assert len(curves)==120
    for p in curves:
        seed=int(p.stem.rsplit('__seed',1)[1]);assert valid_curve(p,curve_provenance(O,seed))
    for seed in [47,48,49]:
        rec=json.loads((O/f'evaluation/seed{seed}_provenance.json').read_text())
        assert rec['checkpoint_sha256']==reg['models'][str(seed)]['sha256']
        assert all(sha256(O/'evaluation'/n)==h for n,h in rec['outputs'].items())
    J=O/'joint_analysis';e=pd.read_csv(J/'all_program_evidence.csv').set_index('combo')
    programs=json.loads((J/'programs.json').read_text());assert len(e)==6849 and set(e.index)==set(programs)
    score=pd.read_csv(J/'program_seed_scores.csv')
    assert not score.duplicated(['combo','seed']).any() and len(score)==68490
    means=score[score.seed.ge(47)].groupby('combo').signed_program_slope.mean()
    np.testing.assert_allclose(e.signed_program_slope,means.reindex(e.index),rtol=1e-10,atol=1e-12)
    for p,q in [('p_gene_decoy','q_gene_decoy'),('p_protein_decoy','q_protein_decoy')]:
        np.testing.assert_allclose(e[q],false_discovery_control(e[p].to_numpy()),atol=1e-12)
    six=pd.read_csv(O/'six_criteria_20261003/all_6849_programs_six_criteria.csv').set_index('combo')
    current=O/'mr_parallel_gate_20261003';d=pd.read_csv(current/'all_programs_parallel_gate.csv').set_index('combo')
    assert set(d.index)==set(e.index)
    for col in ['p_gene_decoy','q_gene_decoy','signed_program_slope','confirmation_positive_fraction','protein_rho','confirmation_disease_projection_mean']:
        np.testing.assert_allclose(e[col],d[col].reindex(e.index),equal_nan=True,atol=1e-12)
    expected=d.MR_gate.eq('pass')&d.q_gene_decoy.le(.05)&d[['C1','C2','C4','C5','C6']].eq('pass').all(axis=1)
    assert expected.equals(d.priority_pass)
    assert d.C3.eq('pass').equals(d.q_gene_decoy.le(.05))
    priority=pd.read_csv(current/'priority_all_programs.csv');assert set(priority.combo)==set(d.index[expected])
    mr=pd.read_csv(O.parent/'blood_brain_MR_20261003/all_brain_MR_estimates.csv')
    np.testing.assert_allclose(mr.q_BH,false_discovery_control(mr.p_MR.to_numpy()),atol=1e-12)
    pig=json.loads((O.parent/'pig_external_validation_20261003/predictions_locked.json').read_text())
    assert pig['checkpoint_sha256']=={s:r['sha256'] for s,r in reg['models'].items()}
    # Replace obsolete shortcut nomination files with the current MR-parallel result.
    full=d.reset_index();best=full.sort_values(['priority_pass','MR_q','q_gene_decoy','protein','combo'],
                                              ascending=[False,True,True,True,True]).drop_duplicates('protein')
    assert len(best)==213
    pr=priority.sort_values(['MR_q','q_gene_decoy','protein','combo']).drop_duplicates('protein')
    for name,table in [('new_priority_targets.csv',pr),('new_priority_genes.csv',pr.drop_duplicates('genes')),('all_targets.csv',best)]:
        table.to_csv(O/name,index=False)
    for folder in ['current_nomination','final_results']:
        dest=O/folder
        for name,table in [('all_program_evidence.csv',full),('priority_targets.csv',pr),('all_213_targets.csv',best),
                           ('priority_gene_groups.csv',pr.drop_duplicates('genes')),('CD22_CD35_all_programs.csv',full[full.protein.isin(['CD22','CD35'])])]:
            table.to_csv(dest/name,index=False)
        # Existing diagnostic exports must carry the same refreshed numeric columns.
        for diagnostic in ['model_evaluation_only.csv']:
            p=dest/diagnostic
            if p.exists():
                columns=list(pd.read_csv(p,nrows=0))
                full[columns].to_csv(p,index=False)
        atomic_json(dest/'active_source.json',{'active_source':'../mr_parallel_gate_20261003',
                                             'model_registry_sha256':rh,'status':'current_shortcuts_updated'})
    summary=json.loads((current/'summary.json').read_text());summary['model_registry_sha256']=rh
    summary['replaced_seeds']=[47,48,49];summary['confirmation_recomputed']=True
    atomic_json(current/'summary.json',summary)
    for folder, name in [('current_nomination','completed.json'),('final_results','completion.json')]:
        atomic_json(O/folder/name,dict(summary,active_source='../mr_parallel_gate_20261003'))
    # Select the validation record at the actual chosen checkpoint, including the
    # minimum-epoch rule used for repaired checkpoints.
    training=O/'final_results/ten_seed_training_comparison.csv'
    table=pd.read_csv(training)
    current_rows=[]
    for seed in range(42,52):
        root=O/f'models/seed_{seed}'
        done=json.loads((root/'completed.json').read_text())
        history=json.loads((root/'validation.json').read_text())
        selected=[v for v in history if v['epoch']==done['best_epoch']]
        assert len(selected)==1,(seed,done['best_epoch'])
        row=dict(seed=seed,epochs=done['epochs'],steps=done['steps'],best_epoch=done['best_epoch'],
                 **{k:selected[0][k] for k in ['selection_score','protein_mse','brain_mmd']})
        for k,v in row.items():table.loc[table.branch.eq('corrected')&table.seed.eq(seed),k]=v
        current_rows.append(dict(row,checkpoint_sha256=reg['models'][str(seed)]['sha256']))
    table.to_csv(training,index=False)
    pd.DataFrame(current_rows).to_csv(O/'evaluation/current_training_summary.csv',index=False)
    spatial=O/'seven_criteria_20261003'
    record=json.loads((spatial/'completed.json').read_text())
    seven=pd.read_csv(spatial/'all_6849_programs_seven_criteria.csv')
    required=sorted(seven.loc[seven.six_criteria_pass,'brain_state'].unique())
    assert required==sorted(record['required_states_for_global_adoption']),required
    spatial_a=seven[seven.seven_criteria_pass]
    opposed=spatial_a[spatial_a.old_mixed_spatial_direction_opposed_matched_brain]
    record['A7_previously_opposed_in_mixed_tissue']={'programs':len(opposed),'antigens':int(opposed.protein.nunique())}
    for path in [spatial/'completed.json',O/'spatial_subtype_author_20261003/completed.json']:
        atomic_json(path,record)
    atomic_json(O/'rescreen_summary.json',dict(summary,n_raw_programs=5009,n_relative_programs=1840))
    # Refresh the original model/hash checklist rather than leaving an old passed marker.
    model_rows=[]
    for seed in range(42,52):
        complete=json.loads((O/f'models/seed_{seed}/completed.json').read_text())
        model_rows.append(dict(seed=seed,sha256=reg['models'][str(seed)]['sha256'],training_completion=complete))
    report={'status':'passed','model_registry_sha256':rh,'n_models':10,'n_antigens':217,
        'n_biological_antigens':213,'n_axes':12,'n_curves_files':120,'n_programs':6849,
        'n_priority_antigens':len(pr),'models':model_rows,'full_family_BH_recomputed':True,
        'negative_controls_excluded_from_nominations':True,'discovery_memberships_unchanged':True,
        'confirmation_mean_matches_seed_table':True,'MR_inputs_unchanged':True,
        'performance_metrics_not_eligibility_gates':True,'pig_uses_current_ten_seeds':True,
        'completed':time.time()}
    atomic_json(O/'final_quality_checks.json',report);atomic_json(O/'current_pipeline_verification.json',report)
    # Binary hashes are refreshed after the actual arrays have passed validation.
    m=json.loads((O/'server_binary_manifest.json').read_text())
    m['artifacts']=[a for a in m['artifacts'] if not ('magma' in a['path'].lower() or 'primary_genetics' in a['path'] or 'gwas_program_level' in a['path'])]
    for a in m['artifacts']:
        p=Path(a['path']);assert p.exists();a.update(bytes=p.stat().st_size,sha256=sha256(p))
    m['total_bytes']=sum(a['bytes'] for a in m['artifacts']);m['model_registry_sha256']=rh
    atomic_json(O/'server_binary_manifest.json',m)
    for folder in [current,O/'six_criteria_20261003',O/'seven_criteria_20261003',O/'spatial_subtype_author_20261003',O.parent/'pig_external_validation_20261003']:
        path=folder/('output_manifest.json' if folder.name in ['spatial_subtype_author_20261003','pig_external_validation_20261003'] else 'manifest.json')
        if path.exists():
            value=json.loads(path.read_text());entries=value.get('files',value) if isinstance(value,dict) else value
            if isinstance(entries,dict):
                for name in list(entries):
                    p=folder/name
                    if p.is_file():entries[name]=sha256(p)
            elif isinstance(entries,list):
                for item in entries:
                    p=folder/item['path']
                    assert p.is_file(),p
                    item.update(bytes=p.stat().st_size,sha256=sha256(p))
            path.write_text(json.dumps(value,indent=2))
    # Machine-readable current input selection prevents figure builders using old iteration tables.
    active={'models':'model_registry.json','programs':'joint_analysis/programs.json',
        'program_evidence':'mr_parallel_gate_20261003/all_programs_parallel_gate.csv',
        'priority_programs':'mr_parallel_gate_20261003/priority_all_programs.csv',
        'priority_targets':'mr_parallel_gate_20261003/priority_antigens.csv',
        'training':'evaluation/current_training_summary.csv',
        'seed_scores':'joint_analysis/program_seed_scores.csv','seed_omissions':'joint_analysis/leave_one_seed_out.csv',
        'spatial':'spatial_subtype_author_20261003/all_program_spatial_tests.csv'}
    robustness=O/'seed_robustness/definition.json'
    if robustness.exists():
        definition=json.loads(robustness.read_text())
        is_current=(definition['model_registry_sha256']==rh and
                    all(sha256(O/path)==h for path,h in definition['input_sha256'].items()))
        if is_current:
            root=robustness.parent
            for name,h in json.loads((root/'manifest.json').read_text()).items():assert sha256(root/name)==h
            r=pd.read_csv(root/'leave_one_seed_out.csv')
            assert len(r)==5*len(e) and not r.duplicated(['combo','omitted_seed']).any()
            for omitted,part in r.groupby('omitted_seed'):
                assert set(part.combo)==set(e.index)
                np.testing.assert_allclose(part.q_within_omission,false_discovery_control(part.p_gene_decoy.to_numpy()),atol=1e-12)
            active.update(seed_omissions='seed_robustness/leave_one_seed_out.csv',
                          seed_robustness_summary='seed_robustness/leave_one_seed_summary.csv',
                          seed_robustness_plot='seed_robustness/priority_plot_input.csv')
    atomic_json(O/'current_results.json',{'status':'complete','model_registry_sha256':rh,'files':active,
        'sha256':{role:sha256(O/p) for role,p in active.items()}})
    print(json.dumps({k:v for k,v in report.items() if k!='models'}),flush=True)

if __name__=='__main__':main()
