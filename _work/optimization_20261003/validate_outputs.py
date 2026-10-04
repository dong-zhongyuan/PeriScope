"""End-to-end contracts for the actual rerun outputs; never substitutes missing evidence."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

O=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
def main():
    registry=json.loads((O/'candidate_registry.json').read_text());programs=json.loads((O/'programs.json').read_text())
    bio={r['target'] for r in registry if not r['is_control']};models=[];axes=None
    for seed in range(42,52):
        root=O/f'models/seed_{seed}';done=json.loads((root/'completed.json').read_text());vals=json.loads((root/'validation.json').read_text())
        assert done['seed']==seed and (root/'evaluated.json').exists() and (root/'screened.json').exists()
        args=json.loads((root/'training_inputs.json').read_text())['args']
        assert args['seed']==seed and args['epochs']==200 and args['steps_per_epoch']==200
        assert args['dsm_mode']==json.loads((O/'frozen_recipe.json').read_text())['dsm_mode']
        assert args['potential_beta']==5 and args['sample_latent_a'] and args['sinkhorn_iters']==20
        assert args['min_epochs']==50 and args['patience']==12 and args['pairing']=='disease' and args['coupling']=='on'
        assert done['epochs']>=50 and done['steps']<=40000
        assert abs(done['best_validation_score']-min(v['selection_score'] for v in vals))<1e-5
        files=sorted((O/'curves').glob(f'*__seed{seed}.npz'));a={f.name.rsplit('__seed',1)[0] for f in files}
        if axes is None:axes=a
        assert a==axes and len(a)>0
        for f in files:
            z=np.load(f);assert z['C'].shape==(len(registry),5,len(z['genes'])) and np.isfinite(z['C']).all()
            assert list(z['proteins'])==[r['target'] for r in registry] and z['clr_zero_sum_error']<.001
            for draw in json.loads(str(z['diagnostics'])):
                assert all(e1<=e0+1e-5 for e0,e1,t in draw['trajectory'])
        ac=pd.read_csv(O/f'evaluation/protein_accuracy_seed{seed}.csv')
        assert set(ac.donor)=={'P6','P7','P8'} and ac.n.min()>=25 and 'citeseq_subtype' in ac
        models.append(dict(seed=seed,steps=done['steps'],best_epoch=done['best_epoch'],validation_score=done['best_validation_score'],sha256=hashlib.sha256((root/'model.pt').read_bytes()).hexdigest()))
    targets=pd.read_csv(O/'all_targets.csv');assert set(targets.protein)==bio and not targets.protein.duplicated().any()
    evidence=pd.read_csv(O/'all_program_evidence.csv');assert set(evidence.combo)==set(programs)
    rand=pd.read_csv(O/'randko_all_programs.csv');assert set(rand.combo)==set(programs)
    if len(rand):
        for p,q in [('p_gene_decoy','q_gene_decoy'),('p_protein_decoy','q_protein_decoy')]:
            assert np.allclose(false_discovery_control(rand[p].to_numpy()),rand[q],atol=1e-10)
        assert (rand.n_confirmation_seeds==5).all() and (rand.n_protein_decoys==len(bio)-1).all()
    priority=pd.read_csv(O/'new_priority_targets.csv');assert priority.tier.isin(['A','B']).all()
    assert set(priority.protein)<=bio
    if len(priority):
        for name in ['technical_pass','reproducible_pass','specificity_pass','brain_pass','independent_direction_pass']:assert priority[name].all()
    coord=json.loads((O/'spatial_coordinate_checks.json').read_text());assert coord['status']=='passed' and coord['n_total']==11859
    program_hash=hashlib.sha256((O/'programs.json').read_bytes()).hexdigest()
    discovery_fingerprint=json.loads((O/'discovery_definition_fingerprint.json').read_text())
    assert all(hashlib.sha256((O/p).read_bytes()).hexdigest()==h for p,h in discovery_fingerprint.items()),'Confirmation pass changed discovery definitions'
    for name in ['brain_competitive','independent_validation','spatial_analysis']:
        cache=json.loads((O/'execution_cache'/(name+'.json')).read_text());assert cache['inputs']['programs']==program_hash
        assert all(hashlib.sha256((O/p).read_bytes()).hexdigest()==h for p,h in cache['outputs'].items())
    ext=pd.read_csv(O/'independent_validation_stats.csv')
    if len(ext):assert np.allclose(false_discovery_control(ext.p_exact_two_sided),ext.q_BH_all_programs)
    for name,p,q in [('signature_size_sensitivity.csv','p_exact_two_sided','q_all_program_sizes'),('spatial_threshold_sensitivity.csv','p_exact_two_sided','q_all_program_thresholds')]:
        frame=pd.read_csv(O/name)
        if len(frame):assert np.allclose(false_discovery_control(frame[p]),frame[q])
    genes=np.load(next((O/'curves').glob('*__seed42.npz')))['genes'].astype(str)
    assert all(set(g)<=set(genes) for g in programs.values())
    qa=dict(status='passed',n_models=len(models),n_antigens=len(registry),n_biological_antigens=len(bio),n_axes=len(axes),n_curves_files=len(axes)*10,
            n_programs=len(programs),n_priority_antigens=len(priority),models=models,full_family_BH_recomputed=True,negative_controls_excluded_from_nominations=True)
    (O/'final_quality_checks.json').write_text(json.dumps(qa,indent=2));print(json.dumps(qa),flush=True)

if __name__=='__main__':main()
