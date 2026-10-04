"""Expanded-panel specificity discovery, defined before corrected-run screening."""

# Current authoritative pipeline dispatch (imports retain helper compatibility).
if __name__ == '__main__':
    from pathlib import Path as _RegistryPath
    import subprocess as _RegistrySubprocess, sys as _RegistrySys
    if _RegistryPath('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/model_registry.json').exists():
        _RegistrySubprocess.run([_RegistrySys.executable, str(_RegistryPath(__file__).with_name('refresh_model_outputs.py'))], check=True)
        _RegistrySubprocess.run([_RegistrySys.executable, str(_RegistryPath(__file__).with_name('refresh_current_pipeline.py'))], check=True)
        raise SystemExit(0)

from pathlib import Path
import os,json,time,subprocess,sys,hashlib
from concurrent.futures import ThreadPoolExecutor
W=Path(__file__).parent;S=W/'target_specific';P=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003');O=P/'target_specific';O.mkdir(exist_ok=True)
protocol=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
    selection='Primary optimized discovery definition; raw-slope discovery retained as the training-only comparator. Definition fixed before corrected models produce any dose curves.',
    rationale='Large absolute response may be common to many antigens. Discover target-relative programs in discovery seeds, then test them using independent confirmation model seeds.',
    discovery_seeds=[42,43,44,45,46],confirmation_seeds=[47,48,49,50,51],
    ranking='Within each seed and gene, subtract the mean slope of the other biological antigens and divide by their SD; average standardized values within each discovery fold. Rank signed contrasts.',
    sign='Require the corresponding fold-mean raw dose slope to have the same sign. Both up and down are tested.',
    other_rules='Keep top-200 maximum, ontology libraries, hypergeometric background, complete overlap genes, fold union, held-out confirmation tests, all biological gates and thresholds unchanged.',
    code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in S.glob('*.py')})
pf=O/'target_specific_protocol.json'
if not pf.exists():
    assert not list((P/'curves').glob('*__seed*.npz')),'Protocol must be frozen before corrected-run screening'
    pf.write_text(json.dumps(protocol,indent=2))
for name in ['candidate_registry.json','dose_definitions.json','all_candidate_input_support.csv','data_checks.json','downstream_input_checks.json',
             'spatial_coordinate_checks.json','all_target_assay_coverage.csv','drug_gene_interactions.csv','frozen_recipe.json','evaluation','models','ablations','drug_annotation_source.json']:
    f=O/name
    if not f.exists() and not f.is_symlink():f.symlink_to(P/name,target_is_directory=name in ['models','evaluation','ablations'])
env=os.environ|dict(PYTHONPATH='/public/home/mengxl/dzy/pd_product/src',OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
def run(script,extra=None):
    with (O/(script.replace('.py','')+'.log')).open('a') as log:
        subprocess.run([sys.executable,'-u',str(S/script)],env=env|(extra or {}),stdout=log,stderr=subprocess.STDOUT,check=True)
    print('TARGET-SPECIFIC COMPLETE',script,flush=True)
def wait_seeds(seeds):
    while not all((P/f'models/seed_{s}/screened.json').exists() for s in seeds):
        if json.loads((P/'status.json').read_text()).get('stage')=='failed':raise RuntimeError('Main run failed')
        time.sleep(20)
    (O/'curves').mkdir(exist_ok=True)
    for s in seeds:
        for source in (P/'curves').glob(f'*__seed{s}.npz'):
            dest=O/'curves'/source.name
            if not dest.exists():dest.symlink_to(source)
wait_seeds(range(42,47))
run('dose_programs.py',dict(CCWM_SEEDS='42,43,44,45,46'))
with ThreadPoolExecutor(4) as pool:
    fs=[pool.submit(run,s) for s in ['brain_competitive.py','independent_validation.py','spatial_analysis.py']]
    for f in fs:f.result()
(O/'discovery_evidence_completed.json').write_text(json.dumps(dict(status='complete')))
wait_seeds(range(42,52))
run('dose_programs.py')
run('randko.py')
while not (P/'final_quality_checks.json').exists():
    if json.loads((P/'status.json').read_text()).get('stage')=='failed':raise RuntimeError('Main run failed')
    time.sleep(20)
for name in ['blood_disease_associations.csv','blood_disease_donor_values.csv','target_gene_expression_by_donor.csv','tissue_axis_associations.csv','blood_disease_covariate_sensitivity.csv']:
    dest=O/name
    if not dest.exists() and not dest.is_symlink():dest.symlink_to(P/name)
run('rank_candidates.py');run('compare_ablations.py')
while not (P/'blood_disease_covariate_sensitivity.csv').exists():time.sleep(15)
run('review_results.py')
(O/'optimization_completed.json').write_text(json.dumps(dict(status='complete')))
