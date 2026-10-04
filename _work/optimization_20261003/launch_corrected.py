"""Freeze a repaired DSM implementation on pilot validation, then rerun all seeds."""
from pathlib import Path
import os,json,sys,subprocess,time,hashlib
from concurrent.futures import ThreadPoolExecutor
W=Path(__file__).parent
O=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
env=os.environ|dict(PYTHONPATH='/public/home/mengxl/dzy/pd_product/src',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')

def execute(name,args=()):
    with (W/(name.replace('.py','')+'.log')).open('a') as f:
        subprocess.run([sys.executable,'-u',str(W/name),*args],stdout=f,stderr=subprocess.STDOUT,env=env,check=True)

while not (O/'pilot_comparison.json').exists():
    if not Path('/proc/89302/cmdline').exists():raise RuntimeError('Pilot controller stopped before its comparison was written')
    time.sleep(10)
rows=json.loads((O/'pilot_comparison.json').read_text())
means={mode:sum(r['selection_score'] for r in rows if r['mode']==mode)/2 for mode in ['joint_legacy','partial_blood','blood_marginal']}
selected=min(['partial_blood','blood_marginal'],key=means.get)
decision=dict(selected_dsm_mode=selected,criterion='minimum two-seed mean of best internal validation score in equal 20-epoch budget',
    mean_validation_scores=means,comparison=rows,
    explanation='Joint legacy derivative is mathematically inconsistent with a blood partial score and is retained only as a reference.',
    downstream_rules='identical to rescreen_20261002; no A-tier, GWAS or external brain evidence used to select the corrected DSM implementation',
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in W.glob('*.py')})
p=O/'training_fix_decision.json'
if p.exists():
    saved=json.loads(p.read_text())
    assert {k:v for k,v in saved.items() if k!='source_sha256'}=={k:v for k,v in decision.items() if k!='source_sha256'}
    # Diagnostic helpers added later do not change the frozen recipe. Check
    # every scientific runtime module, accepting only the recorded input repair.
    amendment=O/'annotation_repair/protocol_amendment.json'
    repair=json.loads(amendment.read_text()) if amendment.exists() else {}
    runtime=['candidate_registry.py','ccwm.py','train_ccwm.py','training_utils.py','evaluate_models.py',
             'screen_candidates.py','run_full.py','run_discovery_evidence.py','run_ablations.py',
             'dose_programs.py','randko.py','rank_candidates.py','brain_competitive.py',
             'independent_validation.py','spatial_analysis.py','input_utils.py',
             'program_effect_sizes.py','validate_outputs.py']
    for name in runtime:
        expected=saved['source_sha256'][name]
        if repair.get('previous_source_sha256',{}).get(name)==expected:
            expected=repair['changed_source_sha256'][name]
        assert hashlib.sha256((W/name).read_bytes()).hexdigest()==expected,('frozen runtime changed',name)
else:p.write_text(json.dumps(decision,indent=2))
execute('candidate_registry.py')
print('FROZEN',selected,means,flush=True)
with ThreadPoolExecutor(3) as pool:
    main=pool.submit(execute,'run_full.py',['--beta','5','--sample-latent-a','--sinkhorn-iters','20','--dsm-mode',selected])
    while not (O/'status.json').exists():
        if main.done():main.result()
        time.sleep(1)
    evidence=pool.submit(execute,'run_discovery_evidence.py')
    controls=pool.submit(execute,'run_ablations.py')
    for task in [main,evidence,controls]:task.result()
for script in ['blood_covariate_sensitivity.py','compare_ablations.py','review_results.py','priority_identity_qc.py']:
    execute(script)
(O/'optimization_completed.json').write_text(json.dumps(dict(status='complete',selected_dsm_mode=selected,finished_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())),indent=2))
print('OPTIMIZATION COMPLETE',flush=True)
