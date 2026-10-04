"""Matched seed-42 training controls, on GPU slots released by seeds 46 and 47."""
from pathlib import Path
import json,os,time,sys
from concurrent.futures import ThreadPoolExecutor
from run_full import command

W=Path(__file__).parent;A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/optimization_20261003'
def run(kind,gpu,wait_seed,extra):
    while not (O/f'models/seed_{wait_seed}/screened.json').exists():
        status=json.loads((O/'status.json').read_text())
        if status.get('stage')=='failed':raise RuntimeError('Main workflow failed; postpone controls')
        time.sleep(15)
    root=O/'ablations'/kind;target=root/'models/seed_42';target.mkdir(parents=True,exist_ok=True)
    for name in ['candidate_registry.json','dose_definitions.json','all_candidate_input_support.csv']:
        p=root/name
        if not p.exists():p.symlink_to(O/name)
    if not (target/'completed.json').exists():
        command('train_ccwm.py',['--data-dir',str(A/'interim/rescreen_20261002'),'--out-dir',str(target),
            '--run-dir',str(target/'run'),'--run-id','rescreen-ablation-'+kind,'--potential-beta','5','--sample-latent-a',
            '--dsm-mode',json.loads((O/'frozen_recipe.json').read_text())['dsm_mode'],'--sinkhorn-iters','20','--seed','42','--epochs','200','--steps-per-epoch','200']+extra,target/'console.log',gpu)
    if not (target/'evaluated.json').exists():
        command('evaluate_models.py',['--seed','42','--output-root',str(root)],target/'evaluation.log',gpu)
        (target/'evaluated.json').write_text(json.dumps(dict(completed=True)))
    if not (target/'screened.json').exists():
        command('screen_candidates.py',['--seed','42','--output-root',str(root)],target/'screen.log',gpu)
        (target/'screened.json').write_text(json.dumps(dict(completed=True)))
    print('CONTROL COMPLETE',kind,flush=True)

if __name__=='__main__':
    with ThreadPoolExecutor(2) as p:
        futures=[p.submit(run,'random_pairing',0,46,['--pairing','random']),p.submit(run,'coupling_off',1,47,['--coupling','off'])]
        for f in futures:f.result()
    (O/'ablations/completed.json').write_text(json.dumps(dict(status='complete',seed=42,controls=['random_pairing','coupling_off'])))
