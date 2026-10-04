"""Resumeable complete rerun after internal validation freezes the recipe."""

# Current authoritative pipeline dispatch (imports retain helper compatibility).
if __name__ == '__main__':
    from pathlib import Path as _RegistryPath
    import subprocess as _RegistrySubprocess, sys as _RegistrySys
    if _RegistryPath('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/model_registry.json').exists():
        _RegistrySubprocess.run([_RegistrySys.executable, str(_RegistryPath(__file__).with_name('refresh_model_outputs.py'))], check=True)
        _RegistrySubprocess.run([_RegistrySys.executable, str(_RegistryPath(__file__).with_name('refresh_current_pipeline.py'))], check=True)
        raise SystemExit(0)

import argparse,json,os,subprocess,sys,time,traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from rank_candidates import RULES

W=Path(__file__).parent;A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/optimization_20261003'

def command(script,args,log,gpu=None):
    env=os.environ.copy();env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
    if gpu is not None:env['CUDA_VISIBLE_DEVICES']=str(gpu)
    with log.open('a') as f:subprocess.run([sys.executable,'-u',str(W/script)]+args,stdout=f,stderr=subprocess.STDOUT,env=env,check=True)

def worker(gpu,seeds,beta,posterior,iters,dsm):
    for seed in seeds:
        target=O/f'models/seed_{seed}';target.mkdir(parents=True,exist_ok=True)
        if not (target/'completed.json').exists() and (target/'adopted_training.json').exists():
            pid=json.loads((target/'adopted_training.json').read_text())['pid']
            while Path(f'/proc/{pid}/cmdline').exists() and str(target).encode() in Path(f'/proc/{pid}/cmdline').read_bytes():time.sleep(10)
            if not (target/'completed.json').exists():raise RuntimeError(f'Adopted training failed for {seed}; inspect console.log')
        if not (target/'completed.json').exists():
            command('train_ccwm.py',['--data-dir',str(A/'interim/rescreen_20261002'),'--out-dir',str(target),
                '--run-dir',str(target/'run'/time.strftime('%Y%m%d_%H%M%S')),'--run-id',f'rescreen-seed-{seed}',
                '--dsm-mode',dsm,'--potential-beta',str(beta),'--sinkhorn-iters',str(iters),'--seed',str(seed),'--epochs','200','--steps-per-epoch','200']+(['--sample-latent-a'] if posterior else []),target/'console.log',gpu)
        if not (target/'evaluated.json').exists():
            command('evaluate_models.py',['--seed',str(seed)],target/'evaluation.log',gpu)
            (target/'evaluated.json').write_text(json.dumps(dict(seed=seed,completed=True)))
        if not (target/'screened.json').exists():
            command('screen_candidates.py',['--seed',str(seed)],target/'screen.log',gpu)
            (target/'screened.json').write_text(json.dumps(dict(seed=seed,completed=True)))
        print('SEED COMPLETE',seed,flush=True)

def main(beta,posterior,iters,dsm):
    assert json.loads((O/'data_checks.json').read_text())['status']=='passed'
    rule=O/'frozen_recipe.json'
    proposed=dict(dsm_mode=dsm,potential_beta=beta,seeds=list(range(42,52)),steps_max=40000,selection='internal donor validation only',
                  inputs=str(A/'interim/rescreen_20261002'),minimum_sampling_group=25,minimum_training_donors_per_stratum=2,screening_axes='all supported PD pairs',sample_latent_a=posterior,sinkhorn_iterations=iters)
    if rule.exists():assert json.loads(rule.read_text())==proposed
    else:rule.write_text(json.dumps(proposed,indent=2))
    (O/'ranking_rules.json').write_text(json.dumps(RULES,indent=2))
    (O/'status.json').write_text(json.dumps(dict(stage='ten_seed_training_and_full_screen',recipe=proposed),indent=2))
    try:
        with ThreadPoolExecutor(6) as pool:
            futures=[pool.submit(worker,slot%2,list(range(42+slot,52,6)),beta,posterior,iters,dsm) for slot in range(6)]
            for f in futures:f.result()
        for script in ['dose_programs.py','randko.py','brain_competitive.py','independent_validation.py','spatial_analysis.py','disease_association.py','rank_candidates.py']:
            (O/'status.json').write_text(json.dumps(dict(stage='evidence_recompute',running=script),indent=2))
            command(script,[],O/(script.replace('.py','')+'.log'))
        summary=json.loads((O/'rescreen_summary.json').read_text())
        (O/'status.json').write_text(json.dumps(dict(stage='analysis_complete',summary=summary),indent=2))
    except Exception:
        (O/'status.json').write_text(json.dumps(dict(stage='failed',traceback=traceback.format_exc()),indent=2));raise

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--dsm-mode',required=True,choices=['partial_blood','blood_marginal']);ap.add_argument('--beta',type=float,required=True);ap.add_argument('--sample-latent-a',action='store_true');ap.add_argument('--sinkhorn-iters',type=int,default=20);args=ap.parse_args();main(args.beta,args.sample_latent_a,args.sinkhorn_iters,args.dsm_mode)
