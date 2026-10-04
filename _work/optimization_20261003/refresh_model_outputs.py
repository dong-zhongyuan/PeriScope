"""Replace repaired seeds in place and refresh their dependent model outputs."""
from pathlib import Path
import os, sys, json, shutil, subprocess, time, fcntl, traceback
from concurrent.futures import ThreadPoolExecutor
from model_registry import ROOT as O, sha256, atomic_json, curve_provenance, stamp_curve

W=Path(__file__).resolve().parent
B=O.parent/'benchmark_20261002'
STATUS=O/'model_refresh_status.json'

def prepare():
    overrides=json.loads((B/'checkpoint_overrides.json').read_text())['checkpoints']
    assert set(overrides)=={'47','48','49'}
    previous=json.loads((O/'model_registry.json').read_text()) if (O/'model_registry.json').exists() else None
    registry={'version':'repaired_ten_seed_20261004','seeds':list(range(42,52)),
              'replaced_seeds':[47,48,49],'discovery_seeds':list(range(42,47)),
              'confirmation_seeds':list(range(47,52)), 'models':{},
              'discovery_files':{str(p.relative_to(O)):sha256(p) for p in [
                  O/'programs.json', O/'target_specific/programs.json',O/'joint_analysis/programs.json',
                  O/'program_folds.json',O/'target_specific/program_folds.json']}}
    for seed in range(42,52):
        dest=O/f'models/seed_{seed}'
        if str(seed) in overrides:
            src=Path(overrides[str(seed)]).parent
            for name in ['model.pt','last_model.pt','train_log.csv','training_inputs.json','validation.json','completed.json']:
                p=src/name
                assert p.exists(),p
                if not (dest/name).exists() or sha256(p)!=sha256(dest/name):
                    tmp=dest/(name+'.pending');shutil.copy2(p,tmp);tmp.replace(dest/name)
            for marker in ['screened.json','evaluated.json']:
                if previous is None or previous['models'].get(str(seed),{}).get('sha256')!=sha256(dest/'model.pt'):
                    (dest/marker).unlink(missing_ok=True)
        registry['models'][str(seed)]={'path':str(dest/'model.pt'),'sha256':sha256(dest/'model.pt')}
    atomic_json(O/'model_registry.json',registry)
    # Only the seven unchanged seeds may inherit existing, hash-verified outputs.
    old_manifest=json.loads((O/'server_binary_manifest.json').read_text())
    known={r['relative_path']:r['sha256'] for r in old_manifest['artifacts']}
    for seed in [42,43,44,45,46,50,51]:
        for p in (O/'curves').glob(f'*__seed{seed}.npz'):
            assert known[str(p.relative_to(O))]==sha256(p),f'Unexpected pre-existing curve change: {p}'
            stamp_curve(p,curve_provenance(O,seed))
    pointer=json.loads((O/'active_nomination.json').read_text())
    pointer['status']='refreshing_model_dependent_results';pointer['model_registry']='model_registry.json'
    atomic_json(O/'active_nomination.json',pointer)
    return registry

def worker(seed,gpu):
    env=os.environ|{'CUDA_VISIBLE_DEVICES':str(gpu),'OMP_NUM_THREADS':'2','MKL_NUM_THREADS':'2',
                       'OPENBLAS_NUM_THREADS':'2','PYTHONPATH':'/public/home/mengxl/dzy/pd_product/src'}
    for script in ['evaluate_models.py','screen_candidates.py']:
        if script=='evaluate_models.py':
            receipt=O/f'evaluation/seed{seed}_provenance.json'
            if receipt.exists():
                r=json.loads(receipt.read_text())
                if r['checkpoint_sha256']==sha256(O/f'models/seed_{seed}/model.pt') and all(
                    (O/'evaluation'/n).exists() and sha256(O/'evaluation'/n)==h for n,h in r['outputs'].items()):
                    print('VERIFIED_CURRENT',seed,script,flush=True)
                    continue
        log=O/f'models/seed_{seed}/refresh_{script[:-3]}.log'
        with log.open('w') as f:
            subprocess.run([sys.executable,'-u',str(W/script),'--seed',str(seed)],env=env,
                           stdout=f,stderr=subprocess.STDOUT,check=True)
        print('UPDATED',seed,script,flush=True)
    for name in ['screened.json','evaluated.json']:
        atomic_json(O/f'models/seed_{seed}'/name,{'seed':seed,'completed':True,
            'checkpoint_sha256':sha256(O/f'models/seed_{seed}/model.pt')})

def main():
    lock=(O/'model_refresh.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    try:
        registry=prepare()
        atomic_json(STATUS,{'stage':'refreshing_predictions_and_36_curves','started':time.time(),
                            'model_registry_sha256':sha256(O/'model_registry.json')})
        # One GPU job per device; seed49 follows seed47 on GPU0.
        def lane(gpu,seeds):
            for seed in seeds:worker(seed,gpu)
        with ThreadPoolExecutor(2) as pool:
            fs=[pool.submit(lane,0,[47,49]),pool.submit(lane,1,[48])]
            for f in fs:f.result()
        for p,h in registry['discovery_files'].items():assert sha256(O/p)==h
        from model_registry import valid_curve
        curves=list((O/'curves').glob('*__seed*.npz'))
        assert len(curves)==120
        for p in curves:
            seed=int(p.stem.rsplit('__seed',1)[1]);assert valid_curve(p,curve_provenance(O,seed)),p
        atomic_json(STATUS,{'stage':'model_outputs_complete','curves':120,'replaced_seeds':[47,48,49],
                            'model_registry_sha256':sha256(O/'model_registry.json'),'finished':time.time()})
    except Exception:
        atomic_json(STATUS,{'stage':'failed','traceback':traceback.format_exc()});raise

if __name__=='__main__':main()
