"""Matched short-budget DSM diagnostics; select on internal validation only."""
from pathlib import Path
import sys,os,json,subprocess
from concurrent.futures import ThreadPoolExecutor
W=Path(__file__).parent;A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/optimization_20261003'
protocol=dict(status='pilot_training',seeds=[42,43],modes=['partial_blood','blood_marginal'],epochs=20,steps_per_epoch=200,
    selection='Internal validation score only; compare equal-budget baseline trajectory through epoch20. No A-tier counts, GWAS, spatial or locked-test outcomes used for pilot selection.',
    basis='Fix diagonal vs partial DSM derivative; compare separate marginal-blood density learning as an ablation.',
    unchanged='Source matrices, splits, genes, doses, neural widths, optimizer, loss weights, beta5, posterior initialization and 20 OT iterations.')
assert json.loads((O/'dsm_derivative_checks.json').read_text())['status']=='passed'
(O/'pilot_protocol.json').write_text(json.dumps(protocol,indent=2))
def run(mode,seed,gpu):
    p=O/'pilots'/mode/f'seed_{seed}';p.mkdir(parents=True,exist_ok=True)
    env=os.environ|dict(CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',PYTHONPATH='/public/home/mengxl/dzy/pd_product/src')
    args=[sys.executable,'-u',str(W/'train_ccwm.py'),'--data-dir',str(A/'interim/rescreen_20261002'),'--out-dir',str(p),'--run-dir',str(p/'run'),'--run-id',f'dsm-{mode}-s{seed}','--dsm-mode',mode,'--seed',str(seed),'--epochs','20','--steps-per-epoch','200','--sample-latent-a','--sinkhorn-iters','20','--potential-beta','5']
    if not (p/'completed.json').exists():
        with (p/'console.log').open('w') as f:subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
    print('PILOT COMPLETE',mode,seed,flush=True)
with ThreadPoolExecutor(4) as pool:
    fs=[pool.submit(run,mode,seed,si) for si,seed in enumerate([42,43]) for mode in protocol['modes']]
    for f in fs:f.result()
rows=[]
for seed in protocol['seeds']:
    for mode in ['joint_legacy']+protocol['modes']:
        p=A/f'results/rescreen_20261002/models/seed_{seed}' if mode=='joint_legacy' else O/'pilots'/mode/f'seed_{seed}'
        vals=[r for r in json.loads((p/'validation.json').read_text()) if r['epoch']<=20]
        best=min(vals,key=lambda r:r['selection_score'])
        rows.append(dict(mode=mode,seed=seed,**{k:best[k] for k in ['epoch','selection_score','protein_mse','brain_mmd']}))
(O/'pilot_comparison.json').write_text(json.dumps(rows,indent=2));print(json.dumps(rows),flush=True)
