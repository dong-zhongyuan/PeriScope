"""Regenerate independent evidence once discovery seeds finish, while confirmation models continue training."""
from pathlib import Path
import subprocess,sys,os,time,json
from concurrent.futures import ThreadPoolExecutor
W=Path(__file__).parent;A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/optimization_20261003'
while not all((O/f'models/seed_{s}/screened.json').exists() for s in range(42,47)):
    if json.loads((O/'status.json').read_text()).get('stage')=='failed':raise RuntimeError('Main run failed')
    time.sleep(15)
env=os.environ.copy();env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2',CCWM_SEEDS='42,43,44,45,46')
def run(name):
    with (O/('early_'+name.replace('.py','')+'.log')).open('w') as f:subprocess.run([sys.executable,'-u',str(W/name)],stdout=f,stderr=subprocess.STDOUT,env=env,check=True)
    print('DISCOVERY EVIDENCE COMPLETE',name,flush=True)
run('dose_programs.py')
with ThreadPoolExecutor(4) as pool:
    futures=[pool.submit(run,s) for s in ['brain_competitive.py','independent_validation.py','spatial_analysis.py']]
    for f in futures:f.result()
(O/'discovery_evidence_completed.json').write_text(json.dumps(dict(status='complete',discovery_seeds=list(range(42,47)),confirmation_models='not used for program definition')))
