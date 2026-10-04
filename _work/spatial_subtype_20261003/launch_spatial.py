"""Run a fixed QC and all ten independent spatial donors; checkpoint each."""
from pathlib import Path
import json,subprocess,os,time,shutil
from concurrent.futures import ThreadPoolExecutor
W=Path(__file__).parent
O=Path(os.environ.get('SPATIAL_OUTPUT_DIR','/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/spatial_subtype_20261003'))
R=os.environ.get('SPATIAL_RSCRIPT','/public/home/mengxl/dzy/envs/rcoloc/bin/Rscript')
if not Path(R).exists():R=shutil.which('Rscript') or R
env=dict(os.environ,LC_ALL='C',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',R_LIBS_USER=str(W/'R_library'))
env['SPATIAL_WORK_DIR']=str(W)
def run(sample):
    if (O/'samples'/sample/'completed.txt').exists():return dict(sample=sample,returncode=0,reused=True)
    (O/'logs').mkdir(exist_ok=True)
    with (O/'logs'/(sample+'.log')).open('w') as log:
        r=subprocess.run([R,str(W/'run_spatial_rctd.R'),sample],env=env,stdout=log,stderr=subprocess.STDOUT)
    return dict(sample=sample,returncode=r.returncode)
samples=['QC']+[f'GSM803136{i}' for i in range(10)]
(O/'run_status.json').write_text(json.dumps(dict(stage='running',samples=samples,started=time.strftime('%Y-%m-%d %H:%M:%S')),indent=2))
results=[]
with ThreadPoolExecutor(max_workers=int(os.environ.get('SPATIAL_PARALLEL_SAMPLES','3'))) as pool:
    for result in pool.map(run,samples):
        results.append(result);print(result,flush=True)
        (O/'run_status.json').write_text(json.dumps(dict(stage='running',completed=results),indent=2))
(O/'run_status.json').write_text(json.dumps(dict(stage='complete' if all(x['returncode']==0 for x in results) else 'failed',completed=results),indent=2))
if any(x['returncode']!=0 for x in results):raise SystemExit(1)
