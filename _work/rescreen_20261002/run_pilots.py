import subprocess,os,json,sys,time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
W=Path(__file__).parent;A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/rescreen_20261002'
def run(beta,gpu):
 out=O/f'pilot_beta{beta}';out.mkdir(exist_ok=True)
 cmd=[sys.executable,'-u',str(W/'train_ccwm.py'),'--data-dir',str(A/'interim/rescreen_20261002'),'--out-dir',str(out),'--run-dir',str(out/'run'),'--run-id',f'rescreen-pilot-beta{beta}','--potential-beta',str(beta),'--seed','42','--epochs','25','--steps-per-epoch','200']
 env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
 with (out/'console.log').open('w') as f:subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 return json.loads((out/'completed.json').read_text())|{'beta':beta}
with ThreadPoolExecutor(2) as pool:
 results=list(pool.map(lambda pair:run(*pair),[(1,0),(5,1)]))
(O/'pilot_comparison.json').write_text(json.dumps(results,indent=2))
(O/'status.json').write_text(json.dumps(dict(stage='pilots_complete',next='inspect_validation_and_freeze_recipe',results=results),indent=2))
print(json.dumps(results),flush=True)
