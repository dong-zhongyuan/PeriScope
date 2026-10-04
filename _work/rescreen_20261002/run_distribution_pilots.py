from pathlib import Path
import subprocess,sys,os,json,time
from concurrent.futures import ThreadPoolExecutor
W=Path(__file__).parent;A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/rescreen_20261002'
while not (O/'pilot_comparison.json').exists():time.sleep(15)
# Check the numerical equivalence before using the faster OT solve.
precision=json.loads((O/'sinkhorn_precision.json').read_text());p20=[r for r in precision if r['iterations']==20]
assert max(r['cost_abs_error'] for r in p20)<1e-5 and max(r['gradient_relative_error'] for r in p20)<1e-3
beta=min(json.loads((O/'pilot_comparison.json').read_text()),key=lambda r:r['best_validation_score'])['beta']
(O/'status.json').write_text(json.dumps(dict(stage='distribution_training_check',beta=beta,variants=['mean_latent','posterior_latent'],sinkhorn_iterations=20)))
def run(posterior,gpu):
 name='posterior' if posterior else 'mean';out=O/('distribution_pilot_'+name);out.mkdir(exist_ok=True)
 cmd=[sys.executable,'-u',str(W/'train_ccwm.py'),'--data-dir',str(A/'interim/rescreen_20261002'),'--out-dir',str(out),'--run-dir',str(out/'run'),'--run-id','distribution-pilot-'+name,'--potential-beta',str(beta),'--sinkhorn-iters','20','--seed','42','--epochs','25','--steps-per-epoch','200']
 if posterior:cmd.append('--sample-latent-a')
 env=os.environ.copy();env.update(CUDA_VISIBLE_DEVICES=str(gpu),OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2',MKL_NUM_THREADS='2')
 with (out/'console.log').open('w') as f:subprocess.run(cmd,env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
 return json.loads((out/'completed.json').read_text())|{'beta':beta,'sample_latent_a':posterior}
with ThreadPoolExecutor(2) as pool:results=list(pool.map(lambda p:run(*p),[(False,0),(True,1)]))
(O/'distribution_pilot_comparison.json').write_text(json.dumps(results,indent=2))
(O/'status.json').write_text(json.dumps(dict(stage='distribution_pilots_complete',results=results),indent=2))
