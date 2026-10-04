"""Recompute missing CellOT artifacts with the recovered original numerical recipe."""
import concurrent.futures,json,subprocess,sys,shutil,hashlib
from pathlib import Path
W=Path(__file__).resolve().parent
R=Path('/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002/cellot_recomputed_20261003')
C=R.with_name('cellot_seed42_cpu_direct_recomputed_20261004')
def run(args):
 subprocess.run([sys.executable]+[str(W/args[0])]+args[1:],check=True)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if __name__=='__main__':
 if not (R/'projection_reconstruction.json').exists():run(['prepare_projection.py'])
 # Worker route: seeds43/44 and pDC seed42. Direct route: original six seed42 axes.
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  pending=[pool.submit(run,['recompute_cellot.py','--workers','18']),pool.submit(run,['cpu_seed42_direct.py'])]
  for job in pending:job.result()
 for group in ['models/seed_42','results/seed_42']:
  for src in (C/group).glob('cellot*'):
   dst=R/group/src.name;dst.parent.mkdir(parents=True,exist_ok=True)
   if dst.exists() and digest(dst)!=digest(src):
    backup=R/'before_verified_recipe'/group/src.name;backup.parent.mkdir(parents=True,exist_ok=True)
    if not backup.exists():shutil.copy2(dst,backup)
   shutil.copy2(src,dst)
 if not (R/'execution_modes.json').exists():
  (R/'execution_modes.json').write_text(json.dumps(dict(seed42_original_six_axes='CPU; torch=2, BLAS=2',all_other_axis_seed_jobs='CPU; torch=1, BLAS=1',PCA='Original six axes=2 BLAS threads; pDC=1',source='Recovered original recipe, confirmed against archived numeric outputs'),indent=2))
 run(['capture_recovery_environment.py']);run(['complete_readouts.py'])
 check=json.loads((R/'recovery_validation.json').read_text())
 assert check['status']=='passed','Archived numeric checks did not pass; no canonical files installed'
 print('Verified CellOT reconstruction complete. Install with install_recovery.py.',flush=True)
