"""Reconstruct a seed42 execution candidate and compare archived numeric values."""
import os
for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[name]='2'
import argparse,subprocess,sys,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import torch
import recompute_cellot as r
G=r.b.O/'cellot_seed42_cpu_direct_recomputed_20261004'
AXES=[a for a in r.b.AXES if not a.startswith('pDC__')]

def run(axis,gpu):
 with (G/(axis+'.log')).open('w') as f:subprocess.run([sys.executable,str(Path(__file__).resolve()),'--axis',axis,'--gpu',str(gpu)],stdout=f,stderr=subprocess.STDOUT,check=True)
 return axis
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--axis');ap.add_argument('--gpu',type=int,default=0);args=ap.parse_args();torch.set_num_threads(2)
 G.mkdir(exist_ok=True)
 if args.axis:
  r.DEVICE='cpu';r.MODEL_DIR=G/'models/seed_42';r.MODEL_DIR.mkdir(parents=True,exist_ok=True)
  r.b.RESULT_DIR=G/'results/seed_42';r.b.RESULT_DIR.mkdir(parents=True,exist_ok=True)
  r.cellot(42,args.axis)
 else:
  done=[]
  with ThreadPoolExecutor(max_workers=6) as pool:
   for f in as_completed([pool.submit(run,a,j%2) for j,a in enumerate(AXES)]):
    done.append(f.result());print('CPU_DIRECT_DONE',done[-1],flush=True)
  (G/'run_status.json').write_text(json.dumps(dict(complete=True,axes=done)))
