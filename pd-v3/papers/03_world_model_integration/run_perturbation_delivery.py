"""Run the frozen perturbation extraction and verification after the live acquisition completes."""
from pathlib import Path
import json,os,time,subprocess,traceback
from datetime import datetime,timezone
R=Path(__file__).resolve().parents[2];O=R/'papers/03_world_model_integration/results/perturbation_readout_20261117';D=R/'data/perturbation_sources_20261117/GSE311359'
def now():return datetime.now(timezone.utc).isoformat()
def write(j):
 p=O/'DELIVERY_RUN.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def main():
 if (O/'DELIVERY_RUN.json').exists():raise RuntimeError('Existing controller')
 run=dict(status='waiting_for_live_acquisition',pid=os.getpid(),complete=False,started_at_utc=now(),stages=[]);write(run)
 env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 try:
  while True:
   a=json.loads((D/'COUNT_ACQUISITION.json').read_text())
   if a['complete']:break
   if a['status']!='running':raise RuntimeError('Acquisition terminal without completion: '+str(a.get('error')))
   pid=a['pid'];proc=Path('/proc')/str(pid)
   if not proc.exists() or (proc/'stat').read_text().split()[2]=='Z':raise RuntimeError('Acquisition process no longer live')
   cmd=(proc/'cmdline').read_bytes();assert b'acquire_perturbation_ranges.py' in cmd
   run.update(acquisition_pid=pid,acquisition_chunks=a.get('chunks_complete',0),last_live_check_at_utc=now());write(run);time.sleep(10)
  for script,receipt in [('extract_perturbation_readouts.py','RUN.json'),('verify_perturbation_readouts.py','INDEPENDENT_READOUT_VERIFICATION.json'),('compare_perturbation_directions.py','COMPARISON_RUN.json')]:
   run.update(status='running',stage=script);write(run)
   with (O/(script.removesuffix('.py')+'.log')).open('ab',buffering=0) as log:
    job=subprocess.Popen(['/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python','-u',str(R/'papers/03_world_model_integration'/script)],cwd=str(R),stdout=log,stderr=subprocess.STDOUT,env=env)
    run['child_pid']=job.pid;write(run);code=job.wait()
   assert code==0,(script,code)
   status=json.loads((O/receipt).read_text());assert status['complete']
   run['stages'].append(dict(script=script,receipt=receipt,exit_code=code,finished_at_utc=now()));write(run)
  run.update(status='complete',complete=True,finished_at_utc=now(),model_validation_complete=False)
 except Exception as e:
  run.update(status='needs_attention',error=repr(e),finished_at_utc=now());traceback.print_exc();raise
 finally:write(run)
if __name__=='__main__':main()
