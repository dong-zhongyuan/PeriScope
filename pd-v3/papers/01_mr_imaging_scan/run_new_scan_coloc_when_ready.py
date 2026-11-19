"""Run the original regional continuation after all current-batch native sources complete."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,subprocess,time
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'papers/01_mr_imaging_scan/results/new_scan_regions_20261117';PY='/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python';RS='/public/home/mengxl/dzy/envs/rcoloc/bin/Rscript'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(j):
 p=OUT/'RUN.json';t=p.with_suffix('.json.tmp');t.write_text(json.dumps(j,indent=2)+'\n');t.replace(p)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 if (OUT/'RUN.json').exists():raise RuntimeError('Existing regional controller; inspect before restarting')
 scripts=['prepare_new_scan_coloc.py','new_scan_coloc.R','verify_new_scan_coloc.py','run_new_scan_coloc_when_ready.py'];sources={str((ROOT/'papers/01_mr_imaging_scan'/s).relative_to(ROOT)):sha(ROOT/'papers/01_mr_imaging_scan'/s) for s in scripts};sources['reference/coloc/claudia.R']=sha(ROOT/'reference/coloc/claudia.R');assert sources['reference/coloc/claudia.R']=='428f68c70d9351668c2eb61b0c870a31877dae731ff34b83d7c3096edd57b2cc'
 run=dict(status='waiting_for_native_imaging_regions',pid=os.getpid(),complete=False,started_at_utc=datetime.now(timezone.utc).isoformat(),source_code_sha256=sources,steps={});dump(run)
 try:
  while True:
   j=json.loads((ROOT/'data/new_scan_regional_20261117/IMAGING_RUN.json').read_text())
   if j['complete']:break
   if j.get('failures'):raise RuntimeError('Native source errors: '+repr(j['failures']))
   time.sleep(20)
  commands=[('prepare',[PY,str(ROOT/'papers/01_mr_imaging_scan/prepare_new_scan_coloc.py')]),('coloc',[RS,'--vanilla',str(ROOT/'papers/01_mr_imaging_scan/new_scan_coloc.R'),str(ROOT)]),('verify',[PY,str(ROOT/'papers/01_mr_imaging_scan/verify_new_scan_coloc.py')])]
  for name,cmd in commands:
   for f,h in sources.items():assert sha(ROOT/f)==h,f
   run['status']='running_'+name;dump(run)
   with (OUT/(name+'.log')).open('ab') as f:p=subprocess.run(cmd,cwd=ROOT,env=dict(os.environ,LC_ALL='C'),stdout=f,stderr=subprocess.STDOUT)
   run['steps'][name]=dict(exit_code=p.returncode,command=cmd);dump(run)
   if p.returncode:raise RuntimeError(name+' failed; inspect preserved log')
  run.update(status='complete',complete=True,finished_at_utc=datetime.now(timezone.utc).isoformat());dump(run)
  files=[p for p in OUT.rglob('*') if p.is_file() and p.name not in ['MANIFEST.json','controller.log','JOB.json']];(OUT/'MANIFEST.json').write_text(json.dumps(dict(complete=True,artifacts_sha256={str(p.relative_to(ROOT)):sha(p) for p in files}),indent=2)+'\n')
 except Exception as error:
  run.update(status='needs_attention',error=repr(error),finished_at_utc=datetime.now(timezone.utc).isoformat());dump(run);raise
if __name__=='__main__':main()
