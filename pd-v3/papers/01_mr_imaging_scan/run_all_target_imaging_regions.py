"""Acquire every planned BIG40 regional source with bounded native downloads."""
from pathlib import Path
import concurrent.futures,fcntl,json,os,subprocess,sys,threading,time
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(Path(__file__).parent))
import native_region_io as native
BASE=ROOT/"data/all_target_followup_20261118/regional_all/imaging"
# Oxford shows the same refusal pattern as EBI under bursty launches; space them out.
STAGGER_SECONDS=0.25
SOURCE_WORKERS=24
def main():
 with (BASE/".big40_controller.lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  path=BASE/"BIG40_CONTROLLER.json"
  if path.exists():raise RuntimeError("Inspect previous controller before explicit recovery")
  for name in ["NATIVE_EXTRACTION_PILOT_VERIFICATION.json","OXFORD_RANGE_TRANSPORT_PROBE.json"]:
   assert json.loads((BASE/name).read_text())["complete"]
  plan=json.loads((BASE/"PLAN.json").read_text())
  order=[plan["phenotypes"][p]["idp"] for p in plan["phenotype_order"] if plan["phenotypes"][p]["family"]=="BIG40"]
  script=Path(__file__).with_name("acquire_fast_imaging_regions.py")
  dependencies=[Path(__file__),script,script.with_name("acquire_all_target_imaging_regions.py"),script.with_name("fast_range_common.py"),script.with_name("fast_range")/"pfz.py",script.with_name("fast_range")/"imgfetch.py",script.with_name("fast_range")/"libpfz.so",Path(native.__file__),script.with_name("s09_acquire_qsm.py"),
   script.with_name("s07_acquire_big40.py"),BASE/"PLAN.json",BASE/"variant_reference/MANIFEST.json"]
  hashes={str(p.relative_to(ROOT)):native.sha(p) for p in dependencies}
  state=dict(pid=os.getpid(),status="running",complete=False,total_native_sources=len(order),source_workers=SOURCE_WORKERS,connections_per_source=4,
   started_at_unix=time.time(),completed=[],failures={},running={},remaining_unscheduled=len(order),consecutive_failures=0,inputs_sha256=hashes)
  guard=threading.Lock();native.dump(path,state)
  def one(idp):
   for p,h in hashes.items():assert native.sha(ROOT/p)==h,"Frozen imaging dependency changed"
   folder=BASE/"BIG40"/idp;folder.mkdir(parents=True,exist_ok=True)
   env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
   with (folder/"worker.log").open("ab") as log:
    p=subprocess.Popen(["/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python",str(script),"--idp",idp,"--connections","4"],
     cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env=env)
    with guard:state["running"][idp]=p.pid;native.dump(path,state)
    code=p.wait()
   if code:raise RuntimeError("Native source process exited "+str(code))
   d=json.loads((folder/"COMPLETE.json").read_text());assert d["complete"]
   assert d["plan_sha256"]==hashes[str((BASE/"PLAN.json").relative_to(ROOT))]
   return idp
  submitted=0;stop=False
  with concurrent.futures.ThreadPoolExecutor(max_workers=SOURCE_WORKERS) as pool:
   active={}
   while active or (submitted<len(order) and not stop):
    while len(active)<SOURCE_WORKERS and submitted<len(order) and not stop:
     idp=order[submitted];active[pool.submit(one,idp)]=idp;submitted+=1
     if len(active)<SOURCE_WORKERS and submitted<len(order) and not (BASE/"BIG40"/idp/"COMPLETE.json").exists():time.sleep(STAGGER_SECONDS)
    finished,_=concurrent.futures.wait(active,timeout=30,return_when=concurrent.futures.FIRST_COMPLETED)
    with guard:
     for f in finished:
      idp=active.pop(f);state["running"].pop(idp,None)
      try:state["completed"].append(f.result());state["consecutive_failures"]=0
      except Exception as exc:
       state["failures"][idp]=repr(exc);state["consecutive_failures"]+=1
       if state["consecutive_failures"]>=4:stop=True
      print(json.dumps({"idp":idp,"complete_sources":len(state["completed"]),"failed_sources":len(state["failures"])}),flush=True)
     state.update(updated_at_unix=time.time(),remaining_unscheduled=len(order)-submitted);native.dump(path,state)
  state.update(status="needs_attention_source_errors" if stop else ("complete" if not state["failures"] else "finished_with_source_failures"),
   complete=not state["failures"] and submitted==len(order),finished_at_unix=time.time())
  native.dump(path,state)
if __name__=="__main__":main()
