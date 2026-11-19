"""Run original-plan regional ABFs as each verified protein/image source becomes available."""
from pathlib import Path
import concurrent.futures,csv,fcntl,json,os,subprocess,time
import run_one_all_target_coloc as worker
ROOT=worker.ROOT;BASE=worker.BASE;OUT=worker.OUT
def main():
 OUT.mkdir(exist_ok=True)
 with (OUT/".controller.lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  status=OUT/"CONTROLLER.json"
  if status.exists():raise RuntimeError("Inspect previous controller before explicit recovery")
  verification=OUT/"RECIPE_VERIFICATION.json";assert json.loads(verification.read_text())["complete"]
  plan=json.loads((BASE/"PLAN.json").read_text());ip=json.loads((BASE/"imaging/PLAN.json").read_text())
  order=[rid for acc in plan["source_order"] for rid in plan["sources"][acc]["regions"]]
  assert len(order)==1375 and len(set(order))==1375
  script=Path(worker.__file__)
  deps=[Path(__file__),script,script.with_name("all_target_coloc.R"),script.with_name("prepare_new_scan_coloc.py"),
   ROOT/"reference/coloc/claudia.R",BASE/"PLAN.json",BASE/"imaging/PLAN.json",verification]
  hashes={str(p.relative_to(ROOT)):worker.sha(p) for p in deps}
  state=dict(pid=os.getpid(),status="running",complete=False,total_regions=len(order),process_workers=12,
   completed=[],failures={},running={},started_at_unix=time.time(),inputs_sha256=hashes,consecutive_failures=0)
  worker.dump(status,state);done=set();active={};stop_reason=None
  env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
  def one(rid):
   folder=OUT/rid;folder.mkdir(exist_ok=True)
   with (folder/"worker.log").open("a") as log:
    p=subprocess.Popen(["/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python",str(script),"--region",rid],
     cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    code=p.wait()
   if code:raise RuntimeError("Regional worker failed: exit="+str(code)+", pid="+str(p.pid))
   d=json.loads((folder/"COMPLETE.json").read_text());assert d["complete"]
   return d["window_status"]
  with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
   while True:
    for fut in list(active):
     if not fut.done():continue
     rid=active.pop(fut);state["running"].pop(rid,None);done.add(rid)
     try:
      statuses=fut.result();state["completed"].append(rid);state["consecutive_failures"]=0
      print(json.dumps(dict(region=rid,window_status=statuses,completed=len(state["completed"]))),flush=True)
     except Exception as exc:
      state["failures"][rid]=repr(exc);state["consecutive_failures"]+=1
      if state["consecutive_failures"]>=4:stop_reason="needs_attention_region_errors"
    coverage=[]
    for rid in order:
     rec=plan["regions"][rid];ph=ip["phenotypes"][rec["phenotype"]]
     pc=BASE/"sources"/rec["accession"]/"COMPLETE.json"
     ic=BASE/"imaging"/ph["family"]/ph["idp"]/"COMPLETE.json"
     ready=pc.exists() and ic.exists()
     if rid in state["completed"]:stage="complete"
     elif rid in state["failures"]:stage="failed"
     elif rid in active.values():stage="running"
     elif ready:stage="ready"
     elif not pc.exists() and not ic.exists():stage="waiting_protein_and_imaging"
     elif not pc.exists():stage="waiting_protein"
     else:stage="waiting_imaging"
     coverage.append(dict(region_id=rid,accession=rec["accession"],gene=rec["gene"],phenotype=rec["phenotype"],
      model_linked=rec["model_linked"],status=stage,reason=state["failures"].get(rid,"")))
     if stage=="ready" and len(active)<12 and not stop_reason:
      for p,h in hashes.items():assert worker.sha(ROOT/p)==h,"Frozen regional analysis dependency changed"
      fut=pool.submit(one,rid);active[fut]=rid;state["running"][rid]="submitted";coverage[-1]["status"]="running"
    temp=OUT/"coverage.csv.pending"
    with temp.open("w",newline="") as f:
     w=csv.DictWriter(f,fieldnames=list(coverage[0]));w.writeheader();w.writerows(coverage)
    temp.replace(OUT/"coverage.csv")
    counts={}
    for row in coverage:counts[row["status"]]=counts.get(row["status"],0)+1
    state.update(coverage_counts=counts,updated_at_unix=time.time());worker.dump(status,state)
    if len(done)==len(order):break
    if stop_reason and not active:break
    if not active and not any(v["status"]=="ready" for v in coverage):
     upstreams=[BASE/"CONTROLLER.json",BASE/"imaging/BIG40_CONTROLLER.json",BASE/"imaging/OxfordQSM2022/RUN.json"]
     upstream=[];alive=False;stale=[]
     for p in upstreams:
      if not p.exists():upstream.append({"path":str(p.relative_to(ROOT)),"status":"absent"});continue
      d=json.loads(p.read_text());s=d["status"];running=s in ["running","extracting","checking_archive_identity"]
      proc=Path("/proc")/str(d["pid"])/"cmdline";live=proc.exists() and bool(proc.read_bytes())
      upstream.append(dict(path=str(p.relative_to(ROOT)),status=s,pid=d["pid"],live=live))
      alive=alive or (running and live)
      if running and not live:stale.append(str(p.relative_to(ROOT)))
     state["upstream_status"]=upstream
     if stale:stop_reason="needs_attention_upstream_not_live"
     elif not alive:stop_reason="finished_with_unavailable_sources"
     if stop_reason:break
    time.sleep(2)
  state.update(complete=len(state["completed"])==len(order) and not state["failures"],
   status=stop_reason or ("complete" if not state["failures"] else "finished_with_region_failures"),finished_at_unix=time.time())
  worker.dump(status,state)
if __name__=="__main__":main()
