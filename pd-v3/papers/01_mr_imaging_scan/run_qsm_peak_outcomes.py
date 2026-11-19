"""Complete QSM outcome supplements as each original source becomes available."""
from pathlib import Path
import concurrent.futures,fcntl,json,os,subprocess,time
import acquire_all_target_protein_regions as original
ROOT=original.ROOT;BASE=original.BASE;SUP=BASE/"reverse_supplement";OUT=SUP/"protein_outcomes"
def main():
 OUT.mkdir(exist_ok=True)
 with (OUT/".controller.lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (OUT/"CONTROLLER.json").exists()
  plan=json.loads((BASE/"PLAN.json").read_text());supp=json.loads((SUP/"PLAN.json").read_text())
  order=[a for a in plan["source_order"] if a in supp["sources"]];assert len(order)==22
  script=Path(__file__).with_name("extract_qsm_peak_outcomes.py")
  deps=[Path(__file__),script,Path(original.__file__),Path(original.bridge.__file__),BASE/"PLAN.json",SUP/"PLAN.json"]
  hashes={str(p.relative_to(ROOT)):original.sha(p) for p in deps}
  state=dict(pid=os.getpid(),status="running",complete=False,total_sources=len(order),completed=[],failures={},running={},
   inputs_sha256=hashes,started_at_unix=time.time())
  path=OUT/"CONTROLLER.json";original.dump(path,state);active={};done=set()
  def one(acc):
   folder=OUT/acc;folder.mkdir(exist_ok=True)
   with (folder/"worker.log").open("a") as log:
    p=subprocess.Popen(["/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python",str(script),"--accession",acc],
     cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    code=p.wait()
   if code:raise RuntimeError("Outcome extraction exit="+str(code)+", pid="+str(p.pid))
   d=json.loads((folder/"COMPLETE.json").read_text());assert d["complete"];return d["outcome_rows"]
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   while len(done)<len(order):
    for f in list(active):
     if not f.done():continue
     acc=active.pop(f);done.add(acc);state["running"].pop(acc,None)
     try:count=f.result();state["completed"].append(acc);print("QSM_SOURCE_COMPLETE",acc,count,flush=True)
     except Exception as exc:state["failures"][acc]=repr(exc)
    for acc in order:
     if len(active)>=2:break
     if acc in done or acc in active.values():continue
     if supp["sources"][acc]["peak_rows"] and not (BASE/"sources"/acc/"COMPLETE.json").exists():continue
     for p,h in hashes.items():assert original.sha(ROOT/p)==h,"Frozen dependency changed"
     active[pool.submit(one,acc)]=acc;state["running"][acc]="submitted"
    state["updated_at_unix"]=time.time();original.dump(path,state)
    if len(done)==len(order):break
    if not active:
     up=json.loads((BASE/"CONTROLLER.json").read_text());proc=Path("/proc")/str(up["pid"])/"cmdline"
     if up["status"]!="running" or not proc.exists() or not proc.read_bytes():
      state.update(status="needs_attention_upstream_unavailable",upstream_status=up["status"]);original.dump(path,state);return
    time.sleep(15)
  state.update(complete=not state["failures"],status="complete" if not state["failures"] else "finished_with_source_errors",finished_at_unix=time.time());original.dump(path,state)
if __name__=="__main__":main()
