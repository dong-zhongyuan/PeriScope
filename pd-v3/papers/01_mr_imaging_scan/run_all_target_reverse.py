"""Run every original planned reverse region after its exact source and reference are verified."""
from pathlib import Path
import concurrent.futures,csv,fcntl,hashlib,json,os,subprocess,time
import reverse_followup_inputs as routing
ROOT=routing.ROOT;BASE=routing.BASE;SUP=routing.SUP;OUT=ROOT/"papers/01_mr_imaging_scan/results/all_target_reverse_20261118"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):
 t=p.with_name(p.name+".pending");t.write_text(json.dumps(d,indent=2)+"\n");t.replace(p)
def main():
 OUT.mkdir(exist_ok=True)
 with (OUT/".controller.lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (OUT/"CONTROLLER.json").exists()
  plan=json.loads((BASE/"PLAN.json").read_text());supp=json.loads((SUP/"PLAN.json").read_text())
  order=[rid for acc in plan["source_order"] for rid in plan["sources"][acc]["regions"]];assert len(order)==1375
  recipe=OUT/"RECIPE_VERIFICATION.json";assert json.loads(recipe.read_text())["complete"]
  script=Path(__file__).with_name("run_one_all_target_reverse.py")
  deps=[Path(__file__),script,Path(routing.__file__),Path(routing.core.__file__),ROOT/"pdv3/mr.py",BASE/"PLAN.json",SUP/"PLAN.json",recipe]
  hashes={str(p.relative_to(ROOT)):sha(p) for p in deps}
  needs={rid:routing.requirements(plan["regions"][rid],plan,supp)[-1] for rid in order}
  counts={rid:len(routing.source_peaks(plan["regions"][rid],plan,supp)) for rid in order}
  state=dict(pid=os.getpid(),status="running",complete=False,total_regions=1375,workers=4,completed=[],failures={},running={},
   started_at_unix=time.time(),inputs_sha256=hashes,consecutive_failures=0)
  target=OUT/"CONTROLLER.json";dump(target,state);active={};done=set();stop=None
  env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",NUMEXPR_NUM_THREADS="1")
  def one(rid):
   folder=OUT/rid;folder.mkdir(exist_ok=True)
   with (folder/"worker.log").open("a") as log:
    p=subprocess.Popen(["/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python",str(script),"--region",rid],
     cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    code=p.wait()
   if code:raise RuntimeError("Reverse region exit="+str(code)+", pid="+str(p.pid))
   d=json.loads((folder/"COMPLETE.json").read_text());assert d["complete"];return d["scientific_status_counts"]
  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
   while True:
    for f in list(active):
     if not f.done():continue
     rid=active.pop(f);done.add(rid);state["running"].pop(rid,None)
     try:
      statuses=f.result();state["completed"].append(rid);state["consecutive_failures"]=0
      print(json.dumps(dict(region=rid,statuses=statuses,completed=len(state["completed"]))),flush=True)
     except Exception as exc:
      state["failures"][rid]=repr(exc);state["consecutive_failures"]+=1
      if state["consecutive_failures"]>=4:stop="needs_attention_region_errors"
    coverage=[]
    for rid in order:
     rec=plan["regions"][rid];missing=[str(p.relative_to(ROOT)) for p in needs[rid] if not p.exists()]
     if rid in state["completed"]:status="complete"
     elif rid in state["failures"]:status="failed"
     elif rid in active.values():status="running"
     elif missing:status="waiting_inputs"
     else:status="ready"
     if status=="ready" and len(active)<4 and not stop:
      for p,h in hashes.items():assert sha(ROOT/p)==h,"Frozen reverse dependency changed"
      active[pool.submit(one,rid)]=rid;state["running"][rid]="submitted";status="running"
     coverage.append(dict(region_id=rid,accession=rec["accession"],gene=rec["gene"],phenotype=rec["phenotype"],
      source_peaks=counts[rid],status=status,missing_inputs=json.dumps(missing),error=state["failures"].get(rid,"")))
    temp=OUT/"coverage.csv.pending"
    with temp.open("w",newline="") as f:
     w=csv.DictWriter(f,fieldnames=list(coverage[0]));w.writeheader();w.writerows(coverage)
    temp.replace(OUT/"coverage.csv");cc={}
    for row in coverage:cc[row["status"]]=cc.get(row["status"],0)+1
    state.update(coverage_counts=cc,updated_at_unix=time.time());dump(target,state)
    if len(done)==len(order):break
    if stop and not active:break
    if not active and not cc.get("ready"):
     live=False;upstates=[]
     for p in [BASE/"CONTROLLER.json",SUP/"protein_outcomes/CONTROLLER.json",BASE/"reference_verification/CONTROLLER.json"]:
      u=json.loads(p.read_text());proc=Path("/proc")/str(u["pid"])/"cmdline"
      live=live or (u["status"]=="running" and proc.exists() and bool(proc.read_bytes()))
      upstates.append(dict(path=str(p.relative_to(ROOT)),status=u["status"]))
     state["upstream_status"]=upstates
     if not live:stop="finished_with_unavailable_inputs";break
    time.sleep(5 if active else 15)
  state.update(status=stop or ("complete" if not state["failures"] else "finished_with_region_errors"),
   complete=len(state["completed"])==len(order) and not state["failures"],finished_at_unix=time.time());dump(target,state)
if __name__=="__main__":main()
