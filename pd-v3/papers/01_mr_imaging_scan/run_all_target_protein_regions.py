"""Bounded source acquisition queue with an operational stop on repeated failures."""
from pathlib import Path
import concurrent.futures,fcntl,hashlib,json,os,subprocess,threading,time
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"data/all_target_followup_20261118/regional_all"
PY="/public/home/mengxl/.conda/envs/state_env/bin/python"
# ftp.ebi.ac.uk refuses bursts of connection attempts; space the source launches out.
STAGGER_SECONDS=0.25
# EBI refused connections at 12x4 plus the fetcher (19 failed sources, all ConnectionRefused).
# Range route: one keep-alive connection per source; full-file exceptions use 2.
SOURCE_WORKERS=14
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def dump(p,d):
 q=p.with_name(p.name+".pending");q.write_text(json.dumps(d,indent=2)+"\n");q.replace(p)
def main():
 with (BASE/".controller.lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (BASE/"CONTROLLER.json").exists():raise RuntimeError("Inspect existing controller before explicit recovery")
  plan=json.loads((BASE/"PLAN.json").read_text());script=Path(__file__).with_name("acquire_fast_protein_regions.py")
  deps=[Path(__file__),script,script.with_name("acquire_all_target_protein_regions.py"),script.with_name("fast_range_common.py"),script.with_name("fast_range")/"pfz.py",script.with_name("fast_range")/"libpfz.so",script.with_name("reconstruct_ages_formatted.py"),script.with_name("s09_acquire_qsm.py"),script.with_name("s07_acquire_big40.py"),BASE/"PLAN.json"]
  hashes={str(p.relative_to(ROOT)):sha(p) for p in deps}
  state=dict(pid=os.getpid(),status="running",complete=False,source_workers=SOURCE_WORKERS,connections_per_source=4,
   started_at_unix=time.time(),inputs_sha256=hashes,completed=[],failures={},running={},total_sources=len(plan["sources"]),
   consecutive_failures=0,remaining_unscheduled=len(plan["sources"]))
  guard=threading.Lock();dump(BASE/"CONTROLLER.json",state)
  def one(acc):
   for f,h in hashes.items():
    if sha(ROOT/f)!=h:raise RuntimeError("Frozen input changed: "+f)
   folder=BASE/"sources"/acc;folder.mkdir(parents=True,exist_ok=True)
   with (folder/"worker.log").open("ab") as log:
    child=subprocess.Popen([PY,str(script),"--accession",acc,"--connections","4"],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    with guard:
     state["running"][acc]=child.pid;dump(BASE/"CONTROLLER.json",state)
    code=child.wait()
   if code:raise RuntimeError("Source process failed with code "+str(code))
   d=json.loads((folder/"COMPLETE.json").read_text())
   assert d["complete"] and d["plan_sha256"]==hashes[str((BASE/"PLAN.json").relative_to(ROOT))]
   return acc
  order=iter(plan["source_order"]);submitted=0;stop_scheduling=False
  with concurrent.futures.ThreadPoolExecutor(max_workers=SOURCE_WORKERS) as pool:
   active={}
   while active or (submitted<len(plan["source_order"]) and not stop_scheduling):
    while len(active)<SOURCE_WORKERS and submitted<len(plan["source_order"]) and not stop_scheduling:
     acc=next(order);active[pool.submit(one,acc)]=acc;submitted+=1
     if len(active)<SOURCE_WORKERS and submitted<len(plan["source_order"]) and not (BASE/"sources"/acc/"COMPLETE.json").exists():time.sleep(STAGGER_SECONDS)
    finished,_=concurrent.futures.wait(active,timeout=30,return_when=concurrent.futures.FIRST_COMPLETED)
    with guard:
     for future in finished:
      acc=active.pop(future);state["running"].pop(acc,None)
      try:
       state["completed"].append(future.result());state["consecutive_failures"]=0
      except Exception as exc:
       state["failures"][acc]=repr(exc);state["consecutive_failures"]+=1
       if state["consecutive_failures"]>=4:stop_scheduling=True
      print(json.dumps({"finished":acc,"complete_sources":len(state["completed"]),"failed_sources":len(state["failures"])}),flush=True)
     state.update(updated_at_unix=time.time(),remaining_unscheduled=len(plan["source_order"])-submitted)
     dump(BASE/"CONTROLLER.json",state)
  state.update(status="needs_attention_source_errors" if stop_scheduling else ("complete" if not state["failures"] else "finished_with_source_failures"),
   complete=not bool(state["failures"]) and submitted==len(plan["source_order"]),finished_at_unix=time.time())
  dump(BASE/"CONTROLLER.json",state)
if __name__=="__main__":main()
