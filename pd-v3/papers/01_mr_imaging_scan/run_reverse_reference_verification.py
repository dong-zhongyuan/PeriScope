"""Verify each completed base or supplemental reference without restarting acquisitions."""
from pathlib import Path
import concurrent.futures,fcntl,hashlib,json,os,subprocess,time
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/"data/all_target_followup_20261118/regional_all"
OUT=BASE/"reference_verification"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):
 t=p.with_name(p.name+".pending");t.write_text(json.dumps(d,indent=2)+"\n");t.replace(p)
def main():
 OUT.mkdir(exist_ok=True)
 with (OUT/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (OUT/"CONTROLLER.json").exists()
  splan=BASE/"reverse_supplement/PLAN.json";supp=json.loads(splan.read_text())
  chroms=sorted({23 if str(r["chromosome"])=="X" else int(r["chromosome"]) for r in supp["reference_requests"]})
  folders=[BASE/"reverse_LD"/f"chr{c}" for c in range(1,23)]+[BASE/"reverse_supplement/reference_LD"/f"chr{c}" for c in chroms]
  script=Path(__file__).with_name("verify_reverse_reference.py");hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),script,splan]}
  state=dict(pid=os.getpid(),status="running",complete=False,total_chromosome_sources=len(folders),completed=[],failures={},running={},inputs_sha256=hashes)
  target=OUT/"CONTROLLER.json";dump(target,state);active={};done=set()
  def one(folder):
   key=str(folder.relative_to(BASE))
   with (OUT/(key.replace("/","_")+".log")).open("a") as log:
    p=subprocess.Popen(["/public/home/mengxl/.conda/envs/CellFM/bin/python",str(script),"--folder",str(folder)],
     cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    code=p.wait()
   if code:raise RuntimeError("Verification exit="+str(code)+", pid="+str(p.pid))
   assert json.loads((folder/"INDEPENDENT_VERIFICATION.json").read_text())["complete"]
   return key
  with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
   while len(done)<len(folders):
    for f in list(active):
     if not f.done():continue
     folder=active.pop(f);key=str(folder.relative_to(BASE));done.add(key);state["running"].pop(key,None)
     try:state["completed"].append(f.result())
     except Exception as exc:state["failures"][key]=repr(exc)
    for folder in folders:
     key=str(folder.relative_to(BASE))
     if len(active)>=2:break
     if key in done or folder in active.values() or not (folder/"MANIFEST.json").exists():continue
     for p,h in hashes.items():assert sha(ROOT/p)==h
     active[pool.submit(one,folder)]=folder;state["running"][key]="submitted"
    state["updated_at_unix"]=time.time();dump(target,state)
    if len(done)==len(folders):break
    if not active:
     live=False;states=[]
     for source in [BASE/"reverse_LD/RUN.json",BASE/"reverse_supplement/reference_LD/RUN.json"]:
      u=json.loads(source.read_text());proc=Path("/proc")/str(u["pid"])/"cmdline"
      live=live or (u["status"]=="running" and proc.exists() and bool(proc.read_bytes()))
      states.append(u["status"])
     if not live:state.update(status="needs_attention_upstream_sources",upstream_status=states);dump(target,state);return
    time.sleep(10)
  state.update(status="complete" if not state["failures"] else "finished_with_verification_errors",complete=not state["failures"]);dump(target,state)
if __name__=="__main__":main()
