"""Acquire one full native BIG40 file and extract all its frozen regional requests."""
from pathlib import Path
import argparse,fcntl,gzip,hashlib,json,os,sys,time,traceback
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(Path(__file__).parent))
import native_region_io as native
import s09_acquire_qsm as downloader
BASE=ROOT/"data/all_target_followup_20261118/regional_all/imaging"
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--idp",required=True);ap.add_argument("--connections",type=int,choices=range(1,9),default=4)
 args=ap.parse_args();plan=json.loads((BASE/"PLAN.json").read_text());phenotype="BIG40:release2:discovery:"+args.idp
 rec=plan["phenotypes"][phenotype];assert rec["family"]=="BIG40"
 folder=BASE/"BIG40"/args.idp;folder.mkdir(parents=True,exist_ok=True)
 with (folder/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  ph=native.sha(BASE/"PLAN.json")
  if (folder/"COMPLETE.json").exists():
   d=json.loads((folder/"COMPLETE.json").read_text());assert d["complete"] and d["plan_sha256"]==ph
   for v in d["windows"].values():assert native.sha(ROOT/v["output"])==v["sha256"]
   print(json.dumps({"cached_verified":args.idp}));return
  if (folder/"FAILED.json").exists():raise RuntimeError("Inspect/archive prior failure before explicit retry")
  rm=BASE/"variant_reference/MANIFEST.json";ref_hash=native.sha(rm);ref=json.loads(rm.read_text());assert ref["complete"] and ref["plan_sha256"]==ph
  deps=[Path(__file__),Path(native.__file__),Path(downloader.__file__),Path(downloader._shared.__file__)]
  code={str(p.relative_to(ROOT)):native.sha(p) for p in deps}
  started=time.time();state=dict(pid=os.getpid(),status="running",complete=False,plan_sha256=ph,code_sha256=code,
   phenotype=phenotype,started_at_unix=started)
  native.dump(folder/"RUN.json",state)
  try:
   known=rec["prior_full_native_source_identities"]
   for x in known:assert native.sha(ROOT/x["receipt"])==x["receipt_sha256"]
   hashes={x["source_sha256"] for x in known};assert len(hashes)<=1
   target=folder/(args.idp+".txt.gz")
   cache=[ROOT/p for p in rec["existing_full_payloads"] if (ROOT/p).exists()]
   transport=None;staged=False
   if cache:
    path=cache[0];digest=native.sha(path);assert hashes and digest in hashes
    transport=dict(route="reuse_existing_validated_native_payload",path=str(path.relative_to(ROOT)),sha256=digest)
   elif target.exists():
    receipt=json.loads((folder/"PAYLOAD.json").read_text());assert native.sha(target)==receipt["sha256"]
    path=target;digest=receipt["sha256"];transport=receipt
   else:
    downloader.URL=rec["native_url"];downloader.BLOCK=4<<20
    identity=downloader.http_identity()
    path,blocks=downloader.download_blocks(target,identity,args.connections,folder/"DOWNLOAD_PROGRESS.json")
    digest=native.sha(path);staged=True
    transport=dict(route="parallel_verified_HTTP_ranges",http=identity,sha256=digest,bytes=path.stat().st_size,
     block_size=downloader.BLOCK,blocks=blocks["completed_blocks"])
   if hashes and digest not in hashes:raise ValueError("Native source changed from prior verified payload")
   windows={wid:plan["windows"][wid] for wid in rec["windows"]}
   def progress(n):
    state.update(status="validating_and_extracting",source_rows_processed=n,elapsed_seconds=time.time()-started)
    native.dump(folder/"RUN.json",state);print(json.dumps({"idp":args.idp,"source_rows":n}),flush=True)
   with gzip.open(path,"rb") as gz:frames,validation=native.extract(gz,rec,windows,progress)
   assert validation["chromosome_rows"]==ref["chromosome_rows"]
   validation["gzip_crc_eof_verified"]=True;validation["source_sha256"]=digest
   outputs={}
   for wid,frame in frames.items():
    rp=ROOT/ref["outputs"][wid]["output"];assert native.sha(rp)==ref["outputs"][wid]["sha256"]
    frequency=pd.read_csv(rp,dtype={"chromosome":str},float_precision="round_trip",keep_default_na=False)
    result=native.attach_frequency(frame,frequency)
    out=folder/(wid+".csv.gz");result.to_csv(out,index=False,compression={"method":"gzip","mtime":0})
    outputs[wid]=dict(output=str(out.relative_to(ROOT)),sha256=native.sha(out),rows=len(result),
     frequency_reference_sha256=native.sha(rp),region=windows[wid],
     frequency_definition="Official full BIG40 catalog a2 frequency, not per-IDP per-SNP sample-specific frequency")
   for f,h in code.items():assert native.sha(ROOT/f)==h
   assert native.sha(BASE/"PLAN.json")==ph and native.sha(rm)==ref_hash
   if staged:
    native.dump(folder/"PAYLOAD.json",dict(**transport,validation=validation))
    path.replace(target);target.with_name(target.name+".part.json").unlink()
   report=dict(**state,**validation,windows=outputs,transport=transport,
    variant_reference_manifest_sha256=native.sha(rm),matched_prior_native_source_identity=bool(hashes),
    matched_existing_scientific_recipe="Full original GRCh37 region; all alleles retained; exact native catalog frequency join")
   report.update(complete=True,status="complete",elapsed_seconds=time.time()-started);native.dump(folder/"COMPLETE.json",report);native.dump(folder/"RUN.json",report)
   print(json.dumps({"complete":args.idp,"regions":len(rec["regions"]),"windows":len(outputs),"source_rows":validation["source_rows"]}),flush=True)
  except BaseException as exc:
   state.update(status="failed",error=repr(exc),traceback=traceback.format_exc());native.dump(folder/"FAILED.json",state);native.dump(folder/"RUN.json",state);raise
if __name__=="__main__":main()
