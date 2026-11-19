"""Extract all frozen QSM regional requests from the already verified local archive."""
from pathlib import Path
import fcntl,gzip,json,os,sys,tarfile,time,traceback
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(Path(__file__).parent))
import native_region_io as native
BASE=ROOT/"data/all_target_followup_20261118/regional_all/imaging"
def main():
 out=BASE/"OxfordQSM2022";out.mkdir(exist_ok=True)
 with (out/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (out/"RUN.json").exists():raise RuntimeError("Inspect previous archive extraction before recovery")
  plan=json.loads((BASE/"PLAN.json").read_text());ph=native.sha(BASE/"PLAN.json")
  selected={r["source_member"]:r for r in plan["phenotypes"].values() if r["family"]=="OxfordQSM2022"}
  source_manifest=ROOT/"data/oxford_qsm/MANIFEST.json";manifest=json.loads(source_manifest.read_text())
  archive=source_manifest.parent/manifest["path"];assert archive.stat().st_size==manifest["size_bytes"]
  code={str(p.relative_to(ROOT)):native.sha(p) for p in [Path(__file__),Path(native.__file__)]}
  started=time.time();state=dict(pid=os.getpid(),status="checking_archive_identity",complete=False,plan_sha256=ph,
   code_sha256=code,started_at_unix=started,completed_members=[],planned_members=sorted(selected))
  native.dump(out/"RUN.json",state)
  try:
   assert native.sha(archive)==manifest["sha256"]
   assert manifest["validation"]["gzip_crc_and_eof_verified"]
   records={};inventory=[];seen=set()
   with gzip.open(archive,"rb") as gz:
    with tarfile.open(fileobj=gz,mode="r|") as tar:
     for member in tar:
      if not member.isfile():continue
      name=Path(member.name).name
      inventory.append(dict(name=member.name,bytes=member.size,selected=name in selected))
      stream=tar.extractfile(member)
      if name not in selected:
       while stream.read(8<<20):pass
       continue
      assert name not in seen;seen.add(name);rec=selected[name];dest=out/rec["idp"];dest.mkdir(exist_ok=True)
      windows={wid:plan["windows"][wid] for wid in rec["windows"]}
      def progress(n):
       state.update(status="extracting",current_member=name,source_rows_processed=n,elapsed_seconds=time.time()-started)
       native.dump(out/"RUN.json",state);print(json.dumps({"member":name,"source_rows":n}),flush=True)
      frames,validation=native.extract(stream,rec,windows,progress)
      outputs={}
      for wid,frame in frames.items():
       p=dest/(wid+".csv.gz");frame.to_csv(p,index=False,compression={"method":"gzip","mtime":0})
       outputs[wid]=dict(output=str(p.relative_to(ROOT)),sha256=native.sha(p),rows=len(frame),region=windows[wid],
        exact_IDP_N_and_EAF="not supplied; retained as missing, no cohort-size substitution")
      records[rec["phenotype"]]=dict(phenotype=rec["phenotype"],member=member.name,windows=outputs,**validation)
      state["completed_members"].append(name);native.dump(out/"RUN.json",state)
      print(json.dumps({"extracted_member":name,"windows":len(outputs)}),flush=True)
    # Consume tar padding to verify this full gzip pass, as in the original recipe.
    while gz.read(8<<20):pass
   assert seen==set(selected)
   for p,h in code.items():assert native.sha(ROOT/p)==h
   assert native.sha(BASE/"PLAN.json")==ph
   for phenotype,rec in records.items():
    idp=phenotype.split(":")[-1]
    report=dict(complete=True,plan_sha256=ph,**rec,source_archive_sha256=manifest["sha256"],
     archive_manifest_sha256=native.sha(source_manifest),gzip_crc_eof_verified=True,code_sha256=code,
     frequency_source=None,imaging_SD=1,units="Gaussianized and deconfounded SD units",
     exact_IDP_N_not_imputed=True,EAF_not_imputed=True,biological_inference_computed=False)
    native.dump(out/idp/"COMPLETE.json",report)
   native.dump(out/"MANIFEST.json",dict(complete=True,plan_sha256=ph,source_archive_sha256=manifest["sha256"],
    source_manifest_sha256=native.sha(source_manifest),gzip_crc_eof_verified=True,records=records,inventory=inventory,
    code_sha256=code,elapsed_seconds=time.time()-started,full_archive_download_repeated=False))
   state.update(status="complete",complete=True,elapsed_seconds=time.time()-started);native.dump(out/"RUN.json",state)
  except BaseException as exc:
   state.update(status="failed",error=repr(exc),traceback=traceback.format_exc());native.dump(out/"RUN.json",state);raise
if __name__=="__main__":main()
