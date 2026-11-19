"""Extract published QSM peak outcomes from already verified complete AGES source files."""
from pathlib import Path
import argparse,copy,csv,fcntl,json,os,sys,time,traceback
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(Path(__file__).parent))
import acquire_all_target_protein_regions as original
BASE=original.BASE;SUP=BASE/"reverse_supplement";OUT=SUP/"protein_outcomes"
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--accession",required=True);a=ap.parse_args()
 plan=json.loads((BASE/"PLAN.json").read_text());supp=json.loads((SUP/"PLAN.json").read_text())
 request=supp["sources"][a.accession];folder=OUT/a.accession;folder.mkdir(parents=True,exist_ok=True)
 with (folder/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (folder/"COMPLETE.json").exists():
   d=json.loads((folder/"COMPLETE.json").read_text());assert d["complete"]
   for p,h in d["input_sha256"].items():assert original.sha(ROOT/p)==h
   for p,h in d["outputs_sha256"].items():assert original.sha(folder/p)==h
   print("CACHED_QSM_OUTCOME",a.accession);return
  assert not (folder/"FAILED.json").exists(),"Inspect prior failure"
  dependencies=[BASE/"PLAN.json",SUP/"PLAN.json",Path(__file__),Path(original.__file__),Path(original.bridge.__file__)]
  source=copy.deepcopy(plan["sources"][a.accession]);source["peak_rows"]=request["peak_rows"];source["regions"]=[]
  started=time.time();state=dict(pid=os.getpid(),status="running",complete=False,accession=a.accession,started_at_unix=started)
  original.dump(folder/"RUN.json",state)
  try:
   metadata=None
   if source["peak_rows"]:
    cp=BASE/"sources"/a.accession/"COMPLETE.json";completed=json.loads(cp.read_text())
    assert completed["complete"] and completed["plan_sha256"]==original.sha(BASE/"PLAN.json")
    path=ROOT/completed["input_path"];assert original.sha(path)==completed["input_sha256"]
    formatted=completed["input_format"]=="GWAS_Catalog_formatted_reconstructed_author"
    dependencies.extend([cp,path])
    _,outcomes,metadata=original.extract(path,formatted,source,{"regions":{}})
    assert metadata["original_md5"]==completed["original_md5"]
    assert metadata["original_sha256"]==completed["original_sha256"] and metadata["source_rows"]==completed["source_rows"]
   else:outcomes=[]
   inputs={str(p.relative_to(ROOT)):original.sha(p) for p in dependencies}
   original.write_table(folder/"imaging_peak_outcomes.csv.gz",outcomes,original.FIELDS)
   fields=list(source["peak_rows"][0]) if source["peak_rows"] else ["phenotype","variant_key","source_variant_id"]
   original.write_table(folder/"all_source_QSM_peaks.csv.gz",source["peak_rows"],fields)
   obtained={original.key(x["chromosome"],x["base_pair_location"],x["effect_allele"],x["other_allele"]) for x in outcomes}
   assert len(obtained)==len(outcomes),"Ambiguous exact protein outcome identity"
   coverage=[]
   for p in source["peak_rows"]:
    k=original.key(p["chromosome"],p["position_grch37"],p["effect_allele"],p["other_allele"])
    coverage.append(dict(phenotype=p["phenotype"],variant_key=p["variant_key"],source_variant_id=p["source_variant_id"],
     protein_outcome_available=k in obtained,protein_outcome_status="available_exact_identity" if k in obtained else "absent_from_complete_author_source"))
   original.write_table(folder/"outcome_coverage.csv.gz",coverage,["phenotype","variant_key","source_variant_id","protein_outcome_available","protein_outcome_status"])
   for p,h in inputs.items():assert original.sha(ROOT/p)==h,"Frozen input changed"
   state.update(status="complete",complete=True,regions=request["regions"],source_peak_rows=len(source["peak_rows"]),
    outcome_rows=len(outcomes),no_entry_in_complete_author_peak_table=not bool(source["peak_rows"]),
    source_file_read=bool(source["peak_rows"]),full_author_identity=metadata,input_sha256=inputs,
    outputs_sha256={p.name:original.sha(p) for p in folder.glob("*.csv.gz")},
    source_download_repeated=False,elapsed_seconds=time.time()-started)
   original.dump(folder/"COMPLETE.json",state);original.dump(folder/"RUN.json",state)
   print("QSM_OUTCOME",a.accession,len(source["peak_rows"]),len(outcomes),flush=True)
  except BaseException as exc:
   state.update(status="failed",error=repr(exc),traceback=traceback.format_exc());original.dump(folder/"FAILED.json",state);original.dump(folder/"RUN.json",state);raise
if __name__=="__main__":main()
