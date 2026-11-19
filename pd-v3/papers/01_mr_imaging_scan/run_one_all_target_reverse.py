"""Original reverse-direction analysis for one frozen all-target region; no additional selection layer."""
from pathlib import Path
import argparse,fcntl,hashlib,json,os,time,traceback
import numpy as np,pandas as pd
import reverse_followup_core as core
import reverse_followup_inputs as routing
ROOT=core.ROOT;BASE=routing.BASE;SUP=routing.SUP;OUT=ROOT/"papers/01_mr_imaging_scan/results/all_target_reverse_20261118"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,d):
 t=p.with_name(p.name+".pending");t.write_text(json.dumps(d,indent=2)+"\n");t.replace(p)
def read(p):return pd.read_csv(p,float_precision="round_trip",keep_default_na=False)
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--region",required=True);a=ap.parse_args()
 plan=json.loads((BASE/"PLAN.json").read_text());supp=json.loads((SUP/"PLAN.json").read_text());rec=plan["regions"][a.region]
 folder=OUT/a.region;folder.mkdir(parents=True,exist_ok=True)
 with (folder/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (folder/"COMPLETE.json").exists():
   d=json.loads((folder/"COMPLETE.json").read_text());assert d["complete"]
   for p,h in d["input_sha256"].items():assert sha(ROOT/p)==h
   for p,h in d["outputs_sha256"].items():assert sha(folder/p)==h
   print("CACHED_REVERSE",a.region);return
  assert not (folder/"FAILED.json").exists(),"Inspect previous failure"
  check=json.loads((OUT/"RECIPE_VERIFICATION.json").read_text());assert check["complete"]
  for p,h in check["code_sha256"].items():assert sha(ROOT/p)==h
  peaks,source,refdirs,paths=routing.requirements(rec,plan,supp);assert all(p.exists() for p in paths),"Required inputs are not yet complete"
  dependencies=[BASE/"PLAN.json",SUP/"PLAN.json",Path(__file__),Path(core.__file__),Path(routing.__file__),ROOT/"pdv3/mr.py",
   OUT/"RECIPE_VERIFICATION.json"]+paths
  x=pd.DataFrame(peaks) if peaks else pd.DataFrame(columns=core.PEAK_COLUMNS)
  # Python float preserves published decimal conversion, including long workbook values.
  for column in ["effect_allele_frequency","discovery_neglog10_p","discovery_beta","discovery_se","replication_neglog10_p","replication_beta","replication_se"]:
   if column in x and len(x):x[column]=x[column].map(float)
  y=pd.DataFrame(columns=core.OUTCOME_COLUMNS);source_receipt=None
  if source:
   source_receipt=json.loads(source.read_text());assert source_receipt["complete"]
   yp=source.parent/"imaging_peak_outcomes.csv.gz"
   if rec["phenotype"].startswith("BIG40:"):
    assert source_receipt["plan_sha256"]==sha(BASE/"PLAN.json")
    assert sha(yp)==source_receipt["outputs_sha256"][str(yp.relative_to(ROOT))]
   else:
    assert sha(yp)==source_receipt["outputs_sha256"][yp.name]
    for p,h in source_receipt["input_sha256"].items():
     # Full original sources were verified by their acquisition and outcome receipts;
     # compare the small immutable control inputs here, not 0.5 GB per repeated readout.
     if not p.endswith((".tsv",".tsv.gz")):assert sha(ROOT/p)==h
   dependencies.append(yp);y=read(yp)
  genotypes={};refcoverage={};reference_samples=None
  for directory in refdirs:
   m=json.loads((directory/"MANIFEST.json").read_text());v=json.loads((directory/"INDEPENDENT_VERIFICATION.json").read_text())
   assert m["complete"] and v["complete"] and v["manifest_sha256"]==sha(directory/"MANIFEST.json")
   for name,h in m["outputs_sha256"].items():
    assert sha(directory/name)==h;dependencies.append(directory/name)
   cv=read(directory/"coverage.csv")
   with np.load(directory/"effect_allele_dosages.npz",allow_pickle=False) as z:
    samples=z["samples"];assert len(samples)==503
    if reference_samples is None:reference_samples=samples.copy()
    else:assert np.array_equal(samples,reference_samples)
    for key,g in zip(z["variant_keys"],z["dosages"]):
     assert key not in genotypes,"Duplicate genotype source";genotypes[key]=g.copy()
   for row in cv.to_dict("records"):
    assert row["variant_key"] not in refcoverage;refcoverage[row["variant_key"]]=row
  for row in peaks:
   if routing.snv(row):
    key=core.key(row["chromosome"],row["position_grch37"],row["effect_allele"],row["other_allele"])
    assert key in refcoverage,"Missing reference coverage must not be silently treated as absent SNP"
  hashes={str(p.relative_to(ROOT)):sha(p) for p in dependencies}
  run=dict(pid=os.getpid(),status="running",complete=False,region_id=a.region,started_at_unix=time.time(),input_sha256=hashes)
  dump(folder/"RUN.json",run)
  try:
   result=core.calculate(rec,x,y,genotypes)
   names={"audit":"all_peaks_exposure_only_LD_audit.csv","pairs":"all_peak_pair_reference_LD.csv",
    "ratios":"all_peak_reverse_ratios_and_coverage.csv","summary":"reverse_direction_estimates.csv"}
   for name,frame in result.items():
    frame["region_id"]=a.region;frame.to_csv(folder/names[name],index=False)
   assert len(result["summary"])==4 and len(result["ratios"])==4*len(peaks)
   assert len(result["audit"])==len(peaks)
   reference_table=[]
   for row in peaks:
    k=core.key(row["chromosome"],row["position_grch37"],row["effect_allele"],row["other_allele"])
    reference_table.append(dict(variant_key=k,reference_status=refcoverage[k]["status"] if k in refcoverage else "unsupported_non_biallelic_SNV"))
   pd.DataFrame(reference_table,columns=["variant_key","reference_status"]).to_csv(folder/"reference_coverage.csv",index=False)
   for p,h in hashes.items():assert sha(ROOT/p)==h,"Frozen input changed during calculation"
   run.update(status="complete",complete=True,accession=rec["accession"],gene=rec["gene"],phenotype=rec["phenotype"],
    source_peaks=len(peaks),ratio_coverage_rows=len(result["ratios"]),summary_rows=4,
    scientific_status_counts=result["summary"].status.value_counts().to_dict(),
    no_published_peak=not bool(peaks),method="Unchanged original published-peak reverse recipe; see RECIPE_VERIFICATION",
    LD_selection_before_outcome_join=True,imaging_replication="Same discovery-selected IVs, replication effect sensitivity; no replication-P selection",
    reference="External 1000 Genomes Phase3 EUR503 GRCh37; not cohort-matched",
    X_reference="Frozen supplemental ploidy convention only where required; original autosomal arithmetic unchanged",
    no_proxy_substitution=True,nominal_P_not_candidate_acceptance=True,new_main_MR_screen=False,
    protein_units="Original AGES units, not assumed equivalent to Olink or ADT CLR",
    limitations=["External small LD reference cannot establish cohort independence","Discovery instrument effects may have winner bias",
     "Horizontal pleiotropy, tissue mediation and therapeutic direction are not established"],
    outputs_sha256={p.name:sha(p) for p in folder.glob("*.csv")},elapsed_seconds=time.time()-run["started_at_unix"])
   dump(folder/"COMPLETE.json",run);dump(folder/"RUN.json",run)
   print(json.dumps(dict(region=a.region,source_peaks=len(peaks),statuses=run["scientific_status_counts"])),flush=True)
  except BaseException as exc:
   run.update(status="failed",error=repr(exc),traceback=traceback.format_exc());dump(folder/"FAILED.json",run);dump(folder/"RUN.json",run);raise
if __name__=="__main__":main()
