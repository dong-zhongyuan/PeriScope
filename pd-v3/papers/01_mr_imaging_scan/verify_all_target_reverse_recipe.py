"""Replay the three verified original analyses and test explicit unavailable-source states."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
import reverse_followup_core as core
ROOT=core.ROOT;DATA=ROOT/"data/new_scan_regional_20261117"
OLD=ROOT/"papers/01_mr_imaging_scan/results/new_scan_reverse_20261117"
OUT=ROOT/"papers/01_mr_imaging_scan/results/all_target_reverse_20261118"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=True);assert not (OUT/"RECIPE_VERIFICATION.json").exists()
 assert json.loads((OLD/"INDEPENDENT_VERIFICATION.json").read_text())["complete"]
 g={}
 for folder in (DATA/"reverse_LD").glob("chr*"):
  m=json.loads((folder/"MANIFEST.json").read_text());assert m["complete"]
  for f,h in m["outputs_sha256"].items():assert sha(folder/f)==h
  with np.load(folder/"effect_allele_dosages.npz",allow_pickle=False) as z:
   for k,d in zip(z["variant_keys"],z["dosages"]):g[k]=d.copy()
 names={"audit":"all_peaks_exposure_only_LD_audit.csv","pairs":"all_peak_pair_reference_LD.csv",
 "ratios":"all_peak_reverse_ratios_and_coverage.csv","summary":"reverse_direction_estimates.csv"}
 originals={k:pd.read_csv(OLD/n,float_precision="round_trip") for k,n in names.items()}
 recs=json.loads((DATA/"FIRST_DOWNLOAD_BATCH.json").read_text())["records"];details=[]
 for rec in recs:
  acc=rec["accession"];rec=dict(rec,phenotype=rec["best_phenotype"])
  x=pd.read_csv(DATA/acc/"all_source_imaging_peaks.csv",float_precision="round_trip")
  y=pd.read_csv(DATA/acc/"imaging_peak_outcomes.csv.gz",float_precision="round_trip")
  result=core.calculate(rec,x,y,g)
  folder=OUT/"recipe_verification"/acc;folder.mkdir(parents=True,exist_ok=True)
  for k,n in names.items():
   result[k].to_csv(folder/n,index=False)
   current=pd.read_csv(folder/n,float_precision="round_trip")
   expected=originals[k][originals[k].accession.eq(acc)].reset_index(drop=True)
   # This deliberately compares original coverage/selection/output fields, not just P values.
   pd.testing.assert_frame_equal(expected,current[expected.columns],check_dtype=False,rtol=1e-12,atol=1e-14)
  empty=core.calculate(rec,x.iloc[:0],y.iloc[:0],g)
  assert len(empty["summary"])==4 and empty["summary"].status.eq("not_estimated_no_published_peak").all()
  missing=core.calculate(rec,x,y.iloc[:0],g)
  assert missing["summary"].status.eq("not_estimated_no_usable_instrument").all()
  # Outcome removal must not alter discovery-only clumping.
  pd.testing.assert_frame_equal(result["audit"][["variant_key","ld_selection_status","ld_blocker"]],
   missing["audit"][["variant_key","ld_selection_status","ld_blocker"]],check_exact=True)
  details.append(dict(accession=acc,source_rows=len(x),pair_LD_rows=len(result["pairs"]),
   all_ratio_rows=len(result["ratios"]),summary_rows=len(result["summary"]),all_original_fields_reproduced=True,
   source_empty_and_outcome_missing_states_explicit=True,LD_selection_invariant_to_outcome_removal=True))
  print("REVERSE_REPRODUCED",acc,len(x),flush=True)
 report=dict(complete=True,records=details,method_parameters_changed=False,
  code_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),Path(core.__file__),ROOT/"pdv3/mr.py"]},
  legacy_verification_sha256=sha(OLD/"INDEPENDENT_VERIFICATION.json"))
 (OUT/"RECIPE_VERIFICATION.json").write_text(json.dumps(report,indent=2)+"\n")
if __name__=="__main__":main()
