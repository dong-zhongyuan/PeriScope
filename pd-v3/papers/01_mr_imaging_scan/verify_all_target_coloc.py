"""Regression validation against all three previously completed regional analyses."""
from pathlib import Path
import json,subprocess,time
import numpy as np,pandas as pd
import run_one_all_target_coloc as worker
ROOT=worker.ROOT;OUT=worker.OUT;OLD=ROOT/"papers/01_mr_imaging_scan/results/new_scan_regions_20261117"
def main():
 plan=json.loads((worker.BASE/"PLAN.json").read_text());records=[]
 source=ROOT/"data/new_scan_regional_20261117"
 for acc in ["GCST90087997","GCST90090368","GCST90089324"]:
  rec=next(v for v in plan["regions"].values() if v["accession"]==acc)
  xp=source/acc/"regional_statistics.csv.gz";yp=source/acc/"imaging_region_GRCh37.csv.gz"
  audit=json.loads((OLD/acc/"INPUT_AUDIT.json").read_text())
  for p in [xp,yp]:assert worker.sha(p)==audit["input_sha256"][str(p.relative_to(ROOT))]
  if acc=="GCST90087997":
   folder=OUT/rec["region_id"];assert json.loads((folder/"COMPLETE.json").read_text())["complete"]
   route="complete all-target source extraction and worker"
  else:
   folder=OUT/"recipe_verification"/acc;folder.mkdir(parents=True,exist_ok=True)
   a=worker.prepare_frames(rec,worker.read_frame(xp),worker.read_frame(yp),folder,plan["sources"][acc]["metadata"]["sample_size"])
   worker.dump(folder/"INPUT_AUDIT.json",a)
   with (folder/"R.log").open("w") as log:
    subprocess.run(["/public/home/mengxl/dzy/envs/rcoloc/bin/Rscript",str(Path(__file__).with_name("all_target_coloc.R")),str(ROOT),str(folder)],stdout=log,stderr=subprocess.STDOUT,check=True)
   worker.validate_outputs(folder);route="same frozen legacy source files, generalized preparation and R"
  old=pd.read_csv(OLD/acc/"shared_snps.csv.gz",float_precision="round_trip").sort_values("variant_key").reset_index(drop=True)
  new=pd.read_csv(folder/"shared_snps.csv.gz",float_precision="round_trip").sort_values("variant_key").reset_index(drop=True)
  pd.testing.assert_frame_equal(old,new[old.columns],check_dtype=False,check_exact=True)
  prior=pd.read_csv(OLD/"coloc_prior_sensitivity.csv",float_precision="round_trip")
  prior=prior.loc[prior.accession.eq(acc)].sort_values(["window_bp","p12_requested"]).reset_index(drop=True)
  current=pd.read_csv(folder/"coloc_prior_sensitivity.csv",float_precision="round_trip").sort_values(["window_bp","p12_requested"]).reset_index(drop=True)
  pd.testing.assert_frame_equal(prior,current[prior.columns],check_dtype=False,rtol=1e-12,atol=1e-14)
  posterr={}
  for w in [500000,1000000]:
   p=pd.read_csv(OLD/acc/f"window_{w}_snp_posteriors.csv.gz",float_precision="round_trip").sort_values("variant_key").reset_index(drop=True)
   q=pd.read_csv(folder/f"window_{w}_snp_posteriors.csv.gz",float_precision="round_trip").sort_values("variant_key").reset_index(drop=True)
   pd.testing.assert_frame_equal(p,q[p.columns],check_dtype=False,rtol=1e-12,atol=1e-14)
   names=["log_abf_protein","log_abf_imaging","snp_pp_H4_conditional","snp_pp_protein_conditional","snp_pp_imaging_conditional"]
   posterr[str(w)]=float(np.max(abs(p[names].to_numpy()-q[names].to_numpy())))
  records.append(dict(accession=acc,gene=rec["gene"],shared_snps=len(new),prior_rows=len(current),route=route,
   full_shared_input_fields_exact=True,all_prior_result_fields_match=True,all_SNP_posterior_fields_match=True,
   posterior_max_absolute_error=posterr,input_sha256={str(p.relative_to(ROOT)):worker.sha(p) for p in [xp,yp]},
   result_sha256=worker.sha(folder/"coloc_prior_sensitivity.csv")))
  print("REPRODUCED",acc,len(new),posterr,flush=True)
 # Exercise genuine absent-data handling using empty copies; these are software tests, never biological results.
 rec=dict(rec);rec["region_id"]="software_empty_window_test"
 folder=OUT/"recipe_verification"/"empty_window";folder.mkdir(exist_ok=True)
 worker.prepare_frames(rec,worker.read_frame(xp).iloc[:0],worker.read_frame(yp).iloc[:0],folder,plan["sources"][rec["accession"]]["metadata"]["sample_size"])
 with (folder/"R.log").open("w") as log:
  subprocess.run(["/public/home/mengxl/dzy/envs/rcoloc/bin/Rscript",str(Path(__file__).with_name("all_target_coloc.R")),str(ROOT),str(folder)],stdout=log,stderr=subprocess.STDOUT,check=True)
 statuses=worker.validate_outputs(folder);assert set(statuses.values())=={"undefined_no_shared_snps"}
 report=dict(complete=True,scientific_regression=records,empty_window_software_test=statuses,
  no_change_to_prior_grid_window_or_effect_prior=True,code_sha256={str(p.relative_to(ROOT)):worker.sha(p) for p in [
   Path(__file__),Path(worker.__file__),Path(__file__).with_name("all_target_coloc.R"),Path(__file__).with_name("prepare_new_scan_coloc.py"),ROOT/"reference/coloc/claudia.R"]})
 worker.dump(OUT/"RECIPE_VERIFICATION.json",report)
if __name__=="__main__":main()
