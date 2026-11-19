"""Execute the existing ABF recipe for one frozen all-target region, with full coverage audit."""
from pathlib import Path
import argparse,fcntl,hashlib,json,os,subprocess,sys,time,traceback
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent))
from prepare_new_scan_coloc import key as autosomal_key,write
from acquire_all_target_protein_regions import normchr
def chrom_label(v):
 if isinstance(v,float) and float(v).is_integer():v=int(v)
 return normchr(v)
def region_key(row):
 """Frozen key for numeric chromosomes; normchr label for chrX."""
 label=chrom_label(row.chromosome)
 if label.isdigit():return autosomal_key(row)
 a,b=row.effect_allele,row.other_allele
 if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b:return None
 return f'{label}:{int(row.base_pair_location)}:'+'_'.join(sorted([a,b]))
BASE=ROOT/"data/all_target_followup_20261118/regional_all"
OUT=ROOT/"papers/01_mr_imaging_scan/results/all_target_coloc_20261118"
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def dump(p,d):
 q=p.with_name(p.name+".pending");q.write_text(json.dumps(d,indent=2)+"\n");q.replace(p)
def read_frame(p):
 return pd.read_csv(p,float_precision="round_trip",keep_default_na=False,
  na_values={k:[""] for k in ["n_y","effect_allele_frequency","imputation_info"]})
def prepare_frames(rec,x,y,folder,protein_n):
 folder.mkdir(parents=True,exist_ok=True);x=x.copy();y=y.copy()
 for frame in [x,y]:
  assert np.isfinite(frame[["beta","standard_error"]].to_numpy(dtype=float)).all()
  assert frame.standard_error.gt(0).all()
  assert frame.chromosome.map(chrom_label).eq(chrom_label(rec["chromosome"])).all()
  assert frame.base_pair_location.between(rec["start_grch37"],rec["end_grch37"]).all()
  frame["variant_key"]=frame.apply(region_key,axis=1) if len(frame) else pd.Series(dtype=object)
  assert not frame.loc[frame.variant_key.notna(),"variant_key"].duplicated().any(),"Ambiguous regional SNP identity"
 assert np.isfinite(x.effect_allele_frequency.to_numpy(dtype=float)).all()
 assert x.effect_allele_frequency.between(0,1).all() and x.n_metadata.eq(protein_n).all()
 if rec["phenotype"].startswith("BIG40:"):
  assert np.isfinite(y.effect_allele_frequency.to_numpy(dtype=float)).all()
  assert y.effect_allele_frequency.between(0,1).all()
 else:
  assert y.effect_allele_frequency.isna().all() and y.n_y.isna().all(),"QSM missing frequency/N must remain undefined"
 cols=["variant_key","chromosome","base_pair_location","variant_id","effect_allele","other_allele","effect_allele_frequency","beta","standard_error","p_value"]
 a=x.loc[x.variant_key.notna(),cols+["n_metadata"]].rename(columns={"n_metadata":"n"})
 b=y.loc[y.variant_key.notna(),cols+["n_y"]].rename(columns={"n_y":"n"})
 a=a.rename(columns={c:c+"_x" for c in a if c!="variant_key"})
 b=b.rename(columns={c:c+"_y" for c in b if c!="variant_key"})
 union=a.merge(b,on="variant_key",how="outer",validate="one_to_one",indicator=True)
 counts={str(k):int(v) for k,v in union["_merge"].value_counts().items()}
 shared=union.loc[union["_merge"].eq("both")].copy().drop(columns="_merge")
 assert shared.chromosome_x.map(chrom_label).eq(shared.chromosome_y.map(chrom_label)).all()
 assert shared.base_pair_location_x.eq(shared.base_pair_location_y).all()
 shared["pos37"]=shared.base_pair_location_x.astype(int)
 for c in ["accession","gene","phenotype"]:shared[c]=rec[c]
 shared["palindromic"]=shared.variant_key.str.endswith(("A_T","C_G")).astype(bool)
 shared["direction_status"]=np.where(shared.palindromic,"not_resolved_for_signed_use","exact_allele_alignment_required_for_signed_use")
 shared=shared.sort_values(["pos37","variant_key"])
 outputs={}
 for name,frame,filename in [
  ("union",union,"variant_identity_audit.csv.gz"),("shared",shared,"shared_snps.csv.gz"),
  ("unsupported_protein",x[x.variant_key.isna()],"unsupported_protein_rows.csv.gz"),
  ("unsupported_imaging",y[y.variant_key.isna()],"unsupported_imaging_rows.csv.gz")]:
  outputs[name]=write(frame,folder/filename)
 pd.DataFrame([dict(region_id=rec["region_id"],accession=rec["accession"],gene=rec["gene"],
  phenotype=rec["phenotype"],center=rec["center_grch37"],protein_N_metadata=protein_n)]).to_csv(folder/"region.csv",index=False)
 audit=dict(region=rec,protein_source_rows=len(x),imaging_source_rows=len(y),match_counts=counts,
  shared_snps=len(shared),unsupported_protein_rows=int(x.variant_key.isna().sum()),
  unsupported_imaging_rows=int(y.variant_key.isna().sum()),outputs=outputs,
  imaging_frequency_missing_rows=int(y.effect_allele_frequency.isna().sum()),imaging_N_missing_rows=int(y.n_y.isna().sum()),
  matching="GRCh37 chromosome, coordinate, exact unordered alleles; biallelic SNVs only",
  P_or_MAF_filter=None,protein_SD="Pinned coloc sdY.est with reported assay N and protein frequency",
  imaging_SD="1: original Gaussianized/deconfounded phenotype SD units; not inferred from missing QSM N/EAF",
  imaging_SD_reference="https://chr1swallace.github.io/coloc/reference/check_dataset.html",
  palindromes="Squared-Z ABF only; no signed-direction conclusion",
  cross_platform_scope="AGES aptamer and UKB Olink share author gene label; assay equivalence not established",
  model="Single causal variant ABF; does not resolve multiple causal signals",
  windows_bp=[500000,1000000],p12="Existing 32-value grid from 1e-6 through 3e-5, including 1e-5",
  decision_threshold=None)
 if not chrom_label(rec["chromosome"]).isdigit():
  audit["chromosome_scope"]=("chrX included on the author-published per-variant statistics as given: "
   "no ploidy recoding, no sex-stratified reanalysis, and sdY.est uses the reported assay N and protein "
   "frequency exactly as for autosomes; the reverse-MR chrX ploidy question remains open")
 return audit
def validate_outputs(folder):
 result=pd.read_csv(folder/"coloc_prior_sensitivity.csv",float_precision="round_trip",keep_default_na=False,na_values=["NA"])
 assert len(result)==64 and sorted(result.window_bp.unique())==[500000,1000000]
 assert result.groupby("window_bp").size().eq(32).all()
 shared=read_frame(folder/"shared_snps.csv.gz");rec=pd.read_csv(folder/"region.csv").iloc[0]
 probs=["PP.H"+str(i)+".abf" for i in range(5)]
 for window,group in result.groupby("window_bp"):
  n=int((abs(shared.pos37-rec.center)<=window).sum());assert group.nsnps.eq(n).all()
  assert group.status.nunique()==1
  if group.iloc[0].status=="defined":
   assert np.isfinite(group[probs].to_numpy()).all()
   np.testing.assert_allclose(group[probs].sum(axis=1),1,rtol=0,atol=1e-10)
   assert group[probs].ge(0).all().all() and group[probs].le(1).all().all()
   post=pd.read_csv(folder/f"window_{int(window)}_snp_posteriors.csv.gz")
   assert len(post)==n
   np.testing.assert_allclose(post[["snp_pp_H4_conditional","snp_pp_protein_conditional","snp_pp_imaging_conditional"]].sum(),1,rtol=0,atol=1e-10)
  else:assert group[probs].isna().all().all()
 return {str(int(w)):g.iloc[0].status for w,g in result.groupby("window_bp")}
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--region",required=True);args=ap.parse_args()
 plan=json.loads((BASE/"PLAN.json").read_text());ip=json.loads((BASE/"imaging/PLAN.json").read_text())
 rec=plan["regions"][args.region];source=plan["sources"][rec["accession"]];ph=ip["phenotypes"][rec["phenotype"]]
 folder=OUT/args.region;folder.mkdir(parents=True,exist_ok=True)
 with (folder/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (folder/"COMPLETE.json").exists():
   d=json.loads((folder/"COMPLETE.json").read_text());assert d["complete"]
   for p,h in d["input_sha256"].items():assert sha(ROOT/p)==h
   for p,h in d["outputs_sha256"].items():assert sha(folder/p)==h
   print("CACHED_VERIFIED",args.region,flush=True);return
  if (folder/"FAILED.json").exists():raise RuntimeError("Inspect previous failed region before explicit recovery")
  pcpath=BASE/"sources"/rec["accession"]/"COMPLETE.json"
  icpath=BASE/"imaging"/ph["family"]/ph["idp"]/"COMPLETE.json"
  pc=json.loads(pcpath.read_text());ic=json.loads(icpath.read_text())
  assert pc["complete"] and ic["complete"]
  assert pc["plan_sha256"]==sha(BASE/"PLAN.json") and ic["plan_sha256"]==sha(BASE/"imaging/PLAN.json")
  xp=pcpath.parent/args.region/"protein_regional_statistics.csv.gz"
  assert sha(xp)==pc["outputs_sha256"][str(xp.relative_to(ROOT))]
  wid=f'chr{rec["chromosome"]}_center{rec["center_grch37"]}'
  wr=ic["windows"][wid];yp=ROOT/wr["output"];assert sha(yp)==wr["sha256"]
  dependencies=[Path(__file__),Path(__file__).with_name("all_target_coloc.R"),Path(__file__).with_name("prepare_new_scan_coloc.py"),
   Path(__file__).with_name("acquire_all_target_protein_regions.py"),ROOT/"reference/coloc/claudia.R",
   BASE/"PLAN.json",BASE/"imaging/PLAN.json",pcpath,icpath,xp,yp]
  inputs={str(p.relative_to(ROOT)):sha(p) for p in dependencies}
  run=dict(pid=os.getpid(),status="running",complete=False,region_id=args.region,started_at_unix=time.time(),input_sha256=inputs)
  dump(folder/"RUN.json",run)
  try:
   audit=prepare_frames(rec,read_frame(xp),read_frame(yp),folder,source["metadata"]["sample_size"])
   audit["input_sha256"]=inputs;dump(folder/"INPUT_AUDIT.json",audit)
   cmd=["/public/home/mengxl/dzy/envs/rcoloc/bin/Rscript",str(Path(__file__).with_name("all_target_coloc.R")),str(ROOT),str(folder)]
   with (folder/"R.log").open("w") as log:subprocess.run(cmd,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
   statuses=validate_outputs(folder)
   for p,h in inputs.items():assert sha(ROOT/p)==h,"Input changed during ABF calculation"
   outputs={p.name:sha(p) for p in folder.iterdir() if p.suffix in [".gz",".csv"] or p.name in ["INPUT_AUDIT.json","R_sessionInfo.txt"]}
   run.update(status="complete",complete=True,window_status=statuses,outputs_sha256=outputs,elapsed_seconds=time.time()-run["started_at_unix"])
   dump(folder/"COMPLETE.json",run);dump(folder/"RUN.json",run)
   print(json.dumps(dict(region=args.region,window_status=statuses,shared_snps=audit["shared_snps"])),flush=True)
  except BaseException as exc:
   run.update(status="failed",error=repr(exc),traceback=traceback.format_exc());dump(folder/"FAILED.json",run);dump(folder/"RUN.json",run);raise
if __name__=="__main__":main()
