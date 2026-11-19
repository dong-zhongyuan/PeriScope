"""Shared public EUR reference retrieval for the frozen all-target published-peak worklist."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import fcntl,hashlib,json,os,re,sys,time,traceback
import requests,numpy as np,pandas as pd,pysam
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"data/all_target_followup_20261118/regional_all"
OUT=BASE/"reverse_LD";OLD=ROOT/"data/new_scan_regional_20261117/reverse_LD"
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
PANEL="integrated_call_samples_v3.20130502.ALL.panel"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,j):
 t=p.with_name(p.name+".pending");t.write_text(json.dumps(j,indent=2)+"\n");t.replace(p)
def chrom(value):return 23 if str(value)=="X" else int(value)
def identity(url):
 response=requests.head(url,timeout=45);response.raise_for_status()
 result={k:response.headers.get(k) for k in ["ETag","Last-Modified","Content-Length"]}
 assert result["Content-Length"] and (result["ETag"] or result["Last-Modified"])
 return result
def main():
 OUT.mkdir(exist_ok=True)
 with (OUT/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  if (OUT/"RUN.json").exists():raise RuntimeError("Inspect previous run before explicit recovery")
  ph=sha(BASE/"PLAN.json");code=sha(__file__);plan=json.loads((BASE/"PLAN.json").read_text())
  rows=[]
  for accession in plan["source_order"]:
   for row in plan["sources"][accession]["peak_rows"]:rows.append(dict(row,accession=accession))
  peaks=pd.DataFrame(rows);peaks["chromosome"]=peaks.chromosome.map(chrom)
  peaks["position_grch37"]=peaks.position_grch37.astype(int)
  peaks["variant_key"]=peaks.apply(lambda r:f"{r.chromosome}:{r.position_grch37}:"+ "_".join(sorted([r.effect_allele,r.other_allele])),axis=1)
  peaks.to_csv(OUT/"all_source_peak_links.csv",index=False)
  unique=peaks.drop_duplicates("variant_key").sort_values(["chromosome","position_grch37","variant_key"]).copy()
  snv=unique.effect_allele.str.fullmatch("[ACGT]")&unique.other_allele.str.fullmatch("[ACGT]")
  unique["reference_request_status"]=np.where(~snv,"unsupported_non_biallelic_SNV",
   np.where(unique.chromosome.eq(23),"pending_X_ploidy_method_review","planned_autosomal_diploid_EUR"))
  unique.to_csv(OUT/"all_unique_peak_coverage.csv",index=False)
  # X records remain visible and pending. Do not silently reinterpret haploid male GTs as diploid.
  wanted=unique[unique.reference_request_status.eq("planned_autosomal_diploid_EUR")]
  assert (wanted.chromosome.between(1,22)).all()
  oldrun=json.loads((OLD/"RUN.json").read_text());assert oldrun["complete"]
  panelpath=OLD/PANEL;panelhash=sha(panelpath)
  panel=pd.read_csv(panelpath,sep="\t",usecols=["sample","pop","super_pop"])
  samples=panel.loc[panel.super_pop.eq("EUR"),"sample"].tolist();assert len(samples)==len(set(samples))==503
  listing=OLD/"source_listing.html";names=re.findall(r'href="([^"]+)"',listing.read_text())
  state=dict(pid=os.getpid(),status="running",complete=False,started_at_unix=time.time(),plan_sha256=ph,
   code_sha256=code,all_peak_links=len(peaks),unique_requested_variants=len(unique),
   autosomal_SNV_requests=len(wanted),unsupported_non_SNV=int((~snv).sum()),
   X_SNV_pending_method=int((snv&unique.chromosome.eq(23)).sum()),
   chromosome_workers=4,completed_chromosomes=[],failures={},
   reference="1000 Genomes Phase3 EUR503 GRCh37",cohort_matched=False,
   LD_selection_before_outcome_join=True,new_peak_discovery=False,
   requested_tables_sha256={p.name:sha(p) for p in [OUT/"all_source_peak_links.csv",OUT/"all_unique_peak_coverage.csv"]},
   source_listing_sha256=sha(listing),panel_sha256=panelhash)
  dump(OUT/"RUN.json",state)
  def chromosome_job(c):
   folder=OUT/f"chr{c}";folder.mkdir(exist_ok=True)
   entries=wanted[wanted.chromosome.eq(c)]
   choices=[n for n in names if n.startswith(f"ALL.chr{c}.") and n.endswith(".genotypes.vcf.gz")]
   assert len(choices)==1;name=choices[0];url=URL+name;initial=identity(url)
   oldfolder=OLD/f"chr{c}";oldmanifest=oldfolder/"MANIFEST.json";cached={}
   index=folder/(name+".tbi")
   if oldmanifest.exists():
    meta=json.loads(oldmanifest.read_text())
    assert meta["complete"] and meta["panel_sha256"]==panelhash and meta["url"]==url
    assert meta["http_identity"]==initial,"Reference changed since prior LD extraction"
    oldindex=oldfolder/index.name;assert sha(oldindex)==meta["index_sha256"]
    index.write_bytes(oldindex.read_bytes())
    for n,h in meta["outputs_sha256"].items():assert sha(oldfolder/n)==h
    table=pd.read_csv(oldfolder/"coverage.csv",keep_default_na=False)
    with np.load(oldfolder/"effect_allele_dosages.npz",allow_pickle=False) as z:
     assert np.array_equal(z["samples"],np.array(samples));g=z["dosages"].copy()
    for row in table.to_dict("records"):
     cached[row["variant_key"]]=(row,g[int(row["dosage_row"])].copy() if row["status"]=="matched" else None)
   else:
    response=requests.get(url+".tbi",timeout=60);response.raise_for_status();index.write_bytes(response.content)
   records=[];dosages=[];reused=0;fetched=0
   with pysam.VariantFile(url,index_filename=str(index)) as source:
    source.subset_samples(samples);assert set(source.header.samples)==set(samples)
    for entry in entries.itertuples():
     base=dict(variant_key=entry.variant_key,chromosome=int(c),position_grch37=int(entry.position_grch37),
      effect_allele=entry.effect_allele,other_allele=entry.other_allele)
     if entry.variant_key in cached:
      row,a=cached[entry.variant_key]
      if a is not None:
       if row["effect_allele"]!=entry.effect_allele:
        a=a.copy();ok=a>=0;a[ok]=2-a[ok]
       base.update(status="matched",ref=row["ref"],alt=row["alt"])
      else:base.update(status="not_matched_biallelic_SNV_in_reference",ref="",alt="")
      reused+=1;route="previous_verified_reference_subset"
     else:
      candidates=[]
      for v in source.fetch(str(c),entry.position_grch37-1,entry.position_grch37):
       if v.pos!=entry.position_grch37 or not v.alts or len(v.alts)!=1:continue
       if len(v.ref)!=1 or len(v.alts[0])!=1:continue
       if sorted([v.ref,v.alts[0]])==sorted([entry.effect_allele,entry.other_allele]):candidates.append(v)
      assert len(candidates)<=1,"Ambiguous reference allele identity"
      a=None
      if candidates:
       v=candidates[0];g=[]
       for sample in samples:
        gt=v.samples[sample]["GT"]
        g.append(-1 if len(gt)!=2 or any(a is None for a in gt) else sum(gt))
       a=np.array(g,dtype=np.int8);ok=a>=0
       if entry.effect_allele==v.ref:a[ok]=2-a[ok]
       base.update(status="matched",ref=v.ref,alt=v.alts[0])
      else:base.update(status="not_matched_biallelic_SNV_in_reference",ref="",alt="")
      fetched+=1;route="indexed_official_remote_VCF"
     base["route"]=route
     if a is None:base.update(dosage_row=-1,missing_calls=None,effect_allele_frequency=None)
     else:
      assert a.shape==(503,) and np.isin(a,[-1,0,1,2]).all()
      ok=a>=0;base.update(dosage_row=len(dosages),missing_calls=int((~ok).sum()),
       effect_allele_frequency=float(a[ok].mean()/2) if ok.any() else None)
      dosages.append(a)
     records.append(base)
     dump(folder/"PROGRESS.json",dict(chromosome=int(c),pid=os.getpid(),requested=len(entries),processed=len(records),reused=reused,fetched=fetched))
   assert identity(url)==initial
   table=pd.DataFrame(records);table.to_csv(folder/"coverage.csv",index=False)
   g=np.stack(dosages) if dosages else np.empty((0,503),dtype=np.int8)
   np.savez_compressed(folder/"effect_allele_dosages.npz",dosages=g,samples=np.array(samples),
    variant_keys=table.loc[table.status.eq("matched"),"variant_key"].to_numpy(str))
   report=dict(complete=True,chromosome=int(c),url=url,http_identity=initial,index_sha256=sha(index),panel_sha256=panelhash,
    plan_sha256=ph,code_sha256=code,requested=len(entries),matched=len(dosages),reused_previous_variants=reused,
    remote_queries=fetched,source_reference="1000 Genomes Phase3 EUR503 GRCh37",cohort_matched=False,
    whole_remote_file_checksum_verified=False,
    outputs_sha256={n:sha(folder/n) for n in ["coverage.csv","effect_allele_dosages.npz"]})
   dump(folder/"MANIFEST.json",report);return int(c)
  with ThreadPoolExecutor(max_workers=4) as pool:
   futures={pool.submit(chromosome_job,int(c)):int(c) for c in sorted(wanted.chromosome.unique())}
   for future in as_completed(futures):
    c=futures[future]
    try:state["completed_chromosomes"].append(future.result());print("REFERENCE_CHR_COMPLETE",c,flush=True)
    except Exception as exc:
     state["failures"][str(c)]=repr(exc);dump(OUT/f"chr{c}"/"FAILED.json",dict(error=repr(exc),traceback=traceback.format_exc()))
    dump(OUT/"RUN.json",state)
  assert sha(BASE/"PLAN.json")==ph and sha(__file__)==code
  state.update(complete=not state["failures"],status="autosomal_reference_complete" if not state["failures"] else "finished_with_chromosome_errors",
   finished_at_unix=time.time(),pysam_version=pysam.__version__,
   full_reverse_MR_complete=False,X_SNV_pending_method_review=state["X_SNV_pending_method"])
  dump(OUT/"RUN.json",state)
if __name__=="__main__":main()
