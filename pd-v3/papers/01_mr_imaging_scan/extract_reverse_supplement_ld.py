"""Acquire only additional QSM and chrX reference SNPs under the frozen ploidy specification."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import fcntl,hashlib,json,os,re,time,traceback
import requests,numpy as np,pandas as pd,pysam
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/"data/all_target_followup_20261118/regional_all"
SUP=BASE/"reverse_supplement";OUT=SUP/"reference_LD";OLD=ROOT/"data/new_scan_regional_20261117/reverse_LD"
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,j):
 t=p.with_name(p.name+".pending");t.write_text(json.dumps(j,indent=2)+"\n");t.replace(p)
def identity(url):
 r=requests.head(url,timeout=45);r.raise_for_status();d={k:r.headers.get(k) for k in ["ETag","Last-Modified","Content-Length"]}
 assert d["Content-Length"] and (d["ETag"] or d["Last-Modified"]);return d
def main():
 OUT.mkdir(exist_ok=True)
 with (OUT/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);assert not (OUT/"RUN.json").exists()
  plan=json.loads((SUP/"PLAN.json").read_text());ph=sha(SUP/"PLAN.json");code=sha(__file__)
  table=pd.DataFrame(plan["reference_requests"]);table["chromosome"]=table.chromosome.map(lambda x:23 if str(x)=="X" else int(x))
  table["position_grch37"]=table.position_grch37.astype(int)
  assert len(table)==79 and table.variant_key.is_unique
  pp=OLD/"integrated_call_samples_v3.20130502.ALL.panel"
  panel=pd.read_csv(pp,sep="\t",usecols=["sample","super_pop","gender"]);panel=panel[panel.super_pop.eq("EUR")]
  samples=panel["sample"].tolist();gender=panel.gender.tolist();assert len(samples)==len(set(samples))==503
  assert set(gender)=={"male","female"}
  names=re.findall(r'href="([^"]+)"',(OLD/"source_listing.html").read_text())
  state=dict(pid=os.getpid(),status="running",complete=False,plan_sha256=ph,code_sha256=code,
   requested_SNVs=len(table),completed_chromosomes=[],failures={},started_at_unix=time.time(),method=plan["method"])
  dump(OUT/"RUN.json",state)
  def one(c):
   folder=OUT/f"chr{c}";folder.mkdir(exist_ok=True);label="X" if c==23 else str(c)
   names_c=[n for n in names if n.startswith(f"ALL.chr{label}.") and n.endswith(".genotypes.vcf.gz")];assert len(names_c)==1
   name=names_c[0];url=URL+name;initial=identity(url);index=folder/(name+".tbi")
   # Reuse exact source indexes when available; published HTTP identity must still match.
   candidates=[OLD/f"chr{c}",BASE/"reverse_LD"/f"chr{c}"];copied=False
   for prior in candidates:
    mp=prior/"MANIFEST.json"
    if mp.exists():
     m=json.loads(mp.read_text());assert m["complete"] and m["url"]==url and m["http_identity"]==initial
     assert sha(prior/index.name)==m["index_sha256"];index.write_bytes((prior/index.name).read_bytes());copied=True;break
   if not copied:
    response=requests.get(url+".tbi",timeout=60);response.raise_for_status();index.write_bytes(response.content)
   rows=[];doses=[];raws=[];ploidies=[]
   wanted=table[table.chromosome.eq(c)].sort_values(["position_grch37","variant_key"])
   with pysam.VariantFile(url,index_filename=str(index)) as source:
    source.subset_samples(samples);assert set(source.header.samples)==set(samples)
    for entry in wanted.itertuples():
     matches=[]
     for v in source.fetch(label,entry.position_grch37-1,entry.position_grch37):
      if v.pos==entry.position_grch37 and v.alts and len(v.alts)==1 and sorted([v.ref,v.alts[0]])==sorted([entry.effect_allele,entry.other_allele]):
       if len(v.ref)==len(v.alts[0])==1:matches.append(v)
     assert len(matches)<=1,"Ambiguous reference identity"
     ispar=c==23 and (entry.position_grch37<=2699520 or entry.position_grch37>=154931044)
     row=dict(variant_key=entry.variant_key,chromosome=int(c),position_grch37=int(entry.position_grch37),effect_allele=entry.effect_allele,
      other_allele=entry.other_allele,is_PAR_GRCh37=ispar,request_reason=entry.request_reason)
     if not matches:row.update(status="not_matched_biallelic_SNV_in_reference",dosage_row=-1,ref="",alt="")
     else:
      v=matches[0];g=[];raw=[];pl=[]
      for sample,sex in zip(samples,gender):
       gt=v.samples[sample]["GT"];assert len(gt) in [1,2]
       raw.append([(-1 if a is None else int(a)) for a in gt]+[-2]*(2-len(gt)));pl.append(len(gt))
       if any(a is None for a in gt):g.append(-1);continue
       assert all(a in [0,1] for a in gt)
       biological=1 if c==23 and not ispar and sex=="male" else 2
       if biological==1:
        assert len(gt)==1 or (len(gt)==2 and gt[0]==gt[1]),"Heterozygous male non-PAR GT requires inspection"
        alt=float(gt[0])
       else:
        assert len(gt)==2,"Unexpected haploid call for autosome/female/PAR requires inspection"
        alt=sum(gt)/2
       g.append(int(round(2*(1-alt if entry.effect_allele==v.ref else alt))))
      a=np.array(g,dtype=np.int8);assert np.isin(a,[-1,0,1,2]).all()
      row.update(status="matched",ref=v.ref,alt=v.alts[0],dosage_row=len(doses),
       missing_calls=int((a<0).sum()),male_samples=int(np.sum(np.array(gender)=="male")),
       female_samples=int(np.sum(np.array(gender)=="female")),observed_haploid_GT=int(np.sum(np.array(pl)==1)),
       observed_diploid_GT=int(np.sum(np.array(pl)==2)))
      doses.append(a);raws.append(np.array(raw,dtype=np.int8));ploidies.append(np.array(pl,dtype=np.int8))
     rows.append(row);dump(folder/"PROGRESS.json",dict(pid=os.getpid(),chromosome=int(c),processed=len(rows),requested=len(wanted)))
   assert identity(url)==initial
   cov=pd.DataFrame(rows);cov.to_csv(folder/"coverage.csv",index=False)
   np.savez_compressed(folder/"effect_allele_dosages.npz",dosages=np.stack(doses) if doses else np.empty((0,503),dtype=np.int8),
    raw_GT=np.stack(raws) if raws else np.empty((0,503,2),dtype=np.int8),
    source_GT_ploidy=np.stack(ploidies) if ploidies else np.empty((0,503),dtype=np.int8),
    samples=np.array(samples),gender=np.array(gender),variant_keys=cov.loc[cov.status.eq("matched"),"variant_key"].to_numpy(str))
   dump(folder/"MANIFEST.json",dict(complete=True,chromosome=int(c),url=url,http_identity=initial,index_sha256=sha(index),
    panel_sha256=sha(pp),plan_sha256=ph,code_sha256=code,requested=len(wanted),matched=len(doses),
    source_reference="1000 Genomes Phase3 EUR503 GRCh37",whole_remote_file_checksum_verified=False,
    cohort_matched=False,raw_GT_retained=True,method=plan["method"],
    outputs_sha256={n:sha(folder/n) for n in ["coverage.csv","effect_allele_dosages.npz"]}))
   print("SUPPLEMENT_CHR",c,len(doses),len(wanted),flush=True);return int(c)
  with ThreadPoolExecutor(max_workers=2) as pool:
   futures={pool.submit(one,int(c)):int(c) for c in sorted(table.chromosome.unique(),reverse=True)}
   for future in as_completed(futures):
    c=futures[future]
    try:state["completed_chromosomes"].append(future.result())
    except Exception as exc:
     state["failures"][str(c)]=repr(exc);dump(OUT/f"chr{c}"/"FAILED.json",dict(error=repr(exc),traceback=traceback.format_exc()))
    dump(OUT/"RUN.json",state)
  assert sha(SUP/"PLAN.json")==ph and sha(__file__)==code
  state.update(status="complete" if not state["failures"] else "finished_with_chromosome_errors",complete=not state["failures"],finished_at_unix=time.time())
  dump(OUT/"RUN.json",state)
if __name__=="__main__":main()
