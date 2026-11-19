"""Make a deduplicated execution plan for the already defined all-assay regional worklist."""
from pathlib import Path
from collections import defaultdict,Counter
import csv,hashlib,json
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"data/all_target_followup_20261118"
OUT=BASE/"regional_all";WORK=ROOT/"papers/01_mr_imaging_scan/results/all_target_coverage_20261118"
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
def read(p):return list(csv.DictReader(Path(p).open()))
def main():
 OUT.mkdir(exist_ok=True)
 if (OUT/"PLAN.json").exists():raise RuntimeError("Inspect existing plan before replacement")
 wp=WORK/"all_assays_first_pair_followup_worklist.csv";work=read(wp)
 lp=ROOT/"papers/03_world_model_integration/results/cross_domain_inputs_20261117/all_model_target_assay_links.csv"
 model={r["exposure_id"] for r in read(lp) if r["exposure_id"]}
 peakpath=ROOT/"data/big40_peak_catalog/all_10889_source_peaks_with_frequencies.csv"
 coveragepath=ROOT/"data/big40_peak_catalog/all_3935_phenotype_coverage.csv"
 peakverification=ROOT/"data/big40_peak_catalog/VERIFICATION.json"
 verification=json.loads(peakverification.read_text())
 assert verification["complete"] and verification["source_rows"]==10889
 for f,h in verification["output_sha256"].items():assert sha(coveragepath.parent/f)==h
 big40covered={r["phenotype"] for r in read(coveragepath)}
 assert len(big40covered)==3935
 peaks=defaultdict(list)
 for row in read(peakpath):peaks[row["phenotype"]].append(row)
 sources={};regions={};links=[]
 for row in work:
  if row["status"] not in ["pending_region_and_reverse","already_completed_first_batch"]:
   links.append({**row,"region_id":"","execution_status":row["status"]});continue
  acc=row["accession"];chrom=str(int(row["chromosome"]));center=int(row["center_grch37"])
  key=[acc,chrom,center,row["phenotype"]]
  rid=acc+"_"+hashlib.sha256(json.dumps(key).encode()).hexdigest()[:12]
  if rid not in regions:
   regions[rid]=dict(region_id=rid,accession=acc,gene=row["gene"],chromosome=chrom,center_grch37=center,
    start_grch37=max(1,center-1000000),end_grch37=center+1000000,phenotype=row["phenotype"],
    phenotype_description=row["phenotype_description"],exposure_ids=[],policies=[],
    original_worklist_status=row["status"],model_linked=False,
    reverse_peak_source_status="official_BIG40_peak_table" if row["phenotype"] in peaks else ("no_entry_in_complete_official_BIG40_peak_table" if row["phenotype"] in big40covered else "published_peak_source_not_yet_resolved"))
  reg=regions[rid]
  if row["exposure_id"] not in reg["exposure_ids"]:reg["exposure_ids"].append(row["exposure_id"])
  if row["policy"] not in reg["policies"]:reg["policies"].append(row["policy"])
  reg["model_linked"] |= row["exposure_id"] in model
  if acc not in sources:
   path=BASE/"sources"/acc/"SOURCE.json";meta=json.loads(path.read_text())
   assert meta["complete"] and meta["original_build"]=="GRCh37"
   for f,h in meta["input_sha256"].items():assert sha(ROOT/f)==h
   sources[acc]=dict(accession=acc,metadata=meta,metadata_sha256=sha(path),
    regions=[],model_linked=False,peak_rows=[],peak_phenotypes=[])
  source=sources[acc]
  if rid not in source["regions"]:source["regions"].append(rid)
  source["model_linked"] |= row["exposure_id"] in model
  phenotype=row["phenotype"]
  if phenotype not in source["peak_phenotypes"]:
   source["peak_rows"].extend(peaks.get(phenotype,[]));source["peak_phenotypes"].append(phenotype)
  links.append({**row,"region_id":rid,"execution_status":"planned"})
 assert len(links)==len(work)
 plan=dict(complete=True,scope="Execution of existing one-best-imaging-pair-per-policy all-assay worklist; no additional result cutoff",
  window_bp=1000000,sources=sources,regions=regions,
  source_order=sorted(sources,key=lambda a:(not sources[a]["model_linked"],a)),
  counts=dict(worklist_rows=len(work),sources=len(sources),regions=len(regions),
   model_linked_sources=sum(s["model_linked"] for s in sources.values()),
   phenotypes=len({r["phenotype"] for r in regions.values()}),
   source_family=dict(Counter(r["phenotype"].split(":")[0] for r in regions.values())),
   unresolved_reverse_peak_regions=sum(r["reverse_peak_source_status"]=="published_peak_source_not_yet_resolved" for r in regions.values()),
   no_official_BIG40_peak_regions=sum(r["reverse_peak_source_status"]=="no_entry_in_complete_official_BIG40_peak_table" for r in regions.values())),
  source_order_reason="Schedule model-covered sources first, then all other planned sources; none excluded by effect or score.",
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [wp,lp,peakpath,coveragepath,peakverification,Path(__file__)]})
 (OUT/"PLAN.json").write_text(json.dumps(plan,indent=2))
 with (OUT/"all_worklist_execution_links.csv").open("w") as f:
  w=csv.DictWriter(f,fieldnames=list(links[0]),lineterminator="\n");w.writeheader();w.writerows(links)
 print(json.dumps(plan["counts"]),flush=True)
if __name__=="__main__":main()
