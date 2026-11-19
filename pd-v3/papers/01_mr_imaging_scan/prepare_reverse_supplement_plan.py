"""Separate immutable source supplement for the existing all-target reverse worklist."""
from pathlib import Path
import csv,hashlib,json,datetime
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/"data/all_target_followup_20261118/regional_all"
QSM=ROOT/"data/all_target_followup_20261118/qsm_published_peaks";OUT=BASE/"reverse_supplement"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rows(p):
 with Path(p).open() as f:return list(csv.DictReader(f))
def snv(r):return all(len(r[k])==1 and r[k] in "ACGT" for k in ["effect_allele","other_allele"])
def main():
 OUT.mkdir(exist_ok=True);path=OUT/"PLAN.json";assert not path.exists()
 plan=json.loads((BASE/"PLAN.json").read_text());qsm=rows(QSM/"all_requested_QSM_author_peaks.csv")
 freq=json.loads((QSM/"FREQUENCY_AUDIT.json").read_text());assert freq["complete"]
 full=rows(BASE/"reverse_LD/all_unique_peak_coverage.csv");old={r["variant_key"]:r for r in full}
 coverage=[];sources={}
 for rec in plan["regions"].values():
  if not rec["phenotype"].startswith("OxfordQSM2022:"):continue
  peaks=[x for x in qsm if x["phenotype"]==rec["phenotype"]]
  coverage.append(dict(region_id=rec["region_id"],accession=rec["accession"],phenotype=rec["phenotype"],
   source_peak_rows=len(peaks),status="published_Table2_available" if peaks else "no_entry_in_complete_author_Table2"))
  acc=rec["accession"]
  if acc not in sources:sources[acc]=dict(accession=acc,regions=[],peak_rows=[])
  sources[acc]["regions"].append(rec["region_id"]);sources[acc]["peak_rows"].extend(peaks)
 for source in sources.values():
  bykey={(p["phenotype"],p["variant_key"]):p for p in source["peak_rows"]};source["peak_rows"]=list(bykey.values())
 requests={}
 for row in full:
  if row["reference_request_status"]=="pending_X_ploidy_method_review":
   requests[row["variant_key"]]=dict(row,request_reason="complete_original_BIG40_X_reference")
 for row in qsm:
  if snv(row) and row["variant_key"] not in old:requests[row["variant_key"]]=dict(row,request_reason="additional_published_QSM_SNV")
 unique_qsm={p["variant_key"]:p for p in qsm}
 qcov=[]
 for key,row in unique_qsm.items():
  route="unsupported_non_biallelic_SNV" if not snv(row) else ("existing_base_reference_queue" if key in old else "supplemental_reference_queue")
  qcov.append(dict(variant_key=key,reference_route=route))
 method=dict(predeclared_before_supplemental_reference_or_outcome_read=True,
  X_dosage="2 * effect-allele count / biological ploidy; males outside PAR 0/2, females and PAR diploid 0/1/2",
  X_mode_reference="https://www.cog-genomics.org/plink/1.9/ld",
  X_mode="PLINK --ld-xchr 2 dosage convention; calculation is NumPy, not claimed to be a PLINK run",
  PAR_GRCh37_boundary_reference="https://www.cog-genomics.org/plink/1.9/data#split_x",
  PAR1_last_position=2699520,PAR2_first_position=154931044,
  sex_source="Official 1000 Genomes Phase3 panel gender field, EUR503",
  raw_GT_retained=True,missing_calls="Unfilled; pairwise complete existing LD recipe",
  unexpected_ploidy="Explicit source exception requiring inspection, not silent genotype reassignment",
  doubled_haploid_representation="Homozygous diploid male non-PAR GTs accepted as the same haploid allele, source ploidy recorded",
  autosomal_dosage="Unchanged 0/1/2 effect-allele count",
  LD_selection="Existing imaging discovery P<=5e-8, greedy discovery-P/key order, r2>.001 within10Mb; no outcome-based choice",
  unchanged_original_autosomal_jobs=True)
 deps=[BASE/"PLAN.json",BASE/"reverse_LD/all_unique_peak_coverage.csv",QSM/"AUTHOR_PEAK_MANIFEST.json",QSM/"FREQUENCY_AUDIT.json",QSM/"all_requested_QSM_author_peaks.csv",Path(__file__)]
 result=dict(complete=True,created_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sources=sources,
  all_QSM_region_coverage=coverage,QSM_unique_variant_routes=qcov,reference_requests=list(requests.values()),method=method,
  counts=dict(QSM_regions=len(coverage),QSM_sources=len(sources),QSM_peak_rows=len(qsm),QSM_unique_variants=len(unique_qsm),
   additional_reference_SNVs=len(requests),additional_QSM_SNVs=sum(r["request_reason"].startswith("additional") for r in requests.values()),
   BIG40_X_SNVs=sum(r["request_reason"].startswith("complete") for r in requests.values())),
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in deps},
  purpose="Fill missing public-source inputs for existing planned regions; no new MR screening layer")
 path.write_text(json.dumps(result,indent=2)+"\n");print(json.dumps(result["counts"]),flush=True)
if __name__=="__main__":main()
