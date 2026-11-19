"""Resolve frozen source/LD dependencies for one existing all-target reverse region."""
from pathlib import Path
import reverse_followup_core as core
ROOT=core.ROOT;BASE=ROOT/"data/all_target_followup_20261118/regional_all";SUP=BASE/"reverse_supplement"
def source_peaks(rec,plan,supp):
 if rec["phenotype"].startswith("BIG40:"):rows=plan["sources"][rec["accession"]]["peak_rows"]
 else:rows=supp["sources"][rec["accession"]]["peak_rows"]
 return [dict(x) for x in rows if x["phenotype"]==rec["phenotype"]]
def snv(row):return all(len(row[k])==1 and row[k] in "ACGT" for k in ["effect_allele","other_allele"])
def requirements(rec,plan,supp):
 peaks=source_peaks(rec,plan,supp);folders=set();source=None
 additional={x["variant_key"] for x in supp["reference_requests"]}
 if peaks:
  source=(BASE/"sources"/rec["accession"]/"COMPLETE.json") if rec["phenotype"].startswith("BIG40:") else (SUP/"protein_outcomes"/rec["accession"]/"COMPLETE.json")
  for row in peaks:
   if not snv(row):continue
   chrom=core.chromosome(row["chromosome"]);key=core.key(chrom,row["position_grch37"],row["effect_allele"],row["other_allele"])
   folder=(SUP/"reference_LD" if key in additional else BASE/"reverse_LD")/f"chr{chrom}";folders.add(folder)
 paths=[] if source is None else [source]
 for folder in sorted(folders):paths.extend([folder/"MANIFEST.json",folder/"INDEPENDENT_VERIFICATION.json"])
 return peaks,source,sorted(folders),paths
