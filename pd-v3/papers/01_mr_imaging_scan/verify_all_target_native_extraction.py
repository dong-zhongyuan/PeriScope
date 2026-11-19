"""Compare batched regional extraction with the three earlier independently acquired regions."""
from pathlib import Path
import hashlib,json
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"data/all_target_followup_20261118/regional_all/imaging"
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def read(p):return pd.read_csv(p,dtype={"chromosome":str},float_precision="round_trip",keep_default_na=False)
def main():
 plan=json.loads((BASE/"PLAN.json").read_text());reference=json.loads((BASE/"variant_reference/MANIFEST.json").read_text())
 assert reference["complete"] and reference["source_rows"]==17103079
 comparisons=[]
 for receipt in (ROOT/"data/new_scan_regional_20261117").glob("GCST*/IMAGING_REGION.json"):
  previous=json.loads(receipt.read_text());assert previous["complete"]
  acq=json.loads((receipt.parent/"ACQUISITION.json").read_text())
  region=acq["download"]["validation"]["region"];wid=f'chr{region["chrom"]}_center{region["center"]}'
  old=receipt.parent/"complete_regional_variant_frequencies.csv.gz";new=ROOT/reference["outputs"][wid]["output"]
  assert sha(new)==reference["outputs"][wid]["sha256"]
  a,b=read(old),read(new)
  pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=True)
  comparisons.append(dict(kind="complete_variant_frequency_window",accession=receipt.parent.name,rows=len(a),
   reference=str(old.relative_to(ROOT)),reference_sha256=sha(old),candidate=str(new.relative_to(ROOT)),candidate_sha256=sha(new),all_fields_exact=True))
  idp=previous["phenotype_record"]["idp"]
  done=BASE/"BIG40"/idp/"COMPLETE.json"
  if done.exists():
   d=json.loads(done.read_text());assert d["complete"] and d["source_sha256"]==previous["source_sha256"]
   old=receipt.parent/"imaging_region_GRCh37.csv.gz";new=ROOT/d["windows"][wid]["output"]
   assert sha(old)==previous["output_sha256"] and sha(new)==d["windows"][wid]["sha256"]
   a,b=read(old),read(new)
   pd.testing.assert_frame_equal(a,b[a.columns],check_dtype=False,check_exact=True)
   comparisons.append(dict(kind="complete_native_imaging_region",accession=receipt.parent.name,phenotype=d["phenotype"],
    rows=len(a),reference=str(old.relative_to(ROOT)),reference_sha256=sha(old),candidate=str(new.relative_to(ROOT)),
    candidate_sha256=sha(new),all_original_fields_exact=True,extra_field="source_data_row only"))
 assert len([x for x in comparisons if x["kind"]=="complete_variant_frequency_window"])==3
 assert any(x["kind"]=="complete_native_imaging_region" and x["phenotype"].endswith(":1441") for x in comparisons)
 result=dict(complete=True,plan_sha256=sha(BASE/"PLAN.json"),reference_manifest_sha256=sha(BASE/"variant_reference/MANIFEST.json"),
  comparisons=comparisons,script_sha256=sha(__file__))
 (BASE/"NATIVE_EXTRACTION_PILOT_VERIFICATION.json").write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=="__main__":main()
