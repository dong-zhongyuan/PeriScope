"""Prepare one native imaging request per phenotype and shared reference windows.

This executes the frozen all-assay regional plan. It does not rescan, rank, or
filter MR results. Multiple protein sources targeting the same IDP share one
download, and all windows for that IDP are extracted in a single pass.
"""
from pathlib import Path
from collections import defaultdict,Counter
import csv,hashlib,json
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"data/all_target_followup_20261118/regional_all"
OUT=BASE/"imaging"
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def main():
 OUT.mkdir(exist_ok=True)
 if (OUT/"PLAN.json").exists():raise RuntimeError("Inspect existing imaging plan")
 plan=json.loads((BASE/"PLAN.json").read_text())
 bp=ROOT/"data/big40_atlas/phenotypes_3935.csv";qp=ROOT/"data/oxford_qsm/qsm_28_phenotypes.csv"
 bm={x["phenotype"]:x for x in read(bp)};qm={x["phenotype"]:x for x in read(qp)}
 phenotypes={};windows={}
 source_rank={acc:i for i,acc in enumerate(plan["source_order"])}
 for rid,reg in plan["regions"].items():
  phenotype=reg["phenotype"];family=phenotype.split(":")[0]
  if phenotype not in phenotypes:
   meta=(bm if family=="BIG40" else qm)[phenotype]
   phenotypes[phenotype]=dict(phenotype=phenotype,family=family,idp=phenotype.split(":")[-1],
    metadata=meta,regions=[],windows=[],model_linked=False,source_priority=len(source_rank))
  rec=phenotypes[phenotype];rec["regions"].append(rid)
  rec["model_linked"]|=reg["model_linked"];rec["source_priority"]=min(rec["source_priority"],source_rank[reg["accession"]])
  wid=f'chr{reg["chromosome"]}_center{reg["center_grch37"]}'
  if wid not in rec["windows"]:rec["windows"].append(wid)
  if wid not in windows:
   windows[wid]={k:reg[k] for k in ["chromosome","center_grch37","start_grch37","end_grch37"]}
   windows[wid].update(window_id=wid,regions=[],BIG40_phenotypes=[],QSM_phenotypes=[])
  windows[wid]["regions"].append(rid)
  key="BIG40_phenotypes" if family=="BIG40" else "QSM_phenotypes"
  if phenotype not in windows[wid][key]:windows[wid][key].append(phenotype)
 for phenotype,rec in phenotypes.items():
  if rec["family"]=="BIG40":
   rec["native_url"]=f'https://open.oxcin.ox.ac.uk/ukbiobank/big40/release2/stats/{rec["idp"]}.txt.gz'
   rec["expected_source_rows"]=17103079;known=[]
   previous=ROOT/"data/big40_atlas/native"/(rec["idp"]+".json")
   if previous.exists():
    p=json.loads(previous.read_text());assert p["complete"]
    if "source_sha256" in p["source_receipt"]:
     known.append(dict(receipt=str(previous.relative_to(ROOT)),receipt_sha256=sha(previous),source_sha256=p["source_receipt"]["source_sha256"]))
   cached_receipt=ROOT/"data/big40/download/stats"/(rec["idp"]+".txt.gz.download.json")
   if cached_receipt.exists():
    p=json.loads(cached_receipt.read_text())
    assert p["validation"]["gzip_crc_and_eof_verified"] and p["validation"]["source_rows"]==17103079
    known.append(dict(receipt=str(cached_receipt.relative_to(ROOT)),receipt_sha256=sha(cached_receipt),source_sha256=p["sha256"]))
   for p in (ROOT/"data/new_scan_regional_20261117").glob("GCST*/IMAGING_REGION.json"):
    d=json.loads(p.read_text())
    if d["phenotype_record"]["phenotype"]==phenotype:
     known.append(dict(receipt=str(p.relative_to(ROOT)),receipt_sha256=sha(p),source_sha256=d["source_sha256"]))
   assert len({x["source_sha256"] for x in known})<=1
   rec["prior_full_native_source_identities"]=known
   rec["existing_full_payloads"]=[str(p.relative_to(ROOT)) for p in (ROOT/"data/big40/download/stats").rglob(rec["idp"]+".txt.gz")]
  else:
   rec["native_url"]=rec["metadata"]["source_y"];rec["expected_source_rows"]=20193843
   rec["source_member"]=rec["idp"]+".txt"
 order=sorted(phenotypes,key=lambda p:(not phenotypes[p]["model_linked"],phenotypes[p]["source_priority"],p))
 data=dict(complete=True,scope="All 1375 frozen regional requests, grouped by native source and window, without additional result cutoff",
  regions=plan["regions"],phenotypes=phenotypes,windows=windows,phenotype_order=order,
  counts=dict(regions=len(plan["regions"]),phenotypes=len(phenotypes),windows=len(windows),
   source_families=dict(Counter(x["family"] for x in phenotypes.values())),
   BIG40_reference_windows=sum(bool(x["BIG40_phenotypes"]) for x in windows.values())),
  source_order_reason="Follow the existing protein-source queue and prioritize model-linked sources; every planned source retained",
  input_sha256={str(p.relative_to(ROOT)):sha(p) for p in [BASE/"PLAN.json",bp,qp,Path(__file__)]})
 (OUT/"PLAN.json").write_text(json.dumps(data,indent=2));print(json.dumps(data["counts"]),flush=True)
if __name__=="__main__":main()
