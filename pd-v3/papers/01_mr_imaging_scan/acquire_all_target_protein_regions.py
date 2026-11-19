"""Acquire exact AGES author summaries and extract the existing all-target worklist.

Source extraction is provisional until the full author MD5 has passed. Formatted
bytes are accepted only through the frozen, reversible original-byte bridge.
"""
from pathlib import Path
from collections import defaultdict
import argparse,csv,fcntl,gzip,hashlib,io,json,math,os,sys,time,traceback
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(Path(__file__).parent))
import reconstruct_ages_formatted as bridge
import s09_acquire_qsm as downloader
BASE=ROOT/"data/all_target_followup_20261118/regional_all"
PLAN=BASE/"PLAN.json"
COLUMNS=[x.decode() for x in bridge.COLUMNS]
FIELDS=COLUMNS+["source_row","gene","assay","accession","build","n_metadata","n_scope"]
class OriginalIdentityMismatch(ValueError):pass
def sha(p,kind="sha256"):return bridge.sha(p,kind)
def dump(p,d):
 p=Path(p);q=p.with_name(p.name+".pending");q.write_text(json.dumps(d,indent=2)+"\n");q.replace(p)
def normchr(v):
 v=str(v).removeprefix("chr").lstrip("0");return "X" if v in ["23","X"] else v
def key(c,p,a,b):return normchr(c),int(p),tuple(sorted([a,b]))
def source_lines(path,formatted):
 if formatted:yield from bridge.canonical_lines(path)
 else:
  with path.open("rb") as f:
   for line in f:yield line,False
def write_table(path,rows,fields):
 with path.open("wb") as raw:
  with gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0) as gz:
   with io.TextIOWrapper(gz,encoding="utf-8",newline="") as out:
    w=csv.DictWriter(out,fieldnames=fields,lineterminator="\n")
    w.writeheader();w.writerows(rows)
def extract(path,formatted,source,plan):
 regions=[plan["regions"][rid] for rid in source["regions"]]
 bychrom=defaultdict(list)
 for reg in regions:bychrom[normchr(reg["chromosome"])].append(reg)
 pk={key(r["chromosome"],r["position_grch37"],r["effect_allele"],r["other_allele"]) for r in source["peak_rows"]}
 wanted={r["region_id"]:[] for r in regions};reverse=[]
 meta=source["metadata"];md5=hashlib.md5();strong=hashlib.sha256()
 size=0;corrected=0;rows=0
 for number,(line,change) in enumerate(source_lines(path,formatted)):
  md5.update(line);strong.update(line);size+=len(line);corrected+=int(change)
  fields=line.rstrip(b"\n").split(b"\t")
  if number==0:
   if fields!=bridge.COLUMNS:raise ValueError("Unexpected author schema")
   continue
  rows=number
  if len(fields)!=9 or b"\r" in line or not line.endswith(b"\n"):raise ValueError("Malformed author row")
  chrom=normchr(fields[2].decode("ascii"));pos=int(fields[3])
  selected=[r["region_id"] for r in bychrom.get(chrom,[]) if r["start_grch37"]<=pos<=r["end_grch37"]]
  ispeak=key(chrom,pos,fields[4].decode(),fields[5].decode()) in pk
  if not selected and not ispeak:continue
  decoded=next(csv.reader([line.decode("utf-8")],delimiter="\t",strict=True))
  if len(decoded)!=9:raise ValueError("Malformed quoted author row")
  row=dict(zip(COLUMNS,decoded));values=[float(row[k]) for k in ["p_value","effect_allele_frequency","beta","standard_error"]]
  if not all(math.isfinite(x) for x in values) or not (0<=values[0]<=1 and 0<=values[1]<=1 and values[3]>0):
   raise ValueError("Nonfinite/invalid selected source statistics")
  row.update(source_row=number,gene=meta["gene"],assay=meta["ages_assay"],accession=source["accession"],build="GRCh37",
   n_metadata=meta["sample_size"],n_scope="reported_assay_sample_size_not_per_variant_N")
  for rid in selected:wanted[rid].append(row)
  if ispeak:reverse.append(row)
 if not rows:raise ValueError("Empty source")
 if md5.hexdigest()!=meta["original_published_md5"]:raise OriginalIdentityMismatch("Full reconstructed author MD5 mismatch: "+md5.hexdigest())
 return wanted,reverse,dict(source_rows=rows,reconstructed_original_bytes=size,original_md5=md5.hexdigest(),
  original_sha256=strong.hexdigest(),original_published_md5_verified=True,source_eof_verified=True,
  formatted_gzip_crc_to_EOF_verified=formatted,quote_escape_corrections=corrected,
  input_format="GWAS_Catalog_formatted_reconstructed_author" if formatted else "original_author_TSV")
def acquire(source,formatted,folder,workers):
 meta=source["metadata"];prefix="formatted" if formatted else "original"
 url=meta[prefix+"_url"]
 if not url:raise ValueError("No "+prefix+" source")
 # Same official EBI endpoint over HTTPS, retaining the author URL in metadata.
 url=url.replace("http://ftp.ebi.ac.uk/","https://ftp.ebi.ac.uk/")
 target=folder/meta[prefix+"_name"];receipt=target.with_name(target.name+".download.json")
 expected=meta[prefix+"_published_md5"]
 if target.exists():
  if not receipt.exists():raise RuntimeError("Existing download lacks receipt")
  d=json.loads(receipt.read_text())
  if sha(target)!=d["sha256"] or sha(target,"md5")!=expected:raise RuntimeError("Cached download identity mismatch")
  return target
 downloader.URL=url;downloader.BLOCK=4<<20
 identity=downloader.http_identity()
 part,state=downloader.download_blocks(target,identity,workers,folder/"DOWNLOAD_PROGRESS.json")
 digest=sha(part,"md5")
 if digest!=expected:raise OriginalIdentityMismatch("Downloaded published "+prefix+" MD5 mismatch")
 d=dict(http=identity,sha256=sha(part),md5=digest,published_md5_verified=True,
  size_bytes=part.stat().st_size,completed_blocks=state["completed_blocks"],block_size=downloader.BLOCK)
 dump(receipt,d);part.replace(target);target.with_name(target.name+".part.json").unlink()
 return target
def verify_complete(folder,plan_hash):
 d=json.loads((folder/"COMPLETE.json").read_text())
 if not d["complete"] or d["plan_sha256"]!=plan_hash:raise RuntimeError("Completed extraction has changed plan")
 for f,h in d["outputs_sha256"].items():
  if sha(ROOT/f)!=h:raise RuntimeError("Completed output hash mismatch")
 return d
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--accession",required=True)
 ap.add_argument("--connections",type=int,choices=range(1,9),default=4)
 ap.add_argument("--input-path",type=Path);ap.add_argument("--input-format",choices=["original","formatted"])
 args=ap.parse_args();plan=json.loads(PLAN.read_text());source=plan["sources"][args.accession]
 folder=BASE/"sources"/args.accession;folder.mkdir(parents=True,exist_ok=True)
 with (folder/".lock").open("a+") as lock:
  fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
  ph=sha(PLAN)
  if (folder/"COMPLETE.json").exists():
   d=verify_complete(folder,ph);print(json.dumps({"cached_verified":args.accession,"source_rows":d["source_rows"]}));return
  if (folder/"FAILED.json").exists():raise RuntimeError("Inspect and explicitly archive prior failure before retry")
  meta=source["metadata"]
  mp=ROOT/"data/all_target_followup_20261118/sources"/args.accession/"SOURCE.json"
  assert sha(mp)==source["metadata_sha256"]
  for f,h in meta["input_sha256"].items():assert sha(ROOT/f)==h
  code=[Path(__file__),Path(bridge.__file__),Path(downloader.__file__),Path(downloader._shared.__file__)]
  hashes={str(p.relative_to(ROOT)):sha(p) for p in code}
  started=time.time();run=dict(pid=os.getpid(),status="running",accession=args.accession,plan_sha256=ph,code_sha256=hashes,started_at_unix=started)
  dump(folder/"RUN.json",run)
  try:
   original=ROOT/"data/new_scan_regional_20261117/download"/meta["original_name"]
   fallback=None
   if args.input_path:
    if not args.input_format:raise ValueError("--input-format required")
    path=args.input_path.resolve();formatted=args.input_format=="formatted"
   elif original.exists():path=original;formatted=False
   else:
    formatted=bool(meta.get("formatted_url"));path=acquire(source,formatted,folder,args.connections)
   if formatted and sha(path,"md5")!=meta["formatted_published_md5"]:raise OriginalIdentityMismatch("Formatted source MD5 mismatch")
   try:wanted,reverse,identity=extract(path,formatted,source,plan)
   except (ValueError,gzip.BadGzipFile) as exc:
    if not formatted or args.input_path:raise
    fallback=dict(reason=repr(exc),formatted_path=str(path.relative_to(ROOT)),formatted_sha256=sha(path))
    dump(folder/"FORMATTED_BRIDGE_REJECTED.json",fallback)
    path=acquire(source,False,folder,args.connections);formatted=False
    wanted,reverse,identity=extract(path,False,source,plan)
   outputs={}
   for rid,rows in wanted.items():
    out=folder/rid;out.mkdir(exist_ok=True);p=out/"protein_regional_statistics.csv.gz"
    write_table(p,rows,FIELDS);outputs[str(p.relative_to(ROOT))]=sha(p)
   p=folder/"imaging_peak_outcomes.csv.gz";write_table(p,reverse,FIELDS);outputs[str(p.relative_to(ROOT))]=sha(p)
   p=folder/"all_source_imaging_peaks.csv.gz"
   peakfields=list(source["peak_rows"][0]) if source["peak_rows"] else list(csv.DictReader((ROOT/"data/big40_peak_catalog/all_10889_source_peaks_with_frequencies.csv").open()).fieldnames)
   write_table(p,source["peak_rows"],peakfields);outputs[str(p.relative_to(ROOT))]=sha(p)
   for f,h in hashes.items():assert sha(ROOT/f)==h,"Code changed during extraction"
   assert sha(PLAN)==ph
   receipt=dict(**run,complete=True,**identity,input_path=str(path.relative_to(ROOT)),input_sha256=sha(path),
    regions={rid:{"rows":len(rows),"coverage":"extracted" if rows else "no_source_variants_in_region"} for rid,rows in wanted.items()},
    peak_outcome_rows=len(reverse),source_peak_rows=len(source["peak_rows"]),outputs_sha256=outputs,
    source_full_file_retained=True,formatted_fallback=fallback,elapsed_seconds=time.time()-started,
    biological_inference_computed=False)
   receipt["status"]="complete";dump(folder/"COMPLETE.json",receipt);dump(folder/"RUN.json",receipt)
   print(json.dumps({"complete":args.accession,"source_rows":identity["source_rows"],"region_rows":{k:len(v) for k,v in wanted.items()},"peak_outcomes":len(reverse),"elapsed_seconds":receipt["elapsed_seconds"]}),flush=True)
  except BaseException as exc:
   run.update(status="failed",error=repr(exc),traceback=traceback.format_exc(),elapsed_seconds=time.time()-started)
   dump(folder/"FAILED.json",run);dump(folder/"RUN.json",run);raise
if __name__=="__main__":main()
