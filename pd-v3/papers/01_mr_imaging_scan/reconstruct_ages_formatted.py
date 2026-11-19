"""Reconstruct exact AGES author bytes from a GWAS Catalog formatted gzip stream.

This narrow reversible text rule is accepted only when the complete reconstructed
stream matches that accession's published ORIGINAL-file MD5. No build is inferred
from formatted filenames or sidecars. Every future file must pass independently.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,time
COLUMNS=b"variant_id p_value chromosome base_pair_location effect_allele other_allele effect_allele_frequency beta standard_error".split()
def sha(p,kind="sha256"):
 h=hashlib.new(kind)
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8*1024*1024),b""):h.update(b)
 return h.hexdigest()
def canonical_lines(formatted):
 with gzip.open(formatted,"rb") as f:
  header=f.readline()
  if not header.endswith(b"\n") or b"\r" in header:raise ValueError("Unexpected line ending")
  fields=header[:-1].split(b"\t")
  if fields[:9]!=COLUMNS or len(fields)!=12:raise ValueError("Unexpected formatted schema")
  yield b"\t".join(fields[:9])+b"\n",False
  for line in f:
   if not line.endswith(b"\n") or b"\r" in line:raise ValueError("Unexpected line ending")
   fields=line[:-1].split(b"\t")
   if len(fields)!=12:raise ValueError("Malformed formatted row")
   changed=fields[0]==b'\"\"\"\"\"\"'
   if changed:fields[0]=b'\"\"'
   yield b"\t".join(fields[:9])+b"\n",changed
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--formatted",required=True)
 ap.add_argument("--author-md5",required=True);ap.add_argument("--formatted-md5")
 ap.add_argument("--original-for-comparison");ap.add_argument("--report",required=True)
 args=ap.parse_args();report=Path(args.report)
 if report.exists():raise RuntimeError("Existing report; inspect before repeating")
 started=time.monotonic();formatted=Path(args.formatted)
 if args.formatted_md5:assert sha(formatted,"md5")==args.formatted_md5
 original=Path(args.original_for_comparison) if args.original_for_comparison else None
 if original:assert sha(original,"md5")==args.author_md5
 md5=hashlib.md5();strong=hashlib.sha256();corrected=0;n=0;different=0;bytecount=0
 ref=original.open("rb") if original else None
 try:
  for n,(line,change) in enumerate(canonical_lines(formatted)):
   md5.update(line);strong.update(line);corrected+=int(change);bytecount+=len(line)
   if ref and ref.readline()!=line:different+=1
  if ref and ref.read(1):raise ValueError("Original has more lines")
 finally:
  if ref:ref.close()
 complete=md5.hexdigest()==args.author_md5 and different==0
 result=dict(complete=complete,accession_source_file=str(formatted),formatted_sha256=sha(formatted),
  formatted_published_md5=args.formatted_md5,formatted_published_md5_verified=bool(args.formatted_md5),
  data_rows=n,reconstructed_bytes=bytecount,quote_escape_corrections=corrected,other_edits=0,
  reconstructed_original_md5=md5.hexdigest(),published_original_md5=args.author_md5,
  original_published_md5_verified=md5.hexdigest()==args.author_md5,
  reconstructed_original_sha256=strong.hexdigest(),direct_original_comparison=bool(original),
  differing_rows_after_reconstruction=different if original else None,
  original_sha256=sha(original) if original else None,
  gzip_crc_to_EOF_verified=True,script_sha256=sha(__file__),elapsed_seconds=time.monotonic()-started,
  rule="First nine fields; only exact six quote bytes in variant_id become two; LF. Each full output must match its own published original-file MD5.",
  build="Use original author metadata only after successful exact original identity; formatted metadata is not authoritative.")
 report.parent.mkdir(parents=True,exist_ok=True);tmp=report.with_name(report.name+".pending")
 tmp.write_text(json.dumps(result,indent=2));tmp.replace(report);print(json.dumps(result),flush=True)
 if not complete:raise SystemExit(1)
if __name__=="__main__":main()
