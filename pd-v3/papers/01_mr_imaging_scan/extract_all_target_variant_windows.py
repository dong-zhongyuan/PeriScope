"""Extract all complete BIG40 variant-reference windows in one verified local pass."""
from pathlib import Path
from collections import defaultdict,Counter
import bisect,csv,gzip,hashlib,io,json,math,os,time
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/"data/all_target_followup_20261118/regional_all/imaging"
SOURCE=ROOT/"data/big40_variants/download/variants.txt.gz"
EXPECTED="152517ada2fb5c572b26959d184f099369e521eeac681c01b3a3ca94b1e7299a"
FIELDS=["chromosome","base_pair_location","rsid","other_allele","effect_allele","effect_allele_frequency","imputation_info"]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def dump(p,d):
 t=p.with_name(p.name+".pending");t.write_text(json.dumps(d,indent=2));t.replace(p)
def main():
 out=BASE/"variant_reference";out.mkdir(exist_ok=True)
 if (out/"MANIFEST.json").exists() or (out/"RUN.json").exists():raise RuntimeError("Inspect previous reference extraction")
 assert sha(SOURCE)==EXPECTED
 plan=json.loads((BASE/"PLAN.json").read_text());ph=sha(BASE/"PLAN.json");code=sha(__file__)
 bychr=defaultdict(list)
 for wid,win in plan["windows"].items():
  if win["BIG40_phenotypes"]:bychr[int(win["chromosome"])].append((win["center_grch37"],wid))
 for c in bychr:bychr[c].sort()
 centers={c:[v[0] for v in items] for c,items in bychr.items()}
 count={wid:0 for items in bychr.values() for _,wid in items}
 state=dict(pid=os.getpid(),status="running",complete=False,plan_sha256=ph,source_sha256=EXPECTED,code_sha256=code)
 dump(out/"RUN.json",state)
 opened=[];writers={};last=(0,0);chrcount=Counter();started=time.time()
 try:
  with gzip.open(SOURCE,"rt",encoding="ascii") as source:
   assert source.readline().split()==["chr","rsid","pos","a1","a2","af","info"]
   for number,line in enumerate(source,1):
    v=line.split();assert len(v)==7
    c=int(v[0].lstrip("0").replace("X","23"));pos=int(v[2]);af=float(v[5]);info=float(v[6])
    assert 1<=c<=23 and pos>0 and (c,pos)>=last and math.isfinite(af) and 0<=af<=1 and math.isfinite(info)
    if c!=last[0]:
     for handle in opened:handle.close()
     opened=[];writers={}
     for _,wid in bychr.get(c,[]):
      p=out/(wid+".csv.gz")
      raw=p.open("wb");zipped=gzip.GzipFile(filename="",mode="wb",fileobj=raw,mtime=0)
      text=io.TextIOWrapper(zipped,encoding="ascii",newline="")
      opened.extend([text,raw]);w=csv.writer(text,lineterminator="\n");w.writerow(FIELDS);writers[wid]=w
    cc=centers.get(c,[])
    for i in range(bisect.bisect_left(cc,pos-1000000),bisect.bisect_right(cc,pos+1000000)):
     wid=bychr[c][i][1];writers[wid].writerow([c,pos,v[1],v[3],v[4],af,info]);count[wid]+=1
    last=(c,pos);chrcount[str(c)]+=1
    if number%1000000==0:
     state.update(source_rows_processed=number,chromosome=c,elapsed_seconds=time.time()-started);dump(out/"RUN.json",state);print(json.dumps(state),flush=True)
  assert number==17103079
 finally:
  for handle in opened:handle.close()
 assert sha(BASE/"PLAN.json")==ph and sha(__file__)==code
 outputs={}
 for wid,n in count.items():
  p=out/(wid+".csv.gz")
  assert p.exists()
  outputs[wid]=dict(rows=n,output=str(p.relative_to(ROOT)),sha256=sha(p),coverage="full_native_region" if n else "no_native_catalog_variants")
 report=dict(complete=True,source_rows=number,chromosome_rows=dict(chrcount),gzip_crc_eof_verified=True,all_rows_order_coordinate_and_frequency_checked=True,
  source_sha256=EXPECTED,plan_sha256=ph,code_sha256=code,outputs=outputs,elapsed_seconds=time.time()-started,
  allele_or_P_filter=None,QSM_frequencies_assigned=False)
 dump(out/"MANIFEST.json",report);state.update(status="complete",complete=True);dump(out/"RUN.json",state)
 print(json.dumps(dict(complete=True,windows=len(outputs),total_regional_rows=sum(count.values()),elapsed_seconds=report["elapsed_seconds"])),flush=True)
if __name__=="__main__":main()
