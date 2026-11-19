"""Independent all-call audits plus deterministic public-source replay for each reference chromosome."""
from pathlib import Path
import argparse,hashlib,json,time
import numpy as np,pandas as pd,pysam,requests
ROOT=Path(__file__).resolve().parents[2];OLD=ROOT/"data/new_scan_regional_20261117/reverse_LD"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def identity(url):
 r=requests.head(url,timeout=45);r.raise_for_status();return {k:r.headers.get(k) for k in ["ETag","Last-Modified","Content-Length"]}
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--folder",type=Path,required=True);a=ap.parse_args();folder=a.folder
 receipt=folder/"MANIFEST.json";m=json.loads(receipt.read_text());assert m["complete"]
 target=folder/"INDEPENDENT_VERIFICATION.json"
 if target.exists():
  d=json.loads(target.read_text());assert d["complete"] and d["manifest_sha256"]==sha(receipt)
  for n,h in m["outputs_sha256"].items():assert sha(folder/n)==h
  print("CACHED_REFERENCE_VERIFIED",folder);return
 started=time.time()
 for n,h in m["outputs_sha256"].items():assert sha(folder/n)==h
 panel=OLD/"integrated_call_samples_v3.20130502.ALL.panel";assert sha(panel)==m["panel_sha256"]
 p=pd.read_csv(panel,sep="\t",usecols=["sample","super_pop","gender"]);p=p[p.super_pop.eq("EUR")]
 samples=p["sample"].to_numpy(str);sex=p.gender.to_numpy(str);assert len(samples)==503
 cov=pd.read_csv(folder/"coverage.csv",keep_default_na=False,float_precision="round_trip")
 with np.load(folder/"effect_allele_dosages.npz",allow_pickle=False) as z:arrays={k:z[k] for k in z.files}
 assert np.array_equal(samples,arrays["samples"])
 match=cov[cov.status.eq("matched")].copy()
 assert match.variant_key.is_unique and cov.variant_key.is_unique and len(cov)==m["requested"]
 assert np.array_equal(match.dosage_row.astype(int),np.arange(len(match)))
 assert np.array_equal(match.variant_key.to_numpy(str),arrays["variant_keys"])
 g=arrays["dosages"];assert g.shape==(len(match),503) and np.isin(g,[-1,0,1,2]).all()
 raw_count=0
 if "raw_GT" in arrays:
  raw=arrays["raw_GT"];pl=arrays["source_GT_ploidy"];assert raw.shape==(len(match),503,2) and pl.shape==g.shape
  assert np.array_equal(arrays["gender"],sex)
  assert np.isin(raw,[-2,-1,0,1]).all() and np.isin(pl,[1,2]).all()
  assert np.array_equal(pl,np.sum(raw!=-2,axis=2))
  for row in match.itertuples():
   i=int(row.dosage_row);v=raw[i];missing=(v==-1).any(axis=1)
   ea_index=0 if row.effect_allele==row.ref else 1
   expected=2*(v==ea_index).sum(axis=1)/pl[i];expected[missing]=-1
   np.testing.assert_array_equal(expected,g[i]);raw_count+=503
   if int(m["chromosome"])==23:
    par=int(row.position_grch37)<=2699520 or int(row.position_grch37)>=154931044
    required_diploid=(sex=="female")|par
    assert (pl[i,required_diploid]==2).all()
    male_nonpar=(sex=="male")&(not par)
    diploid_male=(pl[i]==2)&male_nonpar&~missing
    assert (v[diploid_male,0]==v[diploid_male,1]).all()
 # Replay first/middle/last observed SNP across all 503 samples from the official VCF.
 choices=np.unique(np.linspace(0,len(match)-1,min(3,len(match)),dtype=int)) if len(match) else []
 assert identity(m["url"])==m["http_identity"]
 index=next(folder.glob("*.tbi"));assert sha(index)==m["index_sha256"];replayed=[]
 label="X" if int(m["chromosome"])==23 else str(m["chromosome"])
 with pysam.VariantFile(m["url"],index_filename=str(index)) as source:
  source.subset_samples(samples.tolist())
  for j in choices:
   row=match.iloc[int(j)];found=[v for v in source.fetch(label,int(row.position_grch37)-1,int(row.position_grch37))
    if v.pos==int(row.position_grch37) and v.alts and len(v.alts)==1 and v.ref==row.ref and v.alts[0]==row.alt]
   assert len(found)==1;v=found[0];ea_index=0 if row.effect_allele==v.ref else 1;expected=[]
   for sample in samples:
    gt=v.samples[str(sample)]["GT"]
    if not gt or any(x is None for x in gt):expected.append(-1)
    else:expected.append(2*sum(x==ea_index for x in gt)/len(gt))
   np.testing.assert_array_equal(np.array(expected),g[int(row.dosage_row)])
   replayed.append(row.variant_key)
 assert identity(m["url"])==m["http_identity"]
 for n,h in m["outputs_sha256"].items():assert sha(folder/n)==h
 report=dict(complete=True,chromosome=m["chromosome"],manifest_sha256=sha(receipt),coverage_rows=len(cov),
  matched_rows=len(match),all_reference_calls_shape_and_domain_checked=True,source_sample_order_exact=True,
  retained_raw_GT_calls_independently_reconstructed=raw_count,
  remote_replayed_SNPs=replayed,remote_replayed_calls=len(replayed)*503,all_replayed_genotypes_identical=True,
  whole_remote_file_checksum_verified=False,verifier_sha256=sha(__file__),elapsed_seconds=time.time()-started)
 target.write_text(json.dumps(report,indent=2)+"\n");print("REFERENCE_VERIFIED",str(folder.relative_to(ROOT)),len(match),flush=True)
if __name__=="__main__":main()
