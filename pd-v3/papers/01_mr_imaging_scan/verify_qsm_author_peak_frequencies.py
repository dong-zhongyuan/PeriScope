"""Audit published QSM peak allele/frequency orientation against the original complete BIG40 catalog."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/"data/all_target_followup_20261118/qsm_published_peaks"
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def main():
 source=ROOT/"data/big40_variants/download/variants.txt.gz"
 assert sha(source)=="152517ada2fb5c572b26959d184f099369e521eeac681c01b3a3ca94b1e7299a"
 path=BASE/"all_517_author_peaks.csv";meta=json.loads((BASE/"AUTHOR_PEAK_MANIFEST.json").read_text())
 assert sha(path)==meta["outputs_sha256"][path.name]
 peaks=pd.read_csv(path,dtype={"chromosome":str,"idp":str},keep_default_na=False,float_precision="round_trip")
 ids=set(peaks.source_variant_id);parts=[];total=0
 for frame in pd.read_csv(source,sep=r"\s+",chunksize=500000,dtype={"chr":str,"rsid":str,"a1":str,"a2":str},
  keep_default_na=False,float_precision="round_trip"):
  total+=len(frame);selected=frame[frame.rsid.isin(ids)]
  if len(selected):parts.append(selected.copy())
 assert total==17103079
 reference=pd.concat(parts,ignore_index=True)
 reference["chromosome"]=reference.chr.map(lambda x:"X" if x in ["23","X","0X"] else str(int(x)))
 reference=reference.rename(columns={"rsid":"source_variant_id","pos":"position_grch37","a1":"other_allele","a2":"effect_allele","af":"reference_A2_frequency","info":"reference_info"})
 keys=["chromosome","position_grch37","source_variant_id","other_allele","effect_allele"]
 assert not reference.duplicated(keys).any()
 joined=peaks.merge(reference[keys+["reference_A2_frequency","reference_info"]],on=keys,how="left",validate="many_to_one",indicator=True)
 assert len(joined)==517
 matched=joined["_merge"].eq("both")
 delta=abs(joined.effect_allele_frequency-joined.reference_A2_frequency)
 joined["frequency_absolute_difference"]=delta
 joined["frequency_audit_status"]=np.where(~matched,"no_exact_catalog_identity",np.where(delta.le(1e-12),"identical_to_official_A2_frequency","same_alleles_different_frequency"))
 joined.to_csv(BASE/"all_517_author_peak_frequency_audit.csv",index=False)
 result=dict(complete=True,source_sha256=sha(source),source_rows=total,gzip_crc_eof_verified=True,
  source_peak_rows=len(joined),status_counts=joined.frequency_audit_status.value_counts().to_dict(),
  maximum_frequency_difference_matched=float(delta[matched].max()) if matched.any() else None,
  no_source_frequency_changed_or_imputed=True,reference_frequency_scope="Original complete BIG40 a2 catalog, not asserted per-IDP per-SNP N",
  code_sha256=sha(__file__),input_sha256={str(path.relative_to(ROOT)):sha(path)},
  output_sha256=sha(BASE/"all_517_author_peak_frequency_audit.csv"))
 (BASE/"FREQUENCY_AUDIT.json").write_text(json.dumps(result,indent=2)+"\n")
 print(json.dumps(result),flush=True)
if __name__=="__main__":main()
