"""Native Oxford region extraction shared by BIG40 and the already-local QSM archive."""
from collections import Counter
import hashlib,io,json
from pathlib import Path
import numpy as np,pandas as pd
HEADER=["chr","rsid","pos","a1","a2","beta","se","pval(-log10)"]
FIELDS=["chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error",
 "p_value","neg_log10_p_value","variant_id","rsid","phenotype","source_y","build_y","n_y","n_y_status","source_data_row"]
class SequentialReader(io.RawIOBase):
 def __init__(self,binary):super().__init__();self.binary=binary
 def readable(self):return True
 def seekable(self):return False
 def readinto(self,buffer):
  value=self.binary.read(len(buffer));buffer[:len(value)]=value;return len(value)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(8<<20),b""):h.update(b)
 return h.hexdigest()
def dump(p,d):
 p=Path(p);t=p.with_name(p.name+".pending");t.write_text(json.dumps(d,indent=2));t.replace(p)
def extract(binary,record,windows,progress=None):
 """Validate every native source row, retain every allele in each complete window."""
 rows=0;pieces={wid:[] for wid in windows};chromosomes=Counter()
 with io.TextIOWrapper(io.BufferedReader(SequentialReader(binary)),encoding="ascii") as handle:
  assert handle.readline().split()==HEADER
  for chunk in pd.read_csv(handle,sep=r"\s+",names=HEADER,chunksize=300000,float_precision="round_trip",keep_default_na=False,
   dtype={"chr":str,"rsid":str,"a1":str,"a2":str}):
   chrom=chunk["chr"].str.lstrip("0").replace({"X":"23"})
   errors=dict(bad_coordinate=int((~chrom.isin([str(i) for i in range(1,24)]) | (chunk.pos<=0) | (chunk.pos!=np.floor(chunk.pos))).sum()),
    bad_beta=int((~np.isfinite(chunk.beta)).sum()),bad_se=int((~np.isfinite(chunk.se)|(chunk.se<=0)).sum()),
    bad_logp=int((~np.isfinite(chunk["pval(-log10)"])|(chunk["pval(-log10)"]<0)).sum()),
    bad_alleles=int((chunk.a1.isna()|chunk.a2.isna()|chunk.a1.eq(chunk.a2)).sum()))
   if any(errors.values()):raise ValueError("Invalid native rows: "+str(errors))
   for c,n in chrom.value_counts().items():chromosomes[c]+=int(n)
   for wid,w in windows.items():
    take=chrom.eq(str(w["chromosome"])) & chunk.pos.between(w["start_grch37"],w["end_grch37"])
    sub=chunk.loc[take].copy()
    if len(sub):
     sub["source_data_row"]=np.flatnonzero(take.to_numpy())+rows+1
     sub["chromosome"]=chrom.loc[take].to_numpy();pieces[wid].append(sub)
   rows+=len(chunk)
   if progress is not None and rows%3000000==0:progress(rows)
 if rows!=record["expected_source_rows"]:raise ValueError(f"Row count {rows} != {record['expected_source_rows']}")
 outputs={}
 for wid,parts in pieces.items():
  data=pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=HEADER+["source_data_row","chromosome"])
  frame=pd.DataFrame(dict(chromosome=data.chromosome,base_pair_location=data.pos.astype("int64"),
   effect_allele=data.a2,other_allele=data.a1,beta=data.beta.astype(float),standard_error=data.se.astype(float),
   p_value=[10**(-float(x)) for x in data["pval(-log10)"]],neg_log10_p_value=data["pval(-log10)"].astype(float),
   variant_id=data.rsid,rsid=data.rsid,phenotype=record["phenotype"],source_y=record["native_url"],build_y="GRCh37",
   n_y=np.nan,n_y_status="exact_IDP_sample_size_not_verified",source_data_row=data.source_data_row.astype("int64")))
  if record["family"]=="BIG40":
   auto=frame.chromosome.ne("23")
   n=int(record["metadata"]["n_discovery"]);assert n>0
   frame.loc[auto,"n_y"]=n
   frame["n_y_status"]=np.where(auto,"catalog_autosomal_exact","X_not_assigned")
  else:
   frame["effect_allele_frequency"]=np.nan;frame["imputation_info"]=np.nan
   frame["source_member"]=record["source_member"]
   frame["build_evidence"]="inherited_official_BIG40_convention"
   frame["units_y"]="SD of Gaussianized and deconfounded imaging phenotype"
  outputs[wid]=frame
 return outputs,dict(source_rows=rows,chromosome_rows=dict(chromosomes),all_records_numerically_validated=True,
  region_allele_or_P_filter=None)
def attach_frequency(frame,reference):
 keys=["chromosome","base_pair_location","effect_allele","other_allele","rsid"]
 for f in [frame,reference]:
  f["chromosome"]=f["chromosome"].astype(str)
  assert not f.duplicated(keys).any(),"Duplicate native regional identity"
 merged=frame.merge(reference[keys+["effect_allele_frequency","imputation_info"]],on=keys,how="outer",validate="one_to_one",indicator=True)
 if not merged["_merge"].eq("both").all():
  raise ValueError("Native statistics differ from complete official regional variant catalog: "+str(merged["_merge"].value_counts().to_dict()))
 return merged.drop(columns="_merge")
