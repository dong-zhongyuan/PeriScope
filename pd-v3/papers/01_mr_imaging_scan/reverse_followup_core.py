"""Original published-peak reverse follow-up arithmetic, generalized to each frozen region."""
from pathlib import Path
import itertools,sys
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
from pdv3.mr import harmonize,wald,ivw
NUMERIC_RESULT_COLUMNS="beta,se,p,log_p,ci_lower,ci_upper,se_delta,p_delta,log_p_delta,ci_lower_delta,ci_upper_delta,F,n_instruments,se_random,p_random,log_p_random,ci_lower_random,ci_upper_random,Q,Q_p,df,dispersion,F_min".split(",")
PEAK_COLUMNS="phenotype,source_variant_id,chromosome,position_grch37,effect_allele,other_allele,effect_allele_frequency,discovery_neglog10_p,discovery_beta,discovery_se,replication_neglog10_p,replication_beta,replication_se".split(",")
OUTCOME_COLUMNS="chromosome,base_pair_location,effect_allele,other_allele,beta,standard_error,effect_allele_frequency".split(",")
PAIR_COLUMNS="accession,variant_a,variant_b,n_complete_reference_samples,r,r2,within_10Mb_same_chromosome".split(",")
RATIO_COLUMNS="accession,gene,phenotype,variant_key,source_variant_id,chromosome,position_grch37,ld_selection_status,palindrome_policy,imaging_effect_set,selected_after_LD,protein_outcome_available,beta_imaging,se_imaging,status,protein_beta_source,protein_se,protein_effect_allele,imaging_effect_allele,harmonization_status,protein_beta_aligned,method".split(",")
def chromosome(c):
 value=str(c).lstrip("0")
 return 23 if value=="X" else int(value)
def key(c,p,a,b):return f"{chromosome(c)}:{int(p)}:"+"_".join(sorted([a,b]))
def canonical_inputs(x,y):
 x=x.copy();y=y.copy()
 for c in PEAK_COLUMNS:
  if c not in x:
   assert len(x)==0; x[c]=pd.Series(dtype=object)
 for c in OUTCOME_COLUMNS:
  if c not in y:
   assert len(y)==0;y[c]=pd.Series(dtype=object)
 if len(x):
  x["chromosome"]=x.chromosome.map(chromosome)
  for c in ["position_grch37","effect_allele_frequency","discovery_neglog10_p","discovery_beta","discovery_se","replication_neglog10_p","replication_beta","replication_se"]:
   x[c]=pd.to_numeric(x[c],errors="raise")
 if len(y):
  y["chromosome"]=y.chromosome.map(chromosome)
  for c in ["base_pair_location","beta","standard_error","effect_allele_frequency"]:y[c]=pd.to_numeric(y[c],errors="raise")
 x["variant_key"]=[key(v.chromosome,v.position_grch37,v.effect_allele,v.other_allele) for v in x.itertuples()]
 y["variant_key"]=[key(v.chromosome,v.base_pair_location,v.effect_allele,v.other_allele) for v in y.itertuples()]
 assert x.variant_key.is_unique and y.variant_key.is_unique
 return x,y
def calculate(rec,x,y,genotypes):
 x,y=canonical_inputs(x,y);acc=rec["accession"]
 assert not len(x) or x.phenotype.eq(rec["phenotype"]).all()
 correlations={};pairs=[]
 for a,b in itertools.combinations(x.to_dict("records"),2):
  ka,kb=a["variant_key"],b["variant_key"];within=a["chromosome"]==b["chromosome"] and abs(a["position_grch37"]-b["position_grch37"])<=10000000
  correlation=np.nan;n=0
  if ka in genotypes and kb in genotypes:
   ga,gb=genotypes[ka],genotypes[kb];assert ga.shape==gb.shape==(503,)
   valid=(ga>=0)&(gb>=0);n=int(valid.sum())
   if n>=3 and np.var(ga[valid])>0 and np.var(gb[valid])>0:correlation=float(np.corrcoef(ga[valid],gb[valid])[0,1])
  correlations[frozenset([ka,kb])]=correlation
  pairs.append(dict(accession=acc,variant_a=ka,variant_b=kb,n_complete_reference_samples=n,r=correlation,r2=correlation**2,within_10Mb_same_chromosome=within))
 # Selection depends only on imaging discovery association and the reference.
 kept=[];clump=[]
 for a in x.sort_values(["discovery_neglog10_p","variant_key"],ascending=[False,True]).to_dict("records"):
  k=a["variant_key"];status="retained_by_exposure_only_LD_pruning";blocker="";r2=np.nan
  if not (len(a["effect_allele"])==len(a["other_allele"])==1 and a["effect_allele"] in "ACGT" and a["other_allele"] in "ACGT"):status="unsupported_alleles"
  elif a["discovery_neglog10_p"] < -np.log10(5e-8):status="outside_discovery_instrument_P_threshold"
  elif k not in genotypes:status="missing_reference_variant"
  elif np.all(genotypes[k]<0) or np.var(genotypes[k][genotypes[k]>=0])==0:status="monomorphic_or_no_reference_calls"
  else:
   for b in kept:
    if a["chromosome"]==b["chromosome"] and abs(a["position_grch37"]-b["position_grch37"])<=10000000:
     corr=correlations[frozenset([k,b["variant_key"]])]
     if not np.isfinite(corr):status="LD_not_estimable";blocker=b["variant_key"];break
     if corr**2>0.001:status="LD_pruned";blocker=b["variant_key"];r2=corr**2;break
  if status=="retained_by_exposure_only_LD_pruning":kept.append(a)
  a.update(accession=acc,gene=rec["gene"],ld_selection_status=status,ld_blocker=blocker,ld_blocker_r2=r2);clump.append(a)
 kept_keys={a["variant_key"] for a in kept};cy=y.set_index("variant_key")
 audit=pd.DataFrame(clump) if clump else pd.DataFrame(columns=list(x.columns)+["accession","gene","ld_selection_status","ld_blocker","ld_blocker_r2"])
 audit["protein_outcome_available"]=audit.variant_key.isin(set(y.variant_key))
 estimates=[];summary=[]
 for policy in ["drop","infer"]:
  for effect_set in ["discovery","replication"]:
   valid=[]
   for a in clump:
    row={k:a[k] for k in ["accession","gene","phenotype","variant_key","source_variant_id","chromosome","position_grch37","ld_selection_status"]}
    row.update(palindrome_policy=policy,imaging_effect_set=effect_set,selected_after_LD=a["variant_key"] in kept_keys,protein_outcome_available=a["variant_key"] in cy.index)
    bx=float(a[effect_set+"_beta"]);sx=float(a[effect_set+"_se"]);row.update(beta_imaging=bx,se_imaging=sx)
    if a["variant_key"] not in cy.index:row["status"]="protein_outcome_unavailable";estimates.append(row);continue
    b=cy.loc[a["variant_key"]]
    h=harmonize(a["effect_allele"],a["other_allele"],b.effect_allele,b.other_allele,b.beta,a["effect_allele_frequency"],b.effect_allele_frequency,
     exposure_build="GRCh37",outcome_build="GRCh37",palindrome=policy)
    row.update(protein_beta_source=float(b.beta),protein_se=float(b.standard_error),protein_effect_allele=b.effect_allele,
     imaging_effect_allele=a["effect_allele"],harmonization_status=h["status"],protein_beta_aligned=h["beta_y"])
    if not h["usable"]:row["status"]=h["status"];estimates.append(row);continue
    if not np.isfinite([bx,sx]).all() or bx==0 or sx<=0:row["status"]="invalid_imaging_effect_for_ratio";estimates.append(row);continue
    row.update(wald(bx,sx,h["beta_y"],b.standard_error));row["status"]="estimated";estimates.append(row)
    if row["selected_after_LD"]:valid.append(row)
   common=dict(accession=acc,gene=rec["gene"],phenotype=rec["phenotype"],palindrome_policy=policy,imaging_effect_set=effect_set,
    source_peaks=len(x),LD_selected_before_outcome_join=len(kept),usable_instruments=len(valid),
    direction="imaging_to_AGES_protein",p_scope="nominal_followup_P_not_a_support_decision")
   if not len(x):common.update(status="not_estimated_no_published_peak",method="none",n_instruments=0)
   elif not valid:common.update(status="not_estimated_no_usable_instrument",method="none",n_instruments=0)
   elif len(valid)==1:common.update(status="estimated",**wald(valid[0]["beta_imaging"],valid[0]["se_imaging"],valid[0]["protein_beta_aligned"],valid[0]["protein_se"]))
   else:common.update(status="estimated",**ivw([v["beta_imaging"] for v in valid],[v["se_imaging"] for v in valid],
    [v["protein_beta_aligned"] for v in valid],[v["protein_se"] for v in valid],independent=True))
   summary.append(common)
 ratios=pd.DataFrame(estimates) if estimates else pd.DataFrame(columns=RATIO_COLUMNS)
 for c in RATIO_COLUMNS+NUMERIC_RESULT_COLUMNS:
  if c not in ratios:ratios[c]=np.nan
 combined=pd.DataFrame(summary)
 for c in NUMERIC_RESULT_COLUMNS:
  if c not in combined:combined[c]=np.nan
 return dict(audit=audit,pairs=pd.DataFrame(pairs,columns=PAIR_COLUMNS),ratios=ratios,summary=combined)
