"""Verify actual reverse follow-up outputs without importing the MR producer or kernels."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np,pandas as pd
from scipy.stats import norm,chi2
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'papers/01_mr_imaging_scan/results/new_scan_reverse_20261117';DATA=ROOT/'data/new_scan_regional_20261117'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 assert json.loads((OUT/'RUN.json').read_text())['complete']
 audit=pd.read_csv(OUT/'all_peaks_exposure_only_LD_audit.csv',float_precision='round_trip');pairs=pd.read_csv(OUT/'all_peak_pair_reference_LD.csv',float_precision='round_trip');ratios=pd.read_csv(OUT/'all_peak_reverse_ratios_and_coverage.csv',float_precision='round_trip');summary=pd.read_csv(OUT/'reverse_direction_estimates.csv',float_precision='round_trip');worst=0.;genotypes={}
 for p in (DATA/'reverse_LD').glob('chr*/effect_allele_dosages.npz'):
  with np.load(p,allow_pickle=False) as z:
   for k,g in zip(z['variant_keys'],z['dosages']):genotypes[k]=g.copy()
 def check(actual,expected):
  nonlocal worst
  a,b=np.asarray(actual,dtype=float),np.asarray(expected,dtype=float);assert np.allclose(a,b,rtol=1e-10,atol=1e-12,equal_nan=True),(actual,expected)
  if np.isfinite(a).any():worst=max(worst,float(np.nanmax(np.abs(a-b))))
 for row in pairs.itertuples():
  if row.variant_a not in genotypes or row.variant_b not in genotypes:assert np.isnan(row.r);continue
  a,b=genotypes[row.variant_a],genotypes[row.variant_b];valid=(a>=0)&(b>=0);assert row.n_complete_reference_samples==valid.sum();a=a[valid].astype(float);b=b[valid].astype(float);a-=a.mean();b-=b.mean();den=np.linalg.norm(a)*np.linalg.norm(b)
  expected=np.dot(a,b)/den if den>0 else np.nan;check([row.r,row.r2],[expected,expected**2])
 for acc,group in audit.groupby('accession'):
  # Independently replay selection over the full declared exposure-only order.
  remaining=group.sort_values(['discovery_neglog10_p','variant_key'],ascending=[False,True]);kept=[]
  for row in remaining.itertuples():
   if row.ld_selection_status in ['unsupported_alleles','outside_discovery_instrument_P_threshold','missing_reference_variant','monomorphic_or_no_reference_calls']:continue
   eligible=True
   for lead in kept:
    if row.chromosome==lead.chromosome and abs(row.position_grch37-lead.position_grch37)<=10000000:
     a,b=genotypes[row.variant_key],genotypes[lead.variant_key];valid=(a>=0)&(b>=0);a=a[valid].astype(float);b=b[valid].astype(float);a-=a.mean();b-=b.mean();r=np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b))
     if not np.isfinite(r) or r*r>.001:eligible=False;break
   assert eligible==(row.ld_selection_status=='retained_by_exposure_only_LD_pruning')
   if eligible:kept.append(row)
  raw=pd.read_csv(DATA/acc/'imaging_peak_outcomes.csv.gz',float_precision='round_trip');raw['variant_key']=[f'{v.chromosome}:{v.base_pair_location}:'+ '_'.join(sorted([v.effect_allele,v.other_allele])) for v in raw.itertuples()];raw=raw.set_index('variant_key');x=group.set_index('variant_key')
  for row in ratios[ratios.accession==acc].itertuples():
   if row.status!='estimated':continue
   protein=raw.loc[row.variant_key];exposure=x.loc[row.variant_key];sign=1 if protein.effect_allele==exposure.effect_allele else -1
   if row.harmonization_status.startswith('palindromic'):
    choices=[s for s in [1,-1] if abs(exposure.effect_allele_frequency-(protein.effect_allele_frequency if s==1 else 1-protein.effect_allele_frequency))<=.1];assert len(choices)==1 and max(min(exposure.effect_allele_frequency,1-exposure.effect_allele_frequency),min(protein.effect_allele_frequency,1-protein.effect_allele_frequency))<=.42;sign=choices[0]
   check(row.protein_beta_aligned,protein.beta*sign);check(row.beta_imaging,exposure[row.imaging_effect_set+'_beta']);check(row.se_imaging,exposure[row.imaging_effect_set+'_se'])
   b=row.protein_beta_aligned/row.beta_imaging;s=abs(row.protein_se/row.beta_imaging);sd=np.sqrt(s*s+row.protein_beta_aligned**2*row.se_imaging**2/row.beta_imaging**4);check([row.beta,row.se,row.se_delta,row.p,row.p_delta,row.F],[b,s,sd,2*norm.sf(abs(b/s)),2*norm.sf(abs(b/sd)),(row.beta_imaging/row.se_imaging)**2])
 for row in summary.itertuples():
  selected=ratios[(ratios.accession==row.accession)&(ratios.palindrome_policy==row.palindrome_policy)&(ratios.imaging_effect_set==row.imaging_effect_set)&ratios.selected_after_LD&(ratios.status=='estimated')];assert len(selected)==row.n_instruments
  x=selected.beta_imaging.to_numpy();y=selected.protein_beta_aligned.to_numpy();s=selected.protein_se.to_numpy();design=x/s;target=y/s
  b=float(np.linalg.lstsq(design[:,None],target,rcond=None)[0][0]);se=1/np.linalg.norm(design);check([row.beta,row.se,row.p],[b,se,2*norm.sf(abs(b/se))]);check([row.ci_lower,row.ci_upper],[b-norm.ppf(.975)*se,b+norm.ppf(.975)*se])
  if len(selected)>1:
   q=float(np.linalg.norm(target-design*b)**2);df=len(selected)-1;se_random=se*np.sqrt(max(1,q/df));check([row.Q,row.Q_p,row.se_random,row.p_random],[q,chi2.sf(q,df),se_random,2*norm.sf(abs(b/se_random))])
 assert len(audit)==48 and len(summary)==12 and len(ratios)==192
 report=dict(complete=True,verified_at_utc=datetime.now(timezone.utc).isoformat(),source_peak_rows=48,all_ratio_coverage_rows=192,summary_rows=12,LD_pair_rows=len(pairs),maximum_absolute_difference=worst,methods='Centered genotype dot products, exposure-only pruning replay, raw-source allele/estimate checks, closed-form scalar delta ratios, independent weighted least-squares reconstruction of combined estimates and Q',verifier_sha256=sha(Path(__file__)),input_sha256={str(p.relative_to(ROOT)):sha(p) for p in OUT.glob('*.csv')})
 (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
