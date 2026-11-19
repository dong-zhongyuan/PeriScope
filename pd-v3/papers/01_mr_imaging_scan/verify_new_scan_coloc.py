"""Independent Python reconstruction of the R coloc scales, ABFs and all posterior grids."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np,pandas as pd
from scipy.special import logsumexp
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'papers/01_mr_imaging_scan/results/new_scan_regions_20261117'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 audit=json.loads((OUT/'INPUT_AUDIT.json').read_text());results=pd.read_csv(OUT/'coloc_prior_sensitivity.csv',float_precision='round_trip');checks=[];worst=0.
 for rec in audit['datasets']:
  acc=rec['accession'];d=OUT/acc;shared=pd.read_csv(d/'shared_snps.csv.gz',float_precision='round_trip')
  assert sha(d/'shared_snps.csv.gz')==rec['outputs']['shared']['sha256']
  for window in [500000,1000000]:
   v=shared.loc[(shared.pos37-rec['region']['center']).abs()<=window].copy();variance=v.standard_error_x.to_numpy()**2;freq=v.effect_allele_frequency_x.to_numpy();inverse=1/variance;nvx=2*v.n_x.to_numpy()*freq*(1-freq);sd=np.sqrt(np.dot(inverse,nvx)/np.dot(inverse,inverse));a=[]
   for suffix,scale in [('x',sd),('y',1.)]:
    s=v['standard_error_'+suffix].to_numpy();z=v['beta_'+suffix].to_numpy()/s;W=(.15*scale)**2;weight=W/(s*s+W);a.append(.5*(np.log1p(-weight)+weight*z*z))
   saved=pd.read_csv(d/f'window_{window}_snp_posteriors.csv.gz',float_precision='round_trip');assert saved.variant_key.tolist()==v.variant_key.tolist()
   for col,expected in [('log_abf_protein',a[0]),('log_abf_imaging',a[1]),('snp_pp_H4_conditional',np.exp(a[0]+a[1]-logsumexp(a[0]+a[1]))),('snp_pp_protein_conditional',np.exp(a[0]-logsumexp(a[0]))),('snp_pp_imaging_conditional',np.exp(a[1]-logsumexp(a[1])))]:
    error=float(np.max(np.abs(saved[col].to_numpy()-expected)));worst=max(worst,error);assert np.allclose(saved[col],expected,rtol=1e-10,atol=1e-10),(acc,window,col,error)
   group=results[(results.accession==acc)&(results.window_bp==window)];assert len(group)==32 and group.nsnps.eq(len(v)).all();assert np.allclose(group.sdY_protein_estimated,sd,rtol=1e-12,atol=1e-12)
   l1,l2=map(logsumexp,a);lj=logsumexp(a[0]+a[1]);lp=l1+l2
   assert lp>=lj
   ld=lp+np.log(-np.expm1(lj-lp)) if lp>lj else -np.inf
   for row in group.itertuples(index=False,name=None):
    q=dict(zip(group.columns,row));lp1=np.log(q['p1']);lp2=np.log(q['p2']);logs=np.array([0.,lp1+l1,lp2+l2,lp1+lp2+ld,np.log(q['p12'])+lj]);post=np.exp(logs-logsumexp(logs));actual=np.array([q[f'PP.H{i}.abf'] for i in range(5)]);error=float(np.max(np.abs(actual-post)));worst=max(worst,error);assert np.allclose(actual,post,rtol=1e-10,atol=1e-10),(acc,window,error)
   checks.append(dict(accession=acc,window_bp=window,shared_snps=len(v),protein_SD_estimated=float(sd),prior_rows=len(group)))
 assert len(results)==192
 report=dict(complete=True,verified_at_utc=datetime.now(timezone.utc).isoformat(),independent_implementation='NumPy closed-form through-origin variance estimate and Wakefield ABFs; SciPy logsumexp H0-H4 reconstruction. No production or R function imports.',checks=checks,posterior_rows_verified=len(results),maximum_absolute_difference=worst,verifier_sha256=sha(Path(__file__)),input_audit_sha256=sha(OUT/'INPUT_AUDIT.json'),coloc_results_sha256=sha(OUT/'coloc_prior_sensitivity.csv'))
 (OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(report,indent=2)+'\n');print('COLOC_INDEPENDENT_VERIFIED',len(results),worst,flush=True)
if __name__=='__main__':main()
