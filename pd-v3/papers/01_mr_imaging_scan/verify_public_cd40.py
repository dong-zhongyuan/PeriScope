"""Check pinned code and independently verify every R ABF posterior in Python."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.special import logsumexp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pdv3.tabix_http import save_json


def posterior(bx,sx,by,sy,p1,p2,p12):
    # Independent direct marginal-likelihood expression; no R intermediates.
    def abf(b,s):
        v=s*s; w=.15**2
        return .5*(np.log(v/(v+w)) + b*b*w/(v*(v+w)))
    x,y=abf(bx,sx),abf(by,sy)
    # Stable off-diagonal sum for H3, without subtracting two nearly equal sums.
    left=np.r_[-np.inf,np.logaddexp.accumulate(y)[:-1]]
    right=np.r_[np.logaddexp.accumulate(y[::-1])[-2::-1],-np.inf]
    distinct=logsumexp(x+np.logaddexp(left,right))
    logs=np.array([0,np.log(p1)+logsumexp(x),np.log(p2)+logsumexp(y),
                   np.log(p1*p2)+distinct,np.log(p12)+logsumexp(x+y)])
    return np.exp(logs-logsumexp(logs))


def main():
    out=ROOT/'papers/01_mr_imaging_scan/results/public_cd40_20261116'
    ref=json.loads((ROOT/'data/scallop_cd40/COLOC_REFERENCE.json').read_text())
    for record in ref['files'].values():
        assert hashlib.sha256((ROOT/record['path']).read_bytes()).hexdigest()==record['sha256']
    audit=json.loads((out/'INPUT_AUDIT.json').read_text()); assert audit['complete']
    for d in audit['outcomes'].values():
        record=d['coloc_input'];assert hashlib.sha256((ROOT/record['path']).read_bytes()).hexdigest()==record['sha256']
    subprocess.run(['Rscript',str(Path(__file__).with_name('public_cd40_coloc.R')),str(ROOT)],check=True)
    r=pd.read_csv(out/'coloc_prior_sensitivity.csv',dtype={'idp':str})
    differences=[]
    for idp in ['0017','0018']:
        d=pd.read_csv(out/f'{idp}_shared_snps.csv.gz')
        for row in r[r.idp.eq(idp)].itertuples(index=False,name=None):
            rec=dict(zip(r.columns,row));v=d[abs(d.pos37-44747947)<=rec['window_bp']]
            assert len(v)==rec['nsnps']
            pp=posterior(v.beta_x.to_numpy(),v.se_x.to_numpy(),v.beta_y.to_numpy(),v.se_y.to_numpy(),rec['p1'],rec['p2'],rec['p12'])
            values=np.array([rec[f'PP.H{i}.abf'] for i in range(5)])
            assert np.isfinite(pp).all() and abs(pp.sum()-1)<1e-10
            differences.append(float(abs(pp-values).max()))
    assert len(differences)==len(r)>0 and max(differences)<1e-9
    save_json(out/'NUMERICAL_AUDIT.json',{'complete':True,'rows_checked':len(r),'max_abs_posterior_difference':max(differences),
        'official_coloc_commit':ref['commit'],'official_functions':['approx.bf.estimates','combine.abf'],
        'independent_method':'Python direct ABF expression; H3 via stable off-diagonal log-sum',
        'coloc_model':'At most one causal variant per trait in the supplied region; no LD conditioning',
        'not_completed':['UKB-PPP regional colocalization','multi-signal coloc with matched signed LD','reverse MR'],
        'code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('public_cd40_coloc.R')]},
        'source_data':'SCALLOP CD40 / BIG40 discovery bilateral pallidum; public external protein source',
        'posterior_sha256':hashlib.sha256((out/'coloc_prior_sensitivity.csv').read_bytes()).hexdigest()})
    print(r[np.isclose(r.p12_requested,1e-5,rtol=1e-12)].to_string(index=False))


if __name__=='__main__':main()
