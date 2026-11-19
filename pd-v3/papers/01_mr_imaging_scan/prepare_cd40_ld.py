"""Align public CD40 statistics and external EUR ALT dosages for LD sensitivity."""
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pdv3.mr import harmonize
from pdv3.tabix_http import save_json


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    source=ROOT/'papers/01_mr_imaging_scan/results/public_cd40_20261116'
    ref=ROOT/'data/ld_1000g/extracted/CD40'
    out=ROOT/'papers/01_mr_imaging_scan/results/public_cd40_ld_20261116'
    out.mkdir(exist_ok=True)
    m=json.loads((ref/'MANIFEST.json').read_text());assert m['complete'] and not m['cohort_matched']
    for p,h in m['outputs'].items():assert sha(ROOT/p)==h
    v=pd.read_csv(ref/'variants.csv.gz');g=np.load(ref/'alt_dosages.npz')
    assert v.variant_key.tolist()==g['variant_keys'].tolist()
    v['genotype_row']=np.arange(len(v))
    tables={};audits=[]
    for idp in ['0017','0018']:
        d=pd.read_csv(source/f'{idp}_shared_snps.csv.gz')
        joined=d.merge(v,on='variant_key',how='left',suffixes=('_summary','_reference'),validate='one_to_one')
        records=[]
        for r in joined.to_dict('records'):
            row=dict(variant_key=r['variant_key'],idp=idp,source_variant=r['variant_y'],pos37=r['pos37_summary'],
                     genotype_row=r['genotype_row'],n_x=r['n_x'],n_y=r['n_y'],status='reference_missing')
            if pd.notna(r['genotype_row']):
                freq=r['alt_frequency'];maf=min(freq,1-freq)
                hx=harmonize(r['alt'],r['ref'],r['ea_x'],r['oa_x'],r['beta_x'],eaf_x=freq,eaf_y=r['eaf_x'],exposure_build='GRCh37',outcome_build='GRCh37',palindrome='drop')
                hy=harmonize(r['alt'],r['ref'],r['ea_y'],r['oa_y'],r['beta_y'],exposure_build='GRCh37',outcome_build='GRCh37',palindrome='drop')
                row.update(ref=r['ref'],alt=r['alt'],alt_frequency=freq,maf_reference=maf,
                           beta_x_alt=hx['beta_y'],beta_y_alt=hy['beta_y'],se_x=r['se_x'],se_y=r['se_y'],
                           protein_reference_eaf_difference=hx['eaf_difference'],x_alignment=hx['status'],y_alignment=hy['status'])
                if not hx['usable'] or not hy['usable']:row['status']='allele_alignment_unusable'
                elif r['missing_calls']!=0:row['status']='missing_reference_genotypes'
                elif maf<.01:row['status']='reference_maf_below_0p01'
                elif hx['eaf_difference']>.1:row['status']='protein_reference_frequency_discordance'
                else:row['status']='common_aligned'
            records.append(row)
        audit=pd.DataFrame(records)
        audit['near_complete_sample_size']=audit.n_x.ge(.99*audit.n_x.max())
        audit.to_csv(out/f'{idp}_alignment_audit.csv.gz',index=False,compression={'method':'gzip','mtime':0})
        tables[idp]=audit[audit.status.eq('common_aligned') & audit.near_complete_sample_size].sort_values(['pos37','variant_key']).reset_index(drop=True)
        audits.append(dict(idp=idp,status_counts=audit.status.value_counts().to_dict(),
                           near_complete_model_snps=len(tables[idp]),all_common_aligned=int(audit.status.eq('common_aligned').sum())))
    left,right=tables['0017'],tables['0018']
    assert left.variant_key.tolist()==right.variant_key.tolist()
    assert np.allclose(left.beta_x_alt,right.beta_x_alt) and np.allclose(left.se_x,right.se_x)
    combined=left.rename(columns={'beta_y_alt':'beta_y_left','se_y':'se_y_left'}).copy()
    combined['beta_y_right']=right.beta_y_alt;combined['se_y_right']=right.se_y
    combined.to_csv(out/'model_inputs.csv.gz',index=False,compression={'method':'gzip','mtime':0})
    dosage=g['dosages'][left.genotype_row.to_numpy(dtype=int)].astype(np.float64)
    assert np.isin(dosage,[0,1,2]).all() and (dosage.var(axis=1)>0).all()
    pd.DataFrame(dosage.astype(np.int8),index=left.variant_key,columns=g['samples']).to_csv(out/'alt_dosages.csv.gz',index_label='variant_key',compression={'method':'gzip','mtime':0})
    # Export independently computed checks for R's LD construction, including
    # every retained SNP's signed correlation to the historical lead.
    lead_i=int(np.flatnonzero(left.source_variant.eq('rs4810485'))[0])
    centered=dosage-dosage.mean(axis=1,keepdims=True)
    norms=np.sqrt((centered*centered).sum(axis=1))
    correlations=centered.dot(centered[lead_i])/(norms*norms[lead_i])
    pd.DataFrame(dict(variant_key=left.variant_key,alt_frequency=dosage.mean(axis=1)/2,
                      signed_r_to_rs4810485=correlations)).to_csv(out/'python_ld_checks.csv',index=False)
    inputs=[source/f'{i}_shared_snps.csv.gz' for i in ['0017','0018']]+[ref/'MANIFEST.json',ref/'variants.csv.gz',ref/'alt_dosages.npz']
    save_json(out/'INPUT_AUDIT.json',dict(complete=True,source_inputs={str(p.relative_to(ROOT)):sha(p) for p in inputs},
        source_code_sha256=sha(Path(__file__)),variant_audits=audits,model_snps=len(left),reference_samples=503,
        population_counts=m['population_counts'],external_reference_not_cohort_matched=True,
        filters={'reference_maf_min':.01,'palindromes':'drop','protein_reference_eaf_max_difference':.1,
                 'protein_n_min_fraction_of_max':.99,'outcome_p_filter':None},
        n_x_min=float(left.n_x.min()),n_x_max=float(left.n_x.max()),n_x_model=float(left.n_x.median()),
        n_y_model=int(left.n_y.iloc[0]),
        n_x_caveat='Model uses median N in the near-complete stratum; original per-SNP N retained; meta-analysis residual covariance remains approximate',
        interpretation='Multi-signal sensitivity using an external EUR reference; no new support threshold'))
    print(json.dumps(audits),flush=True)


if __name__=='__main__':main()
