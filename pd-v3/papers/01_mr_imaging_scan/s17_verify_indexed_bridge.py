"""Compare independently retrieved EBI point statistics against full Oxford files."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pdv3.tabix_http import save_json, digest


def orientation(ea,oa,target_ea,target_oa):
    if (ea,oa)==(target_ea,target_oa):return 1
    if (ea,oa)==(target_oa,target_ea):return -1
    if all(a in {'A','C','G','T'} for a in [ea,oa,target_ea,target_oa]):
        complement=str.maketrans('ACGT','TGCA')
        if (ea.translate(complement),oa.translate(complement))==(target_ea,target_oa):return 1
        if (ea.translate(complement),oa.translate(complement))==(target_oa,target_ea):return -1
    return None


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--idps',required=True)
    args=parser.parse_args()
    native_path=ROOT/'data/big40/extracted/phase1_17_outcomes.csv'
    native=pd.read_csv(native_path)
    inst_path=ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv'
    inst=pd.read_csv(inst_path)
    inst=inst.drop_duplicates(['variant','chrom','pos37','pos38','ea_x','oa_x'])
    records=[];inputs={str(native_path.relative_to(ROOT)):digest(native_path.read_bytes()),
                       str(inst_path.relative_to(ROOT)):digest(inst_path.read_bytes())}
    for idp in args.idps.split(','):
        idp=idp.zfill(4);phenotype=f'BIG40:release2:discovery:{idp}'
        path=ROOT/f'data/big40_atlas/extracted/{idp}.csv.gz'
        inputs[str(path.relative_to(ROOT))]=digest(path.read_bytes())
        ebi=pd.read_csv(path,dtype={'chromosome':str})
        old=native[native.phenotype.eq(phenotype)]
        for row in inst.to_dict('records'):
            n=old[old.chrom_y.astype(str).str.replace('.0','',regex=False).replace({'X':'23'}).eq(str(row['chrom'])) & old.pos_y.eq(row['pos37'])]
            e=ebi[ebi.chromosome.replace({'X':'23'}).eq(str(row['chrom'])) & ebi.base_pair_location.eq(row['pos38'])]
            n=n[n.variant_y.eq(row['variant'])]
            e=e[e.variant_id.eq(row['variant']) | e.rsid.eq(row['variant'])]
            n=n[n.apply(lambda r:orientation(r.ea_y,r.oa_y,row['ea_x'],row['oa_x']) is not None,axis=1)] if len(n) else n
            e=e[e.apply(lambda r:orientation(r.effect_allele,r.other_allele,row['ea_x'],row['oa_x']) is not None,axis=1)] if len(e) else e
            rec={'idp':idp,'variant':row['variant'],'native_rows':len(n),'indexed_rows':len(e),
                 'exposure_ea':row['ea_x'],'exposure_oa':row['oa_x'],
                 'status':'not_shared','beta_abs_difference':None,'se_abs_difference':None,
                 'neglog10p_abs_difference':None}
            if len(n)==1 and len(e)==1:
                n,e=n.iloc[0],e.iloc[0]
                sign=orientation(n.ea_y,n.oa_y,e.effect_allele,e.other_allele)
                if sign is None:rec['status']='allele_mismatch'
                else:
                    rec.update(beta_abs_difference=abs(n.beta_y-sign*e.beta),
                               se_abs_difference=abs(n.se_y-e.standard_error),
                               neglog10p_abs_difference=abs(n.neg_log10_p_y-e.neg_log10_p_value))
                    rec['status']='matched' if max(rec[k] for k in ['beta_abs_difference','se_abs_difference','neglog10p_abs_difference'])<1e-10 else 'numerical_mismatch'
            elif len(n)>1 or len(e)>1:rec['status']='ambiguous_variant'
            records.append(rec)
    out=ROOT/'papers/01_mr_imaging_scan/results/full_atlas_20261116/bridge_audit'
    out.mkdir(parents=True,exist_ok=True)
    table=pd.DataFrame(records);table.to_csv(out/'variant_comparison.csv',index=False)
    shared=table[table.native_rows.eq(1)&table.indexed_rows.eq(1)]
    summary={'complete_comparison':True,'idps':args.idps.split(','),'compared_instrument_variants':len(table),
             'status_counts':table.status.value_counts().to_dict(),
             'all_shared_rows_match':bool(len(shared)>0 and shared.status.eq('matched').all()),
             'max_beta_abs_difference':float(shared.beta_abs_difference.max()),
             'max_se_abs_difference':float(shared.se_abs_difference.max()),
             'max_neglog10p_abs_difference':float(shared.neglog10p_abs_difference.max()),
             'inputs':inputs,'output_sha256':digest((out/'variant_comparison.csv').read_bytes()),
             'note':'Not-shared rows remain explicit; this is a source/extraction bridge check, not independent biological replication.'}
    save_json(out/'SUMMARY.json',summary)
    print(json.dumps(summary,indent=2))
    if not summary['all_shared_rows_match']:raise SystemExit(1)


if __name__=='__main__':main()
