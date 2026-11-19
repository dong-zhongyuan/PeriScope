"""Exercise real native and EBI schemas through the full-atlas source loader."""
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import s18_run_full_atlas as atlas
from pdv3.tabix_http import save_json

FIELDS=['chrom_y','pos_y','ea_y','oa_y','beta_y','se_y','p_y','neg_log10_p_y',
        'variant_y','phenotype','source_y','build_y','n_y','n_y_status']
NUMERIC=['beta','se','p','log_p','se_delta','p_delta','log_p_delta','ci_lower',
         'ci_upper','F','beta_y_aligned','eaf_y_aligned']


def chrom(series):
    return series.astype(str).str.removesuffix('.0').replace({'23':'X'})


def main():
    data=ROOT/'data/big40_atlas'
    out=ROOT/'papers/01_mr_imaging_scan/results/full_atlas_20261116/native_route_fix'
    out.mkdir(exist_ok=True)
    original=ROOT/'papers/01_mr_imaging_scan/results/expanded_20261116/combined_45'
    cov_path=ROOT/'data/ukb_ppp/extracted/assay_coverage.csv'
    inst_path=ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv'
    freq_path=ROOT/'data/big40_variants/extracted/cis_instrument_frequencies.csv'
    native_original=ROOT/'data/big40/extracted/phase1_17_outcomes.csv'
    catalog_path=data/'phenotypes_3935.csv'
    cov=pd.read_csv(cov_path);inst=pd.read_csv(inst_path);af=pd.read_csv(freq_path)
    old_y=pd.read_csv(native_original)
    catalog=pd.read_csv(catalog_path,dtype={'idp':str})
    catalog['trait']=catalog.short_name
    meta={r['idp']:r for r in catalog.to_dict('records')}
    input_hashes={str(p.relative_to(ROOT)):atlas.sha(p) for p in (cov_path,inst_path,freq_path,native_original,catalog_path)}
    reference={}
    for policy in ('drop','infer'):
        path=original/f'full_family_palindrome_{policy}.csv.gz'
        reference[policy]=pd.read_csv(path,usecols=['exposure_id','phenotype','status',*NUMERIC])
        input_hashes[str(path.relative_to(ROOT))]=atlas.sha(path)
    native_checks=[];numeric_checks=[]
    for receipt_path in sorted((data/'native').glob('*.json')):
        receipt=json.loads(receipt_path.read_text())
        if not receipt.get('complete'):continue
        key=receipt['phenotype_record']['idp']
        source_path=atlas.checked(data/receipt['output'],receipt['output_sha256'])
        input_hashes[str(source_path.relative_to(ROOT))]=receipt['output_sha256']
        input_hashes[str(receipt_path.relative_to(ROOT))]=atlas.sha(receipt_path)
        raw=pd.read_csv(source_path)
        actual=atlas.load_outcomes(dict(route='native',path=str(source_path)),inst,af)
        # Independent positional reconstruction of the acquisition fields.
        expected=pd.DataFrame(dict(chrom_y=chrom(raw.chromosome),pos_y=raw.base_pair_location,
            ea_y=raw.effect_allele,oa_y=raw.other_allele,beta_y=raw.beta,se_y=raw.standard_error,
            p_y=raw.p_value,neg_log10_p_y=raw.neg_log10_p_value,variant_y=raw.variant_id,
            phenotype=raw.phenotype,source_y=raw.source_y,build_y=raw.build_y,
            n_y=raw.n_y,n_y_status=raw.n_y_status))
        pd.testing.assert_frame_equal(actual[FIELDS],expected[FIELDS],check_dtype=False,check_exact=True)
        previous=old_y.loc[old_y.phenotype.eq(meta[key]['phenotype'])].copy()
        is_reused=len(previous)>0
        if is_reused:
            previous['chrom_y']=chrom(previous.chrom_y)
            pd.testing.assert_frame_equal(actual[FIELDS].reset_index(drop=True),
                previous[FIELDS].reset_index(drop=True),check_dtype=False,check_exact=True)
            expected_af=atlas.attach_frequencies(previous,af)
            pd.testing.assert_frame_equal(actual[['eaf_y','eaf_attachment_status']].reset_index(drop=True),
                expected_af[['eaf_y','eaf_attachment_status']].reset_index(drop=True),check_dtype=False,check_exact=True)
        pairs=atlas.build_pairs(cov,inst,actual,pd.DataFrame([meta[key]]))
        counts={}
        for policy in ('drop','infer'):
            full,_=atlas.expanded_scan(pairs,palindrome=policy)
            assert len(full)==2941 and full.exposure_id.tolist()==cov.exposure_id.tolist()
            counts[policy]=full.status.value_counts().to_dict()
            if is_reused:
                old=reference[policy].loc[reference[policy].phenotype.eq(meta[key]['phenotype'])].set_index('exposure_id').sort_index()
                new=full.set_index('exposure_id').sort_index()
                assert old.index.equals(new.index) and old.status.equals(new.status)
                differences={}
                for field in NUMERIC:
                    assert old[field].isna().equals(new[field].isna())
                    difference=float(abs(old[field]-new[field]).max())
                    assert difference<1e-12,(key,policy,field,difference)
                    differences[field]=difference
                numeric_checks.append(dict(route='native',idp=key,policy=policy,rows=len(full),
                    all_statuses_match=True,max_absolute_differences=differences))
        native_checks.append(dict(idp=key,rows=len(actual),raw_values_preserved=True,
            matched_original_17_source=is_reused,status_counts=counts))
    assert sum(r['matched_original_17_source'] for r in native_checks)==17
    for key in ('0013','1438'):
        path=data/f'extracted/{key}.csv.gz'
        receipt_path=path.with_name(key+'.json')
        receipt=json.loads(receipt_path.read_text());atlas.checked(path,receipt['output_sha256'])
        input_hashes[str(path.relative_to(ROOT))]=receipt['output_sha256']
        input_hashes[str(receipt_path.relative_to(ROOT))]=atlas.sha(receipt_path)
        y=atlas.load_outcomes(dict(route='ebi',path=str(path)),inst,af)
        pairs=atlas.build_pairs(cov,inst,y,pd.DataFrame([meta[key]]))
        for policy in ('drop','infer'):
            full,_=atlas.expanded_scan(pairs,palindrome=policy)
            old=reference[policy].loc[reference[policy].phenotype.eq(meta[key]['phenotype'])].set_index('exposure_id').sort_index()
            new=full.set_index('exposure_id').sort_index()
            assert old.index.equals(new.index) and len(new)==2941
            shared=old.status.eq('estimated') & new.status.eq('estimated')
            differences={}
            for field in NUMERIC:
                difference=float(abs(old.loc[shared,field]-new.loc[shared,field]).max())
                assert difference<1e-12,(key,policy,field,difference)
                differences[field]=difference
            numeric_checks.append(dict(route='ebi',idp=key,policy=policy,planned_rows=len(full),
                shared_estimates=int(shared.sum()),max_absolute_differences=differences,
                nonshared_variants_retained_as_unestimated=True))
    for relative,digest in input_hashes.items():
        assert atlas.sha(ROOT/relative)==digest,relative
    result=dict(complete=True,created_at_utc=datetime.now(timezone.utc).isoformat(),
        issue='Original full-atlas native route passed acquisition columns directly to a normalized-schema frequency join, causing KeyError chrom_y.',
        fix='Explicit GRCh37 field conversion before unchanged exact-allele frequency attachment; no allele or effect reversal.',
        native_source_files_checked=len(native_checks),native_checks=native_checks,
        real_source_numerical_checks=numeric_checks,input_sha256=input_hashes,
        code_sha256={str(p.relative_to(ROOT)):atlas.sha(p) for p in
            (Path(__file__),ROOT/'pdv3/atlas_outcomes.py',ROOT/'papers/01_mr_imaging_scan/s18_run_full_atlas.py')},
        full_scan_complete=False)
    save_json(out/'NATIVE_AND_EBI_VERIFICATION.json',result)
    print(json.dumps(dict(complete=True,native_source_files_checked=len(native_checks),
        original_17_each_policy_all_assay_rows=17*2941,ebi_bridge_traits=2,
        maximum_numerical_difference=max(max(r['max_absolute_differences'].values()) for r in numeric_checks)),indent=2))


if __name__=='__main__':main()
