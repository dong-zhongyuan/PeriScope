"""Compare the bounded-memory pipeline against the audited complete 45-IDP scan."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pdv3.tabix_http import save_json


def main():
    original=ROOT/'papers/01_mr_imaging_scan/results/expanded_20261116/combined_45'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,default=ROOT/'papers/01_mr_imaging_scan/results/full_atlas_20261116/validation45')
    parser.add_argument('--source-route-verification',type=Path)
    args=parser.parse_args()
    new=args.results.resolve()
    columns='exposure_id phenotype status beta se p log_p se_delta p_delta q_bh_global q_by_global q_bh_family q_bh_cis_family q_by_cis_family q_storey_lambda_0p25 q_storey_lambda_0p5 q_storey_lambda_0p75'.split()
    results={}
    for policy in ['drop','infer']:
        name=f'full_family_palindrome_{policy}.csv.gz'
        a=pd.read_csv(original/name,usecols=columns).set_index(['exposure_id','phenotype']).sort_index()
        b=pd.read_csv(new/name,usecols=columns).set_index(['exposure_id','phenotype']).sort_index()
        assert len(a)==len(b)==132345 and a.index.equals(b.index)
        assert a.status.equals(b.status)
        differences={}
        for c in columns[3:]:
            assert a[c].isna().equals(b[c].isna()),c
            difference=abs(a[c]-b[c]).max()
            assert difference<1e-10,(policy,c,difference)
            differences[c]=float(difference)
        results[policy]=dict(rows=len(a),statuses_identical=True,missingness_identical=True,max_absolute_differences=differences)
    state=json.loads((new/'RUN.json').read_text())
    if not state.get('complete'):raise RuntimeError('Validation run is incomplete')
    route_check=None
    if args.source_route_verification:
        path=args.source_route_verification.resolve()
        route=json.loads(path.read_text())
        if not route.get('complete'):raise RuntimeError('Source-route verification is incomplete')
        from pdv3.tabix_http import digest
        for relative,expected in route['code_sha256'].items():
            if digest((ROOT/relative).read_bytes())!=expected:raise RuntimeError(f'Source-route verification code changed: {relative}')
        for relative in ('papers/01_mr_imaging_scan/s18_run_full_atlas.py','pdv3/atlas_outcomes.py'):
            if route['code_sha256'][relative]!=state['source_code_sha256'][relative]:raise RuntimeError('Source-route and numerical validations use different implementations')
        route_check=dict(path=str(path.relative_to(ROOT)),sha256=digest(path.read_bytes()),
                         native_source_files_checked=route['native_source_files_checked'])
    save_json(new/'PIPELINE_VERIFICATION.json',dict(complete=True,checks=results,source_code_sha256=state['source_code_sha256'],
              source_route_verification=route_check,
              original='Audited 45-IDP full scan; independent family and numerical audits already passed',
              result='Batched estimates and whole-family corrections reproduce both original policies'))
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
