"""Independently verify every presentation value against a completed atlas.

Recovery verifier: auxiliary F_x annotation alone permits a measured two-ULP
decimal-parser difference. Every MR value and correction remains exact.
Original failed verifier and its log are retained without edits.

Whole-table pandas grouping is deliberately separate from the production
per-trait Profile accumulator. No production MR/presentation function imported.
"""
import argparse
from datetime import datetime, timezone
import gc, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'papers/01_mr_imaging_scan/results/full_atlas_20261116'
BEST=['phenotype','trait','family','beta','se','ci_lower','ci_upper','p','log_p',
      'q_bh_global','q_by_global','q_bh_family','F']
F_ANNOTATION_AUDIT=[]
ANNOTATION=['exposure_id','variant','chrom','pos37','pos38','pav_consequence',
            'pav_coding_consequence','pav_annotation_scope','F_x']

def require(condition,message):
    if not condition: raise ValueError(message)

def sha(p):
    with p.open('rb') as s:return hashlib.file_digest(s,'sha256').hexdigest()

def read(p,**kwargs):return pd.read_csv(p,float_precision='round_trip',**kwargs)

def dense_profile(raw,cov):
    rows=raw.loc[raw.status.eq('estimated')].copy()
    rows['z']=rows.beta/rows.se
    result=cov[['exposure_id','gene','assay']].copy()
    result['n_planned']=result.exposure_id.map(raw.groupby('exposure_id').size())
    result['n_estimated']=result.exposure_id.map(rows.groupby('exposure_id').size()).fillna(0).astype(int)
    result['n_not_estimated']=result.n_planned-result.n_estimated
    result['z_min']=result.exposure_id.map(rows.groupby('exposure_id').z.min())
    result['z_max']=result.exposure_id.map(rows.groupby('exposure_id').z.max())
    best=rows.sort_values(['exposure_id','log_p','phenotype']).drop_duplicates('exposure_id')
    best=best[['exposure_id',*BEST]].rename(columns={c:'best_'+c for c in BEST})
    return result.merge(best,on='exposure_id',how='left',validate='one_to_one')

def assert_frame(actual,expected,keys,subset=False):
    require(len(actual)==len(expected),'Presentation row count differs')
    require(not actual.duplicated(keys).any() and not expected.duplicated(keys).any(),'Duplicate presentation key')
    if not subset:require(set(actual.columns)==set(expected.columns),'Presentation columns differ')
    a=actual.sort_values(keys).reset_index(drop=True)[expected.columns]
    b=expected.sort_values(keys).reset_index(drop=True)
    if 'F_x' in a.columns:
        # The producer reads the auxiliary instrument annotation with pandas' default
        # parser; this independent verifier uses round_trip. Decimal parsing can
        # differ by up to two binary64 units. MR estimates/q-values remain exact.
        av=a.pop('F_x').to_numpy(float);bv=b.pop('F_x').to_numpy(float)
        require(np.array_equal(np.isnan(av),np.isnan(bv)),'F_x missingness differs')
        finite=np.isfinite(av)&np.isfinite(bv)
        require(np.array_equal(np.isfinite(av),np.isfinite(bv)),'F_x finite mask differs')
        errors=np.abs(av[finite]-bv[finite]);spacing=np.abs(np.spacing(bv[finite]))
        require(np.all(errors<=2*spacing),'F_x differs beyond two binary64 ULPs')
        F_ANNOTATION_AUDIT.append(dict(rows=int(finite.sum()),differing=int((errors>0).sum()),
            maximum_absolute_difference=float(errors.max(initial=0.)),
            maximum_ULP=float(np.max(np.divide(errors,spacing),initial=0.))))
    pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=True)

def complete_run(folder,mode):
    path=folder/'RUN.json'
    require(path.is_file(),'No completed atlas; presentation verification disabled')
    run=json.loads(path.read_text())
    n=3963 if mode=='full' else 45
    require(run.get('complete') is True and run.get('status')=='complete' and run.get('mode')==mode
      and run.get('assays')==2941 and run.get('phenotypes')==n and run.get('completed_traits')==n
      and run.get('planned_pairs')==2941*n,'Atlas does not cover the complete declared scope')
    return run

def metadata_qc(meta,mode):
    require(not meta.phenotype.duplicated().any(),'Duplicate metadata phenotype')
    qc=meta.family.eq('QC')
    for row in meta.to_dict('records'):
        flag=str(row.get('is_qc','')).lower()
        if flag in ('true','false'):
            require((flag=='true')==(row['family']=='QC'),'QC metadata disagreement')
        else:
            qsm=str(row['phenotype']).startswith('OxfordQSM2022:') and row['family'] in ('QSM','T2star')
            require(mode=='validation45' or qsm,'Missing frozen QC designation')
            require(row['family']!='QC','QC designation missing')
    require(int(qc.sum())==(16 if mode=='full' else 0),'Wrong number of frozen QC traits')
    return set(meta.loc[qc,'phenotype'])

def compare_policy(raw,cov,meta,qc_ids,annotation,source,policy,display=None):
    require(len(raw)==len(cov)*len(meta),'Not the complete Cartesian family')
    require(not raw.status.isna().any(),'Missing status')
    require(set(raw.phenotype)==set(meta.phenotype),'Source phenotype coverage differs')
    require(set(raw.exposure_id)==set(cov.exposure_id),'Source assay coverage differs')
    require(raw.groupby('phenotype').size().eq(len(cov)).all(),'Unequal phenotype coverage')
    bio=raw.loc[~raw.phenotype.isin(qc_ids)]
    expected=dense_profile(bio,cov).sort_values(['best_log_p','exposure_id'],na_position='last').reset_index(drop=True)
    ranks=pd.Series(pd.NA,index=expected.index,dtype='Int64');valid=expected.n_estimated.gt(0)
    ranks.loc[valid]=np.arange(1,int(valid.sum())+1);expected.insert(0,'descriptive_rank',ranks)
    expected=expected.merge(annotation,on='exposure_id',how='left',validate='one_to_one',sort=False)
    actual=read(source/f'all_assay_ranking_{policy}.csv')
    require(actual.exposure_id.tolist()==expected.exposure_id.tolist(),'Assay rank order differs')
    assert_frame(actual,expected,['exposure_id'])
    profiles=[]
    for family,rows in raw.groupby('family',sort=False):
        group=dense_profile(rows,cov);group.insert(3,'imaging_family',family);group.insert(4,'is_qc',family=='QC')
        profiles.append(group)
    expected_profiles=pd.concat(profiles,ignore_index=True)
    assert_frame(read(source/f'all_assay_category_profiles_{policy}.csv.gz'),expected_profiles,['imaging_family','exposure_id'])
    coverage=raw.groupby(['phenotype','family','status']).size().reset_index(name='n_assays')
    coverage.insert(2,'is_qc',coverage.phenotype.isin(qc_ids))
    assert_frame(read(source/f'coverage_by_phenotype_{policy}.csv'),coverage,['phenotype','family','status'])
    if policy=='drop':
        display=expected.loc[valid].head(24).copy();display.insert(0,'display_row',range(1,len(display)+1))
        actual_display=read(source/'overview_display_assays.csv')
        require(actual_display.exposure_id.tolist()==display.exposure_id.tolist(),'Primary display selection differs')
        assert_frame(actual_display,display,['exposure_id'])
    require(display is not None,'Primary display must be verified before infer')
    selected=expected_profiles.loc[~expected_profiles.is_qc & expected_profiles.exposure_id.isin(display.exposure_id)].copy()
    selected['display_row']=selected.exposure_id.map(dict(zip(display.exposure_id,display.display_row)))
    assert_frame(read(source/f'overview_category_profiles_{policy}.csv'),selected,['imaging_family','exposure_id'])
    report=dict(planned_rows=len(raw),all_assays=len(actual),ranked_assays=int(valid.sum()),
      all_category_profiles=len(expected_profiles),coverage_rows=len(coverage),display_assays=len(display),
      biological_phenotypes=len(meta)-len(qc_ids),qc_phenotypes=len(qc_ids),
      independent_MR_values_exact=True,all_nonfloating_annotations_exact=True,F_x_annotation_tolerance_binary64_ULP=2,qc_excluded_from_ranking=True,
      qc_category_and_coverage_retained=True,missing_assays_retained=True)
    return report,display

def verify(mode,output):
    folder=BASE/mode;run=complete_run(folder,mode)
    require(not output.exists(),'Report exists; inspect before replacing')
    source=folder/'presentation';mp=source/'MANIFEST.json';manifest=json.loads(mp.read_text())
    require(manifest.get('complete') is True and manifest.get('mode')==mode
      and manifest.get('formal_full_atlas_result')==(mode=='full')
      and manifest.get('historical_identity_used') is False
      and manifest.get('full_original_family_denominators_preserved') is True,'Presentation manifest scope differs')
    evidence={**manifest['inputs'],**manifest['code_sha256'],str(mp.relative_to(ROOT)):sha(mp)}
    expected_outputs={'overview_display_assays.csv'}
    for p in ('drop','infer'):
        expected_outputs|={f'all_assay_ranking_{p}.csv',f'all_assay_category_profiles_{p}.csv.gz',
          f'coverage_by_phenotype_{p}.csv',f'overview_category_profiles_{p}.csv'}
    require(set(manifest['outputs'])==expected_outputs,'Presentation artifact family differs')
    for name,digest in manifest['outputs'].items():evidence[str((source/name).relative_to(ROOT))]=digest
    if mode=='full':
        vp=folder/'COMPLETE_RESULT_VERIFICATION.json';v=json.loads(vp.read_text())
        require(v.get('complete') is True and v.get('formal_full_atlas_verified') is True and v.get('mode')=='full','Full numerical verification missing')
        require(all(v['policies'][p]['rows']==11655183 for p in ('drop','infer')),'Full verification scope differs')
        for rel,digest in manifest['inputs'].items():require(v['inputs_sha256'].get(rel)==digest,'Presentation/result verification inputs differ: '+rel)
        evidence.update(v['inputs_sha256']);evidence[str(vp.relative_to(ROOT))]=sha(vp)
        evidence['papers/01_mr_imaging_scan/verify_complete_atlas.py']=v['verifier_sha256']
    for rel,digest in evidence.items():require(sha(ROOT/rel)==digest,'Input SHA mismatch: '+rel)
    cov_path=ROOT/'data/ukb_ppp/extracted/assay_coverage.csv';inst_path=ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv'
    for p in [cov_path,inst_path,folder/'phenotypes.csv']:
        rel=str(p.relative_to(ROOT));require(rel in evidence,'Missing source provenance: '+rel)
    cov=read(cov_path);inst=read(inst_path);meta=read(folder/'phenotypes.csv')
    require(len(cov)==2941 and not cov.exposure_id.duplicated().any(),'Assay inventory differs')
    require(len(meta)==run['phenotypes'],'Metadata dimensions differ')
    qc_ids=metadata_qc(meta,mode)
    annotation=cov[['exposure_id']].merge(inst[ANNOTATION],on='exposure_id',how='left',validate='one_to_one')
    reports={};display=None;verifier_sha=sha(Path(__file__))
    for policy in ('drop','infer'):
        raw=read(folder/f'full_family_palindrome_{policy}.csv.gz',usecols=['exposure_id','status',*BEST])
        reports[policy],display=compare_policy(raw,cov,meta,qc_ids,annotation,source,policy,display)
        print(f'Verified {policy} presentation: {len(raw):,} source pairs',flush=True)
        del raw;gc.collect()
    for rel,digest in evidence.items():require(sha(ROOT/rel)==digest,'Input changed during verification: '+rel)
    require(sha(Path(__file__))==verifier_sha,'Verifier changed during verification')
    report=dict(complete=True,mode=mode,formal_full_presentation_verified=mode=='full',
      verified_at_utc=datetime.now(timezone.utc).isoformat(),verifier_sha256=verifier_sha,
      source_manifest_sha256=sha(mp),inputs_sha256=evidence,policies=reports,
       auxiliary_F_x_parser_comparisons=F_ANNOTATION_AUDIT,
      production_MR_or_presentation_functions_imported=False,new_MR_or_correction=False,
      target_nomination=False,scope='All ranking, best-pair values, z ranges, annotations, category/coverage rows and displayed subsets for both policies; not causal validation or a visual review.')
    output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--mode',choices=['full','validation45'],default='full');parser.add_argument('--output',type=Path)
    args=parser.parse_args();out=args.output or BASE/args.mode/'presentation/COMPLETE_PRESENTATION_VERIFICATION.json'
    report=verify(args.mode,out.resolve());print(json.dumps(dict(complete=report['complete'],mode=args.mode,output=str(out)),indent=2))
