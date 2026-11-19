"""Bounded-memory complete atlas MR; final corrections require every source.

Per-phenotype computations reuse the audited expanded_scan kernel. Their local
q-values are discarded. Final q-values are calculated only after the complete
declared Cartesian product has been written, including every missing slot.
"""
import argparse
import csv
import gc
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pdv3.atlas_outcomes import normalize_ebi,normalize_native
from pdv3.expanded import build_pairs,expanded_scan
from pdv3.multiplicity import adjust_pvalues,storey_qvalues
from pdv3.tabix_http import save_json
from s11_run_expanded import attach_frequencies

KEEP='exposure_id gene assay phenotype trait family pair_status status usable method has_cis_instrument n_instruments variant chrom pos37 pos38 build_x chrom_x pos_x ea_x oa_x beta_x se_x eaf_x build_y chrom_y pos_y variant_y ea_y oa_y beta_y se_y eaf_y n_y n_y_status beta_y_aligned eaf_y_aligned harmonization_status allele_flipped eaf_difference eaf_attachment_status hm_code orientation_recovery outcome_candidate_count rsid_agreement beta se p log_p ci_lower ci_upper se_delta p_delta log_p_delta ci_lower_delta ci_upper_delta F palindrome_policy'.split()


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def utc():return datetime.now(timezone.utc).isoformat()


def write_csv(df,path):
    temp=path.with_name(path.name+'.tmp')
    df.to_csv(temp,index=False,compression={'method':'gzip','mtime':0} if path.name.endswith('.gz') else None)
    temp.replace(path)
    return sha(path)


def checked(path,expected):
    if sha(path)!=expected:raise ValueError(f'Input hash differs: {path}')
    return path


def load_outcomes(source,instruments,frequencies):
    """Convert each acquisition schema before the common MR pairing step."""
    route=source['route']
    if route=='panel45':return attach_frequencies(source['y'],frequencies)
    raw=pd.read_csv(source['path'])
    if route=='ebi':return normalize_ebi(raw,instruments,frequencies)
    if route=='native':return attach_frequencies(normalize_native(raw),frequencies)
    raise ValueError(f'Unknown outcome source route: {route}')


def source_tables(mode,inst_sha):
    combined=ROOT/'data/imaging_combined'
    manifest=json.loads((combined/'MANIFEST.json').read_text())
    assert manifest['complete'] is True
    yp=checked(combined/'extracted/outcomes_45.csv',manifest['output_sha256'])
    mp=checked(combined/'phenotypes_45.csv',manifest['phenotype_plan_sha256'])
    prior_meta=pd.read_csv(mp); prior_y=pd.read_csv(yp)
    if mode=='validation45':
        return prior_meta,[dict(route='panel45',input_sha256=manifest['output_sha256'],
                    y=prior_y[prior_y.phenotype.eq(p)].copy()) for p in prior_meta.phenotype]
    atlas=ROOT/'data/big40_atlas'
    plan=json.loads((atlas/'PLAN.json').read_text())
    meta_path=checked(atlas/'phenotypes_3935.csv',plan['plan_sha256'])
    progress=json.loads((atlas/'extraction_progress.json').read_text())
    if not progress['complete'] or progress['requested']!=3935 or progress['completed']!=3935 or progress['errors']:
        raise RuntimeError('Full atlas extraction is incomplete; global MR publication is disabled')
    if progress['plan_sha256']!=plan['plan_sha256'] or progress['instrument_sha256']!=inst_sha:
        raise RuntimeError('Extraction belongs to a different plan or instrument set')
    catalog=pd.read_csv(meta_path,dtype={'idp':str})
    with meta_path.open() as f:frozen_catalog={r['idp']:r for r in csv.DictReader(f)}
    catalog['trait']=catalog.short_name
    catalog['modality']=catalog['category']
    sources=[]
    for rec in catalog.to_dict('records'):
        idp=rec['idp']
        path=atlas/f'native/{idp}.json'
        if not path.exists():path=atlas/f'extracted/{idp}.json'
        info=json.loads(path.read_text())
        assert info['complete'] and info['instrument_sha256']==inst_sha
        # Compare against the CSV's string record, preserving leading zeros and
        # all source metadata rather than relying on a phenotype name alone.
        if info['phenotype_record']!=frozen_catalog[idp]:raise ValueError('Source phenotype metadata differs from frozen plan')
        yp=checked(atlas/info['output'],info['output_sha256'])
        sources.append(dict(route='native' if path.parent.name=='native' else 'ebi',path=str(yp),
                            input_sha256=info['output_sha256'],receipt_path=str(path.relative_to(ROOT)),receipt_sha256=sha(path)))
    extra=prior_meta[prior_meta.phenotype.str.startswith('OxfordQSM2022:')].copy()
    assert len(extra)==28
    for p in extra.phenotype:
        sources.append(dict(route='panel45',input_sha256=manifest['output_sha256'],y=prior_y[prior_y.phenotype.eq(p)].copy()))
    meta=pd.concat([catalog,extra],ignore_index=True)
    assert len(meta)==3963 and meta.phenotype.nunique()==3963
    return meta,sources


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',choices=['full','validation45'],default='full')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    output=args.output or ROOT/f'papers/01_mr_imaging_scan/results/full_atlas_20261116/{args.mode}'
    output.mkdir(parents=True,exist_ok=True);chunks=output/'chunks';chunks.mkdir(exist_ok=True)
    binary=output/'arrays';binary.mkdir(exist_ok=True)
    cov_path=ROOT/'data/ukb_ppp/extracted/assay_coverage.csv'
    inst_path=ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv'
    freq_path=ROOT/'data/big40_variants/extracted/cis_instrument_frequencies.csv'
    cov=pd.read_csv(cov_path);inst=pd.read_csv(inst_path);af=pd.read_csv(freq_path)
    code_paths=[Path(__file__),Path(__file__).with_name('s11_run_expanded.py')]+[ROOT/f'pdv3/{p}.py' for p in ['atlas_outcomes','expanded','mr','multiplicity']]
    code={str(p.relative_to(ROOT)):sha(p) for p in code_paths}
    shared_signature={str(p.relative_to(ROOT)):sha(p) for p in [cov_path,inst_path,freq_path]}
    meta,sources=source_tables(args.mode,sha(inst_path))
    signature=hashlib.sha256(json.dumps(dict(code=code,inputs=shared_signature),sort_keys=True).encode()).hexdigest()
    started=utc();n_assays=len(cov);n_traits=len(meta);n=n_assays*n_traits
    write_csv(cov,output/'assay_coverage.csv');write_csv(inst,output/'selected_cis_instruments.csv');meta_sha=write_csv(meta,output/'phenotypes.csv')
    state=dict(complete=False,mode=args.mode,started_at_utc=started,planned_pairs=n,phenotypes=n_traits,assays=n_assays,
               source_code_sha256=code,input_sha256=shared_signature,phenotypes_sha256=meta_sha,completed_traits=0)
    save_json(output/'RUN.json',state)
    records=[]
    for i,(record,source) in enumerate(zip(meta.to_dict('records'),sources)):
        base=chunks/f'{i:04d}';manifest_path=base.with_suffix('.json')
        # Include the exact phenotype metadata in cache identity.
        metadata_json=json.dumps(record,sort_keys=True,default=str)
        identity=hashlib.sha256((signature+source['input_sha256']+metadata_json+source.get('receipt_sha256','')).encode()).hexdigest()
        old=json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        paths={policy:chunks/f'{i:04d}_{policy}.csv.gz' for policy in ['drop','infer']}
        reusable=old.get('identity')==identity and all(p.exists() and sha(p)==old.get('output_sha256',{}).get(policy) for policy,p in paths.items())
        if reusable:entry=old
        else:
            y=load_outcomes(source,inst,af)
            assert set(y.phenotype).issubset({record['phenotype']})
            pairs=build_pairs(cov,inst,y,pd.DataFrame([record]))
            entry=dict(identity=identity,phenotype=record['phenotype'],route=source['route'],input_sha256=source['input_sha256'],
                       source_receipt=source.get('receipt_path'),source_receipt_sha256=source.get('receipt_sha256'),output_sha256={},policies={})
            for policy in ['drop','infer']:
                full,_=expanded_scan(pairs,palindrome=policy)
                full=full.reindex(columns=KEEP)
                assert full.exposure_id.tolist()==cov.exposure_id.tolist() and len(full)==n_assays
                entry['output_sha256'][policy]=write_csv(full,paths[policy])
                entry['policies'][policy]=full.status.value_counts().to_dict()
            save_json(manifest_path,entry)
        records.append(entry)
        state.update(completed_traits=i+1,updated_at_utc=utc(),status='estimating')
        save_json(output/'RUN.json',state)
        if (i+1)%10==0:print(f'Estimates {i+1}/{n_traits}',flush=True)
    save_json(output/'SOURCE_RECEIPTS.json',records)
    summaries={}
    cis_assays=cov.exposure_id.isin(inst.exposure_id).to_numpy()
    cis=np.tile(cis_assays,n_traits)
    family_indices={f:np.concatenate([np.arange(i*n_assays,(i+1)*n_assays) for i in range(n_traits) if meta.iloc[i]['family']==f]) for f in meta.family.unique()}
    for policy in ['drop','infer']:
        state.update(status=f'correcting_{policy}',updated_at_utc=utc());save_json(output/'RUN.json',state)
        p=np.full(n,np.nan)
        counts={}
        for i in range(n_traits):
            df=pd.read_csv(chunks/f'{i:04d}_{policy}.csv.gz',usecols=['p','status'])
            assert len(df)==n_assays
            p[i*n_assays:(i+1)*n_assays]=df.p
            for k,v in df.status.value_counts().items():counts[k]=counts.get(k,0)+int(v)
        def save_array(name,value):
            path=binary/f'{policy}_{name}.npy';np.save(path,value,allow_pickle=False)
            return np.load(path,mmap_mode='r')
        q={};pi0s={}
        for label,method in [('q_bh_global','bh'),('q_by_global','by')]:q[label]=save_array(label,adjust_pvalues(p,method));gc.collect()
        for lam in [.25,.5,.75]:
            label=f'q_storey_lambda_{str(lam).replace(".","p")}';values,pi0=storey_qvalues(p,lam)
            pi0s[label.replace('q_storey','pi0')]=pi0
            q[label]=save_array(label,values);del values
        for label,method in [('q_bh_cis_family','bh'),('q_by_cis_family','by')]:
            values=np.full(n,np.nan);values[cis]=adjust_pvalues(p[cis],method)
            q[label]=save_array(label,values);del values
        values=np.full(n,np.nan)
        family_n=np.zeros(n_traits,dtype=int);family_observed=np.zeros(n_traits,dtype=int)
        for family,indices in family_indices.items():
            values[indices]=adjust_pvalues(p[indices],'bh')
            ix=meta.family.eq(family).to_numpy();family_n[ix]=len(indices);family_observed[ix]=np.isfinite(p[indices]).sum()
        q['q_bh_family']=save_array('q_bh_family',values);del values;gc.collect()
        target=output/f'full_family_palindrome_{policy}.csv.gz';temp=target.with_name(target.name+'.tmp')
        historical=[];coverage=[]
        observed=int(np.isfinite(p).sum());cis_n=int(cis.sum());cis_observed=int(np.isfinite(p[cis]).sum())
        with gzip.open(temp,'wt',compresslevel=6) as f:
            for i in range(n_traits):
                df=pd.read_csv(chunks/f'{i:04d}_{policy}.csv.gz');sl=slice(i*n_assays,(i+1)*n_assays)
                for name,array in q.items():df[name]=array[sl]
                for name,value in pi0s.items():df[name]=value
                df['family_n_global']=n;df['family_n_observed_global']=observed
                df['family_n']=family_n[i];df['family_n_observed']=family_observed[i]
                df['cis_family_n']=cis_n;df['cis_family_n_observed']=cis_observed
                df.to_csv(f,index=False,header=i==0)
                historical.append(df[df.gene.isin(['CD22','CD40','TFRC'])])
                c=df.groupby(['phenotype','family','status'],dropna=False).size().rename('n_assays').reset_index();coverage.append(c)
        temp.replace(target)
        write_csv(pd.concat(historical,ignore_index=True),output/f'three_historical_proteins_{policy}.csv')
        write_csv(pd.concat(coverage,ignore_index=True),output/f'coverage_by_phenotype_{policy}.csv')
        summaries[policy]=dict(status_counts=counts,estimated=counts.get('estimated',0),family_n=n,
                               cis_family_n=int(cis.sum()),output_sha256=sha(target))
        # Arrays are reproducible scratch; final tables and source chunks remain.
        q.clear();del p;gc.collect()
        for path in binary.glob(f'{policy}_*.npy'):path.unlink()
    state.update(complete=True,status='complete',finished_at_utc=utc(),policies=summaries,
                 correction_scope='All declared assays x all declared imaging entries; missing estimates retain slots',
                 category_correction='Official BIG40 catalog categories plus distinct QSM/T2* extension families',
                 interpretation='Exploratory one-lead MR; no automatic causal support decisions',
                 sample_overlap='UKB protein and imaging samples may overlap; covariance not estimated',
                 note='Per-trait chunk files contain unadjusted estimates only; final q values are in full_family tables')
    save_json(output/'RUN.json',state)
    print(json.dumps(state,indent=2),flush=True)


if __name__=='__main__':main()
