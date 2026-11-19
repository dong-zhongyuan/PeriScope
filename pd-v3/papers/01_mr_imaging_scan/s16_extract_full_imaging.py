"""Extract every supplied cis position from each unfiltered full-catalog GWAS.

Discovery study/assembly/name checks precede all queries. All high-P, multiallelic
and absent records survive in an explicit audit. No biological decisions here.
"""
import argparse
import concurrent.futures
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.tabix_http import RemotePointTable, HTTPStore, atomic, save_json, digest
from pdv3.native_imaging import stream_oxford, FIELDS as NATIVE_FIELDS

RUN_CODE={str(p.relative_to(ROOT)):digest(p.read_bytes()) for p in
          [Path(__file__),ROOT/'pdv3/tabix_http.py',ROOT/'pdv3/native_imaging.py']}


def normalize_name(value):
    return ' '.join(value.replace('_', ' ').strip().split())


def validate_metadata(raw, meta):
    text = raw.decode()
    if not re.search(r'^gwas_id:\s*' + re.escape(meta['accession']) + r'\s*$', text, re.M):
        raise ValueError('Wrong study accession in metadata')
    if not re.search(r'^genome_assembly:\s*GRCh38\s*$', text, re.M):
        raise ValueError('Harmonized assembly is not verified as GRCh38')
    if not re.search(r'^is_harmonised:\s*true\s*$', text, re.M):
        raise ValueError('Source is not marked harmonized')
    if not re.search(r'^coordinate_system:\s*1-based\s*$', text, re.M):
        raise ValueError('Source is not marked one-based')
    m = re.search(r'^trait_description:\s*\n\s*-\s*(.+)', text, re.M)
    if not m or normalize_name(m[1].strip('"\'')) != normalize_name(meta['short_name']):
        raise ValueError(f'Trait-name mapping mismatch: {m[1] if m else None} != {meta["short_name"]}')
    # EBI sometimes supplies the maximum study n, not IDP-specific n.
    # Keep exact phenotype n from the official Oxford catalog, and record EBI n separately.
    n = re.findall(r'^\s+sample_size:\s*(\d+)\s*$', text, re.M)
    if not n or any(int(x) > 22138 for x in n):
        raise ValueError('Discovery cohort not compatible with 22,138 maximum sample')
    return {'verified_trait': m[1], 'build': 'GRCh38', 'ebi_sample_sizes': list(map(int,n)),
            'n_y_source': 'Oxford exact IDP catalog; X n remains missing', 'metadata_sha256': digest(raw)}


def csv_bytes(rows, fields):
    handle = io.StringIO(newline='')
    writer = csv.DictWriter(handle, fieldnames=fields, lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return handle.getvalue().encode()


def prune_payloads(cache):
    # Only reproducible cache payloads made by this new reader are evicted.
    # Their URL/range/ETag/SHA receipts and final extracted records remain.
    for receipt in cache.glob('*.json'):
        info=json.loads(receipt.read_text())
        if info.get('ranges') or info.get('url','').endswith('.tbi'):
            receipt.with_suffix('.bin').unlink(missing_ok=True)


def write_native(meta,instruments,dest,instrument_sha,rows,receipt):
    key=meta['idp'];output=dest/'native'/f'{key}.csv.gz'
    data=gzip.compress(csv_bytes(rows,NATIVE_FIELDS),mtime=0);atomic(output,data)
    result={'complete':True,'phenotype_record':meta,'instrument_sha256':instrument_sha,
            'output':str(output.relative_to(dest)),'output_sha256':digest(data),'output_rows':len(rows),
            'build':'GRCh37','source_route':'official_Oxford_native','source_receipt':receipt,
            'source_code_sha256':RUN_CODE,'completed_at_utc':datetime.now(timezone.utc).isoformat(),
            'numerical_filter':'none'}
    save_json(dest/'native'/f'{key}.json',result)
    return result


def reuse_native_inputs(plan,instruments,dest,instrument_sha):
    """Reuse the already verified 17 full source extractions without altering them."""
    source=ROOT/'data/big40/extracted/phase1_17_outcomes.csv'
    manifest=ROOT/'data/big40/PHASE1_17_MANIFEST.json'
    if not source.exists() or not manifest.exists():return {}
    receipt=json.loads(manifest.read_text())
    if (receipt.get('complete') is not True or receipt.get('output_sha256')!=digest(source.read_bytes())
            or receipt.get('exposure_sha256')!=instrument_sha):
        raise ValueError('Existing native extraction identity failed')
    grouped={}
    for r in csv.DictReader(source.open()):grouped.setdefault(r['phenotype'],[]).append(r)
    native={}
    for meta in plan:
        if meta['phenotype'] not in grouped:continue
        rows=[]
        for r in grouped[meta['phenotype']]:
            rows.append(dict(zip(NATIVE_FIELDS,[r['chrom_y'].replace('X','23'),r['pos_y'],r['ea_y'],r['oa_y'],
                r['beta_y'],r['se_y'],r['p_y'],r['neg_log10_p_y'],r['variant_y'],r['variant_y'],
                r['phenotype'],r['source_y'],'GRCh37',r['n_y'],r['n_y_status']])))
        native[meta['idp']]=write_native(meta,instruments,dest,instrument_sha,rows,
            {'existing_source':str(source.relative_to(ROOT)),'existing_source_sha256':digest(source.read_bytes()),
             'acquisition_manifest':str(manifest.relative_to(ROOT)),
             'acquisition_manifest_sha256':digest(manifest.read_bytes()),
             'input_instruments_sha256':instrument_sha,'native_original_files_previously_fully_validated':True})
    return native


def work(meta, instruments, dest, instrument_sha, metadata_only=False, keep_blocks=False):
    started = time.monotonic()
    key = meta['idp']
    native_final=dest/'native'/f'{key}.json'
    if native_final.exists() and not metadata_only:
        old=json.loads(native_final.read_text());native_output=dest/old['output']
        if (old.get('complete') and old.get('instrument_sha256')==instrument_sha and
                old.get('phenotype_record')==meta and native_output.exists() and
                digest(native_output.read_bytes())==old.get('output_sha256')):
            return old
        raise ValueError(f'Existing native extraction identity differs: {key}')
    final = dest/'extracted'/f'{key}.json'
    output = dest/'extracted'/f'{key}.csv.gz'
    if final.exists() and not metadata_only:
        old = json.loads(final.read_text())
        if (old.get('complete') and old.get('instrument_sha256') == instrument_sha and
                old.get('phenotype_record') == meta and output.exists() and
                digest(output.read_bytes()) == old.get('output_sha256')):
            if not keep_blocks:prune_payloads(dest/'download'/key)
            return old
        raise ValueError(f'Existing extraction identity differs: {key}')
    cache = dest/'download'/key
    store = HTTPStore(cache)
    try:
        raw, receipt = store.get(meta['metadata_url'])
    except RuntimeError as error:
        if 'HTTP 404' not in str(error):raise
        url=f'https://open.oxcin.ox.ac.uk/ukbiobank/big40/release2/stats/{key}.txt.gz'
        if metadata_only:
            with store.session.head(url,timeout=(20,30),allow_redirects=True) as response:
                response.raise_for_status()
                info={'url':url,'status':response.status_code,'size_bytes':int(response.headers['Content-Length']),
                      'etag':response.headers.get('ETag'),'last_modified':response.headers.get('Last-Modified')}
            result={'phenotype_record':meta,'metadata_verified':True,'complete':False,
                    'source_route':'official_Oxford_native','indexed_metadata_status':'404_unavailable',
                    'native_source_head':info,'build':'GRCh37',
                    'trait_identity_source':'Exact IDP in frozen official Oxford catalog',
                    'source_code_sha256':RUN_CODE}
            save_json(dest/'metadata'/f'{key}.json',result)
            return result
        print(json.dumps({'idp':key,'route':'native_stream_fallback','reason':'EBI harmonized source unavailable'}),flush=True)
        rows,native_receipt=stream_oxford(url,[(r['chrom'],int(r['pos37'])) for r in instruments],meta)
        return write_native(meta,instruments,dest,instrument_sha,rows,native_receipt)
    checked = validate_metadata(raw, meta)
    if metadata_only:
        result = dict(phenotype_record=meta, metadata_validation=checked, metadata_receipt=receipt,
                      metadata_verified=True, complete=False)
        save_json(dest/'metadata'/f'{key}.json', result)
        return result
    table = RemotePointTable(meta['url'], cache)
    points = {(str(int(r['chrom'])), int(r['pos38'])) for r in instruments}
    # Harmonized EBI uses chromosome 23, not X. Use the index's actual label.
    if 'X' in table.refs and '23' not in table.refs:
        points = {('X' if c == '23' else c, p) for c,p in points}
    rows, query = table.fetch_many(points)
    fields = list(table.header) + ['phenotype', 'source_y', 'build_y', 'n_y', 'n_y_status']
    for row in rows:
        c = row['chromosome']
        row.update(phenotype=meta['phenotype'], source_y=meta['url'], build_y='GRCh38',
                   n_y=meta['n_discovery'] if c not in ('X','23') else '',
                   n_y_status='catalog_autosomal_exact' if c not in ('X','23') else 'X_not_assigned')
    rows.sort(key=lambda r: (int(r['chromosome']) if r['chromosome'] != 'X' else 23,
                             int(r['base_pair_location']), r['effect_allele'], r['other_allele']))
    content = gzip.compress(csv_bytes(rows, fields), mtime=0)
    atomic(output, content)
    found = {(r['chromosome'],int(r['base_pair_location'])) for r in rows}
    result = dict(complete=True, phenotype_record=meta, metadata_validation=checked,
                  metadata_receipt=receipt, instrument_sha256=instrument_sha,
                  output=str(output.relative_to(dest)), output_sha256=digest(content),
                  queried_positions=len(points), positions_with_rows=len(found),
                  missing_positions=sorted([list(p) for p in points - found]),
                  output_rows=len(rows), query=query,
                  numerical_filter='none', variant_selection='all supplied exposure cis positions, irrespective of outcome P',
                  source_checks='retrieved BGZF blocks CRC checked; complete file MD5 NOT verified by sparse reading',
                  source_code_sha256=RUN_CODE,
                  missing_chromosomes=[c for c in sorted({str(int(r['chrom'])) for r in instruments})
                                       if c not in table.refs and not (c=='23' and 'X' in table.refs)],
                  completed_at_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=round(time.monotonic()-started,2))
    save_json(final, result)
    if not keep_blocks:prune_payloads(cache)
    print(json.dumps({'idp':key, 'rows':len(rows), 'seconds':result['elapsed_seconds'],
                      'compressed_bytes':query['compressed_bytes_requested']}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,default=ROOT/'data/big40_atlas/phenotypes_3935.csv')
    parser.add_argument('--instruments',type=Path,default=ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv')
    parser.add_argument('--output',type=Path,default=ROOT/'data/big40_atlas')
    parser.add_argument('--idps',help='Comma-separated exact IDs; omission means the entire frozen catalog')
    parser.add_argument('--workers',type=int,default=3)
    parser.add_argument('--metadata-only',action='store_true')
    parser.add_argument('--keep-blocks',action='store_true',help='Keep sparse source byte payloads after successful extraction; receipts always retained')
    parser.add_argument('--reuse-native',action='store_true',help='Use the 17 previously verified native Oxford extractions when available')
    args=parser.parse_args()
    plan=list(csv.DictReader(args.plan.open()))
    if args.idps:
        requested={x.zfill(4) for x in args.idps.split(',')}
        plan=[r for r in plan if r['idp'] in requested]
        if len(plan)!=len(requested):raise ValueError('Requested IDP missing from frozen plan')
    inst=list(csv.DictReader(args.instruments.open()))
    sha=digest(args.instruments.read_bytes())
    args.output.mkdir(parents=True,exist_ok=True)
    mode='metadata' if args.metadata_only else 'extraction'
    with (args.output/f'.{mode}.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        errors=[];done=[]
        native=reuse_native_inputs(plan,inst,args.output,sha) if args.reuse_native and not args.metadata_only else {}
        done.extend(native)
        state=args.output/f'{mode}_progress.json'
        def snapshot():
            save_json(state,{'complete':len(done)==len(plan) and not errors,
                            'metadata_only':args.metadata_only,'requested':len(plan),'completed':len(done),
                            'errors':errors,'completed_idps':sorted(done),
                            'plan_sha256':digest(args.plan.read_bytes()),'instrument_sha256':sha,
                            'updated_at_utc':datetime.now(timezone.utc).isoformat(),
                            'source_code_sha256':RUN_CODE})
        snapshot()
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            jobs={pool.submit(work,r,inst,args.output,sha,args.metadata_only,args.keep_blocks):r for r in plan if r['idp'] not in native}
            for future in concurrent.futures.as_completed(jobs):
                meta=jobs[future]
                try:future.result();done.append(meta['idp'])
                except Exception as error:
                    errors.append({'idp':meta['idp'],'error':str(error)})
                    print(json.dumps(errors[-1]),flush=True)
                snapshot()
                if len(done)%50==0:print(json.dumps({'completed':len(done),'requested':len(plan),'errors':len(errors)}),flush=True)
        snapshot()
        if errors:raise SystemExit(1)


if __name__=='__main__':main()
