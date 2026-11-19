"""Stage the existing atlas's unadjusted chunks while acquisition continues.

The original s18 implementation and its final completion/correction gate stay
unchanged. Completed chunks can be hard-linked into its normal cache while the
main queue remains far from completion. This never publishes partial q-values.
"""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
from datetime import datetime, timezone

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import s18_run_full_atlas as atlas
from pdv3.tabix_http import save_json

RESULTS = ROOT / 'papers/01_mr_imaging_scan/results/full_atlas_20261116'
DATA = ROOT / 'data/big40_atlas'
VALIDATED = RESULTS / 'validation45_nativefix'
CACHE_VALIDATION = RESULTS / 'cache_validation_nativefix'
CACHE_REUSE = RESULTS / 'cache_reuse_validation_nativefix'
CACHE_FULL = RESULTS / 'cache_precompute_nativefix'


def utc():
    return datetime.now(timezone.utc).isoformat()


def context():
    validation = json.loads((VALIDATED / 'PIPELINE_VERIFICATION.json').read_text())
    if not validation['complete']:
        raise RuntimeError('Original statistical implementation is unverified')
    code = validation['source_code_sha256']
    for relative, expected in code.items():
        if atlas.sha(ROOT / relative) != expected:
            raise RuntimeError(f'Frozen statistical code changed: {relative}')
    files = [ROOT / f'data/{path}' for path in (
        'ukb_ppp/extracted/assay_coverage.csv',
        'ukb_ppp/extracted/cis_primary_selected.csv',
        'big40_variants/extracted/cis_instrument_frequencies.csv')]
    inputs = {str(p.relative_to(ROOT)): atlas.sha(p) for p in files}
    frozen_inputs = json.loads((VALIDATED / 'RUN.json').read_text())['input_sha256']
    if inputs != frozen_inputs:
        raise RuntimeError('Frozen exposure/frequency inputs changed')
    signature = hashlib.sha256(json.dumps(dict(code=code, inputs=inputs), sort_keys=True).encode()).hexdigest()
    cov, inst, af = [pd.read_csv(p) for p in files]
    return cov, inst, af, signature, code, inputs


def ready_sources(inst_sha):
    """Construct the exact full metadata, then take only completed source files.

    Indices refer to the original full catalog, never a compacted ready subset.
    s18 will independently reconstruct this catalog and check every identity.
    """
    plan = json.loads((DATA / 'PLAN.json').read_text())
    path = atlas.checked(DATA / 'phenotypes_3935.csv', plan['plan_sha256'])
    catalog = pd.read_csv(path, dtype={'idp': str})
    import csv
    with path.open() as stream:
        frozen = {r['idp']: r for r in csv.DictReader(stream)}
    catalog['trait'] = catalog.short_name
    catalog['modality'] = catalog['category']
    combined = ROOT / 'data/imaging_combined'
    manifest = json.loads((combined / 'MANIFEST.json').read_text())
    if manifest['complete'] is not True:
        raise RuntimeError('QSM extension source is incomplete')
    prior_meta = pd.read_csv(atlas.checked(combined / 'phenotypes_45.csv', manifest['phenotype_plan_sha256']))
    prior_y = pd.read_csv(atlas.checked(combined / 'extracted/outcomes_45.csv', manifest['output_sha256']))
    extra = prior_meta[prior_meta.phenotype.str.startswith('OxfordQSM2022:')].copy()
    if len(extra) != 28:
        raise RuntimeError('Incorrect QSM extension size')
    meta = pd.concat([catalog, extra], ignore_index=True)
    if len(meta) != 3963 or meta.phenotype.nunique() != 3963:
        raise RuntimeError('Full metadata dimensions differ')
    sources = []
    for i, record in enumerate(meta.to_dict('records')):
        if i >= len(catalog):
            source = dict(route='panel45', input_sha256=manifest['output_sha256'],
                y=prior_y[prior_y.phenotype.eq(record['phenotype'])].copy())
        else:
            key = record['idp']
            receipt_path = DATA / f'native/{key}.json'
            if not receipt_path.exists():
                receipt_path = DATA / f'extracted/{key}.json'
            if not receipt_path.exists():
                continue
            info = json.loads(receipt_path.read_text())
            if not info.get('complete'):
                continue
            if info['instrument_sha256'] != inst_sha or info['phenotype_record'] != frozen[key]:
                raise RuntimeError(f'Completed source identity differs: {key}')
            output = atlas.checked(DATA / info['output'], info['output_sha256'])
            source = dict(route='native' if receipt_path.parent.name == 'native' else 'ebi',
                path=str(output), input_sha256=info['output_sha256'],
                receipt_path=str(receipt_path.relative_to(ROOT)), receipt_sha256=atlas.sha(receipt_path))
        sources.append((i, record, source))
    return meta, sources


def stage_one(i, record, source, folder, cov, inst, af, signature):
    # Identity and numerical operations match the unchanged original loop.
    metadata_json = json.dumps(record, sort_keys=True, default=str)
    identity = hashlib.sha256((signature + source['input_sha256'] + metadata_json
        + source.get('receipt_sha256', '')).encode()).hexdigest()
    path = folder / f'{i:04d}.json'
    paths = {policy: folder / f'{i:04d}_{policy}.csv.gz' for policy in ('drop', 'infer')}
    old = json.loads(path.read_text()) if path.exists() else {}
    if old.get('identity') == identity and all(p.exists() and atlas.sha(p) == old.get('output_sha256', {}).get(policy) for policy, p in paths.items()):
        return old
    if source['route'] != 'panel45':
        atlas.checked(Path(source['path']), source['input_sha256'])
        atlas.checked(ROOT / source['receipt_path'], source['receipt_sha256'])
    y = atlas.load_outcomes(source, inst, af)
    if not set(y.phenotype).issubset({record['phenotype']}):
        raise RuntimeError('Outcome phenotype differs')
    pairs = atlas.build_pairs(cov, inst, y, pd.DataFrame([record]))
    entry = dict(identity=identity, phenotype=record['phenotype'], route=source['route'],
        input_sha256=source['input_sha256'], source_receipt=source.get('receipt_path'),
        source_receipt_sha256=source.get('receipt_sha256'), output_sha256={}, policies={})
    for policy in ('drop', 'infer'):
        full, _ = atlas.expanded_scan(pairs, palindrome=policy)
        full = full.reindex(columns=atlas.KEEP)
        if full.exposure_id.tolist() != cov.exposure_id.tolist() or len(full) != len(cov):
            raise RuntimeError('Incomplete assay product')
        if any(c.startswith('q_') for c in full.columns):
            raise RuntimeError('Partial correction must not be exported')
        entry['output_sha256'][policy] = atlas.write_csv(full, paths[policy])
        entry['policies'][policy] = full.status.value_counts().to_dict()
    save_json(path, entry)
    return entry


def import_guard():
    acquisition = json.loads((DATA / 'RUN.json').read_text())
    progress = json.loads((DATA / 'extraction_progress.json').read_text())
    analysis = json.loads((DATA / 'ANALYSIS_RUN.json').read_text())
    if acquisition['status'] != 'running' or analysis['status'] != 'waiting_for_complete_acquisition':
        raise RuntimeError('Main acquisition/statistics state no longer permits import')
    if progress['requested'] - progress['completed'] <= 128 or progress['complete']:
        raise RuntimeError('Stop cache import before the final 128 source completions')
    if (RESULTS / 'full/RUN.json').exists():
        raise RuntimeError('Main statistics have begun; cache import forbidden')
    for pid in (acquisition['pid'], analysis['pid']):
        os.kill(pid, 0)
    return progress['completed']


def link_cache(i, folder, destination):
    entry = json.loads((folder / f'{i:04d}.json').read_text())
    names = [f'{i:04d}_{policy}.csv.gz' for policy in ('drop', 'infer')] + [f'{i:04d}.json']
    # Completion JSON is linked last. No existing destination is overwritten.
    for name in names:
        source, target = folder / name, destination / name
        if target.exists():
            if atlas.sha(target) != atlas.sha(source):
                raise RuntimeError(f'Different cache already exists: {target}')
        else:
            os.link(source, target)
    for policy in ('drop', 'infer'):
        atlas.checked(destination / f'{i:04d}_{policy}.csv.gz', entry['output_sha256'][policy])


def run(mode, do_import=False):
    if do_import and mode != 'full':
        raise RuntimeError('Validation is isolated; cannot import it into full results')
    cov, inst, af, signature, code, inputs = context()
    if do_import:
        verified = json.loads((CACHE_VALIDATION / 'VERIFICATION.json').read_text())
        if not verified['complete'] or verified['precompute_script_sha256'] != atlas.sha(Path(__file__)):
            raise RuntimeError('Current cache implementation has not passed reuse validation')
        import_guard()
    base = CACHE_VALIDATION if mode == 'validation45' else CACHE_FULL
    base.mkdir(exist_ok=True)
    with (base / '.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        chunks = base / 'chunks'
        chunks.mkdir(exist_ok=True)
        if mode == 'validation45':
            meta, sources = atlas.source_tables(mode, inputs['data/ukb_ppp/extracted/cis_primary_selected.csv'])
            ready = [(i, record, source) for i, (record, source) in enumerate(zip(meta.to_dict('records'), sources))]
        else:
            meta, ready = ready_sources(inputs['data/ukb_ppp/extracted/cis_primary_selected.csv'])
        state = dict(status='running', pid=os.getpid(), started_at_utc=utc(), mode=mode,
            snapshot_ready_traits=len(ready), planned_traits=len(meta), completed_traits=0,
            imported_traits=0, precompute_script_sha256=atlas.sha(Path(__file__)),
            source_code_sha256=code, input_sha256=inputs, no_final_corrections=True,
            no_candidate_selection=True, main_scan_complete=False)
        save_json(base / 'RUN.json', state)
        records = []
        try:
            for i, record, source in ready:
                if shutil.disk_usage(base).free < 2 * 2**30:
                    raise RuntimeError('Less than 2 GiB free; stop staging, retain all completed caches')
                if do_import:
                    import_guard()
                entry = stage_one(i, record, source, chunks, cov, inst, af, signature)
                imported = False
                if do_import:
                    main_completed = import_guard()
                    destination = RESULTS / 'full/chunks'
                    destination.mkdir(parents=True, exist_ok=True)
                    link_cache(i, chunks, destination)
                    imported = True
                    state['main_queue_completed_at_latest_import'] = main_completed
                records.append(dict(index=i, phenotype=record['phenotype'], identity=entry['identity'],
                    output_sha256=entry['output_sha256'], imported=imported))
                state.update(completed_traits=len(records), imported_traits=sum(r['imported'] for r in records), updated_at_utc=utc())
                save_json(base / 'RUN.json', state)
                save_json(base / 'CHUNKS.json', records)
                if len(records) % 25 == 0:
                    print(f"Cached {len(records)}/{len(ready)}; imported {state['imported_traits']}", flush=True)
            context()
            state.update(status='snapshot_complete', finished_at_utc=utc())
        except Exception as error:
            state.update(status='stopped_with_reason', error=str(error), ended_at_utc=utc())
            save_json(base / 'RUN.json', state)
            raise
        save_json(base / 'RUN.json', state)
        print(json.dumps(state, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('validation45', 'full'), default='validation45')
    parser.add_argument('--import-main', action='store_true')
    args = parser.parse_args()
    run(args.mode, args.import_main)
