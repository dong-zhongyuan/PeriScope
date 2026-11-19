"""Prove unchanged s18 can reuse staged chunks and reproduce final tables."""
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import precompute_atlas_cache as precompute
from pdv3.tabix_http import save_json


def text_sha(path):
    digest = hashlib.sha256()
    with gzip.open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def identity(path):
    s = path.stat()
    return dict(inode=s.st_ino, mtime_ns=s.st_mtime_ns, sha256=precompute.atlas.sha(path))


def main():
    results = precompute.RESULTS
    original = precompute.VALIDATED
    stage = precompute.CACHE_VALIDATION
    reuse = precompute.CACHE_REUSE
    state = json.loads((stage / 'RUN.json').read_text())
    if state['status'] != 'snapshot_complete' or state['completed_traits'] != 45 or state['imported_traits'] != 0:
        raise RuntimeError('Isolated 45-trait staging is incomplete')
    script = Path(precompute.__file__)
    if precompute.atlas.sha(script) != state['precompute_script_sha256']:
        raise RuntimeError('Staging code changed after validation run')
    precompute.context()
    chunks = []
    for i in range(45):
        names = [f'{i:04d}.json'] + [f'{i:04d}_{policy}.csv.gz' for policy in ('drop', 'infer')]
        for name in names:
            a, b = original / 'chunks' / name, stage / 'chunks' / name
            if precompute.atlas.sha(a) != precompute.atlas.sha(b):
                raise AssertionError(f'Staged chunk differs from original: {name}')
            chunks.append(dict(name=name, sha256=precompute.atlas.sha(b)))
    reuse.mkdir(exist_ok=True)
    destination = reuse / 'chunks'
    destination.mkdir(exist_ok=True)
    for i in range(45):
        precompute.link_cache(i, stage / 'chunks', destination)
    before = {r['name']: identity(destination / r['name']) for r in chunks}
    with (reuse / 'reuse.log').open('w') as log:
        command = [sys.executable, '-u', str(ROOT / 'papers/01_mr_imaging_scan/s18_run_full_atlas.py'),
                   '--mode', 'validation45', '--output', str(reuse)]
        process = subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    if process.returncode != 0:
        raise RuntimeError(f'Original s18 failed reuse check: {process.returncode}')
    after = {r['name']: identity(destination / r['name']) for r in chunks}
    if before != after:
        raise AssertionError('Original reader rewrote a cache instead of reusing it')
    run = json.loads((reuse / 'RUN.json').read_text())
    if not run['complete'] or run['planned_pairs'] != 132345:
        raise AssertionError('Original reader did not finish the complete validation family')
    full_tables = {}
    for policy in ('drop', 'infer'):
        name = f'full_family_palindrome_{policy}.csv.gz'
        a, b = text_sha(original / name), text_sha(reuse / name)
        if a != b:
            raise AssertionError(f'Uncompressed complete table differs: {policy}')
        full_tables[policy] = dict(rows=132345, uncompressed_csv_sha256=a, exact_content_match=True)
    guards = {}
    with tempfile.TemporaryDirectory(prefix='atlas_cache_guards_') as tmp:
        tmp = Path(tmp)
        dest = tmp / 'different'
        dest.mkdir()
        existing = dest / '0000_drop.csv.gz'
        existing.write_bytes(b'different existing content')
        old = existing.read_bytes()
        try:
            precompute.link_cache(0, stage / 'chunks', dest)
        except RuntimeError:
            guards['existing_different_cache_refused'] = existing.read_bytes() == old and len(list(dest.iterdir())) == 1
        else:
            raise AssertionError('Existing cache was replaced')
        saved_data, saved_results = precompute.DATA, precompute.RESULTS
        precompute.DATA = tmp / 'data'
        precompute.RESULTS = tmp / 'results'
        precompute.DATA.mkdir()
        (precompute.RESULTS / 'full').mkdir(parents=True)
        acq = dict(status='running', pid=os.getpid())
        analysis = dict(status='waiting_for_complete_acquisition', pid=os.getpid())
        progress = dict(requested=3935, completed=3806, complete=False)
        def put():
            for name, value in [('RUN.json', acq), ('ANALYSIS_RUN.json', analysis), ('extraction_progress.json', progress)]:
                save_json(precompute.DATA / name, value)
        try:
            put()
            assert precompute.import_guard() == 3806
            guards['129_remaining_permitted'] = True
            for label, change in [
                ('last_128_refused', lambda: progress.update(completed=3807)),
                ('statistics_started_refused', lambda: analysis.update(status='running_statistics')),
                ('acquisition_terminal_refused', lambda: acq.update(status='acquisition_complete')),
                ('existing_full_run_refused', lambda: (precompute.RESULTS / 'full/RUN.json').write_text('{}')),
            ]:
                progress.update(completed=3806)
                analysis.update(status='waiting_for_complete_acquisition')
                acq.update(status='running')
                change()
                put()
                try:
                    precompute.import_guard()
                except RuntimeError:
                    guards[label] = True
                else:
                    raise AssertionError(f'Import guard failed: {label}')
        finally:
            precompute.DATA, precompute.RESULTS = saved_data, saved_results
    result = dict(complete=True, created_at_utc=datetime.now(timezone.utc).isoformat(),
        precompute_script_sha256=precompute.atlas.sha(script),
        verifier_sha256=precompute.atlas.sha(Path(__file__)),
        original_reader_sha256=precompute.atlas.sha(ROOT / 'papers/01_mr_imaging_scan/s18_run_full_atlas.py'),
        identical_chunk_payloads=90, identical_chunk_records=45,
        all_135_original_reader_cache_files_unchanged=True,
        complete_family_tables=full_tables, guards=guards, chunks=chunks,
        formal_full_scan_completed=False, scope='Isolated 45-trait cache interoperability validation')
    save_json(stage / 'VERIFICATION.json', result)
    print(json.dumps({k: v for k, v in result.items() if k != 'chunks'}, indent=2))


if __name__ == '__main__':
    main()
