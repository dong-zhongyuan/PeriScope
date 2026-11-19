"""Validate isolated extraction outputs before reusing them in the ONE main scan.

No main progress file is edited. The existing main reader will encounter these
files in its original queue and perform its normal cache checks. Existing main
outputs are never replaced. The small pilot is deliberately far ahead of the
main queue; an already-started main query prevents importing that entry.
"""
import argparse
import csv
from datetime import datetime, timezone
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.tabix_http import save_json
from s16_extract_full_imaging import work


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def norm_chrom(c):
    return '23' if c == 'X' else str(int(c))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply', action='store_true')
    args = parser.parse_args()
    pilot = ROOT / 'data/big40_tail_pilot'
    dest = ROOT / 'data/big40_atlas'
    plan_path = pilot / 'PLAN.json'
    plan = json.loads(plan_path.read_text())
    run = json.loads((pilot / 'RUN.json').read_text())
    if run.get('status') != 'complete_verified_reference' or run.get('exit_code') != 0:
        raise RuntimeError('Pilot must finish and reproduce the existing reference before import')
    assert sha(plan_path) == run['plan_sha256']
    for relative, digest in plan['unchanged_source_sha256'].items():
        assert sha(ROOT / relative) == digest
    inst_path = ROOT / 'data/ukb_ppp/extracted/cis_primary_selected.csv'
    instruments = list(csv.DictReader(inst_path.open()))
    inst_sha = sha(inst_path)
    points = {(norm_chrom(r['chrom']), int(r['pos38'])) for r in instruments}
    meta = list(csv.DictReader((dest / 'phenotypes_3935.csv').open()))
    by_id = {r['idp']: r for r in meta}
    positions = {r['idp']: i for i, r in enumerate(meta)}
    assert len(plan['new_idps']) == len(set(plan['new_idps'])) == 8
    assert plan['reference_idp'] not in plan['new_idps']
    assert sha(pilot / 'extracted' / (plan['reference_idp'] + '.csv.gz')) == plan['reference_output_sha256']
    prepared = []
    # Validate every source before mutating any destination.
    for key in plan['new_idps']:
        record = json.loads((pilot / 'extracted' / (key + '.json')).read_text())
        output = pilot / record['output']
        assert record['complete'] and record['phenotype_record'] == by_id[key]
        assert record['instrument_sha256'] == inst_sha and record['source_code_sha256'] == plan['source_code_sha256']
        assert record['output'] == f'extracted/{key}.csv.gz'
        assert sha(output) == record['output_sha256'] and record['numerical_filter'] == 'none'
        with gzip.open(output, 'rt') as stream:
            rows = list(csv.DictReader(stream))
        found = {(norm_chrom(r['chromosome']), int(r['base_pair_location'])) for r in rows}
        missing = {(norm_chrom(c), int(p)) for c, p in record['missing_positions']}
        assert found.isdisjoint(missing) and found | missing == points
        assert len(found) == record['positions_with_rows'] and len(rows) == record['output_rows']
        assert record['queried_positions'] == record['query']['requested_positions'] == len(points)
        assert all(r['phenotype'] == by_id[key]['phenotype'] and r['source_y'] == by_id[key]['url']
                   and r['build_y'] == 'GRCh38' for r in rows)
        prepared.append((key, record))
    checks = dict(utc=datetime.now(timezone.utc).isoformat(), dry_run=not args.apply,
                  source_plan_sha256=sha(plan_path), import_code_sha256=sha(Path(__file__)),
                  all_eight_source_tables_validated=True, unchanged_main_progress_file=True,
                  no_new_MR_statistics=True, records=[])
    if not args.apply:
        print(json.dumps(checks), flush=True)
        return
    if (pilot / 'IMPORT.json').exists():
        raise RuntimeError('Import receipt already exists; inspect rather than overwrite')
    with (dest / '.tail_import.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        # A wide queue separation is required in addition to checking cache
        # activity. It avoids racing a main worker that has just started an IDP.
        progress = json.loads((dest / 'extraction_progress.json').read_text())
        for key, record in prepared:
            target = dest / record['output']
            target_receipt = dest / 'extracted' / (key + '.json')
            if target.exists() or target_receipt.exists():
                if not (target.exists() and target_receipt.exists() and sha(target) == record['output_sha256']):
                    raise RuntimeError(f'Main output already exists with different/incomplete content: {key}')
                continue
            if positions[key] - progress['completed'] < 128:
                raise RuntimeError('Main queue too close for a concurrent cache import')
            for receipt in (dest / 'download' / key).glob('*.json'):
                cached = json.loads(receipt.read_text())
                if cached.get('ranges') or cached.get('url', '').endswith('.tbi'):
                    raise RuntimeError(f'Main query has already started: {key}')
        for key, record in prepared:
            target = dest / record['output']
            target_receipt = dest / 'extracted' / (key + '.json')
            already_present = target_receipt.exists()
            if not already_present:
                # Exclusive atomic links cannot overwrite another process's
                # output; the completion receipt is installed last.
                os.link(pilot / record['output'], target)
                os.link(pilot / 'extracted' / (key + '.json'), target_receipt)
            accepted = work(by_id[key], instruments, dest, inst_sha)
            assert accepted['complete'] and accepted['output_sha256'] == record['output_sha256']
            checks['records'].append(dict(idp=key, imported=not already_present,
                output_sha256=sha(target), source_receipt_sha256=sha(pilot / 'extracted' / (key + '.json')),
                main_receipt_sha256=sha(target_receipt), canonical_main_reader_cache_check=True,
                source_HTTP_receipts=str((pilot / 'download' / key).relative_to(ROOT))))
            save_json(pilot / 'IMPORT_IN_PROGRESS.json', checks)
    checks.update(complete=True, newly_imported=sum(r['imported'] for r in checks['records']),
                  main_progress_updated_by_importer=False)
    save_json(pilot / 'IMPORT.json', checks)
    (pilot / 'IMPORT_IN_PROGRESS.json').unlink()
    print(json.dumps(checks, indent=2), flush=True)


if __name__ == '__main__':
    main()
