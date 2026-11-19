"""Bounded operational pilot: same checked reader/positions, HTTP batching only.

No output is imported into the running main scan by this program. The source
identity and exact canonical extracted bytes must match the completed control.
"""
from datetime import datetime, timezone
import csv
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.tabix_http import RemotePointTable, HTTPStore, digest, save_json
from s16_extract_full_imaging import validate_metadata, csv_bytes, RUN_CODE

BASE = ROOT/'data/big40_batching_probe'


def run(batch, name, expected, meta, instruments):
    cache = BASE/'download'/name
    if cache.exists():
        raise ValueError('Use a fresh pilot directory; cached timing is not a network comparison')
    start = time.monotonic()
    store = HTTPStore(cache)
    metadata, metadata_receipt = store.get(meta['metadata_url'])
    validate_metadata(metadata, meta)
    if digest(metadata) != expected['metadata_validation']['metadata_sha256']:
        raise ValueError('Metadata changed since control')
    table = RemotePointTable(meta['url'], cache)
    if digest(table.index_bytes) != expected['query']['index_sha256']:
        raise ValueError('Tabix index changed since control')
    for field in ['etag', 'last_modified', 'size_bytes', 'sha256']:
        if table.identity[field] != expected['query']['source_identity'][field]:
            raise ValueError(f'Source identity differs: {field}')
    points = {(str(int(r['chrom'])), int(r['pos38'])) for r in instruments}
    if 'X' in table.refs and '23' not in table.refs:
        points = {('X' if c == '23' else c, p) for c, p in points}
    rows, query = table.fetch_many(points, batch_size=batch)
    fields = list(table.header)+['phenotype', 'source_y', 'build_y', 'n_y', 'n_y_status']
    for row in rows:
        c = row['chromosome']
        row.update(phenotype=meta['phenotype'], source_y=meta['url'], build_y='GRCh38',
                   n_y=meta['n_discovery'] if c not in ('X','23') else '',
                   n_y_status='catalog_autosomal_exact' if c not in ('X','23') else 'X_not_assigned')
    rows.sort(key=lambda r: (int(r['chromosome']) if r['chromosome'] != 'X' else 23,
                             int(r['base_pair_location']), r['effect_allele'], r['other_allele']))
    output = gzip.compress(csv_bytes(rows, fields), mtime=0)
    identical = digest(output) == expected['output_sha256']
    if not identical:
        raise ValueError('Canonical extracted bytes differ from completed control')
    receipts = [json.loads(p.read_text()) for p in cache.glob('*.json')]
    result = dict(batch_size=batch, trial=name, elapsed_seconds=time.monotonic()-start,
                  output_rows=len(rows), output_sha256=digest(output), exact_control_match=identical,
                  HTTP_requests=len(receipts), HTTP_payload_bytes=sum(r['size_bytes'] for r in receipts),
                  compressed_genomic_bytes=query['compressed_bytes_requested'], query=query,
                  complete=True, completed_at_utc=datetime.now(timezone.utc).isoformat())
    (BASE/f'{name}.csv.gz').write_bytes(output)
    save_json(BASE/f'{name}.json', result)
    print(json.dumps({k:v for k,v in result.items() if k!='query'}), flush=True)
    return result


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    instrument_path = ROOT/'data/ukb_ppp/extracted/cis_primary_selected.csv'
    expected_path = ROOT/'data/big40_atlas/extracted/0444.json'
    expected = json.loads(expected_path.read_text())
    assert expected['complete'] and expected['instrument_sha256'] == digest(instrument_path.read_bytes())
    for path, fingerprint in expected['source_code_sha256'].items():
        assert digest((ROOT/path).read_bytes()) == fingerprint
    meta = expected['phenotype_record']
    instruments = list(csv.DictReader(instrument_path.open()))
    state = dict(complete=False, idp='0444', main_scan_modified=False,
                 planned_trials=[['batch128',128],['batch64',64]],
                 scope='Operational transfer pilot; same positions, reader and scientific selection',
                 control_receipt_sha256=digest(expected_path.read_bytes()),
                 source_code_sha256={**RUN_CODE, str(Path(__file__).relative_to(ROOT)):digest(Path(__file__).read_bytes())},
                 started_at_utc=datetime.now(timezone.utc).isoformat())
    save_json(BASE/'RUN.json',state)
    try:
        results = [run(batch,name,expected,meta,instruments) for name,batch in state['planned_trials']]
        state.update(complete=True, trials=results, completed_at_utc=datetime.now(timezone.utc).isoformat(),
                     inference_limit='Two sequential network timings under the live main load; not a randomized performance study')
        save_json(BASE/'RUN.json',state)
    except Exception as error:
        state.update(failed=True,error=repr(error))
        save_json(BASE/'RUN.json',state)
        raise


if __name__ == '__main__':
    main()
