"""Freeze a finished same-scan cache snapshot after full payload verification."""
from collections import Counter
from datetime import datetime, timezone
import csv
import hashlib
import gzip
import json
import os
from pathlib import Path

import precompute_atlas_cache as pc


def stable_stat(path):
    s = path.stat()
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns)


def main():
    base = pc.CACHE_FULL
    run_bytes, chunk_bytes = (base/'RUN.json').read_bytes(), (base/'CHUNKS.json').read_bytes()
    run, records = json.loads(run_bytes), json.loads(chunk_bytes)
    assert run['status'] == 'snapshot_complete'
    assert run['completed_traits'] == run['imported_traits'] == run['snapshot_ready_traits'] == len(records)
    try:
        os.kill(run['pid'], 0)
    except ProcessLookupError:
        pass
    else:
        raise RuntimeError('Precompute process is still present; wait before freezing its snapshot')
    before_queue = pc.import_guard()
    stamp = datetime.fromisoformat(run['finished_at_utc']).strftime('%Y%m%dT%H%M%SZ')
    output = base/'snapshots'/stamp
    if output.exists():
        raise FileExistsError(f'Snapshot already recorded: {output}')
    cov, inst, af, signature, code, inputs = pc.context()
    assert code == run['source_code_sha256'] and inputs == run['input_sha256']
    assert pc.atlas.sha(Path(pc.__file__)) == run['precompute_script_sha256']
    meta, ready = pc.ready_sources(inputs['data/ukb_ppp/extracted/cis_primary_selected.csv'])
    sources = {i: (record, source) for i, record, source in ready}
    indices = [record['index'] for record in records]
    assert len(indices) == len(set(indices)) and indices == sorted(indices)
    assays = cov.exposure_id.tolist()
    files, routes, policies, source_checks = [], Counter(), Counter(), {}
    for number, record in enumerate(records, 1):
        i = record['index']
        metadata, source = sources[i]
        metadata_json = json.dumps(metadata, sort_keys=True, default=str)
        identity = hashlib.sha256((signature+source['input_sha256']+metadata_json
                                   + source.get('receipt_sha256', '')).encode()).hexdigest()
        assert record['imported'] and record['identity'] == identity
        entry_path = base/'chunks'/f'{i:04d}.json'
        entry = json.loads(entry_path.read_text())
        assert entry['identity'] == identity and entry['phenotype'] == metadata['phenotype'] == record['phenotype']
        assert entry['output_sha256'] == record['output_sha256'] and entry['route'] == source['route']
        assert entry['input_sha256'] == source['input_sha256']
        assert entry['source_receipt'] == source.get('receipt_path')
        assert entry['source_receipt_sha256'] == source.get('receipt_sha256')
        if source['route'] != 'panel45':
            source_checks[source['path']] = source['input_sha256']
            source_checks[str(pc.ROOT/source['receipt_path'])] = source['receipt_sha256']
        routes[entry['route']] += 1
        for name in [f'{i:04d}.json', f'{i:04d}_drop.csv.gz', f'{i:04d}_infer.csv.gz']:
            path, main_path = base/'chunks'/name, pc.RESULTS/'full/chunks'/name
            before = stable_stat(path)
            assert before[:2] == stable_stat(main_path)[:2], 'Main cache is not the same verified inode'
            digest = pc.atlas.sha(path)
            item = dict(name=name, size_bytes=before[2], sha256=digest,
                        shared_inode=before[1], main_path=str(main_path.relative_to(pc.ROOT)))
            if name.endswith('.csv.gz'):
                policy = name.split('_')[1].split('.')[0]
                assert digest == entry['output_sha256'][policy]
                statuses, rows = Counter(), 0
                with gzip.open(path, 'rt', encoding='utf-8', newline='') as stream:
                    reader = csv.reader(stream)
                    header = next(reader)
                    assert header == pc.atlas.KEEP and not any(c.startswith('q_') for c in header)
                    exposure, phenotype, status = [header.index(c) for c in ('exposure_id', 'phenotype', 'status')]
                    for row in reader:
                        assert len(row) == len(header) and rows < len(assays)
                        assert row[exposure] == assays[rows] and row[phenotype] == entry['phenotype']
                        assert row[status]
                        statuses[row[status]] += 1
                        rows += 1
                assert rows == len(assays) and dict(statuses) == entry['policies'][policy]
                policies[policy] += rows
                item.update(rows=rows, explicit_statuses_verified=True, gzip_crc_eof_verified=True,
                            no_partial_q_columns=True)
            assert stable_stat(path) == before and stable_stat(main_path) == before
            assert pc.atlas.sha(path) == digest
            files.append(item)
        if number % 100 == 0:
            print(f'Verified cache {number}/{len(records)}', flush=True)
    for path, expected in source_checks.items():
        assert pc.atlas.sha(Path(path)) == expected, f'Source changed: {path}'
    pc.context()
    after_queue = pc.import_guard()
    assert (base/'RUN.json').read_bytes() == run_bytes and (base/'CHUNKS.json').read_bytes() == chunk_bytes
    output.mkdir(parents=True)
    (output/'SNAPSHOT_RUN.json').write_bytes(run_bytes)
    (output/'SNAPSHOT_CHUNKS.json').write_bytes(chunk_bytes)
    pc.atlas.save_json(output/'CACHE_FILES.json', files)
    report = dict(complete=True, verified_at_utc=datetime.now(timezone.utc).isoformat(),
                  verifier_sha256=pc.atlas.sha(Path(__file__)), snapshot_traits=len(records),
                  route_counts=dict(routes), actual_rows_read_per_policy=dict(policies),
                  cache_files_checked=len(files), all_files_share_main_cache_inodes=True,
                  complete_gzip_CRC_EOF_and_assay_order_verified=True,
                  source_files_and_receipts_checked=len(source_checks),
                  main_queue_before=before_queue, main_queue_after=after_queue,
                  precompute_pid_absent=True, precompute_exit_code_observed=False,
                  main_statistics_started=False, formal_full_scan_complete=False,
                  final_corrections_computed=False, candidates_selected=False,
                  numerical_estimates_recomputed=False,
                  scope='Completed unadjusted caches from the same primary atlas; byte identity, '
                        'canonical input identities, all payload rows and main-cache links verified. '
                        'Final main reader will check each identity again before reuse.')
    pc.atlas.save_json(output/'VERIFICATION.json', report)
    print(json.dumps(dict(output=str(output), **report), indent=2), flush=True)


if __name__ == '__main__':
    main()
