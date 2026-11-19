"""All-target protein source worker: range route first, frozen full-file route where a full file is required."""
from pathlib import Path
import argparse, fcntl, json, os, sys, time, traceback
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE/'fast_range'))
import acquire_all_target_protein_regions as A
import fast_range_common as C
ROOT = A.ROOT; BASE = A.BASE; PLAN = A.PLAN; ORIGINAL = HERE/'acquire_all_target_protein_regions.py'


def needs_full_file(acc, meta, folder):
    supp = json.loads((BASE/'reverse_supplement/PLAN.json').read_text())['sources']
    if acc in supp and supp[acc]['peak_rows']:
        return 'QSM peak outcome step reads the complete verified source'
    if (ROOT/'data/new_scan_regional_20261117/download'/meta['original_name']).exists():
        return 'complete author file already local'
    if meta.get('formatted_name') and (folder/meta['formatted_name']).exists():
        return 'complete formatted file already local'
    return None


def write_outputs(source, wanted, reverse, folder):
    outputs = {}
    for rid, rows in wanted.items():
        out = folder/rid; out.mkdir(exist_ok=True); p = out/'protein_regional_statistics.csv.gz'
        A.write_table(p, rows, A.FIELDS); outputs[str(p.relative_to(ROOT))] = A.sha(p)
    p = folder/'imaging_peak_outcomes.csv.gz'; A.write_table(p, reverse, A.FIELDS); outputs[str(p.relative_to(ROOT))] = A.sha(p)
    p = folder/'all_source_imaging_peaks.csv.gz'
    import csv
    fields = list(source['peak_rows'][0]) if source['peak_rows'] else list(csv.DictReader((ROOT/'data/big40_peak_catalog/all_10889_source_peaks_with_frequencies.csv').open()).fieldnames)
    A.write_table(p, source['peak_rows'], fields); outputs[str(p.relative_to(ROOT))] = A.sha(p)
    return outputs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--accession', required=True); ap.add_argument('--connections', type=int, default=4); a = ap.parse_args()
    plan = json.loads(PLAN.read_text()); source = plan['sources'][a.accession]; meta = source['metadata']
    folder = BASE/'sources'/a.accession; folder.mkdir(parents=True, exist_ok=True)
    lock = (folder/'.lock').open('a+'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB); ph = A.sha(PLAN)
    if (folder/'COMPLETE.json').exists():
        d = A.verify_complete(folder, ph); print(json.dumps({'cached_verified': a.accession, 'source_rows': d['source_rows']})); return
    if (folder/'FAILED.json').exists():
        raise RuntimeError('Inspect and explicitly archive prior failure before retry')
    why = needs_full_file(a.accession, meta, folder)
    if why:
        print(json.dumps({'full_file_route': a.accession, 'reason': why}), flush=True); lock.close()
        os.execv(sys.executable, [sys.executable, str(ORIGINAL), '--accession', a.accession, '--connections', '2'])
    import pfz
    deps = [Path(__file__), Path(C.__file__), Path(pfz.__file__), Path(pfz.__file__).with_name('libpfz.so'), ORIGINAL]
    started = time.time(); code = C.code_hashes(deps)
    run = dict(pid=os.getpid(), status='running', accession=a.accession, plan_sha256=ph, code_sha256=code, started_at_unix=started, acquisition_route=C.ROUTE_VERSION)
    A.dump(folder/'RUN.json', run); http = None
    try:
        cat = C.ProteinCatalog(); http = C.RangeHTTP(meta['original_url'].replace('http://', 'https://'))
        wanted, reverse, stats = C.protein_extract(A, pfz, http, cat, source, plan)
        outputs = write_outputs(source, wanted, reverse, folder)
        assert C.code_hashes(deps) == code, 'Code changed during extraction'; assert A.sha(PLAN) == ph
        receipt = dict(**run, complete=True, source_rows=cat.N, source_rows_basis='shared AGES variant catalog; not counted from this file',
                       input_format='byte_ranges_of_original_author_TSV', input_path=None, input_sha256=None, http=http.identity(), bytes_fetched=http.bytes,
                       http_requests=http.requests, http_retries=http.retries, catalog=cat.meta, **stats,
                       original_published_md5_verified=False, source_eof_verified=False, whole_file_read=False,
                       identity_basis='HTTP validator pinned across ranges; every fetched row key matched the catalog from an MD5-verified file; TLS transport',
                       regions={rid: {'rows': len(rows), 'coverage': 'extracted' if rows else 'no_source_variants_in_region'} for rid, rows in wanted.items()},
                       peak_outcome_rows=len(reverse), source_peak_rows=len(source['peak_rows']), outputs_sha256=outputs, source_full_file_retained=False,
                       formatted_fallback=None, elapsed_seconds=time.time()-started, biological_inference_computed=False)
        receipt['status'] = 'complete'; A.dump(folder/'COMPLETE.json', receipt); A.dump(folder/'RUN.json', receipt)
        print(json.dumps({'complete': a.accession, 'bytes': http.bytes, 'requests': http.requests, 'retries': http.retries, 'region_rows': {k: len(v) for k, v in wanted.items()}, 'peak_outcomes': len(reverse), 'elapsed_seconds': receipt['elapsed_seconds']}))
    except BaseException as exc:
        run.update(status='failed', error=repr(exc), traceback=traceback.format_exc(), elapsed_seconds=time.time()-started)
        A.dump(folder/'FAILED.json', run); A.dump(folder/'RUN.json', run); raise
    finally:
        if http is not None:
            http.close()


if __name__ == '__main__':
    main()
