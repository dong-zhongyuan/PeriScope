"""All-target BIG40 worker: range route first, frozen full-file route if a window cannot be recovered."""
from pathlib import Path
import argparse, fcntl, io, json, os, sys, time, traceback
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE/'fast_range'))
import pandas as pd
import native_region_io as native
import fast_range_common as C
ROOT = C.ROOT; BASE = ROOT/'data/all_target_followup_20261118/regional_all/imaging'; ORIGINAL = HERE/'acquire_all_target_imaging_regions.py'
HEADER = b'chr rsid pos a1 a2 beta se pval(-log10)\n'; L = 1 << 20; E = 48 << 10
FREQ = 'Official full BIG40 catalog a2 frequency, not per-IDP per-SNP sample-specific frequency'


def chrom_number(w):
    return 23 if str(w['chromosome']) in ('X', '23') else int(w['chromosome'])


def window_frame(F, http, cat, rec, wid, w, ref):
    r0, r1 = cat.rows(chrom_number(w), w['start_grch37'], w['end_grch37'])
    if r1 < r0:
        raise F.Failure('window has no catalog rows')
    view = C.Counted(http); got, info = F.fetch_rows(view, cat, r0, r1, L, E)
    rr = dict(rec); rr['expected_source_rows'] = r1-r0+1
    frames, _ = native.extract(io.BytesIO(HEADER+got), rr, {wid: dict(w)}); x = frames[wid]; x['source_data_row'] = x['source_data_row']+r0
    rp = ROOT/ref['outputs'][wid]['output']; assert native.sha(rp) == ref['outputs'][wid]['sha256']
    freq = pd.read_csv(rp, dtype={'chromosome': str}, float_precision='round_trip', keep_default_na=False)
    return native.attach_frequency(x, freq), rp, dict(rows=r1-r0+1, bytes=view.bytes, requests=view.requests, reads=info['reads'], widen=info['widen'], unknown_symbols=info.get('unknown_symbols', 0), solved_chars=info.get('solved_chars'))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--idp', required=True); ap.add_argument('--connections', type=int, default=4); a = ap.parse_args()
    plan = json.loads((BASE/'PLAN.json').read_text()); phenotype = 'BIG40:release2:discovery:'+a.idp; rec = plan['phenotypes'][phenotype]; assert rec['family'] == 'BIG40'
    folder = BASE/'BIG40'/a.idp; folder.mkdir(parents=True, exist_ok=True)
    lock = (folder/'.lock').open('a+'); fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB); ph = native.sha(BASE/'PLAN.json')
    if (folder/'COMPLETE.json').exists():
        d = json.loads((folder/'COMPLETE.json').read_text()); assert d['complete'] and d['plan_sha256'] == ph
        for v in d['windows'].values():
            assert native.sha(ROOT/v['output']) == v['sha256']
        print(json.dumps({'cached_verified': a.idp})); return
    if (folder/'FAILED.json').exists():
        raise RuntimeError('Inspect/archive prior failure before explicit retry')
    def full(reason):
        print(json.dumps({'full_file_route': a.idp, 'reason': reason}), flush=True); lock.close()
        os.execv(sys.executable, [sys.executable, str(ORIGINAL), '--idp', a.idp, '--connections', '2'])
    if rec['existing_full_payloads'] or rec['prior_full_native_source_identities'] or (folder/(a.idp+'.txt.gz')).exists():
        full('verified local or prior payload exists')
    import pfz, imgfetch as F
    rm = BASE/'variant_reference/MANIFEST.json'; ref_hash = native.sha(rm); ref = json.loads(rm.read_text()); assert ref['complete'] and ref['plan_sha256'] == ph
    deps = [Path(__file__), Path(C.__file__), Path(native.__file__), Path(pfz.__file__), Path(pfz.__file__).with_name('libpfz.so'), Path(F.__file__)]
    started = time.time(); code = C.code_hashes(deps)
    state = dict(pid=os.getpid(), status='running', complete=False, plan_sha256=ph, code_sha256=code, phenotype=phenotype, started_at_unix=started, acquisition_route=C.ROUTE_VERSION)
    native.dump(folder/'RUN.json', state); http = None
    try:
        cat = C.imaging_catalog(F); http = C.RangeHTTP(rec['native_url']); windows = {wid: plan['windows'][wid] for wid in rec['windows']}; outputs = {}; stats = {}; pending = {}
        try:
            for wid, w in windows.items():
                pending[wid] = window_frame(F, http, cat, rec, wid, w, ref)
        except F.Failure as exc:
            http.close(); (folder/'RUN.json').unlink(missing_ok=True); full('range recovery failed: %s' % exc)
        for wid, (frame, rp, st) in pending.items():
            out = folder/(wid+'.csv.gz'); frame.to_csv(out, index=False, compression={'method': 'gzip', 'mtime': 0})
            outputs[wid] = dict(output=str(out.relative_to(ROOT)), sha256=native.sha(out), rows=len(frame), frequency_reference_sha256=native.sha(rp), region=windows[wid], frequency_definition=FREQ); stats[wid] = st
        assert C.code_hashes(deps) == code; assert native.sha(BASE/'PLAN.json') == ph and native.sha(rm) == ref_hash
        report = dict(**state, source_rows=cat.N, source_rows_basis='shared BIG40 variant catalog; not counted from this file', chromosome_rows=None,
                      all_records_numerically_validated=False, requested_window_records_numerically_validated=True, region_allele_or_P_filter=None, gzip_crc_eof_verified=False, whole_file_read=False,
                      windows=outputs, transport=dict(route='deflate_block_entry_by_HTTP_byte_range', http=http.identity(), bytes=http.bytes, http_requests=http.requests, http_retries=http.retries, windows=stats, catalog=cat.meta),
                      identity_basis='HTTP validator pinned across ranges; every row prefix matched the catalog; beta, SE and -log10 P agree within rounding; every variant joined the official regional catalog',
                      variant_reference_manifest_sha256=native.sha(rm), matched_prior_native_source_identity=False,
                      matched_existing_scientific_recipe='Full original GRCh37 region; all alleles retained; exact native catalog frequency join')
        report.update(complete=True, status='complete', elapsed_seconds=time.time()-started); native.dump(folder/'COMPLETE.json', report); native.dump(folder/'RUN.json', report)
        print(json.dumps({'complete': a.idp, 'windows': len(outputs), 'bytes': http.bytes, 'requests': http.requests, 'retries': http.retries, 'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    except BaseException as exc:
        if isinstance(exc, SystemExit):
            raise
        state.update(status='failed', error=repr(exc), traceback=traceback.format_exc()); native.dump(folder/'FAILED.json', state); native.dump(folder/'RUN.json', state); raise
    finally:
        if http is not None:
            http.close()


if __name__ == '__main__':
    main()
