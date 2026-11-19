"""Whole-file audit of a chromosome-stratified random sample of range-route units.

Downloads the complete published file for each sampled unit, runs the frozen full-file extraction
(with its published-MD5 gate for proteins and full numeric/row-count validation for BIG40), and
compares every output hash with the range-route receipt. Writes only under fast_range_20261118/audit.
COMPLETE.json files are never changed.
"""
import gzip, json, multiprocessing as mp, sys, time, traceback
from collections import defaultdict
from pathlib import Path
HERE = Path('/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan'); sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE/'fast_range'))
import numpy as np, pandas as pd
import acquire_all_target_protein_regions as A, native_region_io as native, s09_acquire_qsm as D, fast_range_common as C
ROOT = C.ROOT; B = A.BASE; IMG = B/'imaging'; OUT = C.FAST/'audit'; OUT.mkdir(parents=True, exist_ok=True)


def fetch(url, target, progress):
    if target.exists():
        return target
    D.URL = url.replace('http://', 'https://'); D.BLOCK = 4 << 20; ident = D.http_identity()
    part, _ = D.download_blocks(target, ident, 4, progress); part.replace(target); target.with_name(target.name+'.part.json').unlink(missing_ok=True)
    return target


def protein(acc):
    plan = json.loads(A.PLAN.read_text()); s = plan['sources'][acc]; m = s['metadata']; rc = json.loads((B/'sources'/acc/'COMPLETE.json').read_text())
    d = OUT/('protein_'+acc); d.mkdir(exist_ok=True); fmt = bool(m.get('formatted_url')); name = m['formatted_name'] if fmt else m['original_name']; url = m['formatted_url'] if fmt else m['original_url']
    p = fetch(url, d/name, d/'DOWNLOAD_PROGRESS.json'); fm = m.get('formatted_published_md5')
    fmt_ok = (A.sha(p, 'md5') == fm) if (fmt and fm) else None
    wanted, reverse, ident = A.extract(p, fmt, s, plan); base = 'data/all_target_followup_20261118/regional_all/sources/'+acc+'/'; res = {}
    for rid, rows in wanted.items():
        q = d/(rid+'.csv.gz'); A.write_table(q, rows, A.FIELDS); res[rid] = A.sha(q) == rc['outputs_sha256'][base+rid+'/protein_regional_statistics.csv.gz']
    q = d/'imaging_peak_outcomes.csv.gz'; A.write_table(q, reverse, A.FIELDS); res['peak_outcomes'] = A.sha(q) == rc['outputs_sha256'][base+'imaging_peak_outcomes.csv.gz']
    return dict(kind='protein', unit=acc, outputs_match=all(res.values()), outputs=res, published_full_md5_verified=ident['original_published_md5_verified'], formatted_md5_match=fmt_ok,
                source_rows=ident['source_rows'], catalog_rows=rc['source_rows'], rows_match=ident['source_rows'] == rc['source_rows'], peak_rows=len(reverse), etag_same=None)


def imaging(idp):
    ip = json.loads((IMG/'PLAN.json').read_text()); rec = ip['phenotypes']['BIG40:release2:discovery:'+idp]; ref = json.loads((IMG/'variant_reference/MANIFEST.json').read_text()); rc = json.loads((IMG/'BIG40'/idp/'COMPLETE.json').read_text())
    d = OUT/('imaging_'+idp); d.mkdir(exist_ok=True); p = fetch(rec['native_url'], d/(idp+'.txt.gz'), d/'DOWNLOAD_PROGRESS.json'); windows = {w: ip['windows'][w] for w in rec['windows']}
    with gzip.open(p, 'rb') as gz:
        frames, val = native.extract(gz, rec, windows)
    assert val['chromosome_rows'] == ref['chromosome_rows']; res = {}
    for wid, fr in frames.items():
        rp = ROOT/ref['outputs'][wid]['output']; f = pd.read_csv(rp, dtype={'chromosome': str}, float_precision='round_trip', keep_default_na=False); r = native.attach_frequency(fr, f); q = d/(wid+'.csv.gz')
        r.to_csv(q, index=False, compression={'method': 'gzip', 'mtime': 0}); res[wid] = native.sha(q) == rc['windows'][wid]['sha256']
    return dict(kind='imaging', unit=idp, outputs_match=all(res.values()), outputs=res, gzip_crc_eof_verified=True, full_numeric_validation=True, source_rows=val['source_rows'], catalog_rows=rc['source_rows'], rows_match=val['source_rows'] == rc['source_rows'],
                chromosomes=sorted({str(w['chromosome']) for w in windows.values()}), unknown_symbols=sum(v.get('unknown_symbols') or 0 for v in rc['transport']['windows'].values()))


def one(a):
    kind, u = a; t = time.time()
    try:
        r = (protein if kind == 'protein' else imaging)(u)
    except BaseException as e:
        r = dict(kind=kind, unit=u, outputs_match=False, error=repr(e)[:300], traceback=traceback.format_exc()[-600:])
    r['seconds'] = round(time.time()-t, 1); (OUT/((kind)+'_'+u)/'FULL_VERIFICATION.json').parent.mkdir(exist_ok=True); (OUT/(kind+'_'+u)/'FULL_VERIFICATION.json').write_text(json.dumps(r, indent=2))
    return r


def sample(kind, seed=20261118):
    rng = np.random.default_rng(seed); by = defaultdict(list)
    if kind == 'protein':
        plan = json.loads(A.PLAN.read_text())
        for f in (B/'sources').glob('*/COMPLETE.json'):
            if json.loads(f.read_text()).get('acquisition_route'):
                a = f.parent.name; by[A.normchr(plan['regions'][plan['sources'][a]['regions'][0]]['chromosome'])].append(a)
    else:
        ip = json.loads((IMG/'PLAN.json').read_text())
        for f in (IMG/'BIG40').glob('*/COMPLETE.json'):
            if json.loads(f.read_text()).get('acquisition_route'):
                i = f.parent.name; by[str(ip['windows'][ip['phenotypes']['BIG40:release2:discovery:'+i]['windows'][0]]['chromosome'])].append(i)
    out = [sorted(v)[int(rng.integers(len(v)))] for k, v in sorted(by.items())]
    return out


if __name__ == '__main__':
    kind = sys.argv[1]; workers = int(sys.argv[2]); units = sample(kind); print(kind, len(units), units, flush=True)
    with mp.Pool(workers, maxtasksperchild=1) as p, (OUT/('results_'+kind+'.jsonl')).open('a') as w:
        for r in p.imap_unordered(one, [(kind, u) for u in units]):
            w.write(json.dumps(r)+'\n'); w.flush(); print(json.dumps({k: v for k, v in r.items() if k not in ('outputs', 'traceback')}), flush=True)
    print('AUDIT_DONE', kind, flush=True)
