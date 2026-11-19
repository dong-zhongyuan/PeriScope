"""Real-network check: rebuild already completed units through the range route into a private folder and compare hashes."""
import json, os, sys, tempfile, time
from pathlib import Path
HERE = Path('/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan'); sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE/'fast_range'))
import numpy as np, pfz, imgfetch as F
import acquire_all_target_protein_regions as A, native_region_io as native, fast_range_common as C, acquire_fast_imaging_regions as I
ROOT = C.ROOT; B = A.BASE; OUT = C.FAST/'pilot'; OUT.mkdir(parents=True, exist_ok=True)
plan = json.loads(A.PLAN.read_text()); ip = json.loads((B/'imaging/PLAN.json').read_text()); ref = json.loads((B/'imaging/variant_reference/MANIFEST.json').read_text())
n = int(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 20261118); rows = []
acc = sorted(a for a in plan['sources'] if (B/'sources'/a/'COMPLETE.json').exists()); pc = C.ProteinCatalog()
for a in rng.choice(acc, n, replace=False).tolist():
    rc = json.loads((B/'sources'/a/'COMPLETE.json').read_text()); s = plan['sources'][a]; t = time.time()
    if rc.get('input_format', '').startswith('byte_ranges'): continue
    h = C.RangeHTTP(s['metadata']['original_url'].replace('http://', 'https://'))
    wanted, reverse, st = C.protein_extract(A, pfz, h, pc, s, plan); d = OUT/('p_'+a); d.mkdir(exist_ok=True); ok = tot = 0
    base = 'data/all_target_followup_20261118/regional_all/sources/'+a+'/'
    for rid, r in wanted.items():
        p = d/(rid+'.csv.gz'); A.write_table(p, r, A.FIELDS); tot += 1; ok += A.sha(p) == rc['outputs_sha256'][base+rid+'/protein_regional_statistics.csv.gz']
    p = d/'peaks.csv.gz'; A.write_table(p, reverse, A.FIELDS); tot += 1; ok += A.sha(p) == rc['outputs_sha256'][base+'imaging_peak_outcomes.csv.gz']
    print(json.dumps(dict(kind='protein', unit=a, match=ok, total=tot, bytes=h.bytes, requests=h.requests, retries=h.retries, secs=round(time.time()-t, 1))), flush=True); h.close()
idp = sorted(p.parent.name for p in (B/'imaging/BIG40').glob('*/COMPLETE.json')); ic = C.imaging_catalog(F)
for i in rng.choice(idp, n, replace=False).tolist():
    rc = json.loads((B/'imaging/BIG40'/i/'COMPLETE.json').read_text()); rec = ip['phenotypes']['BIG40:release2:discovery:'+i]; t = time.time()
    if rc.get('acquisition_route'): continue
    h = C.RangeHTTP(rec['native_url']); d = OUT/('i_'+i); d.mkdir(exist_ok=True); ok = tot = 0; fails = []
    for wid in rec['windows']:
        try:
            fr, rp, st = I.window_frame(F, h, ic, rec, wid, ip['windows'][wid], ref); p = d/(wid+'.csv.gz'); fr.to_csv(p, index=False, compression={'method': 'gzip', 'mtime': 0}); tot += 1; ok += native.sha(p) == rc['windows'][wid]['sha256']
        except F.Failure as e:
            tot += 1; fails.append(str(e)[:100])
    print(json.dumps(dict(kind='imaging', unit=i, match=ok, total=tot, bytes=h.bytes, requests=h.requests, retries=h.retries, failures=fails, secs=round(time.time()-t, 1))), flush=True); h.close()
print('PILOT_DONE', flush=True)
