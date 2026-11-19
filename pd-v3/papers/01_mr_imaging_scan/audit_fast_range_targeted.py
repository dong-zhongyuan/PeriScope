"""Whole-file audit of every unaudited range-route unit that supports a region with coloc H4 >= 0.8.

Reuses the audit functions in audit_fast_range.py unchanged. Results go to a new directory.
COMPLETE.json files, receipts and earlier audit results are never changed.

Usage: audit_fast_range_targeted.py protein|imaging WORKERS
"""
import hashlib, json, multiprocessing as mp, sys, time
from pathlib import Path
HERE = Path('/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan'); sys.path.insert(0, str(HERE))
import pandas as pd
import audit_fast_range as AU

RES = Path('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/research_20261118/verify-partial_fetch-PF6-a')
OUT = AU.C.FAST/'audit_targeted_H4_20261119'; OUT.mkdir(parents=True, exist_ok=True); AU.OUT = OUT
MAX_PASSES = 3


def units():
    hits = pd.read_csv(RES/'coloc_hits_H4_ge_0.8_by_route.csv', dtype={'idp': str, 'accession': str})
    led = pd.read_csv(RES/'route_ledger_units.csv', dtype={'unit': str}, low_memory=False)
    rng = led[led.route.str.startswith('fast_range') & ~led.audited.astype(bool)]
    p = sorted(set(hits.loc[hits.protein_range.astype(bool) & ~hits.protein_audited.astype(bool), 'accession']) & set(rng.loc[rng.kind == 'protein', 'unit']))
    i = sorted(set(hits.loc[hits.imaging_range.astype(bool) & ~hits.imaging_audited.astype(bool) & (hits.family == 'BIG40'), 'idp']) & set(rng.loc[rng.kind == 'imaging', 'unit']))
    return {'protein': p, 'imaging': i}


def done(path):
    ok = set()
    if path.exists():
        for line in path.read_text().splitlines():
            r = json.loads(line)
            if r.get('outputs_match') and r.get('rows_match'):
                ok.add(r['unit'])
    return ok


if __name__ == '__main__':
    kind, workers = sys.argv[1], int(sys.argv[2]); all_units = units(); todo = all_units[kind]
    m = OUT/'MANIFEST.json'
    if not m.exists():
        body = dict(created_at_unix=int(time.time()), selection='unaudited range-route units supporting coloc H4>=0.8 regions', counts={k: len(v) for k, v in all_units.items()}, units=all_units,
                    source_tables=[str(RES/'coloc_hits_H4_ge_0.8_by_route.csv'), str(RES/'route_ledger_units.csv')], audit_code_sha256=hashlib.sha256((HERE/'audit_fast_range.py').read_bytes()).hexdigest())
        m.write_text(json.dumps(body, indent=2))
    res = OUT/('results_'+kind+'.jsonl'); print(kind, len(todo), flush=True)
    for attempt in range(1, MAX_PASSES+1):
        left = [u for u in todo if u not in done(res)]
        print('PASS', attempt, 'remaining', len(left), flush=True)
        if not left:
            break
        with mp.get_context('fork').Pool(workers, maxtasksperchild=1) as p, res.open('a') as w:
            for r in p.imap_unordered(AU.one, [(kind, u) for u in left]):
                r['pass'] = attempt; w.write(json.dumps(r)+'\n'); w.flush()
                print('RESULT', json.dumps({k: v for k, v in r.items() if k not in ('outputs', 'traceback')}), flush=True)
    left = [u for u in todo if u not in done(res)]
    print('TARGETED_AUDIT_DONE', kind, 'matched', len(todo)-len(left), 'unresolved', left, flush=True)
