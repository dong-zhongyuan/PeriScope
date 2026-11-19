"""Build the shared variant catalogs used by the byte-range route from one verified retained file per family."""
import json, sys, time
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE/'fast_range'))
import numpy as np, pfz, imgfetch as F
import fast_range_common as C
import acquire_all_target_protein_regions as A
ROOT = C.ROOT; FAST = C.FAST; FAST.mkdir(parents=True, exist_ok=True)
B = ROOT/'data/all_target_followup_20261118/regional_all'; plan = json.loads((B/'PLAN.json').read_text())
t = time.time()
acc = sorted(a for a in plan['sources'] if (B/'sources'/a/'COMPLETE.json').exists() and (B/'sources'/a/plan['sources'][a]['metadata']['formatted_name']).exists())
ref = acc[0]; gz = B/'sources'/ref/plan['sources'][ref]['metadata']['formatted_name']
arr, _ = pfz.tsvize(pfz.gunzip(gz)); k = pfz.keyscan(arr, 0, True, 9000000)
key = (k['chrom'].astype(np.uint64) << np.uint64(32)) | k['pos'].astype(np.uint64)
assert np.all(key[1:] >= key[:-1]), 'protein catalog not sorted'
np.save(FAST/'protein_start.npy', k['start']); np.save(FAST/'protein_hash.npy', k['hash']); np.save(FAST/'protein_key.npy', key)
(FAST/'protein_meta.json').write_text(json.dumps(dict(reference_accession=ref, reference_file=str(gz.relative_to(ROOT)), reference_sha256=C.sha_file(gz), rows=int(k['rows']), tsv_bytes=int(len(arr))), indent=2))
print('protein catalog', k['rows'], len(arr), round(time.time()-t), flush=True); del arr, k, key
files = sorted(p for p in (B/'imaging/BIG40').glob('*/*.txt.gz') if (p.parent/'COMPLETE.json').exists()); f = files[0]
txt = pfz.gunzip(f); k = pfz.keyscan(txt, 1, True, 18000000); idx = pfz.index_file(f); c = F.Catalog(txt, k, idx)
(FAST/'imaging_text.bin').write_bytes(txt)
np.save(FAST/'imaging_start.npy', np.asarray(c.start)); np.save(FAST/'imaging_key.npy', np.asarray(c.key)); np.save(FAST/'imaging_pend.npy', np.asarray(c.pend))
np.save(FAST/'imaging_nls.npy', np.asarray(c.nls)); np.save(FAST/'imaging_cbyte.npy', np.asarray(c.cbyte))
(FAST/'imaging_meta.json').write_text(json.dumps(dict(reference_idp=f.parent.name, reference_file=str(f.relative_to(ROOT)), reference_sha256=C.sha_file(f), rows=int(c.N), gz_bytes=int(c.size), text_bytes=len(txt)), indent=2))
print('imaging catalog', c.N, c.size, round(time.time()-t), flush=True)
