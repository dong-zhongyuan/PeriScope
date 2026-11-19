"""Byte-range acquisition of the rows the all-target analysis actually uses.

Protein sources: the author TSV is uncompressed, so the needed rows are read directly by
byte range. BIG40 sources: gzip is entered at a deflate block boundary (fast_range/imgfetch.py).
Every AGES file shares one variant order and every BIG40 file shares another, so a catalog built
from one verified file locates rows in all of them, and every fetched row is checked against it.

What this does NOT do: it never sees the whole remote file, so it cannot check the published
whole-file MD5, count all source rows or verify the gzip CRC. Receipts say so explicitly.
"""
import hashlib, http.client, json, os, random, socket, ssl, sys, time
from pathlib import Path
from urllib.parse import urlsplit
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE/'fast_range'))
ROOT = HERE.parents[1]
FAST = ROOT/'data/all_target_followup_20261118/fast_range_20261118'
UA = 'pd-v3-BIG40-acquisition/1.0'
ROUTE_VERSION = 'fast_range_v1_20261118'


class IdentityChanged(RuntimeError):
    pass


def sha_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


class RangeHTTP:
    """One keep-alive connection, range reads pinned to the validator seen at HEAD."""

    def __init__(self, url, attempts=30, timeout=60):
        self.attempts, self.timeout = attempts, timeout
        self.conn = None; self.requests = 0; self.bytes = 0; self.retries = 0
        self._set(url); self.url = url
        for hop in range(4):
            st, h, _ = self._retry('HEAD', {})
            if st in (301, 302, 303, 307, 308) and h.get('Location'):
                self.close(); self._set(h['Location']); continue
            break
        if st != 200:
            raise IdentityChanged('HEAD status %s for %s' % (st, url))
        self.size = int(h['Content-Length']); self.etag = h.get('ETag'); self.last_modified = h.get('Last-Modified')
        if not (self.etag or self.last_modified):
            raise IdentityChanged('no HTTP validator for ' + url)
        if 'bytes' not in (h.get('Accept-Ranges') or ''):
            raise IdentityChanged('server does not offer byte ranges for ' + url)

    def _set(self, url):
        u = urlsplit(url); self.host = u.netloc; self.path = u.path + (('?' + u.query) if u.query else ''); self.resolved = url

    def close(self):
        if self.conn is not None:
            try:
                self.conn.close()
            except Exception:
                pass
        self.conn = None

    def _once(self, method, headers):
        if self.conn is None:
            self.conn = http.client.HTTPSConnection(self.host, timeout=self.timeout)
        self.conn.request(method, self.path, headers={'Accept-Encoding': 'identity', 'User-Agent': UA, **headers})
        r = self.conn.getresponse(); body = r.read(); st = r.status; h = r.headers
        if r.will_close:
            self.close()
        self.requests += 1
        return st, h, body

    def _retry(self, method, headers, check=None):
        last = None
        for n in range(self.attempts):
            try:
                st, h, body = self._once(method, headers)
                if st in (408, 425, 429) or st >= 500:
                    raise ConnectionError('HTTP %d from %s' % (st, self.host))
                if check is not None:
                    check(st, h, body)
                return st, h, body
            except IdentityChanged:
                self.close(); raise
            except (OSError, http.client.HTTPException, EOFError) as e:
                last = e; self.close(); self.retries += 1
                if n == self.attempts-1:
                    break
                time.sleep(min(2 ** n, 60) + random.uniform(0, 2))
        raise last

    def read(self, a, b):
        a = max(0, int(a)); b = min(self.size, int(b))
        if b <= a:
            return b''
        want = 'bytes %d-%d/%d' % (a, b-1, self.size)
        def check(st, h, body):
            if st == 200:
                raise IdentityChanged('range ignored or validator changed for ' + self.resolved)
            if st != 206:
                raise IdentityChanged('unexpected status %d for %s' % (st, self.resolved))
            if h.get('Content-Range') != want:
                raise IdentityChanged('Content-Range %r != %r' % (h.get('Content-Range'), want))
            for key, val in (('ETag', self.etag), ('Last-Modified', self.last_modified)):
                if val and h.get(key) != val:
                    raise IdentityChanged('validator changed during transfer: ' + key)
            if len(body) != b-a:
                raise EOFError('short range body %d of %d' % (len(body), b-a))
        _, _, body = self._retry('GET', {'Range': 'bytes=%d-%d' % (a, b-1), 'If-Range': self.etag or self.last_modified}, check)
        self.bytes += len(body)
        return body

    def identity(self):
        return dict(url=self.url, resolved_url=self.resolved, http_size_bytes=self.size, etag=self.etag, last_modified=self.last_modified)


class Counted:
    """Per-window byte budget view of a RangeHTTP."""

    def __init__(self, http):
        self.h = http; self.size = http.size; self.bytes = 0; self.requests = 0

    def read(self, a, b):
        d = self.h.read(a, b); self.bytes += len(d); self.requests += 1
        return d


# ---------------------------------------------------------------- protein
class ProteinCatalog:
    def __init__(self):
        m = json.loads((FAST/'protein_meta.json').read_text())
        self.meta = m; self.N = int(m['rows']); self.size = float(m['tsv_bytes'])
        self.start = np.load(FAST/'protein_start.npy', mmap_mode='r')
        self.hash = np.load(FAST/'protein_hash.npy', mmap_mode='r')
        self.key = np.load(FAST/'protein_key.npy', mmap_mode='r')

    def rows(self, c, a, b):
        return int(np.searchsorted(self.key, (int(c) << 32) | int(a), 'left')), int(np.searchsorted(self.key, (int(c) << 32) | int(b), 'right'))-1

    def off(self, r, size):
        return (float(self.start[r]) if r < self.N else self.size)*size/self.size


def cnum(normalized):
    return 23 if normalized == 'X' else int(normalized)


def protein_fetch(pfz, http, cat, r0, r1, E=48 << 10, step=1 << 18):
    a = max(0, int(cat.off(r0, http.size))-E); b = min(http.size, int(cat.off(r1+1, http.size))+E)
    D = http.read(a, b); want = np.asarray(cat.hash[r0:r1+1]); m = r1-r0+1
    for attempt in range(12):
        buf = np.frombuffer(D, np.uint8); nl = np.flatnonzero(buf == 10)
        first = int(nl[0])+1
        seg = buf[first:]; k = pfz.keyscan(seg, 0, False, len(seg)//20+16); h = k['hash']
        found = None
        if len(h):
            for i in np.flatnonzero(h == want[0]).tolist():
                if i+m <= len(h) and np.array_equal(h[i:i+m], want):
                    found = i; break
        if found is not None:
            st = k['start']; out = {}
            for j in range(m):
                p = int(st[found+j])
                if found+j+1 < k['rows']:
                    q = int(st[found+j+1])
                else:
                    tail = np.flatnonzero(seg[p:] == 10)
                    if not len(tail):
                        found = None; break
                    q = p+int(tail[0])+1
                out[r0+j] = seg[p:q].tobytes()
            if found is not None:
                return out
        back = (not len(h)) or ((int(k['chrom'][0]) << 32) | int(k['pos'][0])) > int(cat.key[r0])
        if back and a > 0:
            na = max(0, a-step); D = http.read(na, a)+D; a = na
        elif b < http.size:
            n = http.read(b, b+step); b += len(n); D += n
        else:
            raise RuntimeError('rows not found in source')
    raise RuntimeError('row search exhausted')


def protein_extract(A, pfz, http, cat, source, plan):
    """Same selection, validation and row shape as the frozen extract(), from fetched rows only."""
    from collections import defaultdict
    import csv, math
    regions = [plan['regions'][r] for r in source['regions']]; lines = {}; checked = 0
    for reg in regions:
        r0, r1 = cat.rows(cnum(A.normchr(reg['chromosome'])), reg['start_grch37'], reg['end_grch37'])
        if r1 >= r0:
            lines.update(protein_fetch(pfz, http, cat, r0, r1)); checked += r1-r0+1
    positions = sorted({(cnum(A.normchr(x['chromosome'])), int(x['position_grch37'])) for x in source['peak_rows']}); absent = 0
    for c, p in positions:
        r0, r1 = cat.rows(c, p, p)
        if r1 < r0:
            absent += 1; continue
        if not all(r in lines for r in range(r0, r1+1)):
            lines.update(protein_fetch(pfz, http, cat, r0, r1)); checked += r1-r0+1
    bychrom = defaultdict(list)
    for reg in regions:
        bychrom[A.normchr(reg['chromosome'])].append(reg)
    pk = {A.key(r['chromosome'], r['position_grch37'], r['effect_allele'], r['other_allele']) for r in source['peak_rows']}
    wanted = {r['region_id']: [] for r in regions}; reverse = []; meta = source['metadata']
    for r in sorted(lines):
        line = lines[r]; number = r+1; fields = line.rstrip(b'\n').split(b'\t')
        if len(fields) != 9 or b'\r' in line or not line.endswith(b'\n'):
            raise ValueError('Malformed author row')
        chrom = A.normchr(fields[2].decode('ascii')); pos = int(fields[3])
        sel = [x['region_id'] for x in bychrom.get(chrom, []) if x['start_grch37'] <= pos <= x['end_grch37']]
        ispeak = A.key(chrom, pos, fields[4].decode(), fields[5].decode()) in pk
        if not sel and not ispeak:
            continue
        dec = next(csv.reader([line.decode('utf-8')], delimiter='\t', strict=True))
        if len(dec) != 9:
            raise ValueError('Malformed quoted author row')
        row = dict(zip(A.COLUMNS, dec)); vals = [float(row[k]) for k in ['p_value', 'effect_allele_frequency', 'beta', 'standard_error']]
        if not all(math.isfinite(x) for x in vals) or not (0 <= vals[0] <= 1 and 0 <= vals[1] <= 1 and vals[3] > 0):
            raise ValueError('Nonfinite/invalid selected source statistics')
        row.update(source_row=number, gene=meta['gene'], assay=meta['ages_assay'], accession=source['accession'], build='GRCh37',
                   n_metadata=meta['sample_size'], n_scope='reported_assay_sample_size_not_per_variant_N')
        for rid in sel:
            wanted[rid].append(row)
        if ispeak:
            reverse.append(row)
    return wanted, reverse, dict(rows_checked_against_catalog=checked, peak_positions=len(positions), peak_positions_absent_from_catalog=absent)


# ---------------------------------------------------------------- imaging
def imaging_catalog(F):
    class ImgCat(F.Catalog):
        def __init__(self):
            m = json.loads((FAST/'imaging_meta.json').read_text()); self.meta = m
            self.T = np.memmap(FAST/'imaging_text.bin', dtype=np.uint8, mode='r'); self.N = int(m['rows']); self.size = float(m['gz_bytes'])
            self.start = np.load(FAST/'imaging_start.npy', mmap_mode='r'); self.key = np.load(FAST/'imaging_key.npy', mmap_mode='r')
            self.pend = np.load(FAST/'imaging_pend.npy', mmap_mode='r'); self.nls = np.load(FAST/'imaging_nls.npy'); self.cbyte = np.load(FAST/'imaging_cbyte.npy')
    return ImgCat()


def code_hashes(paths):
    return {str(Path(p).resolve().relative_to(ROOT)): sha_file(p) for p in paths}
