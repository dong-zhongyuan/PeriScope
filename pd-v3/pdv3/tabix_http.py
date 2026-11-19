"""Checked, cached HTTP range reads of unfiltered BGZF/Tabix point tables.

No P-value selection. Index and downloaded byte ranges retain SHA256 receipts.
Only generic, one-based point tables are accepted. See HTS Tabix specification.
"""
import email.parser
import email.policy
import gzip
import hashlib
import json
from pathlib import Path
import re
import struct
import threading
import time
import zlib

import requests

_sessions = threading.local()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def atomic(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_bytes(data)
    tmp.replace(path)


def save_json(path, value):
    atomic(path, (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode())


def bgzf_blocks(data, base=0):
    """Yield complete BGZF blocks with verified gzip CRC; ignore only a final partial block."""
    offset = 0
    while offset + 18 <= len(data):
        if data[offset:offset + 4] != b'\x1f\x8b\x08\x04':
            raise ValueError(f'Invalid BGZF magic at {base + offset}')
        xlen = struct.unpack_from('<H', data, offset + 10)[0]
        if offset + 12 + xlen > len(data):
            break
        extra, size = offset + 12, None
        while extra < offset + 12 + xlen:
            tag, length = data[extra:extra + 2], struct.unpack_from('<H', data, extra + 2)[0]
            if tag == b'BC' and length == 2:
                size = struct.unpack_from('<H', data, extra + 4)[0] + 1
            extra += 4 + length
        if size is None or not 26 <= size <= 65536:
            raise ValueError('Missing/invalid BGZF BC field')
        if offset + size > len(data):
            break
        yield base + offset, zlib.decompress(data[offset:offset + size], 31)
        offset += size


def parse_tbi(data):
    raw = gzip.decompress(data)
    if raw[:4] != b'TBI\x01':
        raise ValueError('Not a Tabix index')
    cursor = 4

    def read(fmt):
        nonlocal cursor
        values = struct.unpack_from(fmt, raw, cursor)
        cursor += struct.calcsize(fmt)
        return values

    nref, fmt, seq, beg, end, meta, skip, length = read('<8i')
    if fmt != 0 or seq != 1 or beg != 2 or end not in (0, 2):
        raise ValueError(f'Only generic one-based chromosome/position tables supported: {(fmt,seq,beg,end)}')
    names = raw[cursor:cursor + length].decode().rstrip('\0').split('\0')
    cursor += length
    if len(names) != nref or len(set(names)) != nref:
        raise ValueError('Invalid reference names')
    refs = {}
    for name in names:
        bins = {}
        for _ in range(read('<i')[0]):
            number, count = read('<Ii')
            bins[number] = [read('<QQ') for _ in range(count)]
        count = read('<i')[0]
        linear = read('<' + 'Q' * count) if count else ()
        refs[name] = (bins, linear)
    if len(raw) - cursor not in (0, 8):
        raise ValueError('Unexpected trailing TBI data')
    return refs


def merge_intervals(intervals, gap=0):
    merged = []
    for start, end in sorted(set(intervals)):
        if start >= end:
            continue
        if merged and start <= merged[-1][1] + gap:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def point_chunks(refs, chrom, pos):
    if pos <= 0:
        raise ValueError('Positions must be positive one-based integers')
    if str(chrom) not in refs:
        return []
    bins, linear = refs[str(chrom)]
    p = pos - 1
    ids = [0, 1 + (p >> 26), 9 + (p >> 23), 73 + (p >> 20),
           585 + (p >> 17), 4681 + (p >> 14)]
    minimum = linear[p >> 14] if p >> 14 < len(linear) else 0
    return merge_intervals((max(a, minimum), b) for i in ids
                           for a, b in bins.get(i, []) if b > minimum)


class HTTPStore:
    def __init__(self, root):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        if not hasattr(_sessions, 'session'):
            _sessions.session = requests.Session()
            _sessions.session.headers.update({'User-Agent': 'pd-v3-unfiltered-summary-reader/1.0',
                                              'Accept-Encoding': 'identity'})
        self.session = _sessions.session

    def get(self, url, ranges=None, expected_identity=None):
        key = digest(json.dumps([url, ranges], separators=(',', ':')).encode())
        path = self.root / (key + '.bin')
        receipt = self.root / (key + '.json')
        if path.exists() and receipt.exists():
            data, info = path.read_bytes(), json.loads(receipt.read_text())
            if digest(data) != info['sha256']:
                raise ValueError('Cached bytes differ from receipt')
            self.check_identity(info, expected_identity)
            return data, info
        headers = {}
        if ranges:
            headers['Range'] = 'bytes=' + ','.join(f'{a}-{b}' for a, b in ranges)
        for attempt in range(4):
            try:
                with self.session.get(url, headers=headers, timeout=(20, 90), stream=True) as response:
                    if response.status_code != (206 if ranges else 200):
                        raise RuntimeError(f'HTTP {response.status_code}: {url}')
                    info = {'url': url, 'resolved_url': response.url, 'ranges': ranges,
                            'status': response.status_code, 'etag': response.headers.get('ETag'),
                            'last_modified': response.headers.get('Last-Modified'),
                            'content_type': response.headers.get('Content-Type'),
                            'content_range': response.headers.get('Content-Range')}
                    self.check_identity(info, expected_identity)
                    limit = sum(b - a + 1 for a, b in ranges) + 1024 * len(ranges) if ranges else 16 << 20
                    parts, total = [], 0
                    for block in response.iter_content(1 << 18):
                        total += len(block)
                        if total > limit:
                            raise ValueError('Response exceeds requested range/metadata size cap')
                        parts.append(block)
                    data = b''.join(parts)
                info.update(sha256=digest(data), size_bytes=len(data))
                if ranges:
                    split_ranges(data, info, ranges)  # validate before committing cache
                atomic(path, data)
                save_json(receipt, info)
                return data, info
            except (requests.RequestException, RuntimeError) as error:
                if attempt == 3 or 'HTTP 404' in str(error):
                    raise
                time.sleep(1 + attempt)

    @staticmethod
    def check_identity(info, expected):
        if expected:
            for field in ('etag', 'last_modified'):
                if expected.get(field) and info.get(field) != expected[field]:
                    raise ValueError(f'Upstream identity changed: {field}')


def split_ranges(data, info, requested):
    ctype = info.get('content_type') or ''
    if ctype.lower().startswith('multipart/byteranges'):
        message = email.parser.BytesParser(policy=email.policy.default).parsebytes(
            ('Content-Type: ' + ctype + '\r\nMIME-Version: 1.0\r\n\r\n').encode() + data)
        parts = [(part.get('Content-Range'), part.get_payload(decode=True)) for part in message.iter_parts()]
    else:
        parts = [(info.get('content_range'), data)]
    result, total = {}, None
    for header, payload in parts:
        match = re.fullmatch(r'bytes (\d+)-(\d+)/(\d+)', header or '')
        if not match:
            raise ValueError(f'Invalid Content-Range: {header}')
        start, end, length = map(int, match.groups())
        if total is not None and total != length:
            raise ValueError('Inconsistent source lengths')
        total = length
        if (start, end) in result or len(payload) != end - start + 1:
            raise ValueError('Duplicate or incomplete byte range')
        result[start, end] = payload
    if set(result) != set(tuple(r) for r in requested):
        raise ValueError('Server did not return exactly the requested ranges')
    return result, total


class RemotePointTable:
    def __init__(self, url, cache):
        self.url, self.store = url, HTTPStore(cache)
        self.index_bytes, self.index_receipt = self.store.get(url + '.tbi')
        self.refs = parse_tbi(self.index_bytes)
        raw, self.identity = self.store.get(url, [(0, 65535)])
        _, self.total = split_ranges(raw, self.identity, [(0, 65535)])
        head = b''.join(b for _, b in bgzf_blocks(raw)).decode().splitlines()[0]
        self.header = head.lstrip('#').split('\t')
        if self.header[:2] != ['chromosome', 'base_pair_location']:
            raise ValueError(f'Unexpected table columns: {self.header}')

    def fetch_many(self, points, batch_size=64):
        """Return all exact-position rows, including other alleles and high P values."""
        wanted = {(str(chrom), int(pos)) for chrom, pos in points}
        chunks = merge_intervals(chunk for c, p in wanted for chunk in point_chunks(self.refs, c, p))
        # Read the last block's 18-byte BGZF header first. Appending 64 KiB
        # to every sparse chunk wastes most transferred bytes for these loci.
        ends = sorted({b >> 16 for _, b in chunks if b & 65535})
        sizes = {}
        for offset in range(0, len(ends), batch_size):
            group = [(p,p+17) for p in ends[offset:offset+batch_size]]
            raw, info = self.store.get(self.url, group, self.identity)
            headers, total = split_ranges(raw, info, group)
            if total != self.total:raise ValueError('Source length changed')
            for (p,_), header in headers.items():
                if header[:4]!=b'\x1f\x8b\x08\x04' or header[12:16]!=b'BC\x02\x00':
                    raise ValueError('Unexpected BGZF last-block header')
                sizes[p]=struct.unpack_from('<H',header,16)[0]+1
        ranges = merge_intervals(((a >> 16), (b >> 16) + (sizes[b >> 16] if b & 65535 else 0)) for a, b in chunks)
        inclusive = [(a, b - 1) for a, b in ranges]
        decoded = {}
        for offset in range(0, len(inclusive), batch_size):
            group = inclusive[offset:offset + batch_size]
            raw, info = self.store.get(self.url, group, self.identity)
            payloads, total = split_ranges(raw, info, group)
            if total != self.total:
                raise ValueError('Source length changed')
            for (start, end), payload in payloads.items():
                for pos, block in bgzf_blocks(payload, start):
                    if pos in decoded and decoded[pos] != block:
                        raise ValueError('Inconsistent overlapping BGZF blocks')
                    decoded[pos] = block
        rows = {}
        ordered = sorted(decoded)
        import bisect
        for start, end in chunks:
            first, last = start >> 16, end >> 16
            offsets = ordered[bisect.bisect_left(ordered, first):bisect.bisect_right(ordered, last)]
            if not offsets or offsets[0] != first or ((end & 65535) and offsets[-1] != last):
                raise ValueError('Missing a required BGZF block')
            parts = [decoded[o][(start & 65535) if o == first else 0:
                                 (end & 65535) if o == last else len(decoded[o])] for o in offsets]
            for line in b''.join(parts).decode().splitlines():
                values = line.split('\t')
                if len(values) != len(self.header):
                    raise ValueError('Truncated/invalid tabix record')
                if (values[0], int(values[1])) in wanted:
                    rows[tuple(values)] = dict(zip(self.header, values))
        return list(rows.values()), {'requested_positions': len(wanted), 'returned_rows': len(rows),
                                    'virtual_chunks': len(chunks), 'http_ranges': len(inclusive),
                                    'header_ranges': len(ends), 'header_bytes_requested': 18*len(ends),
                                    'compressed_bytes_requested': sum(b-a+1 for a,b in inclusive),
                                    'source_bytes': self.total, 'index_sha256': digest(self.index_bytes),
                                    'header': self.header, 'unfiltered_point_query': True,
                                    'source_identity': self.identity}
