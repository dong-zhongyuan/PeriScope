"""Bounded, identity-checked HTTP ZIP access, including stored nested ZIPs."""
import hashlib
import io
import json
from pathlib import Path
import time
import zipfile

import requests


class RangeSource:
    def __init__(self, url, size, cache, published_identity=None):
        self.url, self.size = url, size
        self.cache = Path(cache)
        self.cache.mkdir(parents=True, exist_ok=True)
        self.session = requests.Session()
        self.identity = None
        # Some public archives omit HTTP validators. The caller must then pin
        # and recheck the published version/file id, size and checksum. This is
        # metadata identity, NOT verification of the whole archive checksum.
        self.published_identity = published_identity
        self.last_request = 0.0

    def get(self, start, length):
        if not 0 <= start <= start+length <= self.size:
            raise ValueError('Range outside published source size')
        if length == 0:
            return b''
        if length > 16*1024*1024:
            raise ValueError('Unbounded ZIP read refused')
        end = start+length-1
        key = hashlib.sha256(json.dumps([self.url,start,end,self.published_identity],sort_keys=True).encode()).hexdigest()
        payload, receipt = self.cache/(key+'.bin'), self.cache/(key+'.json')
        if payload.exists() and receipt.exists():
            data, record = payload.read_bytes(), json.loads(receipt.read_text())
            if record['url'] != self.url or record['start'] != start or record['end'] != end or record.get('published_identity') != self.published_identity or len(data) != length or hashlib.sha256(data).hexdigest() != record['sha256']:
                raise ValueError('ZIP range cache checksum mismatch')
            self.check_identity(record['identity'])
            return data
        for attempt in range(4):
            try:
                time.sleep(max(0, 0.6-(time.monotonic()-self.last_request)))
                self.last_request = time.monotonic()
                with self.session.get(self.url, headers={'Accept-Encoding':'identity', 'Range':f'bytes={start}-{end}'}, timeout=(20,60), stream=True) as r:
                    r.raise_for_status()
                    if r.status_code != 206 or r.headers.get('Content-Range') != f'bytes {start}-{end}/{self.size}':
                        raise ValueError('Server did not return the exact ZIP range')
                    identity = {k:r.headers.get(k) for k in ['ETag','Last-Modified']}
                    if not any(identity.values()) and not self.published_identity:
                        raise ValueError('No upstream identity validator or pinned published identity')
                    self.check_identity(identity)
                    data = r.raw.read(length+1, decode_content=False)
                    if len(data) != length:
                        raise ValueError('Range payload length mismatch')
                record = dict(url=self.url,start=start,end=end,size=length,identity=identity,published_identity=self.published_identity,sha256=hashlib.sha256(data).hexdigest())
                payload.write_bytes(data)
                receipt.write_text(json.dumps(record,indent=2)+'\n')
                return data
            except requests.RequestException as exc:
                if attempt == 3:
                    raise
                response = getattr(exc, 'response', None)
                retry = response.headers.get('Retry-After', '') if response is not None else ''
                time.sleep(min(60, int(retry)) if retry.isdigit() else 1+attempt)

    def check_identity(self, identity):
        if self.identity is None:
            self.identity = identity
        if self.identity != identity:
            raise ValueError('ZIP source identity changed between reads')


class RangeReader(io.RawIOBase):
    def __init__(self, source, base=0, size=None):
        super().__init__()
        self.source, self.base = source, base
        self.size = source.size if size is None else size
        self.pos = 0
        if not 0 <= base <= base+self.size <= source.size:
            raise ValueError('Invalid nested ZIP extent')

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos

    def seek(self, offset, whence=0):
        if whence not in [0,1,2]: raise ValueError('Invalid whence')
        value = offset + (0 if whence == 0 else self.pos if whence == 1 else self.size)
        if not 0 <= value <= self.size: raise ValueError('ZIP seek outside extent')
        self.pos = value
        return self.pos

    def read(self, n=-1):
        length = min(self.size-self.pos, n if n >= 0 else self.size-self.pos)
        data = self.source.get(self.base+self.pos, length)
        self.pos += len(data)
        return data


def nested_reader(archive, parent_reader, info):
    """Find a stored child using zipfile's validated local-header offset."""
    if info.compress_type != zipfile.ZIP_STORED or info.compress_size != info.file_size:
        raise ValueError('Random nested access requires an uncompressed child ZIP')
    with archive.open(info) as member:
        # zipfile validates the local-header filename and sets the exact data
        # start on the supplied file-like reader before reading member content.
        start = parent_reader.tell()
        if member.read(4) != b'PK\x03\x04':
            raise ValueError('Nested member does not begin with a ZIP local header')
    return RangeReader(parent_reader.source, parent_reader.base+start, info.file_size)
