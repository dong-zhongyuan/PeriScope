"""Losslessly compact/restore the five already verified historical AGES files.

This changes local storage only. Source identities and scientific outputs stay
unchanged. A complete gzip round trip must match both the source SHA256 and the
published MD5 before its uncompressed local cache can be removed.
"""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'data/ages_public'
RECEIPT = BASE / 'CACHE_STORAGE.json'
BLOCK = 8 * 1024 * 1024


def utc():
    return datetime.now(timezone.utc).isoformat()


def save(value):
    temp = RECEIPT.with_suffix('.json.tmp')
    temp.write_text(json.dumps(value, indent=2) + '\n')
    temp.replace(RECEIPT)


def hashes(stream, destination=None):
    sha, md5, size = hashlib.sha256(), hashlib.md5(), 0
    while block := stream.read(BLOCK):
        sha.update(block)
        md5.update(block)
        size += len(block)
        if destination is not None:
            destination.write(block)
    return dict(bytes=size, sha256=sha.hexdigest(), md5=md5.hexdigest())


def expected(record):
    download = record['download']
    return dict(bytes=download['size_bytes'], sha256=download['sha256'],
                md5=download['validation']['md5'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--restore', action='store_true')
    args = parser.parse_args()
    source_manifest = BASE / 'MANIFEST.json'
    manifest = json.loads(source_manifest.read_text())
    assert manifest['complete'] and len(manifest['datasets']) == 5
    manifest_sha = hashlib.sha256(source_manifest.read_bytes()).hexdigest()
    state = json.loads(RECEIPT.read_text()) if RECEIPT.exists() else dict(
        source_manifest_sha256=manifest_sha, files={}, scientific_outputs_changed=False)
    assert state['source_manifest_sha256'] == manifest_sha
    state.update(status='restoring' if args.restore else 'compacting', updated_at_utc=utc(),
                 code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    save(state)
    for record in manifest['datasets']:
        relative = Path(record['source_file'])
        assert relative.parts[:3] == ('data', 'ages_public', 'download')
        assert len(relative.parts) == 4 and relative.suffix == '.tsv'
        source = ROOT / relative
        cache = source.with_suffix('.tsv.cache.gz')
        truth = expected(record)
        if args.restore:
            if source.exists():
                with source.open('rb') as stream:
                    assert hashes(stream) == truth
            else:
                temp = source.with_suffix('.tsv.restore.part')
                with gzip.open(cache, 'rb') as stream, temp.open('xb') as out:
                    assert hashes(stream, out) == truth
                    out.flush()
                    os.fsync(out.fileno())
                with temp.open('rb') as stream:
                    assert hashes(stream) == truth
                temp.replace(source)
            entry = state['files'].setdefault(str(relative), {})
            entry.update(status='restored_original', original=truth, restored_at_utc=utc())
        else:
            stat_before = source.stat() if source.exists() else None
            if not cache.exists():
                if stat_before is None:
                    raise FileNotFoundError(source)
                temp = cache.with_suffix('.gz.part')
                with source.open('rb') as stream, temp.open('xb') as raw:
                    with gzip.GzipFile(filename='', fileobj=raw, mode='wb', compresslevel=6, mtime=0) as out:
                        assert hashes(stream, out) == truth
                    raw.flush()
                    os.fsync(raw.fileno())
                # Read the actual compressed output all the way through its CRC.
                with gzip.open(temp, 'rb') as stream:
                    assert hashes(stream) == truth
                temp.replace(cache)
            else:
                with gzip.open(cache, 'rb') as stream:
                    assert hashes(stream) == truth
            with cache.open('rb') as stream:
                compressed = hashes(stream)
            entry = dict(status='verified_compressed_original_retained', original=truth,
                         cache=str(cache.relative_to(ROOT)), compressed=compressed,
                         gzip_crc_eof_verified=True, roundtrip_sha256_and_md5_verified=True,
                         verified_at_utc=utc())
            state['files'][str(relative)] = entry
            save(state)
            if source.exists():
                # Reject a source changed by another process during compaction.
                stat_after = source.stat()
                assert stat_before is not None
                assert (stat_after.st_ino, stat_after.st_size, stat_after.st_mtime_ns) == (
                    stat_before.st_ino, stat_before.st_size, stat_before.st_mtime_ns)
                source.unlink()
            entry.update(status='compressed_roundtrip_verified', original_cache_removed_at_utc=utc())
        state.update(updated_at_utc=utc())
        save(state)
        print(json.dumps(dict(source=str(relative), status=entry['status'])), flush=True)
    state.update(status='restored' if args.restore else 'complete', updated_at_utc=utc())
    save(state)


if __name__ == '__main__':
    main()
