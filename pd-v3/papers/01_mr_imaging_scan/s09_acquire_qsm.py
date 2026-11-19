"""Acquire the Oxford QSM/T2* discovery archive without unpacking it.

At most eight 64-MiB HTTP ranges write directly to one sparse .part file.
Only complete blocks enter the atomic resume receipt. Their hashes are verified
on restart. The final file appears only after size/hash, tar and gzip validation.
"""
import argparse
import collections
import concurrent.futures
import csv
import fcntl
import gzip
import hashlib
import html
from http.client import HTTPException, IncompleteRead
import importlib.util
import json
import os
from pathlib import Path
import random
import re
import shutil
import sys
import tarfile
import threading
import time
import urllib.error

sys.dont_write_bytecode = True
_spec = importlib.util.spec_from_file_location("big40_acquire", Path(__file__).with_name("s07_acquire_big40.py"))
_shared = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_shared)
ROOT, utc, sha256, atomic_json, request = _shared.ROOT, _shared.utc, _shared.sha256, _shared.atomic_json, _shared.request
URL = "https://www.fmrib.ox.ac.uk/ukbiobank/gwas_resources/QSMandT2s_IDPs.tar.gz"
INDEX_URL = "https://www.fmrib.ox.ac.uk/ukbiobank/gwas_resources/index.html"
FILENAME = "QSMandT2s_IDPs.tar.gz"
BLOCK = 64 << 20
RESERVE = 3 << 30
BLOCK_ATTEMPTS = 24


def index_mapping(path):
    body = path.read_text(encoding="utf-8")
    matches = re.findall(r"((?:QSM|T2\*)_[A-Za-z0-9_]+)\s*\((\d{4}\.txt)\):\s*([^<\r\n]+)", body)
    result = {name: {"member_basename": name, "short_name": short,
                    "description": html.unescape(desc).strip(),
                    "family": "QSM" if short.startswith("QSM") else "T2star",
                    "n_y": None, "mapping_source": INDEX_URL}
              for short, name, desc in matches}
    if len(matches) != 28 or len(result) != 28:
        raise ValueError(f"Expected 28 unambiguous published mappings; found {len(matches)}/{len(result)}")
    return result


def http_identity(attempts=8):
    """Resolve size and validator, retrying transport refusals.

    Measured 20261118 against ftp.ebi.ac.uk: a burst of connection attempts is refused
    outright (ConnectionRefused in 0.2 s) while 16 steady connections are served with
    16/16 handshakes and TLS at 0.18 s. Without this retry one refused HEAD failed the
    whole source, and four in a row stopped the queue. HTTP status errors below 500 are
    source errors, not transport, and are raised at once.
    """
    for attempt in range(attempts):
        try:
            return _http_identity_once()
        except urllib.error.HTTPError as exc:
            if exc.code < 500 or attempt == attempts-1: raise
        except (OSError, HTTPException) as exc:
            if attempt == attempts-1: raise
        delay = min(90.0, 5.0 * 2 ** attempt) + random.uniform(0, 5)
        print(json.dumps({"http_identity_retry": attempt+1, "sleep_s": round(delay, 1)}), flush=True)
        time.sleep(delay)


def _http_identity_once():
    with request(URL, method="HEAD") as r:
        identity = {"url": URL, "resolved_url": r.url,
                    "http_size_bytes": int(r.headers["Content-Length"]),
                    "etag": r.headers.get("ETag"), "last_modified": r.headers.get("Last-Modified"),
                    "content_type": r.headers.get("Content-Type")}
    if not identity["etag"] and not identity["last_modified"]:
        raise RuntimeError("A stable HTTP validator is required for range downloads")
    return identity


def hash_region(fd, start, length):
    h = hashlib.sha256()
    consumed = 0
    while consumed < length:
        data = os.pread(fd, min(4 << 20, length-consumed), start+consumed)
        if not data:
            raise EOFError("Unexpected EOF while validating a completed block")
        h.update(data)
        consumed += len(data)
    return h.hexdigest()


def pwrite_all(fd, data, offset):
    view = memoryview(data)
    while view:
        written = os.pwrite(fd, view, offset)
        if not written:
            raise OSError("Zero-length pwrite")
        offset += written
        view = view[written:]


def download_blocks(target, identity, workers, progress_path):
    part = target.with_name(target.name + ".part")
    resume_path = target.with_name(target.name + ".part.json")
    total = identity["http_size_bytes"]
    if part.exists():
        if not resume_path.exists():
            raise RuntimeError("Existing staged archive has no resume receipt")
        state = json.loads(resume_path.read_text())
        if state["http"] != identity or state["block_size"] != BLOCK or part.stat().st_size != total:
            raise RuntimeError("Staged archive size, block size, or HTTP identity changed")
    else:
        if resume_path.exists():
            raise RuntimeError("Resume receipt exists without its staged archive")
        state = {"http": identity, "block_size": BLOCK, "started_at_utc": utc(), "completed_blocks": {}}
        atomic_json(resume_path, state)
    fd = os.open(part, os.O_RDWR | os.O_CREAT, 0o644)
    abort, guard = threading.Event(), threading.Lock()
    inflight = {}
    try:
        if os.fstat(fd).st_size == 0:
            os.ftruncate(fd, total)  # Sparse allocation; never create a second archive.
        for key, item in list(state["completed_blocks"].items()):
            start = int(key) * BLOCK
            expected = min(BLOCK, total-start)
            if item["size_bytes"] != expected or hash_region(fd, start, expected) != item["sha256"]:
                print(f"Block {key} failed cache verification; scheduling that block again", flush=True)
                del state["completed_blocks"][key]
        atomic_json(resume_path, state)
        initial = sum(v["size_bytes"] for v in state["completed_blocks"].values())
        remaining = total - initial
        if shutil.disk_usage(target.parent).free < remaining + RESERVE:
            raise RuntimeError("Insufficient free disk for remaining archive plus 3 GiB reserve")
        blocks = (total + BLOCK - 1) // BLOCK
        pending_ids = [i for i in range(blocks) if str(i) not in state["completed_blocks"]]
        started = time.monotonic()

        def worker(number):
            start, end = number * BLOCK, min(total, (number+1)*BLOCK)-1
            length = end-start+1
            for attempt in range(BLOCK_ATTEMPTS):
                if abort.is_set(): raise RuntimeError("Download stopped by another worker")
                with guard: inflight[number] = 0
                try:
                    headers = {"Range": f"bytes={start}-{end}",
                               "If-Range": identity["etag"] or identity["last_modified"]}
                    with request(URL, headers=headers) as r:
                        if r.status != 206 or r.headers.get("Content-Range") != f"bytes {start}-{end}/{total}":
                            raise RuntimeError(f"Invalid range response for block {number}")
                        if int(r.headers.get("Content-Length", -1)) != length:
                            raise RuntimeError(f"Invalid content length for block {number}")
                        for key, header in (("etag", "ETag"), ("last_modified", "Last-Modified")):
                            if identity[key] and r.headers.get(header) != identity[key]:
                                raise RuntimeError("Upstream identity changed during download")
                        h, received = hashlib.sha256(), 0
                        while True:
                            if abort.is_set(): raise RuntimeError("Download stopped by another worker")
                            chunk = r.read(4 << 20)
                            if not chunk: break
                            if received + len(chunk) > length:
                                raise RuntimeError("Range response exceeded requested size")
                            if shutil.disk_usage(target.parent).free < RESERVE + len(chunk):
                                abort.set()
                                raise RuntimeError("Free space reached 3 GiB reserve; paused safely")
                            pwrite_all(fd, chunk, start+received)
                            h.update(chunk)
                            received += len(chunk)
                            with guard: inflight[number] = received
                    if received != length:
                        raise EOFError(f"Short block {number}: {received}/{length}")
                    os.fsync(fd)
                    return number, {"size_bytes": received, "sha256": h.hexdigest(), "completed_at_utc": utc()}
                except (OSError, EOFError, HTTPException) as exc:
                    if attempt == BLOCK_ATTEMPTS-1: raise
                    print(f"Block {number} attempt {attempt+1}: {type(exc).__name__}; retrying", flush=True)
                    abort.wait(min(2 ** attempt, 120) + random.uniform(0, 3))

        pool = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
        futures = {pool.submit(worker, number) for number in pending_ids}
        try:
            last_print = 0
            while futures:
                done, futures = concurrent.futures.wait(futures, timeout=5, return_when=concurrent.futures.FIRST_COMPLETED)
                for future in done:
                    number, receipt = future.result()
                    state["completed_blocks"][str(number)] = receipt
                    with guard: inflight.pop(number, None)
                    state["updated_at_utc"] = utc()
                    atomic_json(resume_path, state)
                committed = sum(x["size_bytes"] for x in state["completed_blocks"].values())
                with guard: downloaded = committed + sum(inflight.values())
                elapsed = time.monotonic() - started
                status = {"phase": "downloading", "updated_at_utc": utc(), "total_bytes": total,
                          "committed_bytes": committed, "including_inflight_bytes": downloaded,
                          "completed_blocks": len(state["completed_blocks"]), "total_blocks": blocks,
                          "elapsed_seconds": round(elapsed, 1),
                          "mean_MiB_per_second": round((downloaded-initial)/max(elapsed, .001)/(1<<20), 3),
                          "free_disk_bytes": shutil.disk_usage(target.parent).free}
                atomic_json(progress_path, status)
                if elapsed - last_print >= 15 or not futures:
                    print(json.dumps(status), flush=True)
                    last_print = elapsed
        except BaseException:
            abort.set()
            for future in futures: future.cancel()
            raise
        finally:
            pool.shutdown(wait=True, cancel_futures=True)
        if len(state["completed_blocks"]) != blocks:
            raise RuntimeError("Not all blocks are committed")
        return part, state
    finally:
        os.close(fd)


def validate_archive(path, mapping, progress_path):
    members, names, total_bytes = [], set(), 0
    with gzip.open(path, "rb") as gz:
        with tarfile.open(fileobj=gz, mode="r|") as archive:
            for member in archive:
                if member.name in names:
                    raise ValueError(f"Duplicate tar member name: {member.name}")
                names.add(member.name)
                item = {"name": member.name, "size_bytes": member.size,
                        "type": "file" if member.isfile() else "other"}
                basename = Path(member.name).name
                if member.isfile():
                    handle = archive.extractfile(member)
                    first = handle.readline(1 << 20)
                    count, byte_count, tail = first.count(b"\n"), len(first), first[-1:]
                    header = first.decode("ascii", errors="strict").split()
                    if basename.endswith(".txt") and header != _shared.SOURCE_HEADER:
                        raise ValueError(f"Unexpected summary header in {member.name}: {header}")
                    for chunk in iter(lambda: handle.read(8 << 20), b""):
                        count += chunk.count(b"\n")
                        byte_count += len(chunk)
                        tail = chunk[-1:]
                    if byte_count != member.size:
                        raise EOFError(f"Member size mismatch for {member.name}")
                    if tail and tail != b"\n": count += 1
                    item.update({"header": header, "data_rows": count-1, "member_bytes_verified": True})
                    total_bytes += byte_count
                item["mapping_status"] = "official_index" if basename in mapping else "UNMAPPED_DO_NOT_INFER"
                item["official_mapping"] = mapping.get(basename)
                item["n_y"] = None
                members.append(item)
                status = {"phase": "validating_tar", "updated_at_utc": utc(), "members_completed": len(members),
                          "uncompressed_member_bytes": total_bytes, "last_member": member.name,
                          "free_disk_bytes": shutil.disk_usage(path.parent).free}
                atomic_json(progress_path, status)
                print(json.dumps(status), flush=True)
        # Tar ends at its zero records. Read beyond them to force gzip trailer/CRC validation.
        trailer_bytes = 0
        for chunk in iter(lambda: gz.read(8 << 20), b""): trailer_bytes += len(chunk)
    found = {Path(x["name"]).name for x in members if x["type"] == "file"}
    missing = sorted(set(mapping) - found)
    if missing: raise ValueError(f"Officially listed members missing from archive: {missing}")
    return {"gzip_crc_and_eof_verified": True, "tar_member_headers_and_sizes_verified": True,
            "regular_file_count": sum(x["type"] == "file" for x in members),
            "uncompressed_member_bytes": total_bytes, "remaining_gzip_payload_bytes": trailer_bytes,
            "official_mapped_files": len(mapping), "unmapped_files": sorted(found-set(mapping)), "members": members,
            "limitations": ["No per-variant statistical validation or MR performed.",
                            "No individual-IDP sample size assigned; n_y remains null."]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, choices=range(1, 9), default=4)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/oxford_qsm")
    args = parser.parse_args()
    dest = args.data_dir.resolve()
    download, extracted = dest / "download", dest / "extracted"
    download.mkdir(parents=True, exist_ok=True)
    extracted.mkdir(parents=True, exist_ok=True)
    progress_path = dest / "progress.json"
    with (dest / ".acquire.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        index = download / "index.html"
        index_receipt = _shared.acquire(INDEX_URL, index, lambda p: {"official_mappings": len(index_mapping(p))})
        mapping = index_mapping(index)
        rows = [mapping[k] for k in sorted(mapping)]
        _shared.write_csv(extracted / "official_28_idp_mapping.csv", rows, list(rows[0]))
        target = download / FILENAME
        final_receipt = target.with_name(target.name + ".download.json")
        identity = http_identity()
        if target.exists():
            if not final_receipt.exists(): raise RuntimeError("Final archive has no validation receipt")
            receipt = json.loads(final_receipt.read_text())
            if receipt["http"] != identity or target.stat().st_size != identity["http_size_bytes"] or sha256(target) != receipt["sha256"]:
                raise RuntimeError("Final archive or upstream identity changed")
            print("Verified complete cached archive", flush=True)
        else:
            part, state = download_blocks(target, identity, args.workers, progress_path)
            atomic_json(progress_path, {"phase": "hashing_archive", "updated_at_utc": utc(), "total_bytes": part.stat().st_size})
            digest = sha256(part)
            validation = validate_archive(part, mapping, progress_path)
            receipt = {"http": identity, "size_bytes": part.stat().st_size, "sha256": digest,
                       "started_at_utc": state["started_at_utc"], "completed_at_utc": utc(),
                       "block_size": BLOCK, "completed_blocks": state["completed_blocks"], "validation": validation}
            atomic_json(final_receipt, receipt)
            part.replace(target)
            target.with_name(target.name + ".part.json").unlink()
        atomic_json(extracted / "archive_validation.json", receipt["validation"])
        atomic_json(dest / "MANIFEST.json", {
            "written_at_utc": utc(), "path": str(target.relative_to(dest)), **receipt,
            "index_path": str(index.relative_to(dest)), "index_receipt": index_receipt,
            "official_mapping_path": "extracted/official_28_idp_mapping.csv",
            "official_mapping_sha256": sha256(extracted / "official_28_idp_mapping.csv"),
            "source": INDEX_URL, "doi": "10.1038/s41593-022-01074-w", "cohort": "discovery",
            "literature_discovery_n": 19720, "literature_n_source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9174052/",
            "n_y": None, "n_y_note": "Per-IDP n not independently verified; literature cohort n is not substituted.",
            "license": "Public summary download; explicit external data-file license UNVERIFIED",
            "script": str(Path(__file__).relative_to(ROOT)), "script_sha256": sha256(__file__),
            "shared_script_sha256": sha256(Path(__file__).with_name("s07_acquire_big40.py")),
            "command": sys.argv, "python": sys.version})
        atomic_json(progress_path, {"phase": "complete", "updated_at_utc": utc(),
                                   "size_bytes": receipt["size_bytes"], "sha256": receipt["sha256"]})
        print(json.dumps({"complete": True, "path": str(target), "sha256": receipt["sha256"]}), flush=True)


if __name__ == "__main__":
    main()
