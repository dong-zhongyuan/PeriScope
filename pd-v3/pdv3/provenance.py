"""Provenance, not gating.

`LOCK.json` here records hashes so that any number can be traced back to the exact code and
inputs that produced it. It is not a threshold and it does not decide anything. Each result
file gets a sibling `<name>_provenance.json`.
"""
import hashlib
import json
import subprocess
import sys
import platform
import importlib.metadata
from datetime import datetime, timezone
from pathlib import Path


def sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def manifest(paths, root):
    root = Path(root).resolve()
    out = {}
    for p in sorted(map(Path, paths)):
        p = Path(p).resolve()
        key = str(p.relative_to(root)) if root in p.parents or p.parent == root else str(p)
        out[key] = sha256(p)
    return out


def git_rev(root):
    try:
        return subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=10).stdout.strip() or None
    except Exception:
        return None


def record(result_path, code_files=(), input_files=(), root=None, extra=None):
    """Write <result_path>_provenance.json beside a result file. Returns the record."""
    result_path = Path(result_path)
    root = Path(root) if root else result_path.parent
    rec = {
        "written_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "argv": sys.argv,
        "git_rev": git_rev(Path(__file__).resolve().parent.parent),
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {name: importlib.metadata.version(name) for name in
                     ('numpy', 'scipy', 'pandas', 'statsmodels', 'matplotlib')},
        "code": manifest(code_files, root) if code_files else {},
        "inputs": manifest(input_files, root) if input_files else {},
        "result_sha256": sha256(result_path) if result_path.exists() else None,
        "extra": extra or {},
    }
    out = result_path.with_name(result_path.name + "_provenance.json")
    out.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    return rec


def write_lock(lock_path, files, root, note=""):
    """Write a hash manifest once, for traceability. Refuses to overwrite."""
    lock_path = Path(lock_path)
    if lock_path.exists():
        raise FileExistsError(f"{lock_path} exists; record changes as a dated amendment")
    rec = {"written_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
           "note": note, "files": manifest(files, root)}
    lock_path.write_text(json.dumps(rec, indent=1, ensure_ascii=False))
    return rec


def verify_lock(lock_path, root):
    rec = json.loads(Path(lock_path).read_text())
    bad = [rel for rel, h in rec["files"].items()
           if not (Path(root) / rel).exists() or sha256(Path(root) / rel) != h]
    if bad:
        raise RuntimeError(f"files changed or missing since the hash manifest was written: {bad}")
    return rec
