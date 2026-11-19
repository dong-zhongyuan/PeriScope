"""Raise the acquisition concurrency to just under each endpoint's measured ceiling.

Measured on A6000 on 2026-11-18, with both queues already running their 4x4 recipe:

  per completed protein source (n=113): total 375.8 s median, of which download
  339.2 s and extraction 36.0 s, so 90% of the wall clock is transport;
  per completed BIG40 IDP (n=74): total 491.7 s, download 453.0 s, 342 MB.

  added parallel range downloads of distinct files, measured over 20 s each:
    ftp.ebi.ac.uk     N=4 0.39  N=8 1.09  N=16 1.99  N=32 2.28 MiB/s (15/32 refused)
    open.oxcin.ox.ac.uk  N=8 1.34  N=16 2.47  N=32 4.07 MiB/s (0 refused)

  Per connection the rate is flat at about 0.127 MiB/s on both hosts regardless of
  N, so aggregate throughput is essentially 0.127 MiB/s x connections. The limit is
  not endpoint bandwidth: EBI starts refusing above roughly 32 concurrent
  connections from this host, while Oxford accepted 48 with no refusal, and the lab
  link carried about 9 MiB/s in total during the Oxford test.

So the queues were using about half the usable capacity. Only the number of worker
slots changes here, in the two operational controllers. Connections per source stay
at 4, and the acquisition scripts, PLAN.json, every identity check and the
consecutive-failure stop are untouched. Extraction is 36-38 s of CPU per item on a
48-core, 1 TB host, so it is not the constraint at these worker counts.
"""
from pathlib import Path
import hashlib, json, shutil, sys
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
SCRIPTS = ROOT / "papers/01_mr_imaging_scan"
STATE = ROOT / "data/all_target_followup_20261118/transport_capacity_20261118"

MEASUREMENT = {
    "measured_at_utc": None,
    "protein_source_seconds": {"n": 113, "total_median": 375.8, "download_median": 339.2,
                               "extract_median": 36.0, "bytes_median_MB": 177.4},
    "imaging_idp_seconds": {"n": 74, "total_median": 491.7, "download_median": 453.0,
                            "bytes_median_MB": 342.0},
    "aggregate_MiB_per_second_by_connections": {
        "ftp.ebi.ac.uk": {"4": 0.39, "8": 1.09, "16": 1.99, "32": 2.28},
        "open.oxcin.ox.ac.uk": {"8": 1.34, "16": 2.47, "32": 4.07},
    },
    "refusals": {"ftp.ebi.ac.uk": "15 of 32 non-2xx while the queue held 16 more",
                 "open.oxcin.ox.ac.uk": "0 of 32 while the queue held 16 more"},
    "per_connection_MiB_per_second": 0.127,
    "note": "Tests ran alongside the live queues, so they measure spare capacity, "
            "not an idle link. Each test connection fetched a distinct file.",
}

EDITS = {
    "run_all_target_protein_regions.py": [
        ('status="running",complete=False,source_workers=4,connections_per_source=4,',
         'status="running",complete=False,source_workers=6,connections_per_source=4,'),
        ("  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:\n",
         "  with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:\n"),
        ('    while len(active)<4 and submitted<len(plan["source_order"]) and not stop_scheduling:\n',
         '    while len(active)<6 and submitted<len(plan["source_order"]) and not stop_scheduling:\n'),
    ],
    "run_all_target_imaging_regions.py": [
        ("total_native_sources=len(order),source_workers=4,connections_per_source=4,",
         "total_native_sources=len(order),source_workers=10,connections_per_source=4,"),
        ("  with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:\n",
         "  with concurrent.futures.ThreadPoolExecutor(max_workers=10) as pool:\n"),
        ("    while len(active)<4 and submitted<len(order) and not stop:\n",
         "    while len(active)<10 and submitted<len(order) and not stop:\n"),
    ],
}

TARGET_CONNECTIONS = {"run_all_target_protein_regions.py": "6 sources x 4 = 24 to ftp.ebi.ac.uk",
                      "run_all_target_imaging_regions.py": "10 IDPs x 4 = 40 to open.oxcin.ox.ac.uk"}


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dump(path, payload):
    tmp = Path(path).with_name(Path(path).name + ".pending")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    tmp.replace(path)


def alive(name):
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            raw = (entry / "cmdline").read_bytes()
        except Exception:
            continue
        if raw and name.encode() in raw:
            return int(entry.name)
    return None


def main():
    blockers = {}
    for name in list(EDITS) + ["acquire_all_target_protein_regions.py",
                               "acquire_all_target_imaging_regions.py"]:
        pid = alive(name)
        if pid:
            blockers[name] = pid
    if blockers:
        print("still running, stop these first: %s" % json.dumps(blockers))
        return 1

    STATE.mkdir(parents=True, exist_ok=True)
    measurement = dict(MEASUREMENT, measured_at_utc=utc())
    dump(STATE / "MEASUREMENT.json", measurement)

    report = dict(applied_at_utc=utc(), measurement=measurement, targets=TARGET_CONNECTIONS,
                  unchanged=["papers/01_mr_imaging_scan/acquire_all_target_protein_regions.py",
                             "papers/01_mr_imaging_scan/acquire_all_target_imaging_regions.py",
                             "connections_per_source=4 (--connections 4)",
                             "consecutive_failures>=4 operational stop",
                             "every published-MD5 / SHA / gzip-CRC identity check",
                             "PLAN.json and both imaging/QSM plans"],
                  files={})
    for name, edits in EDITS.items():
        path = SCRIPTS / name
        before = sha(path)
        text = path.read_text()
        for old, new in edits:
            if text.count(old) != 1:
                print("anchor missing or not unique in %s: %r" % (name, old[:90]))
                return 2
        copy = STATE / (name + ".pre_concurrency_raise")
        shutil.copy2(str(path), str(copy))
        for old, new in edits:
            text = text.replace(old, new)
        path.write_text(text)
        report["files"][name] = dict(sha256_before=before, sha256_after=sha(path),
                                     archived_copy=str(copy.relative_to(ROOT)),
                                     archived_copy_sha256=sha(copy),
                                     edits=[old.strip()[:90] for old, _ in edits])
    dump(STATE / "CONCURRENCY_CHANGE.json", report)
    print(json.dumps({k: {"before": v["sha256_before"][:12], "after": v["sha256_after"][:12]}
                      for k, v in report["files"].items()}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
