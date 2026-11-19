"""Append a progress sample every 15 minutes so the download has a time series.

Counts COMPLETE.json receipts on disk, which is monotonic, rather than the
controllers' in-memory lists, which restart from zero and re-verify. Also records
host interface throughput so a stall can be told apart from a slowdown.
"""
from pathlib import Path
import csv, json, os, sys, time
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
BASE = ROOT / "data/all_target_followup_20261118"
REG = BASE / "regional_all"
RESULTS = ROOT / "papers/01_mr_imaging_scan/results"
OUT = BASE / "progress_log_20261118"
INTERVAL = 900
DEADLINE_HOURS = 72

COUNTS = {
    "protein_sources": (REG / "sources", 1372),
    "big40_idps": (REG / "imaging/BIG40", 915),
    "qsm_peak_outcomes": (REG / "reverse_supplement/protein_outcomes", 22),
    "coloc_regions": (RESULTS / "all_target_coloc_20261118", 1375),
    "reverse_regions": (RESULTS / "all_target_reverse_20261118", 1375),
}


def rx():
    total = 0
    for line in Path("/proc/net/dev").read_text().splitlines()[2:]:
        name, rest = line.split(":", 1)
        name = name.strip()
        if name.startswith(("lo", "veth", "docker", "br-")):
            continue
        total += int(rest.split()[0])
    return total


def done(folder):
    if not Path(folder).exists():
        return 0
    return sum(1 for child in Path(folder).iterdir()
               if child.is_dir() and (child / "COMPLETE.json").is_file())


def originals():
    d = ROOT / "data/new_scan_regional_20261117/download"
    if not d.exists():
        return 0
    return sum(1 for p in d.iterdir() if p.name.endswith("_buildGRCh37.tsv"))


def alive(needle):
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            raw = (entry / "cmdline").read_bytes()
        except Exception:
            continue
        if raw and needle.encode() in raw and b"python" in raw:
            return True
    return False


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "progress.csv"
    fields = (["utc", "MiB_per_s", "supervisor_up", "fetcher_up", "originals_placed"]
              + list(COUNTS) + [k + "_pct" for k in COUNTS])
    if not path.exists():
        with path.open("w", newline="") as handle:
            csv.DictWriter(handle, fieldnames=fields).writeheader()
    (OUT / "LOGGER.json").write_text(json.dumps(
        dict(pid=os.getpid(), interval_seconds=INTERVAL,
             started_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
             note="COMPLETE.json counts on disk, monotonic across controller restarts"),
        indent=2) + "\n")
    deadline = time.time() + DEADLINE_HOURS * 3600
    before = rx()
    mark = time.time()
    while time.time() < deadline:
        time.sleep(INTERVAL)
        now = rx()
        elapsed = max(1e-9, time.time() - mark)
        row = {
            "utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "MiB_per_s": round((now - before) / 1048576 / elapsed, 2),
            "supervisor_up": int(alive("supervise_all_target_queues.py")),
            "fetcher_up": int(alive("fetch_no_md5_original_sources.py")),
            "originals_placed": originals(),
        }
        before, mark = now, time.time()
        for key, (folder, total) in COUNTS.items():
            n = done(folder)
            row[key] = n
            row[key + "_pct"] = round(100.0 * n / total, 1)
        with path.open("a", newline="") as handle:
            csv.DictWriter(handle, fieldnames=fields).writerow(row)
        print(json.dumps(row), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
