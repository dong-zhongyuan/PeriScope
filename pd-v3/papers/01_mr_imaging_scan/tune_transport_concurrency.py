"""Set the worker-slot count of the two operational source controllers.

Usage: tune_transport_concurrency.py <protein_workers> <imaging_workers>

Only the number of concurrent sources changes. connections_per_source stays 4, and
the frozen acquisition scripts, PLAN.json, every published-MD5 / SHA / gzip-CRC
identity check and the consecutive-failure stop are untouched. Each call archives
the pre-edit controller with its sha256 and writes a receipt, so a step back is just
another call with the previous numbers.
"""
from pathlib import Path
import hashlib, json, py_compile, re, shutil, sys, tempfile
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
SCRIPTS = ROOT / "papers/01_mr_imaging_scan"
STATE = ROOT / "data/all_target_followup_20261118/transport_capacity_20261118"

TARGETS = {
    "run_all_target_protein_regions.py": dict(
        endpoint="ftp.ebi.ac.uk",
        patterns=[
            (r'(status="running",complete=False,source_workers=)(\d+)', 1),
            (r"(with concurrent\.futures\.ThreadPoolExecutor\(max_workers=)(\d+)\) as pool:", 1),
            (r'(while len\(active\)<)(\d+)( and submitted<len\(plan\["source_order"\]\))', 1),
        ],
    ),
    "run_all_target_imaging_regions.py": dict(
        endpoint="open.oxcin.ox.ac.uk",
        patterns=[
            (r"(total_native_sources=len\(order\),source_workers=)(\d+)", 1),
            (r"(with concurrent\.futures\.ThreadPoolExecutor\(max_workers=)(\d+)\) as pool:", 1),
            (r"(while len\(active\)<)(\d+)( and submitted<len\(order\))", 1),
        ],
    ),
}


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def alive(name):
    """Pid of a python process running this script, or None.

    A plain substring test over /proc also matches the shell that invoked this
    tool, because that command line mentions the script names, so match the
    script as an actual argument of a python interpreter instead.
    """
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            parts = [p for p in (entry / "cmdline").read_bytes().split(b"\0") if p]
        except Exception:
            continue
        if not parts or b"python" not in Path(parts[0].decode("utf-8", "replace")).name.encode():
            continue
        if any(Path(p.decode("utf-8", "replace")).name == name for p in parts[1:]):
            return int(entry.name)
    return None


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    wanted = {"run_all_target_protein_regions.py": int(sys.argv[1]),
              "run_all_target_imaging_regions.py": int(sys.argv[2])}
    for value in wanted.values():
        if not 1 <= value <= 24:
            print("worker count out of the sane 1-24 range")
            return 2
    blockers = {}
    for name in list(TARGETS) + ["acquire_all_target_protein_regions.py",
                                 "acquire_all_target_imaging_regions.py"]:
        pid = alive(name)
        if pid:
            blockers[name] = pid
    if blockers:
        print("still running, stop these first: %s" % json.dumps(blockers))
        return 1

    STATE.mkdir(parents=True, exist_ok=True)
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    report = dict(applied_at_utc=utc(), requested=wanted, files={},
                  unchanged=["connections_per_source=4 (--connections 4)",
                             "both acquire_* scripts",
                             "consecutive_failures>=4 operational stop",
                             "PLAN.json and every identity check"])
    for name, spec in TARGETS.items():
        path = SCRIPTS / name
        text = path.read_text()
        found = {}
        for pattern, group in spec["patterns"]:
            matches = re.findall(pattern, text)
            if len(matches) != 1:
                print("pattern not unique in %s: %s" % (name, pattern))
                return 3
            found[pattern] = int(matches[0][group])
        if len(set(found.values())) != 1:
            print("inconsistent current worker counts in %s: %s" % (name, found))
            return 3
        current = next(iter(found.values()))
        target = wanted[name]
        entry = dict(endpoint=spec["endpoint"], workers_before=current, workers_after=target,
                     connections_before=current * 4, connections_after=target * 4,
                     sha256_before=sha(path))
        if current == target:
            entry["skipped"] = "already at the requested count"
            report["files"][name] = entry
            continue
        copy = STATE / ("%s.workers%d_%s" % (name, current, stamp))
        shutil.copy2(str(path), str(copy))
        for pattern, group in spec["patterns"]:
            # Replace only the number inside the full match; rebuilding the match
            # from its capture groups silently drops any uncaptured tail.
            text = re.sub(
                pattern,
                lambda m: m.group(0).replace(m.group(1) + m.group(2), m.group(1) + str(target), 1),
                text,
            )
        original = path.read_text()
        path.write_text(text)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                py_compile.compile(str(path), cfile=str(Path(tmp) / "out.pyc"), doraise=True)
        except Exception as exc:
            path.write_text(original)
            print("edited %s does not compile, reverted: %r" % (name, exc))
            return 4
        entry.update(sha256_after=sha(path), archived_copy=str(copy.relative_to(ROOT)),
                     archived_copy_sha256=sha(copy))
        report["files"][name] = entry
    (STATE / ("CONCURRENCY_TUNE_%s.json" % stamp)).write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: {"workers": "%s->%s" % (v["workers_before"], v["workers_after"]),
                          "connections": "%s->%s" % (v["connections_before"], v["connections_after"])}
                      for k, v in report["files"].items()}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
