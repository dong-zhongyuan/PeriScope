"""Repair the two source controllers after a bad regex edit, then set 7/13 correctly.

tune_transport_concurrency.py rebuilt each match from its capture groups and so
dropped the literal ") as pool:" tail, leaving both controllers with a SyntaxError
and both source queues exiting with code 1. This restores the last known-good
copies archived by that same call, applies the intended worker counts by exact
literal replacement, and refuses to keep any file that does not compile.
"""
from pathlib import Path
import hashlib, json, py_compile, shutil, sys, tempfile
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
SCRIPTS = ROOT / "papers/01_mr_imaging_scan"
STATE = ROOT / "data/all_target_followup_20261118/transport_capacity_20261118"

JOBS = {
    "run_all_target_protein_regions.py": dict(
        restore_from="run_all_target_protein_regions.py.workers6_20261118T010358Z",
        edits=[
            ("source_workers=6,connections_per_source=4,", "source_workers=7,connections_per_source=4,"),
            ("ThreadPoolExecutor(max_workers=6) as pool:", "ThreadPoolExecutor(max_workers=7) as pool:"),
            ('while len(active)<6 and submitted<len(plan["source_order"])',
             'while len(active)<7 and submitted<len(plan["source_order"])'),
        ],
        expect="7 sources x 4 = 28 connections to ftp.ebi.ac.uk",
    ),
    "run_all_target_imaging_regions.py": dict(
        restore_from="run_all_target_imaging_regions.py.workers10_20261118T010358Z",
        edits=[
            ("source_workers=10,connections_per_source=4,", "source_workers=13,connections_per_source=4,"),
            ("ThreadPoolExecutor(max_workers=10) as pool:", "ThreadPoolExecutor(max_workers=13) as pool:"),
            ("while len(active)<10 and submitted<len(order)",
             "while len(active)<13 and submitted<len(order)"),
        ],
        expect="13 IDPs x 4 = 52 connections to open.oxcin.ox.ac.uk",
    ),
}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


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


def compiles(path):
    with tempfile.TemporaryDirectory() as tmp:
        try:
            py_compile.compile(str(path), cfile=str(Path(tmp) / "out.pyc"), doraise=True)
            return True, None
        except Exception as exc:
            return False, repr(exc)[:300]


def main():
    blockers = {name: alive(name) for name in list(JOBS) +
                ["acquire_all_target_protein_regions.py", "acquire_all_target_imaging_regions.py"]
                if alive(name)}
    if blockers:
        print("still running, stop these first: %s" % json.dumps(blockers))
        return 1
    report = dict(repaired_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  cause="tune_transport_concurrency.py rebuilt matches from capture groups and "
                        "dropped the ') as pool:' tail, so both controllers raised SyntaxError",
                  files={})
    for name, job in JOBS.items():
        path = SCRIPTS / name
        archive = STATE / job["restore_from"]
        if not archive.is_file():
            print("missing archive for %s: %s" % (name, archive))
            return 2
        broken = STATE / (name + ".syntaxerror_20261118")
        shutil.copy2(str(path), str(broken))
        shutil.copy2(str(archive), str(path))
        ok, err = compiles(path)
        if not ok:
            print("restored copy of %s does not compile: %s" % (name, err))
            return 3
        text = path.read_text()
        for old, new in job["edits"]:
            if text.count(old) != 1:
                print("anchor missing or not unique in %s: %r" % (name, old))
                shutil.copy2(str(archive), str(path))
                return 4
            text = text.replace(old, new)
        path.write_text(text)
        ok, err = compiles(path)
        if not ok:
            shutil.copy2(str(archive), str(path))
            print("edited %s does not compile, reverted to archive: %s" % (name, err))
            return 5
        report["files"][name] = dict(
            restored_from=str(archive.relative_to(ROOT)),
            restored_sha256=sha(archive),
            broken_copy_kept=str(broken.relative_to(ROOT)),
            sha256_after=sha(path),
            compiles=True,
            target=job["expect"],
        )
    (STATE / "CONCURRENCY_REPAIR.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: {"sha": v["sha256_after"][:12], "target": v["target"]}
                      for k, v in report["files"].items()}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
