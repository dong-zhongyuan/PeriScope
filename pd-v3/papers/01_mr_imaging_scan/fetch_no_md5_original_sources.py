"""Fetch the published original author TSV for every planned source that has no
published formatted MD5.

Decision of 2026-11-18: route (a). These 22 EBI sources publish md5sum.txt only for
the harmonised files, so the frozen acquisition script can never pass its early
compressed-MD5 comparison. Instead of editing that frozen script, the original
author TSV (published original MD5 available for all 1372 sources) is placed in the
directory the frozen script already prefers, so its original branch runs unchanged
and extract() still verifies the full reconstructed author MD5.

Operational only: no PLAN.json, recipe, or scientific parameter is touched.
"""
from pathlib import Path
import fcntl, hashlib, json, os, shutil, sys, time
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
SCRIPTS = ROOT / "papers/01_mr_imaging_scan"
BASE = ROOT / "data/all_target_followup_20261118/regional_all"
DOWNLOAD = ROOT / "data/new_scan_regional_20261117/download"
STATE_DIR = ROOT / "data/all_target_followup_20261118/no_formatted_md5_originals_20261118"
# 4 connections, matching one source worker of the protein queue. The protein
# controller is set to 5 workers while this runs, so the total against EBI stays
# at the 24 concurrent connections measured as stable; above roughly 32 that
# endpoint starts refusing. A single connection was starved by the queue and
# averaged 0.062 MiB/s, losing its TLS handshake repeatedly. Each attempt resumes
# the committed blocks, so retries are cheap and worth repeating many times.
CONNECTIONS = 4
MAX_ATTEMPTS = 200
BACKOFF = 120
BACKOFF_MAX = 900

sys.path.insert(0, str(SCRIPTS))
import s09_acquire_qsm as downloader


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def digest(path, kind="sha256"):
    h = hashlib.md5() if kind == "md5" else hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dump(path, payload):
    tmp = Path(path).with_name(Path(path).name + ".pending")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    tmp.replace(path)


def log_line(text):
    line = "%s %s" % (utc(), text)
    with (STATE_DIR / "fetch.log").open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def targets():
    plan = json.loads((BASE / "PLAN.json").read_text())
    picked = []
    for position, accession in enumerate(plan["source_order"]):
        meta = plan["sources"][accession]["metadata"]
        if meta.get("formatted_url") and not meta.get("formatted_published_md5"):
            picked.append(
                dict(
                    source_order_position=position,
                    accession=accession,
                    original_name=meta["original_name"],
                    original_url=meta["original_url"],
                    original_published_md5=meta["original_published_md5"],
                    formatted_checksum_status=meta.get("formatted_checksum_status"),
                )
            )
    return plan, picked


def fetch_one(record, staging):
    final = DOWNLOAD / record["original_name"]
    receipt_path = final.with_name(final.name + ".download.json")
    expected = record["original_published_md5"]
    if final.exists():
        md5 = digest(final, "md5")
        if md5 == expected:
            return dict(status="already_present_verified", md5=md5, sha256=digest(final))
        return dict(status="present_but_md5_mismatch", md5=md5, expected_md5=expected)

    url = record["original_url"].replace("http://ftp.ebi.ac.uk/", "https://ftp.ebi.ac.uk/")
    folder = staging / record["accession"]
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / record["original_name"]
    started = utc()
    downloader.URL = url
    downloader.BLOCK = 4 << 20
    identity = downloader.http_identity()
    part, state = downloader.download_blocks(
        target, identity, CONNECTIONS, folder / "DOWNLOAD_PROGRESS.json"
    )
    md5 = digest(part, "md5")
    if md5 != expected:
        return dict(status="downloaded_md5_mismatch", md5=md5, expected_md5=expected,
                    http=identity, staged_path=str(part.relative_to(ROOT)))
    receipt = dict(
        http=identity,
        download_started_at_utc=started,
        completed_at_utc=utc(),
        size_bytes=part.stat().st_size,
        sha256=digest(part),
        md5=md5,
        published_md5_verified=True,
        completed_blocks=state["completed_blocks"],
        block_size=downloader.BLOCK,
        connections=CONNECTIONS,
        supplement_for="all_target_followup_20261118 sources with no published formatted MD5",
        note="Original author TSV only. Regional extraction and the full reconstructed "
             "author MD5 check are done later by the frozen acquisition script.",
    )
    shutil.move(str(part), str(final))
    dump(receipt_path, receipt)
    progress = folder / "DOWNLOAD_PROGRESS.json"
    if progress.exists():
        shutil.move(str(progress), str(STATE_DIR / ("DOWNLOAD_PROGRESS_%s.json" % record["accession"])))
    stale = target.with_name(target.name + ".part.json")
    if stale.exists():
        shutil.move(str(stale), str(STATE_DIR / ("PART_STATE_%s.json" % record["accession"])))
    return dict(status="downloaded_verified", md5=md5, sha256=receipt["sha256"],
                size_bytes=receipt["size_bytes"], receipt=str(receipt_path.relative_to(ROOT)))


def archive_md5_failure(accession, cycle_name):
    folder = BASE / "sources" / accession
    failed = folder / "FAILED.json"
    if not failed.is_file():
        return None
    data = json.loads(failed.read_text())
    if "OriginalIdentityMismatch" not in str(data.get("error", "")):
        return "left_in_place_not_an_md5_failure"
    target = folder / cycle_name
    target.mkdir(exist_ok=True)
    moved = {}
    for name in ("FAILED.json", "RUN.json", "worker.log", "DOWNLOAD_PROGRESS.json"):
        source = folder / name
        if source.is_file():
            destination = target / name
            shutil.move(str(source), str(destination))
            moved[name] = digest(destination)
    dump(target / "ARCHIVE.json", dict(
        archived_at_utc=utc(),
        reason="published formatted MD5 absent upstream; original author TSV supplied "
               "in data/new_scan_regional_20261117/download so the frozen original "
               "branch runs unchanged",
        previous_error=str(data.get("error")),
        files=moved,
    ))
    return str(target.relative_to(ROOT))


def main():
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    lock = (STATE_DIR / ".fetch.lock").open("a+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print("another fetcher holds the lock", flush=True)
        return 1
    staging = STATE_DIR / "staging"
    staging.mkdir(exist_ok=True)
    plan, records = targets()
    cycle_name = "attempt_no_published_formatted_md5_20261118"
    state = dict(
        pid=os.getpid(),
        started_at_utc=utc(),
        decision="route_a_download_published_original_author_TSV_no_code_change",
        plan_sha256=digest(BASE / "PLAN.json"),
        download_dir=str(DOWNLOAD.relative_to(ROOT)),
        connections=CONNECTIONS,
        total=len(records),
        results={},
    )
    dump(STATE_DIR / "FETCH.json", state)
    log_line("fetching %d originals with %d connections each" % (len(records), CONNECTIONS))
    for record in records:
        accession = record["accession"]
        entry = dict(record)
        for attempt in range(1, MAX_ATTEMPTS + 1):
            try:
                outcome = fetch_one(record, staging)
                entry.update(outcome, attempts=attempt)
                break
            except Exception as exc:
                entry.update(status="error", error=repr(exc), attempts=attempt)
                log_line("%s attempt %d failed: %r" % (accession, attempt, exc))
                if attempt < MAX_ATTEMPTS:
                    time.sleep(min(BACKOFF_MAX, BACKOFF * attempt))
        if entry.get("status") in ("downloaded_verified", "already_present_verified"):
            entry["failure_receipt_archived_to"] = archive_md5_failure(accession, cycle_name)
        state["results"][accession] = entry
        state["updated_at_utc"] = utc()
        dump(STATE_DIR / "FETCH.json", state)
        log_line("%s %s" % (accession, entry.get("status")))
    ok = [a for a, e in state["results"].items()
          if e.get("status") in ("downloaded_verified", "already_present_verified")]
    state.update(
        finished_at_utc=utc(),
        verified=len(ok),
        complete=len(ok) == len(records),
        note="Sources whose FAILED.json was archived are retried by the queue on its "
             "next controller pass and then use the local original author TSV.",
    )
    dump(STATE_DIR / "FETCH.json", state)
    log_line("done: %d/%d verified" % (len(ok), len(records)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
