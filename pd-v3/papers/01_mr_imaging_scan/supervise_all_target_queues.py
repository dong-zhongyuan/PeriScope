"""Operational supervisor for the five frozen all-target queues.

Scope is transport/run-protection only. It never edits PLAN.json, the frozen
acquisition/coloc/reverse recipes, or any scientific parameter. It only does what
the project's own recovery convention already does by hand: confirm a controller
is really dead, re-probe the public endpoint, archive the stopped controller state
and the network-failure receipts with sha256, then restart the same runner.

Failures that are not transport failures (for example the 22 EBI sources with no
published formatted MD5, or the chrX coloc cast error) are left untouched as
preserved failures and are never retried or bypassed.
"""
from pathlib import Path
import fcntl, hashlib, json, os, re, shutil, subprocess, sys, time, urllib.request
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
SCRIPTS = ROOT / "papers/01_mr_imaging_scan"
BASE = ROOT / "data/all_target_followup_20261118"
REG = BASE / "regional_all"
RESULTS = ROOT / "papers/01_mr_imaging_scan/results"
PY = "/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python"
STATE_DIR = BASE / "supervisor_20261118"

EBI = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
OXFORD = "https://open.win.ox.ac.uk/ukbiobank/big40/"

NET_PATTERN = re.compile(
    r"handshake operation timed out|ConnectionRefusedError|TimeoutError|"
    r"IncompleteRead|RemoteDisconnected|ConnectionResetError|Connection reset|"
    r"BadStatusLine|temporarily unavailable|HTTPError\('?5\d\d",
    re.I,
)

# Every status the five controllers can write, split by whether a restart is the
# right response. Transport, upstream-availability and crash states are resumable;
# region-level scientific failures are not, and are left for a human.
RESUMABLE = {
    "running",                              # killed or crashed mid-run, pid already proven dead
    "controller_json_missing",               # state archived, nothing holds the queue lock
    "ready",
    "waiting_inputs",
    "needs_attention_source_errors",
    "needs_attention_upstream_unavailable",
    "needs_attention_upstream_not_live",
    "finished_with_source_errors",
    "finished_with_source_failures",
    "finished_with_unavailable_inputs",
    "finished_with_unavailable_sources",
}
HUMAN = {
    "needs_attention_region_errors",
    "finished_with_region_errors",
    "finished_with_region_failures",
    "failed",
}
ARCHIVE_NAMES = ("FAILED.json", "RUN.json", "worker.log", "DOWNLOAD_PROGRESS.json")

POLL_SECONDS = 60
PROBE_TRIES = 2
BACKOFF_START = 120
BACKOFF_MAX = 1800
NO_PROGRESS_LIMIT = 3
CONSUMER_UPSTREAM_DELTA = 15
# Circuit breaker. The endpoint probe can succeed while the actual block downloads
# still time out, so a plain probe-gated restart can thrash: restart, trip the
# four-consecutive-failure stop within minutes, restart again. If a queue stops
# this soon after being restarted, wait out the bad window instead.
QUICK_STOP_SECONDS = 600
COOLDOWN_START = 300
COOLDOWN_MAX = 1800

QUEUES = [
    dict(
        key="protein_sources",
        controller=REG / "CONTROLLER.json",
        log=REG / "controller.log",
        runner=SCRIPTS / "run_all_target_protein_regions.py",
        kind="source",
        failure_root=REG / "sources",
        probes=[EBI],
        upstream=[],
    ),
    dict(
        key="big40_imaging",
        controller=REG / "imaging/BIG40_CONTROLLER.json",
        log=REG / "imaging/big40_controller.log",
        runner=SCRIPTS / "run_all_target_imaging_regions.py",
        kind="source",
        failure_root=REG / "imaging/BIG40",
        probes=[OXFORD],
        upstream=[],
    ),
    dict(
        key="qsm_peak_outcomes",
        controller=REG / "reverse_supplement/protein_outcomes/CONTROLLER.json",
        log=REG / "reverse_supplement/protein_outcomes/controller.log",
        runner=SCRIPTS / "run_qsm_peak_outcomes.py",
        kind="consumer",
        failure_root=None,
        probes=[],
        upstream=["protein_sources"],
    ),
    dict(
        key="all_target_coloc",
        controller=RESULTS / "all_target_coloc_20261118/CONTROLLER.json",
        log=RESULTS / "all_target_coloc_20261118/controller.log",
        runner=SCRIPTS / "run_all_target_coloc.py",
        kind="consumer",
        failure_root=None,
        probes=[],
        upstream=["protein_sources", "big40_imaging"],
    ),
    dict(
        key="all_target_reverse_mr",
        controller=RESULTS / "all_target_reverse_20261118/CONTROLLER.json",
        log=RESULTS / "all_target_reverse_20261118/controller.log",
        runner=SCRIPTS / "run_all_target_reverse.py",
        kind="consumer",
        failure_root=None,
        probes=[],
        upstream=["protein_sources", "big40_imaging"],
    ),
]


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


def log_line(text):
    line = "%s %s" % (utc(), text)
    with (STATE_DIR / "supervisor.log").open("a") as fh:
        fh.write(line + "\n")
    print(line, flush=True)


def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return None


def alive(pid, runner_name):
    if not pid:
        return False
    try:
        raw = Path("/proc/%d/cmdline" % int(pid)).read_bytes()
    except Exception:
        return False
    if not raw:
        return False
    return runner_name.encode() in raw


def probe(url):
    """Is the endpoint answering right now?

    The TLS handshake timeouts seen on these endpoints are intermittent and clear
    on their own while the endpoint keeps answering other clients, so a single
    failed handshake is not evidence that the source is down. curl is tried as a
    second opinion before the queue is held back.
    """
    last = "no attempt"
    for method in ("HEAD", "GET"):
        for _ in range(PROBE_TRIES):
            try:
                request = urllib.request.Request(url, method=method)
                if method == "GET":
                    request.add_header("Range", "bytes=0-1023")
                with urllib.request.urlopen(request, timeout=15) as response:
                    status = response.status
                    if method == "GET":
                        response.read(1024)
                if status < 400:
                    return True, "%s=%s" % (method, status)
                last = "%s=%s" % (method, status)
            except Exception as exc:
                last = "%s:%s" % (method, repr(exc))
            time.sleep(3)
    for _ in range(PROBE_TRIES):
        try:
            result = subprocess.run(
                ["curl", "-sS", "-o", "/dev/null", "-m", "20", "-w", "%{http_code}", url],
                capture_output=True, text=True, timeout=30,
            )
            code = result.stdout.strip()
            if code.isdigit() and int(code) < 400:
                return True, "curl=%s" % code
            last = "curl=%s %s" % (code, result.stderr.strip()[:80])
        except Exception as exc:
            last = "curl:%s" % repr(exc)
        time.sleep(3)
    return False, last


def completed_count(controller):
    data = read_json(controller)
    if not data:
        return None
    value = data.get("completed")
    if isinstance(value, list):
        return len(value)
    if isinstance(value, int):
        return value
    return None


def network_failures(failure_root):
    hits = []
    if not failure_root or not Path(failure_root).exists():
        return hits
    for folder in sorted(Path(failure_root).iterdir()):
        failed = folder / "FAILED.json"
        if not failed.is_file():
            continue
        data = read_json(failed)
        error = "" if not data else str(data.get("error", ""))
        if NET_PATTERN.search(error):
            hits.append((folder, error))
    return hits


def missing_receipts(failure_root):
    """Scheduled members that carry neither a COMPLETE.json nor a FAILED.json.

    This is how a source becomes retryable again after its receipt was archived
    on purpose, for example the sources whose original author TSV was supplied
    because upstream publishes no formatted MD5.
    """
    pending = []
    if not failure_root or not Path(failure_root).exists():
        return pending
    for folder in sorted(Path(failure_root).iterdir()):
        if not folder.is_dir():
            continue
        if (folder / "COMPLETE.json").is_file() or (folder / "FAILED.json").is_file():
            continue
        pending.append(folder.name)
    return pending


def preserved_failures(failure_root):
    kept = {}
    if not failure_root or not Path(failure_root).exists():
        return kept
    for folder in sorted(Path(failure_root).iterdir()):
        failed = folder / "FAILED.json"
        if not failed.is_file():
            continue
        data = read_json(failed)
        error = "" if not data else str(data.get("error", ""))
        if not NET_PATTERN.search(error):
            kept[folder.name] = error
    return kept


def archive_controller(queue, cycle):
    controller = Path(queue["controller"])
    target = controller.parent / ("controller_supervisor_restart_%03d_20261118" % cycle)
    target.mkdir(parents=True, exist_ok=True)
    receipt = {"archived_at_utc": utc(), "files": {}}
    for source in (controller, Path(queue["log"])):
        if source.is_file():
            destination = target / source.name
            shutil.move(str(source), str(destination))
            receipt["files"][source.name] = sha(destination)
    runner = Path(queue["runner"])
    runner_copy = target / runner.name
    shutil.copy2(str(runner), str(runner_copy))
    receipt["files"][runner.name] = sha(runner_copy)
    receipt["runner_sha256_live"] = sha(runner)
    dump(target / "ARCHIVE.json", receipt)
    return str(target.relative_to(ROOT)), receipt


def archive_network_failures(queue, cycle, hits):
    moved = {}
    for folder, error in hits:
        target = folder / ("attempt_supervisor_restart_%03d_20261118" % cycle)
        target.mkdir(parents=True, exist_ok=True)
        entry = {"error": error, "files": {}}
        for name in ARCHIVE_NAMES:
            source = folder / name
            if source.is_file():
                destination = target / name
                shutil.move(str(source), str(destination))
                entry["files"][name] = sha(destination)
        entry["partial_payloads_kept"] = sorted(
            p.name for p in folder.iterdir() if ".part" in p.name
        )
        dump(target / "ARCHIVE.json", {"archived_at_utc": utc(), **entry})
        moved[folder.name] = entry
    return moved


def launch(queue):
    log = Path(queue["log"])
    log.parent.mkdir(parents=True, exist_ok=True)
    handle = log.open("ab")
    try:
        process = subprocess.Popen(
            [PY, str(queue["runner"])],
            cwd=str(ROOT),
            stdout=handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
    finally:
        handle.close()
    deadline = time.time() + 60
    while time.time() < deadline:
        if process.poll() is not None:
            return None, "runner exited early with code %s" % process.returncode
        data = read_json(queue["controller"])
        if data and data.get("pid") == process.pid:
            return process.pid, None
        time.sleep(2)
    return process.pid, "controller json did not confirm pid within 60s"


def main():
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    lock = (STATE_DIR / ".supervisor.lock").open("a+")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        print("another supervisor holds the lock", flush=True)
        return 1

    state = {
        "pid": os.getpid(),
        "started_at_utc": utc(),
        "python": PY,
        "poll_seconds": POLL_SECONDS,
        "scope": "transport recovery only; frozen plans and recipes unchanged",
        "queues": {},
    }
    for queue in QUEUES:
        state["queues"][queue["key"]] = {
            "kind": queue["kind"],
            "controller": str(Path(queue["controller"]).relative_to(ROOT)),
            "runner": str(Path(queue["runner"]).relative_to(ROOT)),
            "runner_sha256": sha(queue["runner"]),
            "restarts": 0,
            "terminal": None,
            "last_seen_status": None,
            "last_completed": completed_count(queue["controller"]),
            "completed_at_last_restart": None,
            "no_progress_restarts": 0,
            "last_restart_unix": 0,
            "quick_stops": 0,
            "cooldown_seconds": 0,
            "backoff_seconds": 0,
            "next_attempt_unix": 0,
            "upstream_completed_at_last_start": None,
            "receipts": [],
        }
    previous = STATE_DIR / "SUPERVISOR.json"
    if previous.is_file():
        keep = STATE_DIR / ("SUPERVISOR_previous_%s.json" % time.strftime("%Y%m%dT%H%M%SZ", time.gmtime()))
        shutil.copy2(str(previous), str(keep))
    dump(STATE_DIR / "SUPERVISOR.json", state)
    log_line("supervisor started pid=%d" % os.getpid())

    cycle = 0
    while True:
        for queue in QUEUES:
            entry = state["queues"][queue["key"]]
            if entry["terminal"]:
                continue
            data = read_json(queue["controller"])
            if data is None:
                # The controller state was archived and no process holds the queue's
                # flock, so the queue is genuinely down and is started below.
                data = {}
                status = "controller_json_missing"
            else:
                status = data.get("status")
            entry["last_seen_status"] = status
            count = completed_count(queue["controller"]) if data else None
            if count is not None:
                entry["last_completed"] = count
            if alive(data.get("pid"), Path(queue["runner"]).name):
                entry["backoff_seconds"] = 0
                continue
            if status == "complete" and data.get("complete"):
                entry["terminal"] = "complete"
                log_line("%s reached complete" % queue["key"])
                continue
            if status in HUMAN:
                entry["terminal"] = "needs_human:%s" % status
                log_line("%s needs a human decision (%s)" % (queue["key"], status))
                continue
            if status not in RESUMABLE:
                entry["terminal"] = "needs_human:unexpected_status:%s" % status
                log_line("%s stopped with unexpected status %s" % (queue["key"], status))
                continue

            hits = network_failures(queue["failure_root"])
            pending = []
            if queue["kind"] == "source":
                remaining = data.get("remaining_unscheduled")
                pending = missing_receipts(queue["failure_root"])
                entry["members_without_receipt"] = pending
                if remaining == 0 and not hits and not pending:
                    entry["terminal"] = "finished_with_preserved_failures"
                    entry["preserved_failures"] = preserved_failures(queue["failure_root"])
                    log_line(
                        "%s submitted every source; only non-transport failures remain"
                        % queue["key"]
                    )
                    continue
            else:
                upstream_now = 0
                for key in queue["upstream"]:
                    value = state["queues"][key]["last_completed"]
                    upstream_now += value or 0
                upstream_done = all(
                    state["queues"][key]["terminal"] for key in queue["upstream"]
                )
                previous = entry["upstream_completed_at_last_start"]
                if (
                    previous is not None
                    and upstream_now - previous < CONSUMER_UPSTREAM_DELTA
                    and not upstream_done
                ):
                    continue

            # Did this queue stop almost immediately after we restarted it?
            since = time.time() - (entry.get("last_restart_unix") or 0)
            if entry.get("last_restart_unix") and since < QUICK_STOP_SECONDS:
                entry["quick_stops"] = entry.get("quick_stops", 0) + 1
                entry["cooldown_seconds"] = min(
                    COOLDOWN_MAX, COOLDOWN_START * (2 ** (entry["quick_stops"] - 1)))
                entry["next_attempt_unix"] = time.time() + entry["cooldown_seconds"]
                log_line("%s stopped %.0fs after restart (%d in a row); cooling down %ds"
                         % (queue["key"], since, entry["quick_stops"], entry["cooldown_seconds"]))
                dump(STATE_DIR / "SUPERVISOR.json", state)
                continue
            if entry.get("last_restart_unix") and since >= QUICK_STOP_SECONDS:
                entry["quick_stops"] = 0
                entry["cooldown_seconds"] = 0

            if time.time() < entry["next_attempt_unix"]:
                continue

            if entry["completed_at_last_restart"] is not None and (
                entry["last_completed"] or 0
            ) <= entry["completed_at_last_restart"]:
                entry["no_progress_restarts"] += 1
            else:
                entry["no_progress_restarts"] = 0
            if entry["no_progress_restarts"] >= NO_PROGRESS_LIMIT:
                entry["terminal"] = "needs_human:no_progress_after_%d_restarts" % (
                    entry["no_progress_restarts"]
                )
                log_line("%s made no progress across restarts; stopping" % queue["key"])
                continue

            probe_results = {}
            blocked = False
            for url in queue["probes"]:
                ok, detail = probe(url)
                probe_results[url] = detail if not ok else "ok:%s" % detail
                if not ok:
                    blocked = True
            if blocked:
                entry["backoff_seconds"] = min(
                    BACKOFF_MAX, max(BACKOFF_START, entry["backoff_seconds"] * 2)
                )
                entry["next_attempt_unix"] = time.time() + entry["backoff_seconds"]
                log_line(
                    "%s endpoint still failing, waiting %ds: %s"
                    % (queue["key"], entry["backoff_seconds"], json.dumps(probe_results))
                )
                dump(STATE_DIR / "SUPERVISOR.json", state)
                continue

            cycle += 1
            archive_path, archive_receipt = archive_controller(queue, cycle)
            moved = archive_network_failures(queue, cycle, hits)
            entry["completed_at_last_restart"] = entry["last_completed"] or 0
            if queue["kind"] == "consumer":
                entry["upstream_completed_at_last_start"] = sum(
                    (state["queues"][key]["last_completed"] or 0)
                    for key in queue["upstream"]
                )
            pid, problem = launch(queue)
            receipt = {
                "cycle": cycle,
                "restarted_at_utc": utc(),
                "queue": queue["key"],
                "previous_status": status,
                "previous_completed": entry["completed_at_last_restart"],
                "controller_archive": archive_path,
                "controller_archive_sha256": archive_receipt["files"],
                "transport_failures_archived": moved,
                "members_without_receipt": pending,
                "endpoint_probes": probe_results,
                "new_pid": pid,
                "launch_problem": problem,
            }
            entry["receipts"].append(receipt)
            entry["restarts"] += 1
            entry["last_restart_unix"] = time.time()
            dump(
                STATE_DIR / ("RESTART_%03d_%s.json" % (cycle, queue["key"])),
                receipt,
            )
            if problem:
                entry["backoff_seconds"] = min(
                    BACKOFF_MAX, max(BACKOFF_START, entry["backoff_seconds"] * 2)
                )
                entry["next_attempt_unix"] = time.time() + entry["backoff_seconds"]
                log_line("%s restart problem: %s" % (queue["key"], problem))
            else:
                entry["backoff_seconds"] = 0
                entry["next_attempt_unix"] = 0
                log_line(
                    "%s restarted pid=%s after %s (archived %d transport failures)"
                    % (queue["key"], pid, status, len(moved))
                )

        state["updated_at_utc"] = utc()
        dump(STATE_DIR / "SUPERVISOR.json", state)
        if all(state["queues"][q["key"]]["terminal"] for q in QUEUES):
            state["finished_at_utc"] = utc()
            dump(STATE_DIR / "SUPERVISOR.json", state)
            log_line("every queue reached a terminal state; supervisor exiting")
            return 0
        time.sleep(POLL_SECONDS)


if __name__ == "__main__":
    sys.exit(main())
