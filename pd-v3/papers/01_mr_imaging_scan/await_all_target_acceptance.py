"""Run the acceptance check as soon as the five all-target queues are all finished.

So the 20-hour run ends with a written verdict instead of needing someone to
remember to check. Polls the supervisor's own state plus each controller receipt,
runs verify_all_target_completion.py once everything is terminal, and records the
result. Read-only apart from its own report.
"""
from pathlib import Path
import json, subprocess, sys, time
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
SCRIPTS = ROOT / "papers/01_mr_imaging_scan"
BASE = ROOT / "data/all_target_followup_20261118"
OUT = BASE / "acceptance_20261118"
PY = "/public/home/mengxl/dzy/pd-v3_runtime/bin/pdv3-python"
POLL = 300
DEADLINE_HOURS = 48

CONTROLLERS = {
    "protein_sources": BASE / "regional_all/CONTROLLER.json",
    "big40_imaging": BASE / "regional_all/imaging/BIG40_CONTROLLER.json",
    "qsm_peak_outcomes": BASE / "regional_all/reverse_supplement/protein_outcomes/CONTROLLER.json",
    "all_target_coloc": ROOT / "papers/01_mr_imaging_scan/results/all_target_coloc_20261118/CONTROLLER.json",
    "all_target_reverse_mr": ROOT / "papers/01_mr_imaging_scan/results/all_target_reverse_20261118/CONTROLLER.json",
}
LIVE_STATES = {"running", "ready", "waiting_inputs"}


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def read(path):
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return None


def pid_alive(pid, needle):
    if not pid:
        return False
    try:
        raw = Path("/proc/%d/cmdline" % int(pid)).read_bytes()
    except Exception:
        return False
    return bool(raw) and needle.encode() in raw


def supervisor_running():
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            raw = (entry / "cmdline").read_bytes()
        except Exception:
            continue
        if raw and b"supervise_all_target_queues.py" in raw and b"python" in raw:
            return True
    return False


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    state = dict(pid=__import__("os").getpid(), started_at_utc=utc(),
                 waiting_for="all five all-target queues to stop being live")
    (OUT / "AWAIT.json").write_text(json.dumps(state, indent=2) + "\n")
    deadline = time.time() + DEADLINE_HOURS * 3600
    while time.time() < deadline:
        snapshot = {}
        live = False
        for name, path in CONTROLLERS.items():
            data = read(path)
            if data is None:
                snapshot[name] = "controller_json_missing"
                live = True
                continue
            status = data.get("status")
            snapshot[name] = status
            if status in LIVE_STATES or pid_alive(data.get("pid"), Path(path).name):
                live = True
        # The supervisor restarts stopped queues, so it being up means work continues.
        if supervisor_running():
            live = True
        state.update(updated_at_utc=utc(), snapshot=snapshot,
                     supervisor_running=supervisor_running(), live=live)
        (OUT / "AWAIT.json").write_text(json.dumps(state, indent=2) + "\n")
        if not live:
            result = subprocess.run(
                [PY, str(SCRIPTS / "verify_all_target_completion.py")],
                cwd=str(ROOT), capture_output=True, text=True)
            report = read(OUT / "ACCEPTANCE.json") or {}
            state.update(
                finished_at_utc=utc(),
                acceptance_exit_code=result.returncode,
                acceptance_stdout_tail=(result.stdout or "")[-2000:],
                acceptance_stderr_tail=(result.stderr or "")[-1000:],
                fully_complete=report.get("fully_complete"),
                unexplained_gaps=len(report.get("unexplained_gaps", [])),
                reasons_needing_review=report.get("reasons_needing_review"),
                verdict=("all planned regions complete" if report.get("fully_complete")
                         else "stopped with regions still incomplete; see ACCEPTANCE.json blocked_by"),
            )
            (OUT / "AWAIT.json").write_text(json.dumps(state, indent=2) + "\n")
            print(json.dumps({k: state[k] for k in
                              ("verdict", "fully_complete", "unexplained_gaps")}, indent=2))
            return 0
        time.sleep(POLL)
    state.update(finished_at_utc=utc(), verdict="deadline reached while queues were still live")
    (OUT / "AWAIT.json").write_text(json.dumps(state, indent=2) + "\n")
    return 1


if __name__ == "__main__":
    sys.exit(main())
