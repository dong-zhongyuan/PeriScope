"""Run a bounded additional data acquisition with the frozen main reader.

Outputs are isolated. No main process is stopped and no MR is calculated.
"""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from pdv3.tabix_http import save_json


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    base = ROOT / 'data/big40_tail_pilot'
    plan_path = base / 'PLAN.json'
    plan = json.loads(plan_path.read_text())
    if (base / 'RUN.json').exists():
        raise RuntimeError('Existing pilot run must be inspected before rerunning')
    for relative, digest in plan['unchanged_source_sha256'].items():
        assert sha(ROOT / relative) == digest
    command = [sys.executable, '-u', str(Path(__file__).with_name('s16_extract_full_imaging.py')),
               '--workers', str(plan['workers']), '--idps',
               ','.join([plan['reference_idp']] + plan['new_idps']), '--output', str(base)]
    started = time.monotonic()
    state = dict(status='running', pid=os.getpid(), command=command, plan_sha256=sha(plan_path),
                 wrapper_sha256=sha(Path(__file__)), started_at_utc=datetime.now(timezone.utc).isoformat())
    with (base / 'acquisition.log').open('x') as log:
        proc = subprocess.Popen(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
        state['child_pid'] = proc.pid
        save_json(base / 'RUN.json', state)
        print(json.dumps(state), flush=True)
        code = proc.wait()
    state.update(exit_code=code, elapsed_seconds=time.monotonic()-started,
                 finished_at_utc=datetime.now(timezone.utc).isoformat(), status='failed')
    for relative, digest in plan['unchanged_source_sha256'].items():
        assert sha(ROOT / relative) == digest
    if code == 0:
        progress = json.loads((base / 'extraction_progress.json').read_text())
        assert progress['complete'] and progress['requested'] == progress['completed'] == 9
        ref = base / 'extracted' / (plan['reference_idp'] + '.csv.gz')
        assert sha(ref) == plan['reference_output_sha256']
        state.update(status='complete_verified_reference', reference_SHA_exact_match=True,
                     all_nine_complete=True, main_data_imported=False)
    save_json(base / 'RUN.json', state)
    print(json.dumps(state), flush=True)
    raise SystemExit(code)


if __name__ == '__main__':
    main()
