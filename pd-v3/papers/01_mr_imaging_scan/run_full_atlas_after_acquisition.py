"""Continue the authorized full scan as soon as all source extractions succeed."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pdv3.tabix_http import save_json


def main():
    atlas=ROOT/'data/big40_atlas'
    output=ROOT/'papers/01_mr_imaging_scan/results/full_atlas_20261116/full'
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verification',type=Path,default=output.parent/'validation45/PIPELINE_VERIFICATION.json')
    verification_path=parser.parse_args().verification.resolve()
    verified=json.loads(verification_path.read_text())
    assert verified['complete']
    # Fail closed if the tested statistical implementation is edited while the
    # source job is running; never silently execute unverified later changes.
    def check_code():
        for relative,expected in verified['source_code_sha256'].items():
            if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=expected:
                raise RuntimeError(f'Validated analysis code changed: {relative}')
    check_code()
    with (atlas/'.analysis_runner.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        state={'pid':os.getpid(),'started_at_utc':datetime.now(timezone.utc).isoformat(),
               'status':'waiting_for_complete_acquisition','scope':'2941 assays x (3935 BIG40 entries + 28 Oxford QSM/T2* entries)',
               'validated_code':verified['source_code_sha256'],
               'verification_path':str(verification_path.relative_to(ROOT)),
               'verification_sha256':hashlib.sha256(verification_path.read_bytes()).hexdigest()}
        save_json(atlas/'ANALYSIS_RUN.json',state)
        while True:
            acquisition=json.loads((atlas/'RUN.json').read_text())
            if acquisition['status']=='acquisition_complete':break
            if acquisition['status']=='needs_attention':
                state['status']='blocked_by_acquisition_errors';save_json(atlas/'ANALYSIS_RUN.json',state);raise SystemExit(1)
            try:os.kill(acquisition['pid'],0)
            except ProcessLookupError:
                state['status']='acquisition_process_missing';save_json(atlas/'ANALYSIS_RUN.json',state);raise SystemExit(1)
            time.sleep(45)
        try:check_code()
        except RuntimeError as error:
            state.update(status='validated_code_changed',error=str(error));save_json(atlas/'ANALYSIS_RUN.json',state);raise
        state['status']='running_statistics';save_json(atlas/'ANALYSIS_RUN.json',state)
        code=subprocess.call([sys.executable,'-u',str(Path(__file__).with_name('s18_run_full_atlas.py'))],cwd=ROOT)
        state.update(status='complete' if code==0 else 'analysis_failed',exit_code=code,
                     ended_at_utc=datetime.now(timezone.utc).isoformat())
        if code==0:
            result=json.loads((output/'RUN.json').read_text())
            if not result['complete'] or result['planned_pairs']!=11655183:raise RuntimeError('Full scan final verification failed')
            state['result']=str((output/'RUN.json').relative_to(ROOT))
        save_json(atlas/'ANALYSIS_RUN.json',state)
        raise SystemExit(code)


if __name__=='__main__':main()
