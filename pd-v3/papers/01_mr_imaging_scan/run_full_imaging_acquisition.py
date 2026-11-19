"""Durable original-plan acquisition; restartable, with bounded retry and live state."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from pdv3.tabix_http import save_json,digest


def main():
    dest=ROOT/'data/big40_atlas'
    with (dest/'.runner.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        state={'pid':os.getpid(),'started_at_utc':datetime.now(timezone.utc).isoformat(),
               'scope':'all 3935 official BIG40 entries, no outcome-based subset',
               'status':'running','attempts':[],
               'plan_sha256':digest((dest/'phenotypes_3935.csv').read_bytes()),
               'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}
        save_json(dest/'RUN.json',state)
        command=[sys.executable,'-u',str(Path(__file__).with_name('s16_extract_full_imaging.py')),
                 '--reuse-native','--workers','8']
        for attempt in range(1,4):
            code=subprocess.call(command,cwd=ROOT)
            state['attempts'].append({'attempt':attempt,'exit_code':code,
                                      'ended_at_utc':datetime.now(timezone.utc).isoformat()})
            if code==0:
                progress=json.loads((dest/'extraction_progress.json').read_text())
                if not progress['complete'] or progress['requested']!=3935:
                    raise RuntimeError('Successful exit without complete original catalog')
                state['status']='acquisition_complete'
                state['completed_at_utc']=datetime.now(timezone.utc).isoformat()
                save_json(dest/'RUN.json',state)
                return
            save_json(dest/'RUN.json',state)
        state['status']='needs_attention';save_json(dest/'RUN.json',state)
        raise SystemExit(1)


if __name__=='__main__':main()
