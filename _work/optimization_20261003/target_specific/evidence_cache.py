"""Reuse only evidence already regenerated in this run with identical program definitions and code."""
from pathlib import Path
import hashlib,json,fcntl
_START={}
O=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/target_specific')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def fingerprint(script):
    return dict(programs=sha(O/'programs.json'),script=sha(script),input_checks=sha(O/'downstream_input_checks.json'),
                symbols=sha(Path(script).parent/'input_utils.py'),reference=sha(Path(script).parent/'reference/hgnc_complete_set.txt'))
def acquire(script):
    p=O/'execution_cache'/'locks'/(Path(script).stem+'.lock');p.parent.mkdir(parents=True,exist_ok=True)
    f=p.open('a');fcntl.flock(f.fileno(),fcntl.LOCK_EX);return f

def hit(script):
    _START[str(script)]=fingerprint(script)
    p=O/'execution_cache'/(Path(script).stem+'.json')
    if not p.exists():return False
    old=json.loads(p.read_text())
    if old['inputs']!=fingerprint(script):return False
    return all((O/name).exists() and sha(O/name)==digest for name,digest in old['outputs'].items())
def finish(script,outputs):
    p=O/'execution_cache'/(Path(script).stem+'.json');p.parent.mkdir(exist_ok=True)
    assert _START[str(script)]==fingerprint(script),'Inputs changed while computing evidence; rerun this stage'
    p.write_text(json.dumps(dict(inputs=_START[str(script)],outputs={name:sha(O/name) for name in outputs}),indent=2))
