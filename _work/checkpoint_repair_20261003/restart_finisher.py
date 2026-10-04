import json,os,signal,subprocess,sys
from pathlib import Path
W=Path('/public/home/mengxl/dzy/pd_product/_work/checkpoint_repair_20261003')
R=Path('/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002/seed_repair')
assert json.loads((R/'run_state.json').read_text())['status']=='training'
for p in Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:args=(p/'cmdline').read_bytes().split(b'\0')
 except (FileNotFoundError,PermissionError):continue
 if str(W/'finish_repair.py').encode() in args:
  os.kill(int(p.name),signal.SIGTERM);print('RESTART_WAITING_FINISHER',p.name)
with (W/'finish_repair.log').open('a') as f:
 p=subprocess.Popen([sys.executable,'-u',str(W/'finish_repair.py')],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
 print('NEW_FINISHER',p.pid)
