from pathlib import Path
import json,time,sys,subprocess,traceback
W=Path(__file__).parent;O=Path('/public/home/mengxl/dzy/pd_product_assets/results/rescreen_20261002')
while json.loads((O/'status.json').read_text()).get('stage')!='inputs_built':time.sleep(15)
try:
 for script in ['data_checks.py','candidate_registry.py','engineering_checks.py']:
  with (O/(script.replace('.py','')+'.log')).open('w') as f:subprocess.run([sys.executable,'-u',str(W/script)],stdout=f,stderr=subprocess.STDOUT,check=True)
 (O/'status.json').write_text(json.dumps(dict(stage='training_pilots',beta=[1,5],seed=42,steps_each=5000)))
 subprocess.run([sys.executable,'-u',str(W/'run_pilots.py')],check=True)
except Exception as e:
 (O/'status.json').write_text(json.dumps(dict(stage='failed',error=str(e)),indent=2));raise
