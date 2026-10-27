import subprocess,sys,json,time,os
from pathlib import Path
W=Path(__file__).parent
if not set(sys.argv[1:]).issubset({'--build-only'}):
 sys.exit(subprocess.call([sys.executable,str(W.parents[1]/'reproduce/build_inputs.py'),*sys.argv[1:]]))
from input_utils import ROOT
O=ROOT/'preparation_status'
O.mkdir(parents=True,exist_ok=True)
stages=['run_build.py'] if '--build-only' in sys.argv else ['input_utils.py','preprocess_citeseq_hao.py','run_build.py']
for name in stages:
 t=time.time();(O/'status.json').write_text(json.dumps(dict(stage='input_rebuild',running=name,started=time.strftime('%Y-%m-%d %H:%M:%S'))))
 print('START',name,flush=True)
 subprocess.run([sys.executable,'-u',str(W/name)],check=True)
 print('DONE',name,'seconds',time.time()-t,flush=True)
(O/'status.json').write_text(json.dumps(dict(stage='inputs_built',ready_for='data_validation_and_pilot')))
