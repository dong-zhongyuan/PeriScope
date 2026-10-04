from pathlib import Path
import subprocess,json,sys,time
A=Path('/public/home/mengxl/dzy/pd_product_assets');W=Path(__file__).parent
r=json.loads((A/'interim/ccwm_v23/build_report.json').read_text())
cmd=[sys.executable,'-u',str(W/'build_training_data.py')]
for arg,key in [('brain-h5ad','brain'),('blood-h5ad','blood223'),('blood2-h5ad','bloodmb'),('split-json','split')]:cmd+=['--'+arg,r['inputs'][key]['path']]
cmd+=['--purify-dir',r['purify_dir'],'--cite-npz',str(A/'interim/rescreen_20261002/bridge_data.npz'),'--out-dir',str(A/'interim/rescreen_20261002'),'--run-dir',str(A/'results/rescreen_20261002/build_run'/time.strftime('%Y%m%d_%H%M%S'))]
subprocess.run(cmd,check=True)
