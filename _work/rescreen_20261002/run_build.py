from pathlib import Path
import subprocess,json,sys,time
from input_utils import REPORT as r, ROOT
W=Path(__file__).parent
cmd=[sys.executable,'-u',str(W/'build_training_data.py')]
for arg,key in [('brain-h5ad','brain'),('blood-h5ad','blood223'),('blood2-h5ad','bloodmb'),('split-json','split')]:cmd+=['--'+arg,r['inputs'][key]['path']]
cmd+=['--purify-dir',r['purify_dir'],'--cite-npz',str(ROOT/'bridge_data.npz'),'--out-dir',str(ROOT),'--run-dir',str(ROOT/'build_run'/time.strftime('%Y%m%d_%H%M%S'))]
subprocess.run(cmd,check=True)
