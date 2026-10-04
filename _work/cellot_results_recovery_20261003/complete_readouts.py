"""Finish CellOT readouts after verified seed42 weights have been selected."""
import subprocess,sys
from pathlib import Path
W=Path(__file__).resolve().parent
for name in ['recompute_cellot_readout.py','recompute_cellot_pig.py','validate_recovery.py']:
 with (W/(name+'.log')).open('w') as log:subprocess.run([sys.executable,str(W/name)],stdout=log,stderr=subprocess.STDOUT,check=True)
 print('FINISHED',name,flush=True)
