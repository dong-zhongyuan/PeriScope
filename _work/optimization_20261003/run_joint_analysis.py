"""Compatibility entry point for the maintained current-model evidence pipeline."""
from pathlib import Path
import subprocess,sys
if __name__ == '__main__':
    W=Path(__file__).resolve().parent
    subprocess.run([sys.executable,str(W/'refresh_model_outputs.py')],check=True)
    subprocess.run([sys.executable,str(W/'refresh_current_pipeline.py')],check=True)
