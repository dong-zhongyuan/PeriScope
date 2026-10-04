"""Current result publishing and shared evidence-table helpers."""
from pathlib import Path
import hashlib,subprocess,sys
import pandas as pd
GATES = ['technical_pass', 'reproducible_pass', 'specificity_pass',
         'brain_pass', 'independent_direction_pass']
SORT = ['tier', 'confirmation_positive_fraction', 'q_protein_decoy',
        'q_gene_decoy', 'brain_meta_q', 'protein_rho']
ASC = [True, False, True, True, True, False]


def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def normalized(root,primary=False):
    e=pd.read_csv(Path(root)/'all_program_evidence.csv')
    return e.sort_values(SORT,ascending=ASC).reset_index(drop=True)
if __name__ == '__main__':
    W=Path(__file__).resolve().parent
    subprocess.run([sys.executable,str(W/'refresh_model_outputs.py')],check=True)
    subprocess.run([sys.executable,str(W/'refresh_current_pipeline.py')],check=True)
