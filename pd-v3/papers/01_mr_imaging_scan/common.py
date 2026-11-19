import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
PAPER = Path(__file__).resolve().parent
SNAPSHOT = ROOT / 'reference/periscope_mr_20261115'
RESULTS = PAPER / 'results'
FIGURES = PAPER / 'figures'


def provenance(path, inputs, extra=None):
    from pdv3.provenance import record
    code = list((ROOT / 'pdv3').glob('*.py')) + list(PAPER.glob('*.py'))
    return record(path, code_files=code, input_files=inputs, root=ROOT, extra=extra)


def save_csv(frame, path, inputs, extra=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False)
    provenance(path, inputs, extra)
