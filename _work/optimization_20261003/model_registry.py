"""One authoritative ten-seed model set, shared by screening and evaluation."""
from pathlib import Path
import hashlib, json

ROOT = Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')

def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def atomic_json(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + '.pending')
    tmp.write_text(json.dumps(value, indent=2))
    tmp.replace(path)

def checkpoint(root, seed):
    root = Path(root)
    path = root / f'models/seed_{seed}/model.pt'
    registry = ROOT / 'model_registry.json'
    # Ablation models retain their own weights and are never substituted.
    if registry.exists() and root.resolve() in (ROOT.resolve(), (ROOT/'target_specific').resolve()):
        rec = json.loads(registry.read_text())['models'][str(seed)]
        if sha256(path) != rec['sha256']:
            raise RuntimeError(f'Active checkpoint mismatch: seed {seed}')
    return path

def curve_provenance(root, seed):
    root = Path(root)
    return {'seed': seed, 'checkpoint_sha256': sha256(checkpoint(root, seed)),
            'candidate_registry_sha256': sha256(root/'candidate_registry.json'),
            'dose_definitions_sha256': sha256(root/'dose_definitions.json')}

def valid_curve(path, provenance):
    path = Path(path)
    sidecar = path.with_suffix('.provenance.json')
    if not path.exists() or not sidecar.exists():
        return False
    rec = json.loads(sidecar.read_text())
    return all(rec.get(k) == v for k,v in provenance.items()) and rec.get('output_sha256') == sha256(path)

def stamp_curve(path, provenance):
    atomic_json(Path(path).with_suffix('.provenance.json'), dict(provenance, output_sha256=sha256(path)))
