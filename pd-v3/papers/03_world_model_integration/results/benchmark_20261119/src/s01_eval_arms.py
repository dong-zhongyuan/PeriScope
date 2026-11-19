"""Score every ablation checkpoint through the frozen benchmark harness on CPU, one code path for all arms.
No training. Reads checkpoints read-only; writes only under benchmark_20261119/src/."""
import os, sys, json, time, hashlib, shutil
sys.dont_write_bytecode = True
for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ[v] = '12'
os.environ['CUDA_VISIBLE_DEVICES'] = ''
os.environ['BENCHMARK_DEVICE'] = 'cpu'
os.environ['BENCHMARK_SEED'] = '42'
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
A = Path('/public/home/mengxl/dzy/pd_product_assets')
BENCH = A / 'results/benchmark_20261002'
OPT = A / 'results/optimization_20261003'
HARNESS = Path('/public/home/mengxl/dzy/pd_product/_work/benchmark_shared_20261003')
RES = Path('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/research_20261118')
ABL = RES / 'ablation'
RE = ABL / 'retrain'
ROOT = HERE / 'harness_root'
ROOT.mkdir(exist_ok=True)
for n in ('cache', 'representations', 'prepared_eight.json'):
    if not (ROOT / n).exists():
        (ROOT / n).symlink_to(BENCH / n)
if not (ROOT / 'checkpoint_overrides.json').exists():
    shutil.copy(BENCH / 'checkpoint_overrides.json', ROOT / 'checkpoint_overrides.json')
sys.path[:0] = [str(ABL / 'train'), str(HARNESS), str(HARNESS / 'source')]
import ccwm  # patched copy with cfg.ablate; identical otherwise
import torch
import run_benchmark as rb
rb.O = ROOT
import current_benchmark as cb
assert cb.O == ROOT
import pandas as pd
from scipy.stats import rankdata

torch.set_num_threads(12)
assert 'research_20261118/ablation/train' in ccwm.__file__, ccwm.__file__
TMP = HERE / 'assess_scratch'
TMP.mkdir(exist_ok=True)
cb.RESULT_DIR = TMP
np.savez_compressed = lambda *a, **k: None  # the harness also dumps every predicted cell; not needed here
EV = HERE / 'eval'
EV.mkdir(exist_ok=True)
(EV / 'deltas').mkdir(exist_ok=True)
CSV = EV / 'per_axis.csv'
META = EV / 'ckpt_meta.jsonl'
over = json.loads((BENCH / 'checkpoint_overrides.json').read_text())['checkpoints']

jobs = []
for s in range(42, 52):
    jobs.append(('full', s, Path(over.get(str(s), OPT / f'models/seed_{s}/model.pt'))))
for arm in ('random_pairing', 'coupling_off'):
    for s in (42, 43, 44):
        jobs.append((f'{arm}__published_ckpt', s, OPT / f'ablations/{arm}/models/seed_{s}/model.pt'))
for s in (42, 43, 44):
    jobs.append(('random_pairing__retrained', s, RE / f'random_pairing/seed_{s}/model.pt'))
jobs.append(('coupling_off__retrained', 44, RE / 'coupling_off/seed_44/model.pt'))
for arm, seeds in (('no_blood', (42, 43)), ('no_cond2', (42, 43, 44)), ('no_cond2_random', (42, 43, 44)), ('label_shuffle', (42,))):
    for s in seeds:
        jobs.append((arm, s, RE / f'{arm}/seed_{s}/model.pt'))


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load(p):
    d = torch.load(p, map_location='cpu', weights_only=False)
    m = ccwm.CCWM(ccwm.CCWMConfig(**d['cfg'])).to('cpu')
    m.load_state_dict(d['state_dict'])
    m.eval()
    for q in m.parameters():
        q.requires_grad_(False)
    a = d.get('args', {}) or {}
    meta = dict(best_epoch=d.get('best_epoch', d.get('epoch')), best_validation_score=d.get('best_validation_score'), step=d.get('step'),
                beta_kl=d['cfg'].get('beta_kl'), kl_anneal_epochs=d['cfg'].get('kl_anneal_epochs'), ablate=d['cfg'].get('ablate', 'none'),
                coupling_cfg=d['cfg'].get('coupling'), pairing=a.get('pairing'), coupling_arg=a.get('coupling'), min_epochs=a.get('min_epochs'),
                patience=a.get('patience'), label_shuffle=a.get('label_shuffle'), epochs_arg=a.get('epochs'))
    return m, meta


def own(a, b):
    return float(np.corrcoef(rankdata(a), rankdata(b))[0, 1])


done = set()
if CSV.exists():
    o = pd.read_csv(CSV)
    done = set(zip(o.arm, o.seed))
t0 = time.time()
for arm, s, p in jobs:
    if (arm, s) in done:
        continue
    m, meta = load(p)
    nb = 'no_blood' in (meta['ablate'] or '')
    rows, deltas = [], {}
    for axis in cb.AXES:
        k = dict(np.load(ROOT / 'cache' / f'{axis}.npz'))
        for bc in (0, 1):
            z, mu, u = cb.encode(m, k[f'test_x{bc}'], seed=42 + 100 * bc)
            if nb:
                z = np.random.default_rng(5000 + bc).normal(size=z.shape).astype(np.float32)
                u = np.zeros_like(u)
            k[f'test_z{bc}'] = z
            k[f'test_u{bc}'] = u
        pred = {c: {bc: cb.full_predict(m, k, c, bc) for bc in (0, 1)} for c in (0, 1)}
        cb.assess('x', pred, k, s, None)
        r = json.loads((TMP / f'x__{k["blood_state"]}__{k["brain_state"]}.json').read_text())
        da, di = r['distribution_average'], r['direction']
        obs = k['profiles1'].mean(0) - k['profiles0'].mean(0)
        ca = pred[1][1].mean(0) - pred[0][0].mean(0)
        deltas[axis] = ca.astype(np.float32)
        ca_d = di['condition_aware']
        rows.append(dict(arm=arm, seed=s, epoch=meta['best_epoch'], axis=axis, blood_state=str(k['blood_state']), brain_state=str(k['brain_state']),
                         rho=ca_d['rho'], own_rho=own(ca, obs), fixed_control_rho=di['fixed_control']['rho'], fixed_PD_rho=di['fixed_PD']['rho'],
                         delta_mse=ca_d['delta_mse'], pred_rms=ca_d['predicted_rms'], obs_rms=ca_d['observed_rms'], n_PD=ca_d['n_brain_PD'], n_ctrl=ca_d['n_brain_control'],
                         mmd=da['mmd'], mean_mse=da['mean_mse'], energy=da['energy_per_sqrt_gene'], variance_ratio=da['variance_ratio']))
    pd.DataFrame(rows).to_csv(CSV, mode='a', header=not CSV.exists(), index=False)
    np.savez(EV / 'deltas' / f'{arm}__s{s}.npz', **deltas)
    meta.update(arm=arm, seed=s, path=str(p), sha256=sha(p))
    with open(META, 'a') as f:
        f.write(json.dumps(meta) + '\n')
    print(arm, s, 'epoch', meta['best_epoch'], 'rho', round(np.mean([x['rho'] for x in rows]), 4), 'mmd', round(np.mean([x['mmd'] for x in rows]), 4), f'{time.time() - t0:.0f}s', flush=True)
shutil.rmtree(TMP, ignore_errors=True)
print('DONE', flush=True)
