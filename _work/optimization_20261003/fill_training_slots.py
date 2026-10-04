"""Start queued seeds in released GPU slots, retaining the frozen recipe.

Only start while the seed's original predecessor is still training, so the
existing scheduler cannot concurrently start the same seed. Its existing
adopted_training contract takes over evaluation and screening afterward.
"""
from pathlib import Path
import json
import os
import subprocess
import sys
import time

W = Path(__file__).parent
A = Path('/public/home/mengxl/dzy/pd_product_assets')
O = A / 'results/optimization_20261003'
pending = [(50, 0, 42, 44), (51, 1, 43, 45)]
recipe = json.loads((O / 'frozen_recipe.json').read_text())
while pending:
    assert json.loads((O / 'status.json').read_text()).get('stage') != 'failed'
    free = [int(s.strip()) for s in subprocess.check_output(
        ['nvidia-smi', '--query-gpu=memory.free', '--format=csv,noheader,nounits'], text=True).splitlines()]
    for item in pending.copy():
        seed, gpu, released_by, original_predecessor = item
        target = O / f'models/seed_{seed}'
        if (target / 'training_inputs.json').exists() or (target / 'adopted_training.json').exists():
            pending.remove(item)
            continue
        if (O / f'models/seed_{original_predecessor}/completed.json').exists():
            print('LEAVE TO EXISTING SCHEDULER', seed, flush=True)
            pending.remove(item)
            continue
        released = (O / f'models/seed_{released_by}/screened.json').exists()
        # An adopted training job can release its GPU before its original
        # worker reaches evaluation; use that idle slot as well.
        if seed == 51:
            released = released or (O / 'models/seed_49/completed.json').exists()
        if not released or free[gpu] < 16000:
            continue
        target.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, '-u', str(W / 'train_ccwm.py'), '--data-dir', str(A / 'interim/rescreen_20261002'),
                   '--out-dir', str(target), '--run-dir', str(target / 'run' / time.strftime('%Y%m%d_%H%M%S')),
                   '--run-id', f'rescreen-seed-{seed}', '--dsm-mode', recipe['dsm_mode'], '--potential-beta', '5',
                   '--sinkhorn-iters', '20', '--seed', str(seed), '--epochs', '200', '--steps-per-epoch', '200',
                   '--sample-latent-a']
        env = os.environ | dict(PYTHONPATH='/public/home/mengxl/dzy/pd_product/src', CUDA_VISIBLE_DEVICES=str(gpu),
                                OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2')
        with (target / 'console.log').open('a') as log:
            p = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                 start_new_session=True, env=env)
        record = dict(seed=seed, gpu=gpu, pid=p.pid, command=command,
                      note='Same frozen recipe; started in a released GPU slot; existing worker adopts and evaluates.')
        (target / 'adopted_training.json').write_text(json.dumps(record, indent=2))
        path = O / 'concurrency_extension.json'
        scheduling = json.loads(path.read_text())
        scheduling['jobs'].append(record)
        path.write_text(json.dumps(scheduling, indent=2))
        print('EARLY START', seed, p.pid, flush=True)
        pending.remove(item)
    if pending:
        time.sleep(30)
