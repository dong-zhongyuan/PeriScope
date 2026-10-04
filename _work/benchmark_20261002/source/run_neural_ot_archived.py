"""p2 benchmark arms — neural OT arm (CellOT route).

Prereg: docs/BENCHMARK_ARMS_PREREGISTRATION.md (frozen 2026-09-05, commit 42f8f0a).
Upstream code: reference/cellot (github.com/bunnech/cellot @522d2b95, BSD-3-Clause) —
ICNN potentials, compute_loss_f/compute_loss_g, and the update order
(10 inner g-updates + fnorm penalty, 1 f-update, f.clamp_w) used verbatim;
donor-aware data plumbing and the frozen evaluations live here.

Leg A (组织-rev): per condition c, blood myeloid cells (state {0,1,3,4}) -> brain
  microglia cells (state {6,7} pooled); Sinkhorn divergence to locked_test brain
  cells of the same condition, vs naive-marginal and identity references.
Leg B (疾病响应): per brain state s in {5,6,7}, Ctrl->PD transport trained on
  brain_train donors only; Spearman(pred_diff, obs_diff) on locked_test;
  20-draw train-donor label-permutation null; one seed-43 stability run.

Fixed budget (chosen before the first result run, no post-result tuning):
n_iters=5000 (upstream cellot.yaml default is 100000; reduced for the 44-model
main+null+stability schedule — recorded), batch=256 sampled with replacement
(uniform handling of small strata), hidden [64,64,64,64], Adam lr=1e-4
betas=(0.5,0.9), g fnorm_penalty=1, uniform kernel init b=0.1 — all other values
are cellot.yaml defaults.
"""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HIDDEN = [64, 64, 64, 64]
N_INNER = 10
BATCH = 256
N_ITERS = 5000
MYELOID_STATES = (0, 1, 3, 4)      # classical_mono, nonclassical_mono, cDC, pDC
MICROGLIA_STATES = (6, 7)          # homeostatic + mhc2 (leg A pools them)
LEGB_STATES = {5: "astro", 6: "microglia_homeostatic", 7: "microglia_mhc2"}
STATE_NAMES = {0: "classical_mono", 1: "nonclassical_mono", 3: "cDC", 4: "pDC",
               5: "astro", 6: "microglia_homeostatic", 7: "microglia_mhc2",
               9: "neural_other", 10: "none"}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_fg(input_dim: int, seed: int, device: str):
    import torch
    from cellot.networks.icnns import ICNN

    torch.manual_seed(seed)
    init = lambda w: torch.nn.init.uniform_(w, 0.0, 0.1)  # cellot.yaml kernel_init_fxn
    f = ICNN(input_dim, HIDDEN, kernel_init_fxn=init).to(device)
    g = ICNN(input_dim, HIDDEN, fnorm_penalty=1, kernel_init_fxn=init).to(device)
    opt_f = torch.optim.Adam(f.parameters(), lr=1e-4, betas=(0.5, 0.9))
    opt_g = torch.optim.Adam(g.parameters(), lr=1e-4, betas=(0.5, 0.9))
    return f, g, opt_f, opt_g


def train_transport(src: np.ndarray, tgt: np.ndarray, seed: int, device: str,
                    n_iters: int = N_ITERS):
    """Upstream CellOT minimax update order, verbatim losses."""
    import torch
    from cellot.models.cellot import compute_loss_f, compute_loss_g

    f, g, opt_f, opt_g = build_fg(src.shape[1], seed, device)
    src_t = torch.from_numpy(np.asarray(src, dtype=np.float32)).to(device)
    tgt_t = torch.from_numpy(np.asarray(tgt, dtype=np.float32)).to(device)
    rng = np.random.default_rng(seed)
    for step in range(n_iters):
        tb = tgt_t[rng.integers(0, tgt_t.shape[0], BATCH)]
        for _ in range(N_INNER):
            sb = src_t[rng.integers(0, src_t.shape[0], BATCH)].detach().requires_grad_(True)
            opt_g.zero_grad()
            gl = compute_loss_g(f, g, sb).mean()
            if not g.softplus_W_kernels and g.fnorm_penalty > 0:
                gl = gl + g.penalize_w()
            gl.backward()
            opt_g.step()
        sb = src_t[rng.integers(0, src_t.shape[0], BATCH)].detach().requires_grad_(True)
        opt_f.zero_grad()
        fl = compute_loss_f(f, g, sb, tb).mean()
        fl.backward()
        opt_f.step()
        f.clamp_w()
        if step % 5000 == 0:
            print(f"    iter {step}/{n_iters} gloss={float(gl):.4f} floss={float(fl):.4f}",
                  flush=True)
    return g


def transport(g, X: np.ndarray, device: str, chunk: int = 2048) -> np.ndarray:
    import torch

    outs = []
    for lo in range(0, X.shape[0], chunk):
        x = torch.from_numpy(np.asarray(X[lo:lo + chunk], dtype=np.float32)).to(device)
        x.requires_grad_(True)
        outs.append(g.transport(x).detach().cpu().numpy())
    return np.concatenate(outs).astype(np.float64)


def diff_spearman(pred_diff: np.ndarray, obs_diff: np.ndarray) -> float:
    from scipy.stats import spearmanr

    if np.ptp(pred_diff) == 0 or np.ptp(obs_diff) == 0:
        return float("nan")
    return float(spearmanr(pred_diff, obs_diff).statistic)


def subsample(X: np.ndarray, n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    idx = rng.choice(X.shape[0], size=min(n, X.shape[0]), replace=False)
    return X[idx]


def sinkhorn(x: np.ndarray, y: np.ndarray, device: str, eps: float = 0.1) -> float:
    import torch
    from pdproduct.simulators.ccwm import sinkhorn_divergence

    xt = torch.from_numpy(np.asarray(x, dtype=np.float32)).to(device)
    yt = torch.from_numpy(np.asarray(y, dtype=np.float32)).to(device)
    return float(sinkhorn_divergence(xt, yt, eps=eps))


def leg_a(args, blood, brain_train, brain_test, device: str) -> dict:
    out = {}
    bm = np.isin(blood["state"], list(MYELOID_STATES))
    gm_tr = np.isin(brain_train["state"], list(MICROGLIA_STATES))
    gm_te = np.isin(brain_test["state"], list(MICROGLIA_STATES))
    for c, cname in ((0, "Ctrl"), (1, "PD")):
        src = blood["X"][bm & (blood["cond2"] == c)]
        tgt = brain_train["X"][gm_tr & (brain_train["cond2"] == c)]
        obs = brain_test["X"][gm_te & (brain_test["cond2"] == c)]
        print(f"  [legA {cname}] src {src.shape} tgt {tgt.shape} obs {obs.shape}", flush=True)
        g = train_transport(src, tgt, args.seed, device, args.n_iters)
        trans = transport(g, subsample(src, 2000, args.seed), device)
        obs_s = subsample(obs, 2000, args.seed)
        out[cname] = {
            "n_src": int(src.shape[0]), "n_tgt_train": int(tgt.shape[0]),
            "n_obs_test": int(obs.shape[0]),
            "sinkhorn_cellot": sinkhorn(trans, obs_s, device),
            "sinkhorn_naive_marginal": sinkhorn(subsample(tgt, 2000, args.seed), obs_s, device),
            "sinkhorn_identity_floor": sinkhorn(subsample(src, 2000, args.seed + 1), obs_s, device),
        }
        print(f"  [legA {cname}] {out[cname]}", flush=True)
    return out


def leg_b_one(brain_train, brain_test, state: int, seed: int, device: str,
              n_iters: int, label_perm: np.ndarray | None = None):
    """One Ctrl->PD transport. label_perm: optional dict donor->permuted cond2."""
    tr, te = brain_train, brain_test
    if label_perm is None:
        cond_tr = tr["cond2"]
    else:
        cond_tr = np.array([label_perm[int(d)] for d in tr["donor"]], dtype=np.int64)
    src = tr["X"][(tr["state"] == state) & (cond_tr == 0)]
    tgt = tr["X"][(tr["state"] == state) & (cond_tr == 1)]
    x_ctrl = te["X"][(te["state"] == state) & (te["cond2"] == 0)]
    x_pd = te["X"][(te["state"] == state) & (te["cond2"] == 1)]
    if src.shape[0] < 50 or tgt.shape[0] < 50 or x_ctrl.shape[0] < 50 or x_pd.shape[0] < 50:
        return {"status": "insufficient_cells", "n": [int(src.shape[0]), int(tgt.shape[0]),
                                                      int(x_ctrl.shape[0]), int(x_pd.shape[0])]}
    g = train_transport(src, tgt, seed, device, n_iters)
    pred = transport(g, x_ctrl, device)
    pred_diff = pred.mean(0) - x_ctrl.mean(0)
    obs_diff = x_pd.mean(0) - x_ctrl.mean(0)
    rho = diff_spearman(pred_diff, obs_diff)
    per_donor = {}
    for d in sorted(set(te["donor"][(te["state"] == state) & (te["cond2"] == 0)].tolist())):
        m = (te["donor"] == d) & (te["state"] == state) & (te["cond2"] == 0)
        if m.sum() < 50:
            continue
        pd_d = transport(g, te["X"][m], device)
        dd = pd_d.mean(0) - te["X"][m].mean(0)
        per_donor[str(d)] = diff_spearman(dd, obs_diff)
    return {"status": "ok", "rho": rho, "per_donor_ctrl": per_donor,
            "n_train_src": int(src.shape[0]), "n_train_tgt": int(tgt.shape[0]),
            "n_test_ctrl": int(x_ctrl.shape[0]), "n_test_pd": int(x_pd.shape[0])}


def leg_b(args, brain_train, brain_test, device: str) -> dict:
    rng = np.random.default_rng(args.seed)
    donors_tr = sorted(set(brain_train["donor"].tolist()))
    base_labels = {int(d): int(brain_train["cond2"][brain_train["donor"] == d][0]) for d in donors_tr}
    out = {}
    for s, sname in LEGB_STATES.items():
        print(f"  [legB {sname}] main seed={args.seed}", flush=True)
        main = leg_b_one(brain_train, brain_test, s, args.seed, device, args.n_iters)
        print(f"  [legB {sname}] main rho={main.get('rho')}", flush=True)
        stab = leg_b_one(brain_train, brain_test, s, 43, device, args.n_iters)
        print(f"  [legB {sname}] seed43 rho={stab.get('rho')}", flush=True)
        nulls = []
        if main.get("status") == "ok":
            for k in range(args.n_null):
                pd_donors = [d for d in donors_tr if base_labels[d] == 1]
                perm_pd = set(rng.choice(donors_tr, size=len(pd_donors), replace=False).tolist())
                perm = {int(d): (1 if d in perm_pd else 0) for d in donors_tr}
                r = leg_b_one(brain_train, brain_test, s, args.seed, device, args.n_iters,
                              label_perm=perm)
                nulls.append(r.get("rho"))
                print(f"  [legB {sname}] null {k + 1}/{args.n_null} rho={r.get('rho')}",
                      flush=True)
        nulls = [x for x in nulls if x is not None and not np.isnan(x)]
        emp_p = float((1 + sum(1 for x in nulls if x >= main["rho"])) / (1 + len(nulls))) \
            if main.get("status") == "ok" else None
        out[sname] = {"main_seed42": main, "seed43_stability": stab,
                      "null_rhos": nulls, "empirical_p": emp_p}
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--blood-npz", required=True)
    ap.add_argument("--brain-train-npz", required=True)
    ap.add_argument("--brain-test-npz", required=True)
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--leg", choices=["A", "B", "AB"], default="AB")
    ap.add_argument("--n-iters", type=int, default=N_ITERS)
    ap.add_argument("--n-null", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    import torch
    from pdproduct.provenance import finalize_run_manifest, new_run_manifest, write_run_manifest

    device = "cuda" if torch.cuda.is_available() else "cpu"
    manifest = new_run_manifest(
        f"neural-ot-{args.leg.lower()}-v1",
        approved_scope="p2 benchmark arms (prereg docs/BENCHMARK_ARMS_PREREGISTRATION.md, commit 42f8f0a)",
        command=f"python scripts/p2_benchmark_arms/run_neural_ot.py --leg {args.leg}",
        seed=args.seed,
    )
    write_run_manifest(manifest, args.run_dir)

    blood = np.load(args.blood_npz)
    brain_train = np.load(args.brain_train_npz)
    brain_test = np.load(args.brain_test_npz)
    t0 = time.time()
    results = {
        "upstream": "github.com/bunnech/cellot @522d2b953da8ad244fcf36f64521487fdf763788 (BSD-3-Clause); ICNN/losses verbatim, donor-aware plumbing here",
        "input_sha256": {"blood": sha256(args.blood_npz),
                          "brain_train": sha256(args.brain_train_npz),
                          "brain_test": sha256(args.brain_test_npz)},
        "budget": {"n_iters": args.n_iters, "batch": BATCH, "hidden": HIDDEN,
                   "lr": 1e-4, "betas": [0.5, 0.9], "n_inner": N_INNER,
                   "g_fnorm_penalty": 1, "init": "uniform[0,0.1)", "seed": args.seed},
        "states": {str(k): v for k, v in STATE_NAMES.items()},
    }
    if args.leg in ("A", "AB"):
        results["leg_A"] = leg_a(args, blood, brain_train, brain_test, device)
    if args.leg in ("B", "AB"):
        results["leg_B"] = leg_b(args, brain_train, brain_test, device)
    results["runtime_min"] = round((time.time() - t0) / 60, 1)

    out = Path(args.out_json)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(results, fh, indent=2, ensure_ascii=False)
    finalize_run_manifest(manifest, status="completed", exit_code=0, output_paths=[out],
                          run_dir=args.run_dir)
    print(json.dumps({"done": True, "runtime_min": results["runtime_min"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
