"""Donor-resolved blood-only readout with an exact donor-label permutation test.
Six held-out blood donors (three control, three PD; all disjoint from the training donors) give C(6,3)=20
label assignments. The whole statistic - equal-route aggregate contrast, Spearman against the observed
held-out brain disease contrast - is recomputed under each assignment, so the test is exact.
usage: python bloodonly_exact.py OUTDIR seed1,seed2,...
"""
import sys, os, itertools, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
sys.path.insert(0, "/public/home/mengxl/dzy/pd_product/_work/benchmark_shared_20261003")
sys.path.insert(0, "/public/home/mengxl/dzy/pd_product/_work/optimization_20261003")
import torch
from ccwm import CCWM, CCWMConfig
A = "/public/home/mengxl/dzy/pd_product_assets"
OUT = sys.argv[1]; SEEDS = [int(x) for x in sys.argv[2].split(",")]
DEV = "cuda" if torch.cuda.is_available() else "cpu"
pd.set_option("display.width", 250)
BS = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]; RS = ["astro", "microglia_mhc2"]
b = dict(np.load(A + "/interim/rescreen_20261002/blood.npz"))
states = list(b["state_vocab"].astype(str)); dn = b["donor_names"].astype(str)
held = b["is_val"] | b["is_test"]
cache = {rc: np.load(A + f"/results/benchmark_20261002/cache/cDC__{rc}.npz", allow_pickle=True) for rc in RS}
obs = {rc: cache[rc]["profiles1"].mean(0) - cache[rc]["profiles0"].mean(0) for rc in RS}

CK = OUT + "/bloodonly_profiles.npz"
if os.path.exists(CK):
    z = np.load(CK, allow_pickle=True); prof = {tuple(k.split("|")): z[k] for k in z.files}
    prof = {(int(a), c, d, int(e), int(f)): v for (a, c, d, e, f), v in prof.items()}
    print("loaded cached profiles:", len(prof))
else:
    def load(seed):
        ck = torch.load(A + f"/results/optimization_20261003/models/seed_{seed}/model.pt", map_location=DEV, weights_only=False)
        m = CCWM(CCWMConfig(**ck["cfg"])).to(DEV); m.load_state_dict(ck["state_dict"]); m.eval(); return m
    def encode(model, x, seed):
        mu, lv, u = [], [], []
        with torch.no_grad():
            for lo in range(0, len(x), 256):
                a, bb = model.enc_blood(torch.as_tensor(x[lo:lo + 256], device=DEV))
                mu.append(a.cpu().numpy()); lv.append(bb.cpu().numpy()); u.append(model.protein_head(a).cpu().numpy())
        mu = np.concatenate(mu); lv = np.concatenate(lv)
        return (mu + np.exp(.5 * lv) * np.random.default_rng(seed).normal(size=mu.shape)).astype(np.float32), np.concatenate(u)
    def predict(model, z, u, cond):
        zt = torch.as_tensor(z, device=DEV); ut = torch.as_tensor(u, device=DEV)
        cc = torch.as_tensor(np.repeat(cond[None], len(z), 0), device=DEV)
        with torch.enable_grad(): _, q = model.solve_equilibrium(zt, zt.clone(), ut, cc)
        qq = q.detach().cpu().numpy()
        with torch.no_grad():
            return np.concatenate([model.dec_brain(torch.as_tensor(qq[lo:lo + 256], device=DEV, dtype=torch.float32)).cpu().numpy() for lo in range(0, len(qq), 256)]).mean(0)
    prof = {}
    for seed in SEEDS:
        model = load(seed); rng = np.random.default_rng(seed)
        for bs in BS:
            si = states.index(bs)
            for c_blood in (0, 1):
                for d in np.unique(b["donor"][held & (b["state"] == si) & (b["cond2"] == c_blood)]):
                    ids = np.flatnonzero(held & (b["state"] == si) & (b["cond2"] == c_blood) & (b["donor"] == d))
                    if len(ids) < 25: continue
                    ids = rng.choice(ids, 128, replace=len(ids) < 128)
                    z, u = encode(model, b["X"][ids], seed=seed + 100 * c_blood)
                    for rs in RS:
                        for c_brain in (0, 1):
                            prof[(seed, bs, rs, c_brain, int(d))] = predict(model, z, u, cache[rs][f"cond{c_brain}"])
        del model; torch.cuda.empty_cache(); print("seed", seed, "done", flush=True)
    np.savez_compressed(CK, **{"|".join(map(str, k)): v for k, v in prof.items()})

ctrl_d = sorted({int(d) for d in np.unique(b["donor"][held & (b["cond2"] == 0)])})
pd_d = sorted({int(d) for d in np.unique(b["donor"][held & (b["cond2"] == 1)])})
ALL = ctrl_d + pd_d
print(f"\nheld-out blood donors: control {[dn[d] for d in ctrl_d]} | PD {[dn[d] for d in pd_d]}")
labelings = [set(c) for c in itertools.combinations(ALL, len(pd_d))]
print(f"exact donor-label assignments: {len(labelings)}")

def contrast(rs, c_brain, pdset, routes):
    agg = []
    for bs in routes:
        p1 = [prof[(s, bs, rs, c_brain, d)] for s in SEEDS for d in ALL if d in pdset and (s, bs, rs, c_brain, d) in prof]
        p0 = [prof[(s, bs, rs, c_brain, d)] for s in SEEDS for d in ALL if d not in pdset and (s, bs, rs, c_brain, d) in prof]
        if not p1 or not p0: continue
        agg.append(np.array(p1).mean(0) - np.array(p0).mean(0))
    return np.mean(agg, 0) if agg else None

rows = []
for routes, rlabel in [(BS, "all four routes"), (["cDC", "classical_mono", "nonclassical_mono"], "three routes with 3 vs 3 donors")]:
    print(f"\n=== {rlabel} | {len(SEEDS)} seeds ===")
    for c_brain, lbl in [(0, "fixed control"), (1, "fixed PD")]:
        pooled_obs, pooled_pred = [], []
        for rs in RS:
            d_true = contrast(rs, c_brain, set(pd_d), routes)
            r = stats.spearmanr(d_true, obs[rs]).correlation
            null = np.array([stats.spearmanr(contrast(rs, c_brain, s, routes), obs[rs]).correlation for s in labelings])
            p = float(np.mean(null >= r - 1e-12))
            print(f"  {lbl:14s} {rs:16s} rho {r:+.3f} | exact P {p:.3f} | null mean {null.mean():+.3f} sd {null.std():.3f}")
            rows.append(dict(routes=rlabel, readout=lbl, receiver=rs, rho=r, exact_p=p, n_seeds=len(SEEDS)))
            pooled_pred.append(d_true); pooled_obs.append(obs[rs])
        rp = stats.spearmanr(np.concatenate(pooled_pred), np.concatenate(pooled_obs)).correlation
        nullp = np.array([stats.spearmanr(np.concatenate([contrast(rs, c_brain, s, routes) for rs in RS]), np.concatenate(pooled_obs)).correlation for s in labelings])
        pp = float(np.mean(nullp >= rp - 1e-12))
        print(f"  {lbl:14s} {'BOTH RECEIVERS':16s} rho {rp:+.3f} | exact P {pp:.3f}")
        rows.append(dict(routes=rlabel, readout=lbl, receiver="pooled", rho=rp, exact_p=pp, n_seeds=len(SEEDS)))
pd.DataFrame(rows).to_csv(OUT + "/bloodonly_exact.csv", index=False)
print("\nwrote", OUT + "/bloodonly_exact.csv")
