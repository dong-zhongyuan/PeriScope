"""Donor-resolved blood-only readout, using the frozen checkpoints and no retraining.

The released blood-only readout contrasts ONE PD blood donor with ONE control blood donor
('test_blood_per_condition': 1 in the benchmark manifest), so disease and donor identity are not separable
in it. This recomputes the same quantity per blood donor for every held-out blood donor available
(validation: 2 control / 2 PD; test: 1 / 1) and asks:
  (1) how large is between-donor variation of the predicted brain profile within a condition, relative to the
      between-condition difference;
  (2) does the blood-only direction correlation change when the contrast uses 3 vs 3 donors instead of 1 vs 1.
usage: python bloodonly_donors.py OUTDIR [seeds]
"""
import sys, os, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
sys.path.insert(0, "/public/home/mengxl/dzy/pd_product/_work/benchmark_shared_20261003")
sys.path.insert(0, "/public/home/mengxl/dzy/pd_product/_work/optimization_20261003")
import torch
from ccwm import CCWM, CCWMConfig

A = "/public/home/mengxl/dzy/pd_product_assets"
OUT = sys.argv[1]; SEEDS = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2 else ["42", "43", "44"])]
DEV = "cuda" if torch.cuda.is_available() else "cpu"
pd.set_option("display.width", 250)
BS = ["cDC", "classical_mono", "nonclassical_mono", "pDC"]; RS = ["astro", "microglia_mhc2"]

b = dict(np.load(A + "/interim/rescreen_20261002/blood.npz"))
states = list(b["state_vocab"].astype(str)); dn = b["donor_names"].astype(str)
held = b["is_val"] | b["is_test"]

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
    z = (mu + np.exp(.5 * lv) * np.random.default_rng(seed).normal(size=mu.shape)).astype(np.float32)
    return z, np.concatenate(u)

def decode(model, z):
    with torch.no_grad():
        return np.concatenate([model.dec_brain(torch.as_tensor(z[lo:lo + 256], device=DEV, dtype=torch.float32)).cpu().numpy() for lo in range(0, len(z), 256)])

def predict(model, z, u, cond):
    zt = torch.as_tensor(z, device=DEV); ut = torch.as_tensor(u, device=DEV)
    cc = torch.as_tensor(np.repeat(cond[None], len(z), 0), device=DEV)
    with torch.enable_grad():
        _, q = model.solve_equilibrium(zt, zt.clone(), ut, cc)
    return decode(model, q.detach().cpu().numpy()).mean(0)

cache = {rc: np.load(A + f"/results/benchmark_20261002/cache/cDC__{rc}.npz", allow_pickle=True) for rc in RS}
obs = {rc: cache[rc]["profiles1"].mean(0) - cache[rc]["profiles0"].mean(0) for rc in RS}

rows = []
prof = {}   # (seed, bs, rs, brain_cond, donor) -> predicted brain profile
for seed in SEEDS:
    model = load(seed); rng = np.random.default_rng(seed)
    for bs in BS:
        si = states.index(bs)
        donors = []
        for c in (0, 1):
            for d in np.unique(b["donor"][held & (b["state"] == si) & (b["cond2"] == c)]):
                ids = np.flatnonzero(held & (b["state"] == si) & (b["cond2"] == c) & (b["donor"] == d))
                if len(ids) < 25: continue
                donors.append((int(d), c, rng.choice(ids, 128, replace=len(ids) < 128)))
        for d, c_blood, ids in donors:
            z, u = encode(model, b["X"][ids], seed=seed + 100 * c_blood)
            for rs in RS:
                for c_brain in (0, 1):
                    p = predict(model, z, u, cache[rs][f"cond{c_brain}"])
                    prof[(seed, bs, rs, c_brain, d)] = p
                    rows.append(dict(seed=seed, blood_state=bs, brain_state=rs, brain_condition=c_brain,
                                     donor=dn[d], blood_condition=c_blood, split="test" if b["is_test"][ids[0]] else "val"))
    del model; torch.cuda.empty_cache()
    print("seed", seed, "done", flush=True)

meta = pd.DataFrame(rows).drop_duplicates(); meta.to_csv(OUT + "/bloodonly_donor_index.csv", index=False)
print("\n=== held-out blood donors per state ===")
print(meta[meta.seed == SEEDS[0]].groupby(["blood_state", "blood_condition"]).donor.nunique().unstack().to_string())

print("\n=== (1) between-donor variation within a condition vs the between-condition difference ===")
print("    (predicted brain profiles, brain condition held at control)")
for rs in RS:
    for bs in BS:
        ctrl = [prof[(s, bs, rs, 0, d)] for s in SEEDS for d in {int(x) for x in np.unique(b["donor"][held & (b["state"] == states.index(bs)) & (b["cond2"] == 0)])} if (s, bs, rs, 0, d) in prof]
        pdd = [prof[(s, bs, rs, 0, d)] for s in SEEDS for d in {int(x) for x in np.unique(b["donor"][held & (b["state"] == states.index(bs)) & (b["cond2"] == 1)])} if (s, bs, rs, 0, d) in prof]
        if len(ctrl) < 2 or len(pdd) < 2: print(f"  {rs:16s} {bs:20s} only {len(ctrl)}/{len(pdd)} donor-seed profiles, skipped"); continue
        ctrl = np.array(ctrl); pdd = np.array(pdd)
        between = np.sqrt(((pdd.mean(0) - ctrl.mean(0)) ** 2).mean())
        within = np.sqrt(np.mean([((x - ctrl.mean(0)) ** 2).mean() for x in ctrl] + [((x - pdd.mean(0)) ** 2).mean() for x in pdd]))
        rho_dd = np.mean([stats.spearmanr(ctrl[i] - ctrl.mean(0), ctrl[j] - ctrl.mean(0)).correlation for i in range(len(ctrl)) for j in range(i + 1, len(ctrl))])
        print(f"  {rs:16s} {bs:20s} between-condition RMS {between:.4f} | within-condition donor RMS {within:.4f} | ratio {between/within:.2f}")

print("\n=== (2) blood-only direction correlation: 1 vs 1 (released) against all held-out donors ===")
out = []
for rs in RS:
    for c_brain, lbl in [(0, "fixed_control"), (1, "fixed_PD")]:
        for scope in ["test only (1 vs 1)", "all held-out donors"]:
            agg = []
            for bs in BS:
                si = states.index(bs)
                sel = lambda c: [int(d) for d in np.unique(b["donor"][held & (b["state"] == si) & (b["cond2"] == c)])
                                 if (scope == "all held-out donors" or b["is_test"][np.flatnonzero((b["donor"] == d) & (b["state"] == si) & (b["cond2"] == c))[0]])]
                c0 = [prof[(s, bs, rs, c_brain, d)] for s in SEEDS for d in sel(0) if (s, bs, rs, c_brain, d) in prof]
                c1 = [prof[(s, bs, rs, c_brain, d)] for s in SEEDS for d in sel(1) if (s, bs, rs, c_brain, d) in prof]
                if not c0 or not c1: continue
                agg.append(np.array(c1).mean(0) - np.array(c0).mean(0))
            if not agg: continue
            d_agg = np.mean(agg, 0)
            r = stats.spearmanr(d_agg, obs[rs]).correlation
            print(f"  {rs:16s} {lbl:14s} {scope:22s} rho {r:+.3f} | contrast RMS {np.sqrt((d_agg**2).mean()):.4f} (observed {np.sqrt((obs[rs]**2).mean()):.4f})")
            out.append(dict(receiver=rs, readout=lbl, scope=scope, rho=r, contrast_rms=float(np.sqrt((d_agg ** 2).mean())), n_seeds=len(SEEDS)))
pd.DataFrame(out).to_csv(OUT + "/bloodonly_donor_resolved.csv", index=False)
print("\nwrote", OUT + "/bloodonly_donor_resolved.csv")
