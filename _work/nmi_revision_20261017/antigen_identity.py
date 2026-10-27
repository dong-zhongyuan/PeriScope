"""Is antigen identity reproducibly encoded in the predicted brain response at all?
Independent of any program-score statistic: for each axis, take each antigen's predicted gene-response profile
in one confirmation model and ask whether it matches the SAME antigen in a different confirmation model better
than it matches the other 212 antigens. Exact permutation of the antigen labels gives the null.
usage: python antigen_identity.py OUTDIR
"""
import sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; O = A + "/results/optimization_20261003"; OUT = sys.argv[1]
pd.set_option("display.width", 250)
rng = np.random.default_rng(7)
rows, det = [], []
for f in sorted(Path(O + "/curves").glob("*__slopes.npz")):
    axis = f.name.replace("__slopes.npz", ""); z = np.load(f)
    D = z["D"]; seeds = list(z["seeds"]); bg = np.where(~z["is_control"])[0]; ps = z["proteins"].astype(str)
    conf = [i for i, s in enumerate(seeds) if s >= 47]
    if len(conf) < 2: continue
    half1 = D[conf[:len(conf) // 2]].mean(0)[bg]      # antigens x genes
    half2 = D[conf[len(conf) // 2:]].mean(0)[bg]
    # centre each half across antigens so a shared response axis cannot create the match
    h1 = half1 - half1.mean(0); h2 = half2 - half2.mean(0)
    h1 = (h1 - h1.mean(1, keepdims=True)) / (h1.std(1, keepdims=True) + 1e-12)
    h2 = (h2 - h2.mean(1, keepdims=True)) / (h2.std(1, keepdims=True) + 1e-12)
    M = h1 @ h2.T / h1.shape[1]
    n = M.shape[0]
    diag = np.diag(M); off = M[~np.eye(n, dtype=bool)]
    top1 = float(np.mean(M.argmax(1) == np.arange(n)))
    null = np.array([np.mean(M[np.arange(n), rng.permutation(n)]) for _ in range(20000)])
    p = float((1 + (null >= diag.mean()).sum()) / (1 + len(null)))
    rows.append(dict(axis=axis, n_antigens=n, same_antigen_r=float(diag.mean()), different_antigen_r=float(off.mean()),
                     top1_accuracy=top1, chance=1.0 / n, perm_p=p))
    print(f"  {axis:36s} n={n:3d} | same-antigen r {diag.mean():+.3f} vs different {off.mean():+.3f} | top-1 match {100*top1:.1f}% (chance {100/n:.1f}%) | perm P {p:.2e}")
    for i, nm in enumerate(ps[bg]):
        if nm in ("CD22", "CD40", "CD71"):
            det.append(dict(axis=axis, antigen=nm, same_r=float(M[i, i]), rank=int((M[i] >= M[i, i]).sum()), n=n))
t = pd.DataFrame(rows)
print(f"\nacross the eight axes: same-antigen r {t.same_antigen_r.min():.3f}-{t.same_antigen_r.max():.3f}, "
      f"top-1 identification {100*t.top1_accuracy.min():.0f}-{100*t.top1_accuracy.max():.0f}% against chance {100*t.chance.mean():.1f}%")
print("\n=== the three nominated antigens, self-match rank among 213 ===")
d = pd.DataFrame(det)
print(d.to_string(index=False))
print(f"  self-match is the top match in {int((d['rank']==1).sum())}/{len(d)} antigen-axis pairs")
t.to_csv(OUT + "/antigen_identity.csv", index=False); d.to_csv(OUT + "/antigen_identity_targets.csv", index=False)
print("\nwrote", OUT + "/antigen_identity.csv")
