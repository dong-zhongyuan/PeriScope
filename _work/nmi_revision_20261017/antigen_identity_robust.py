"""Robustness of the antigen-identity result, and the same test restricted to each program's own genes.
(1) average the self-match over all balanced splits of the ten models, not one arbitrary split;
(2) repeat on the gene set of each priority program, which is what the nomination actually uses.
usage: python antigen_identity_robust.py OUTDIR
"""
import sys, json, itertools, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; O = A + "/results/optimization_20261003"; OUT = sys.argv[1]
pd.set_option("display.width", 250)
programs = json.loads(open(O + "/programs.json").read())
pri = pd.read_csv("/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/figure_source_data/priority.csv")
want = {c.split("::")[-1] for c in pri.combo}

def self_match(D, idx, cols=None):
    """mean over all balanced splits of idx of the per-antigen self-correlation, after removing the
    across-antigen mean profile and standardising each antigen."""
    out = []
    half = len(idx) // 2
    for sel in itertools.combinations(range(len(idx)), half):
        a = [idx[i] for i in sel]; b = [idx[i] for i in range(len(idx)) if i not in sel]
        h1 = D[a].mean(0); h2 = D[b].mean(0)
        if cols is not None: h1 = h1[:, cols]; h2 = h2[:, cols]
        h1 = h1 - h1.mean(0); h2 = h2 - h2.mean(0)
        h1 = (h1 - h1.mean(1, keepdims=True)) / (h1.std(1, keepdims=True) + 1e-12)
        h2 = (h2 - h2.mean(1, keepdims=True)) / (h2.std(1, keepdims=True) + 1e-12)
        out.append((h1 * h2).mean(1))
    return np.mean(out, 0)

rows, prog_rows = [], []
for f in sorted(Path(O + "/curves").glob("*__slopes.npz")):
    axis = f.name.replace("__slopes.npz", ""); z = np.load(f)
    D = z["D"]; bg = np.where(~z["is_control"])[0]; ps = z["proteins"].astype(str); genes = z["genes"].astype(str).tolist()
    gidx = {g: i for i, g in enumerate(genes)}
    Db = D[:, bg]
    for lbl, idx in [("confirmation 47-51", [i for i, s in enumerate(z["seeds"]) if s >= 47]),
                     ("all ten models", list(range(len(z["seeds"]))))]:
        sm = self_match(Db, idx)
        names = ps[bg]
        for t in ["CD22", "CD40", "CD71"]:
            if t not in names: continue
            i = list(names).index(t)
            rows.append(dict(axis=axis, split=lbl, antigen=t, self_match=float(sm[i]),
                             rank=int((sm >= sm[i]).sum()), n=len(sm), panel_median=float(np.median(sm))))
    # restricted to each priority program's own genes
    for combo, gs in programs.items():
        if combo not in want: continue
        pname, rest = combo.split("__", 1); paxis, _ = rest.rsplit("__", 1)
        if paxis != axis: continue
        cols = np.array([gidx[g] for g in gs])
        sm = self_match(Db, list(range(len(z["seeds"]))), cols)
        names = list(ps[bg]); i = names.index(pname)
        prog_rows.append(dict(program=combo, antigen=pname, n_genes=len(gs), self_match=float(sm[i]),
                              rank=int((sm >= sm[i]).sum()), n=len(sm), panel_median=float(np.median(sm))))

t = pd.DataFrame(rows)
print("=== self-match averaged over all balanced splits ===")
for lbl in ["confirmation 47-51", "all ten models"]:
    g = t[t.split == lbl]
    print(f"\n  -- {lbl} (panel median across 213 antigens shown per axis)")
    print(g.pivot_table(index="axis", columns="antigen", values="self_match").round(3).to_string())
    print(f"     panel median range {g.panel_median.min():.3f}-{g.panel_median.max():.3f}")
    for a in ["CD22", "CD40", "CD71"]:
        s = g[g.antigen == a]
        print(f"     {a}: self-match {s.self_match.min():+.3f} to {s.self_match.max():+.3f} | rank {s['rank'].min()}-{s['rank'].max()} of 213 | above panel median in {int((s.self_match>s.panel_median).sum())}/{len(s)} axes")

p = pd.DataFrame(prog_rows)
print("\n=== restricted to each priority program's own gene set (all ten models) ===")
print(p.round(3).to_string(index=False))
print(f"\n  programs whose own antigen is above the panel median on its own genes: {int((p.self_match>p.panel_median).sum())}/{len(p)}")
print(f"  programs whose own antigen ranks in the top 10% on its own genes: {int((p['rank']<=0.1*p.n).sum())}/{len(p)}")
t.to_csv(OUT + "/antigen_identity_robust.csv", index=False); p.to_csv(OUT + "/antigen_identity_programs.csv", index=False)
print("\nwrote", OUT + "/antigen_identity_robust.csv")
