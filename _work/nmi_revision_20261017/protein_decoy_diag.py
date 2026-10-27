"""Why is the own-protein selectivity of the priority programs not an outlier among the 212 other antigens?
(a) how correlated are different antigens' program scores; (b) does measured dynamic range drive the score;
(c) a decoy pool matched for dynamic range in the same source state, with the empirical p recomputed.
usage: python protein_decoy_diag.py OUTDIR
"""
import sys, json, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; O = A + "/results/optimization_20261003"; OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
QR = np.arange(5, dtype=float)
programs = json.loads(open(O + "/programs.json").read())
pri = pd.read_csv("/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/figure_source_data/priority.csv")
cg = np.load(A + "/interim/rescreen_20261002/citeseq.npz")
Y = cg["Y"]; proteins_all = cg["proteins"].astype(str)
registry = json.load(open(O + "/candidate_registry.json"))
rng_ch = {r["target"]: r["channels"] for r in registry}
spread = {r["target"]: float(np.std(Y[:, r["channels"]].mean(1))) for r in registry}

from pathlib import Path
rows, corr_rows = [], []
for f in sorted(Path(O + "/curves").glob("*__slopes.npz")):
    axis = f.name.replace("__slopes.npz", ""); z = np.load(f)
    D = z["D"]; genes = z["genes"].astype(str).tolist(); ps = z["proteins"].astype(str).tolist()
    background = np.where(~z["is_control"])[0]; gidx = {g: i for i, g in enumerate(genes)}
    mu = D[:, background].mean(1, keepdims=True); sd = D[:, background].std(1, keepdims=True); Z = (D - mu) / (sd + 1e-9)
    conf = [i for i, s in enumerate(z["seeds"]) if s >= 47]
    for combo, gs in programs.items():
        pname, rest = combo.split("__", 1); paxis, direction = rest.rsplit("__", 1)
        if paxis != axis: continue
        full = ("raw::" if "::" not in combo else "") + combo
        match = pri[pri.combo.str.endswith(combo)]
        if not len(match): continue
        ri = ps.index(pname); si = np.array([gidx[g] for g in gs]); sgn = 1 if direction == "up" else -1
        pool = np.array([i for i in background if i != ri])
        vset = sgn * Z[:, :, si].mean(2); Iscore = vset - sgn * Z.mean(2)
        panel = (Iscore - Iscore[:, background].mean(1, keepdims=True)) / (Iscore[:, background].std(1, keepdims=True) + 1e-12)
        mz = panel[conf].mean(0)
        # (a) correlation among antigens: how many decoys exceed, and how concentrated is the pool
        own = mz[ri]; pool_v = mz[pool]
        p_all = (1 + (pool_v >= own).sum()) / (1 + len(pool_v))
        # (c) decoys matched for measured dynamic range of the antigen in this source state
        sp_own = spread.get(pname, np.nan)
        sp_pool = np.array([spread.get(ps[i], np.nan) for i in pool])
        ok = np.isfinite(sp_pool) & np.isfinite(sp_own)
        if ok.sum() > 30:
            d = np.abs(np.log10(sp_pool[ok] + 1e-6) - np.log10(sp_own + 1e-6))
            sel = np.argsort(d)[:50]
            matched = pool_v[ok][sel]
            p_matched = (1 + (matched >= own).sum()) / (1 + len(matched))
        else:
            p_matched = np.nan
        rho_sp = stats.spearmanr(pool_v, sp_pool[ok] if ok.sum() == len(pool_v) else sp_pool)[0] if ok.sum() > 10 else np.nan
        rows.append(dict(combo=combo, axis=axis, protein=pname, own_selectivity=float(own), p_all_decoys=float(p_all),
                         p_range_matched=float(p_matched), own_dynamic_range=sp_own,
                         rank_of_own=int((pool_v >= own).sum()) + 1, n_decoys=len(pool_v),
                         decoy_sd=float(pool_v.std()), corr_selectivity_vs_range=float(rho_sp) if np.isfinite(rho_sp) else np.nan,
                         priority=bool(len(match))))
        # (a) pairwise correlation of antigen response profiles in this axis
        if combo == list(programs)[0] or len(corr_rows) < 8:
            Rm = np.corrcoef(sgn * Z[conf].mean(0)[background][:, si])
            iu = np.triu_indices_from(Rm, 1)
            corr_rows.append(dict(axis=axis, combo=combo, median_pairwise_r=float(np.median(Rm[iu])), frac_r_above_0p5=float(np.mean(Rm[iu] > .5))))
t = pd.DataFrame(rows)
t = t[t.combo.isin([c.split("::")[-1] for c in pri.combo])]
print("=== the 16 priority programs: own-protein selectivity against all decoys vs range-matched decoys ===")
print(t[["protein", "axis", "own_selectivity", "rank_of_own", "n_decoys", "p_all_decoys", "p_range_matched", "own_dynamic_range"]].round(4).to_string(index=False))
print(f"\n  significant at 0.05 against all decoys: {int((t.p_all_decoys<=0.05).sum())}/{len(t)}")
print(f"  significant at 0.05 against range-matched decoys: {int((t.p_range_matched<=0.05).sum())}/{len(t)}")
print("\n=== how correlated are different antigens' program scores (confirmation ensemble) ===")
print(pd.DataFrame(corr_rows).round(3).to_string(index=False))
print("\n=== does measured dynamic range drive the panel score? ===")
print(f"  median Spearman(decoy selectivity, decoy dynamic range) across programs: {t.corr_selectivity_vs_range.median():+.3f}")
print(f"  Spearman(own selectivity, own dynamic range) across the 16: {stats.spearmanr(t.own_selectivity, t.own_dynamic_range).correlation:+.3f}")
t.to_csv(OUT + "/protein_decoy_diag.csv", index=False); print("\nwrote", OUT + "/protein_decoy_diag.csv")
