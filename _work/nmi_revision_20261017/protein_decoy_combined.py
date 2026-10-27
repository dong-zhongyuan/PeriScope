"""Higher-powered protein specificity: combine a protein's programs. For every antigen, average its panel score
across exactly the axes where the nominated protein has priority programs, then ask where the nominated protein
ranks among the 212 antigens on that combined score. Same data, same null, one test per protein.
usage: python protein_decoy_combined.py OUTDIR
"""
import sys, json, warnings
from pathlib import Path
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; O = A + "/results/optimization_20261003"; OUT = sys.argv[1]
pd.set_option("display.width", 250)
programs = json.loads(open(O + "/programs.json").read())
pri = pd.read_csv("/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/figure_source_data/priority.csv")
want = {c.split("::")[-1] for c in pri.combo}

panels = {}   # combo -> (protein_names, mz vector over all proteins)
for f in sorted(Path(O + "/curves").glob("*__slopes.npz")):
    axis = f.name.replace("__slopes.npz", ""); z = np.load(f)
    D = z["D"]; genes = z["genes"].astype(str).tolist(); ps = z["proteins"].astype(str).tolist()
    background = np.where(~z["is_control"])[0]; gidx = {g: i for i, g in enumerate(genes)}
    mu = D[:, background].mean(1, keepdims=True); sd = D[:, background].std(1, keepdims=True); Z = (D - mu) / (sd + 1e-9)
    conf = [i for i, s in enumerate(z["seeds"]) if s >= 47]
    for combo, gs in programs.items():
        if combo not in want: continue
        pname, rest = combo.split("__", 1); paxis, direction = rest.rsplit("__", 1)
        if paxis != axis: continue
        si = np.array([gidx[g] for g in gs]); sgn = 1 if direction == "up" else -1
        vset = sgn * Z[:, :, si].mean(2); Iscore = vset - sgn * Z.mean(2)
        panel = (Iscore - Iscore[:, background].mean(1, keepdims=True)) / (Iscore[:, background].std(1, keepdims=True) + 1e-12)
        panels[combo] = (ps, panel[conf].mean(0), background, pname)

print(f"panels built for {len(panels)} priority programs")
rows = []
for target in ["CD22", "CD40", "CD71"]:
    combos = [c for c in panels if panels[c][3] == target]
    for scope, keep in [("all programs", combos),
                        ("MG-AP programs", [c for c in combos if "microglia_mhc2" in c]),
                        ("Astro programs", [c for c in combos if "astro" in c])]:
        if len(keep) < 2: continue
        ps0 = panels[keep[0]][0]; bg = panels[keep[0]][2]
        common = [p for p in ps0 if all(p in panels[c][0] for c in keep)]
        M = np.array([[panels[c][1][panels[c][0].index(p)] for p in common] for c in keep])
        comb = M.mean(0)
        names = np.array(common)
        bg_mask = np.array([i for i, p in enumerate(common) if p in {ps0[j] for j in bg}])
        own = float(comb[list(names).index(target)])
        pool = np.array([comb[i] for i in bg_mask if names[i] != target])
        p = (1 + (pool >= own).sum()) / (1 + len(pool))
        rank = int((pool >= own).sum()) + 1
        print(f"  {target:5s} {scope:16s} n_programs={len(keep)} | combined score {own:+.3f} | rank {rank}/{len(pool)+1} | empirical P {p:.4f}")
        rows.append(dict(target=target, scope=scope, n_programs=len(keep), combined_score=own, rank=rank, n_antigens=len(pool) + 1, p_empirical=p))
t = pd.DataFrame(rows); t.to_csv(OUT + "/protein_decoy_combined.csv", index=False)
print(f"\n  combined tests reaching P<=0.05: {int((t.p_empirical<=0.05).sum())}/{len(t)}")
print("wrote", OUT + "/protein_decoy_combined.csv")
