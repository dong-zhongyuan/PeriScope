"""Protein-bridge diagnostics.
(1) Within-subtype prediction accuracy vs measured dynamic range for every biological antigen in the held-out
    donors (P6-P8), from the released per-antigen rho table and the frozen CITE-seq arrays.
(2) Profile-level accuracy: predicted vs observed ADT profiles per held-out donor x subtype across all antigens,
    and per-antigen ordering of subtypes, using the frozen seed-42 checkpoint (no retraining).
usage: python protein_bridge_diag.py OUTDIR
"""
import sys, json, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; I = A + "/interim/rescreen_20261002"; O = A + "/results/optimization_20261003"
F = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/figure_source_data"; OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
cg = np.load(I + "/citeseq.npz")
Y = cg["Y"]; counts = cg["adt_counts"]; is_test = cg["is_test"].astype(bool); subtype = cg["subtype"].astype(str); donor = cg["donor"]; dnames = cg["donor_names"].astype(str); proteins = cg["proteins"].astype(str)
registry = json.load(open(O + "/candidate_registry.json"))
print("cells:", Y.shape, "| test cells:", int(is_test.sum()), "| channels:", len(proteins), "| registry antigens:", len(registry))
pa = pd.read_csv(F + "/protein_accuracy.csv")
# ---- (1) dynamic range per antigen x subtype in held-out cells
rows = []
for r in registry:
    ch = r["channels"]
    for st in ["CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC"]:
        m = is_test & (subtype == st)
        if m.sum() < 25: continue
        y = Y[m][:, ch].mean(1); c = counts[m][:, ch].mean(1)
        rr = pa[(pa.target == r["target"]) & (pa.citeseq_subtype == st)].rho.mean()
        rows.append(dict(target=r["target"], subtype=st, is_control=r["is_control"], n=int(m.sum()), clr_sd=float(y.std()), clr_iqr=float(np.subtract(*np.percentile(y, [75, 25]))), mean_count=float(c.mean()), frac_detected=float((c > 0).mean()), rho=rr))
d = pd.DataFrame(rows); d = d[~d.is_control]; d.to_csv(OUT + "/protein_dynamic_range.csv", index=False)
print("\n=== within-subtype rho vs measured dynamic range (biological antigens, held-out cells) ===")
for st, g in d.groupby("subtype"):
    g = g.dropna(subset=["rho"]); print(f"{st:10s} n_antigens={len(g):3d} | Spearman(rho, CLR SD)={stats.spearmanr(g.rho, g.clr_sd).correlation:+.2f} | Spearman(rho, mean count)={stats.spearmanr(g.rho, g.mean_count).correlation:+.2f} | median rho={g.rho.median():.2f}")
for t in ["CD22", "CD40", "CD71"]:
    print(f"\n{t}:"); print(d[d.target == t][["subtype", "n", "clr_sd", "clr_iqr", "mean_count", "frac_detected", "rho"]].round(3).to_string(index=False))
    for st in ["cDC2", "CD14 Mono", "pDC"]:
        g = d[d.subtype == st].dropna(subset=["rho"]); x = d[(d.target == t) & (d.subtype == st)]
        if len(x): print(f"   {st}: CLR SD percentile among antigens = {int((g.clr_sd < x.clr_sd.iloc[0]).mean() * 100)}% | rho percentile = {int((g.rho < x.rho.iloc[0]).mean() * 100)}%")
# ---- (2) profile-level accuracy with the frozen seed-42 model
sys.path.insert(0, "/public/home/mengxl/dzy/pd_product/_work/optimization_20261003")
import torch
from ccwm import CCWM, CCWMConfig
from model_registry import checkpoint
dev = "cuda" if torch.cuda.is_available() else "cpu"; torch.set_num_threads(8)
X = cg["X"]; ci = np.where(is_test)[0]
bio = [r for r in registry if not r["is_control"]]
def predict(seed):
    ck = torch.load(checkpoint(O, seed), map_location="cpu", weights_only=False); m = CCWM(CCWMConfig(**ck["cfg"])); m.load_state_dict(ck["state_dict"]); m.eval().to(dev); yp = []
    with torch.no_grad():
        for s in range(0, len(ci), 1024):
            mu, _ = m.enc_blood(torch.from_numpy(X[ci[s:s + 1024]]).to(dev)); yp.append(m.protein_head(mu).cpu().numpy())
    return np.concatenate(yp)
ph = predict(42); Yt = Y[ci]; st_t = subtype[ci]; dn_t = dnames[donor[ci]]
print("\n=== profile-level accuracy (seed 42): per donor x subtype, Pearson r across", len(bio), "biological antigens ===")
rows = []
for dnm in sorted(set(dn_t)):
    for st in ["CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC"]:
        m = (dn_t == dnm) & (st_t == st)
        if m.sum() < 25: continue
        po = np.array([ph[m][:, r["channels"]].mean() for r in bio]); oo = np.array([Yt[m][:, r["channels"]].mean() for r in bio])
        rows.append(dict(donor=dnm, subtype=st, n=int(m.sum()), pearson=float(np.corrcoef(po, oo)[0, 1]), spearman=float(stats.spearmanr(po, oo).correlation)))
pr = pd.DataFrame(rows); pr.to_csv(OUT + "/protein_profile_accuracy.csv", index=False); print(pr.round(3).to_string(index=False))
print("profile r: median %.3f, range %.3f-%.3f" % (pr.pearson.median(), pr.pearson.min(), pr.pearson.max()))
# per-antigen ordering of subtypes (5 subtypes x 3 donors = up to 15 groups)
rows = []
for r in bio:
    po, oo = [], []
    for dnm in sorted(set(dn_t)):
        for st in ["CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC"]:
            m = (dn_t == dnm) & (st_t == st)
            if m.sum() < 25: continue
            po.append(ph[m][:, r["channels"]].mean()); oo.append(Yt[m][:, r["channels"]].mean())
    rows.append(dict(target=r["target"], groups=len(po), spearman_across_groups=float(stats.spearmanr(po, oo).correlation) if len(po) > 3 else np.nan))
ag = pd.DataFrame(rows); ag.to_csv(OUT + "/protein_across_group_accuracy.csv", index=False)
print("\nper-antigen Spearman across donor x subtype groups: median %.3f | fraction > 0.5: %.2f | CD22 %.2f CD40 %.2f CD71 %.2f" % (ag.spearman_across_groups.median(), (ag.spearman_across_groups > 0.5).mean(), *[ag[ag.target == t].spearman_across_groups.iloc[0] for t in ["CD22", "CD40", "CD71"]]))
# variance decomposition of observed CLR in held-out cells: between-subtype vs within
tot = np.array([Yt[:, r["channels"]].mean(1).var() for r in bio]); within = np.array([np.mean([Yt[st_t == st][:, r["channels"]].mean(1).var() for st in ["CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC"] if (st_t == st).sum() >= 25]) for r in bio])
print("observed ADT variance within subtypes as a fraction of total (median across antigens): %.2f" % np.median(within / tot))
