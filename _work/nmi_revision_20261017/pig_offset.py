"""Pig transfer, third pass: is the dominant between-animal axis a per-sample global offset in the
log-ratio response, and does removing it (label-free, applied identically to prediction and observation)
recover the group-level agreement? Exact 210-labeling test throughout.
usage: python pig_offset.py OUTDIR
"""
import sys, itertools, warnings, json
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
P = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/pig_external_validation_20261003"
OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)

obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0)
animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals])
O = obs.values
methods = ["periscope", "cellot", "conditional_icnn", "conditional_mlp", "ot_ridge", "shared_label_mean", "shared_label_marginal"]
preds = {}
for m in methods:
    try:
        preds[m] = pd.read_csv(P + f"/predicted_response_{m}.csv", index_col=0).reindex(obs.index)[animals].values
    except Exception:
        pass
Pm = preds["periscope"]
labelings = [np.isin(np.arange(10), list(s)) for s in itertools.combinations(range(10), 4)]

# ---- (1) is PC1 the per-animal global offset?
Oc = O - O.mean(1, keepdims=True); U, S, Vt = np.linalg.svd(Oc, full_matrices=False)
ve = S ** 2 / (S ** 2).sum(); pc1_scores = Vt[0] * S[0]; off = O.mean(0)
print("=== PC1 identity check ===")
print(f"variance explained PC1 {ve[0]:.3f} PC2 {ve[1]:.3f} PC3 {ve[2]:.3f}")
print(f"Pearson(PC1 animal score, animal mean log2 response) = {np.corrcoef(pc1_scores, off)[0,1]:+.4f}")
print(f"PC1 gene loading: mean {(U[:,0]*S[0]).mean():+.3f}, sd {(U[:,0]*S[0]).std():.3f}, fraction same sign {max((U[:,0]>0).mean(),(U[:,0]<0).mean()):.3f}")
print("animal mean log2 response:", dict(zip(animals, np.round(off, 3))))
print(f"sd of animal offsets {off.std():.3f} | sd of observed group effect {(O[:,lps].mean(1)-O[:,~lps].mean(1)).std():.3f}")
print(f"control-group mean offset {off[~lps].mean():+.3f} | LPS-group mean offset {off[lps].mean():+.3f}")
print(f"fraction of total response variance that is the per-animal constant: {np.var(np.tile(off,(O.shape[0],1)))/np.var(O):.3f}")

# ---- (2) label-free per-animal centering, identical on both sides
def center(M, how):
    if how == "none": return M
    if how == "mean": return M - M.mean(0, keepdims=True)
    if how == "median": return M - np.median(M, axis=0, keepdims=True)
    if how == "trimmed": return M - stats.trim_mean(M, 0.25, axis=0)[None, :]
    raise ValueError(how)

def exact_test(pred, obsm):
    pe = pred[:, lps].mean(1) - pred[:, ~lps].mean(1)
    rho = stats.spearmanr(pe, obsm[:, lps].mean(1) - obsm[:, ~lps].mean(1)).correlation
    null = np.array([stats.spearmanr(pe, obsm[:, m].mean(1) - obsm[:, ~m].mean(1)).correlation for m in labelings])
    return rho, float(np.mean(null >= rho - 1e-12)), null

print("\n=== primary endpoint under label-free per-animal centering (prediction and observation treated identically) ===")
rows = []
for how in ["none", "mean", "median", "trimmed"]:
    Oh = center(O, how); Ph = center(Pm, how)
    rho, p, null = exact_test(Ph, Oh)
    pear = np.corrcoef(Ph[:, lps].mean(1) - Ph[:, ~lps].mean(1), Oh[:, lps].mean(1) - Oh[:, ~lps].mean(1))[0, 1]
    sgn = np.mean(np.sign(Ph[:, lps].mean(1) - Ph[:, ~lps].mean(1)) == np.sign(Oh[:, lps].mean(1) - Oh[:, ~lps].mean(1)))
    print(f"  centering={how:8s} rho={rho:+.3f} Pearson={pear:+.3f} sign agreement={sgn:.3f} exact P={p:.4f} (null mean {null.mean():+.3f}, sd {null.std():.3f}, 95th {np.quantile(null,.95):+.3f})")
    rows.append(dict(centering=how, spearman=rho, pearson=pear, sign_agreement=sgn, exact_p=p, null_mean=null.mean(), null_sd=null.std()))

# centering only the observation, or only the prediction (asymmetric; diagnostic only)
for lbl, Ph, Oh in [("observation only", Pm, center(O, "median")), ("prediction only", center(Pm, "median"), O)]:
    rho, p, _ = exact_test(Ph, Oh); print(f"  [diagnostic] {lbl:17s} rho={rho:+.3f} exact P={p:.4f}")

# ---- (3) all methods under the same centering
print("\n=== all frozen comparators, median-centered (same label-free operation for every method) ===")
mrows = []
for m, Pv in preds.items():
    rho0, p0, _ = exact_test(Pv, O); rho1, p1, _ = exact_test(center(Pv, "median"), center(O, "median"))
    print(f"  {m:24s} original rho {rho0:+.3f} (P {p0:.4f})   centered rho {rho1:+.3f} (P {p1:.4f})")
    mrows.append(dict(method=m, rho_original=rho0, p_original=p0, rho_centered=rho1, p_centered=p1))

# ---- (4) ten-seed sensitivity under centering
z = np.load(P + "/frozen_ten_seed_predictions.npz", allow_pickle=True)
print("\n=== frozen ten-seed file keys ===", list(z.keys())[:10])
np.save(OUT + "/pig_offset_U.npy", U[:, :3])
pd.DataFrame(rows).to_csv(OUT + "/pig_centering.csv", index=False)
pd.DataFrame(mrows).to_csv(OUT + "/pig_centering_methods.csv", index=False)
print("\nwrote", OUT + "/pig_centering.csv")
