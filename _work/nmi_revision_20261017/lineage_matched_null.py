"""Does the PD decrease of the programs exceed what random gene sets with the SAME lineage composition and
expression level show? Null sets are matched gene-by-gene on (dominant lineage from the Ma reference, expression bin).
Runs for (a) spatial recipient-enriched spots (top-25% design) and (b) GSE157783 matched fine states.
usage: python lineage_matched_null.py DIAGDIR
"""
import json, sys, glob, warnings
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/optimization_20261003"
D = sys.argv[1]; rng = np.random.RandomState(7); B = 4000
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
programs = json.load(open(R + "/joint_analysis/programs.json")); pr = pd.read_csv(R + "/mr_parallel_gate_20261003/priority_all_programs.csv")
six = ["raw::CD294__classical_mono__x__astro__up", "raw::CD294__pDC__x__astro__up", "raw::CD30__classical_mono__x__astro__up",
       "raw::CD80__nonclassical_mono__x__astro__up", "raw::CD80__pDC__x__astro__up", "relative::CD274__pDC__x__astro__up"]
rp = pd.read_csv(A + "/results/optimization_20261003/spatial_subtype_author_20261003/inputs/reference_profiles.csv.gz", index_col=0)
lin = {"astro_fibrous": "astro", "astro_protoplasmic": "astro", "microglia_homeostatic": "microglia", "microglia_mhc2": "microglia", "neuron_GABAergic": "neuron", "neuron_MEG3": "neuron", "neuron_dopaminergic": "neuron", "Oligo": "oligo", "OPC": "opc", "Endothelial": "endo", "BAM": "other", "Tcells": "other"}
share = rp.T.groupby(lambda c: lin[c]).sum().T; share = share.div(share.sum(1) + 1e-12, axis=0); dom = share.idxmax(1)
def collapse_matrix(X, names):
    uniq = sorted(set(names)); gidx = {g: i for i, g in enumerate(uniq)}
    C = sparse.csr_matrix((np.ones(len(names)), (np.arange(len(names)), [gidx[g] for g in names])), shape=(len(names), len(uniq)))
    return X @ C, uniq
def welch_from(M, is_pd):
    a_, b_ = M[is_pd], M[~is_pd]
    t = (a_.mean(0) - b_.mean(0)) / np.sqrt(a_.var(0, ddof=1) / len(a_) + b_.var(0, ddof=1) / len(b_) + 1e-12)
    t[(M > 0).sum(0) < 3] = np.nan; return t
# ---------------- (a) spatial recipient-enriched spots
sp = ad.read_h5ad(A + "/processed/gse253975/v0.1/GSE253975_geomx.h5ad"); Xs, gs_sp = collapse_matrix(sp.X, np.array(sp.var_names.astype(str)))
obs = sp.obs.copy(); obs["donor"] = obs.donor.astype(str); obs["condition"] = obs.condition.astype(str)
fr = pd.concat([pd.read_csv(f, index_col=0) for f in sorted(glob.glob(R + "/spatial_subtype_author_20261003/samples/GSM*/cell_state_RNA_fractions.csv"))]).reindex(obs.index)
fr["astro_total"] = fr.astro_fibrous + fr.astro_protoplasmic
donors = sorted(set(obs.donor)); is_pd = np.array([obs.condition[obs.donor == d].iloc[0] == "PD" for d in donors])
def spatial_design(col, q=0.75):
    M = []
    for d in donors:
        m = (obs.donor == d).values; thr = np.nanquantile(fr.loc[m, col], q); pick = m & (fr[col].values >= thr)
        s = np.asarray(Xs[pick].sum(0)).ravel(); M.append(np.log1p(s / max(s.sum(), 1) * 1e4))
    M = np.stack(M); return pd.Series(welch_from(M, is_pd), index=gs_sp), pd.Series(M.mean(0), index=gs_sp)
DES = {"spatial_astro_top25": spatial_design("astro_total"), "spatial_mhc2_top25": spatial_design("microglia_mhc2")}
# ---------------- (b) GSE157783 fine states
a = ad.read_h5ad(A + "/processed/gse157783/v0.1/gse157783_qc.h5ad"); Xg, gs_g = collapse_matrix(a.X, a.var["gene_symbol"].astype(str).values)
ob = a.obs.copy()
for c in ["donor", "condition"]: ob[c] = ob[c].astype(str)
fine = pd.Series("other", index=ob.index)
for fam in ["brain_astro", "brain_microglia"]:
    m = pd.read_csv(A + f"/interim/v0.1/purification/c157_{fam}.csv", index_col=0); m = m[m.kept]; common = fine.index.intersection(m.index); fine.loc[common] = m.loc[common, "state_pure"].astype(str)
gd = sorted(set(ob.donor)); g_pd = np.array([ob.condition[ob.donor == d].iloc[0] == "PD" for d in gd])
def sn_design(state):
    M = []
    for d in gd:
        mm = ((fine == state) & (ob.donor == d)).values; s = np.asarray(Xg[mm].sum(0)).ravel(); M.append(np.log1p(s / max(s.sum(), 1) * 1e4))
    M = np.stack(M); return pd.Series(welch_from(M, g_pd), index=gs_g), pd.Series(M.mean(0), index=gs_g)
DES["snRNA_astro"] = sn_design("astro"); DES["snRNA_mhc2"] = sn_design("microglia_mhc2")
# ---------------- tests
def test(design, gset):
    t, e = design; t = t.dropna(); e = e.reindex(t.index); bins = pd.qcut(e.rank(method="first"), 20, labels=False).values
    lab = dom.reindex(t.index).fillna("other").values
    gi = [g for g in gset if g in t.index]
    if len(gi) < 4: return dict(stat=np.nan, p_expr=np.nan, p_lineage=np.nan, n=len(gi), neuron_frac=np.nan)
    tv = t.values; pos = {g: i for i, g in enumerate(t.index)}; gix = np.array([pos[g] for g in gi]); stat = tv[gix].mean()
    pools_e = {b: np.where(bins == b)[0] for b in np.unique(bins)}
    strata = {}
    for i in range(len(tv)): strata.setdefault((lab[i], bins[i]), []).append(i)
    strata = {k: np.array(v) for k, v in strata.items()}
    null_e = np.array([np.mean([tv[rng.choice(pools_e[bins[i]])] for i in gix]) for _ in range(B)])
    def draw(i):
        k = (lab[i], bins[i]); pool = strata[k] if len(strata.get(k, [])) >= 5 else pools_e[bins[i]]; return tv[rng.choice(pool)]
    null_l = np.array([np.mean([draw(i) for i in gix]) for _ in range(B)])
    return dict(stat=round(float(stat), 3), p_expr=round((1 + np.sum(null_e <= stat)) / (B + 1), 4), p_lineage=round((1 + np.sum(null_l <= stat)) / (B + 1), 4), n=len(gi), neuron_frac=round(float(np.mean(lab[gix] == "neuron")), 2))
rows = []
for combo in list(pr.combo) + six:
    gset = programs[combo]; bs = combo.split("__x__")[1].split("__")[0]
    r = dict(program=combo.replace("__x__", " x ").replace("raw::", "").replace("relative::", "rel::"), priority=combo in set(pr.combo))
    for name in (["spatial_astro_top25", "snRNA_astro"] if bs == "astro" else ["spatial_mhc2_top25", "snRNA_mhc2"]):
        res = test(DES[name], gset); tag = "sp" if name.startswith("spatial") else "sn"
        r[tag + "_stat"] = res["stat"]; r[tag + "_p_expr"] = res["p_expr"]; r[tag + "_p_lineage"] = res["p_lineage"]; r["neuron_frac"] = res["neuron_frac"]
    rows.append(r)
df = pd.DataFrame(rows); df.to_csv(D + "/lineage_matched_null.csv", index=False)
print("=== PD-decrease p-values: expression-matched null vs lineage+expression-matched null ===")
print(df.to_string())
print("\npriority programs p_lineage<=0.05: spatial =", int((df.loc[df.priority, "sp_p_lineage"] <= 0.05).sum()), "| snRNA =", int((df.loc[df.priority, "sn_p_lineage"] <= 0.05).sum()))
# genome-wide reference: mean PD t of neuron-dominant vs other genes in each design
for name, (t, e) in DES.items():
    t = t.dropna(); lab = dom.reindex(t.index).fillna("other")
    print(f"{name:22s} mean t: neuron-dominant {t[lab == 'neuron'].mean():+.3f} (n={int((lab == 'neuron').sum())}) | astro {t[lab == 'astro'].mean():+.3f} | microglia {t[lab == 'microglia'].mean():+.3f} | oligo {t[lab == 'oligo'].mean():+.3f}")
