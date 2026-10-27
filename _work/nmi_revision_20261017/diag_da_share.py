"""Does the glial PD decrease of program genes track dopaminergic-neuron expression share (DA-neuron loss ambient signal)?
GSE157783 astro / MG-AP pseudobulk Welch t vs Ma reference neuron_dopaminergic share. usage: python diag_da_share.py OUTDIR
"""
import json, sys, warnings
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse, stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"; R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/optimization_20261003"
programs = json.load(open(R + "/joint_analysis/programs.json")); pr = pd.read_csv(R + "/mr_parallel_gate_20261003/priority_all_programs.csv")
a = ad.read_h5ad(A + "/processed/gse157783/v0.1/gse157783_qc.h5ad")
sym = a.var["gene_symbol"].astype(str).values; uniq = sorted(set(sym)); gidx = {g: i for i, g in enumerate(uniq)}
X = a.X @ sparse.csr_matrix((np.ones(len(sym)), (np.arange(len(sym)), [gidx[g] for g in sym])), shape=(len(sym), len(uniq)))
obs = a.obs.copy()
for c in ["donor", "condition", "cell_type"]: obs[c] = obs[c].astype(str)
fine = pd.Series("other", index=obs.index)
for fam in ["brain_astro", "brain_microglia"]:
    m = pd.read_csv(A + f"/interim/v0.1/purification/c157_{fam}.csv", index_col=0); m = m[m.kept]; common = fine.index.intersection(m.index); fine.loc[common] = m.loc[common, "state_pure"].astype(str)
obs["fine"] = fine.values
donors = sorted(set(obs.donor)); is_pd = np.array([obs.condition[obs.donor == d].iloc[0] == "PD" for d in donors])
def tseries(mask):
    M = []
    for d in donors:
        mm = (mask & (obs.donor == d)).values; s = np.asarray(X[mm].sum(0)).ravel(); M.append(np.log1p(s / max(s.sum(), 1) * 1e4))
    M = np.stack(M); a_, b_ = M[is_pd], M[~is_pd]
    t = (a_.mean(0) - b_.mean(0)) / np.sqrt(a_.var(0, ddof=1) / len(a_) + b_.var(0, ddof=1) / len(b_) + 1e-12)
    t[(M > 0).sum(0) < 3] = np.nan; return pd.Series(t, index=uniq), pd.Series(M.mean(0), index=uniq)
T = {}; E = {}
for name, mask in [("astro", obs.fine == "astro"), ("microglia_mhc2", obs.fine == "microglia_mhc2"), ("oligo", obs.cell_type == "Oligodendrocytes"), ("opc", obs.cell_type == "OPCs")]:
    T[name], E[name] = tseries(mask)
rp = pd.read_csv(A + "/results/optimization_20261003/spatial_subtype_author_20261003/inputs/reference_profiles.csv.gz", index_col=0)
sh = rp.div(rp.sum(1) + 1e-12, axis=0)
da = sh["neuron_dopaminergic"]; nonda_neuron = sh["neuron_GABAergic"] + sh["neuron_MEG3"]
common = sh.index.intersection(T["astro"].dropna().index)
print("=== genome-wide: Spearman(PD t in cell type, DA-neuron share) ===")
for name in T:
    cc = sh.index.intersection(T[name].dropna().index)
    print(f"{name:15s} DA share rho = {stats.spearmanr(da[cc], T[name][cc]).correlation:+.3f} | non-DA neuron share rho = {stats.spearmanr(nonda_neuron[cc], T[name][cc]).correlation:+.3f} | n={len(cc)}")
# top DA marker genes (highest DA share, expressed): their PD t in glia
top = da[common].sort_values(ascending=False).head(200).index
print("\nTop-200 DA-share genes: mean PD t in astro =", round(float(T['astro'][top].mean()), 3), "| MG-AP =", round(float(T['microglia_mhc2'][top].mean()), 3), "| oligo =", round(float(T['oligo'][top].mean()), 3), "| OPC =", round(float(T['opc'][top].mean()), 3))
print("examples:", list(top[:12]))
print("\n=== per priority program: DA share of genes and correlation of glial t with DA share ===")
rows = []
for combo in pr.combo:
    gs = [g for g in programs[combo] if g in common]; bs = combo.split("__x__")[1].split("__")[0]
    t = T[bs][gs]; d = da[gs]
    hi = d >= 0.3
    rows.append(dict(program=combo.replace("__x__", " x ")[:46], n=len(gs), mean_DA_share=round(float(d.mean()), 3), n_DA_dominant=int(hi.sum()),
                     t_DA_dominant=round(float(t[hi].mean()), 2) if hi.sum() else np.nan, t_other=round(float(t[~hi].mean()), 2), rho_t_vs_DAshare=round(float(stats.spearmanr(d, t).correlation), 2)))
df = pd.DataFrame(rows); df.to_csv(OUT + "/diag_da_share.csv" if (OUT := sys.argv[1]) else None, index=False); print(df.to_string())
