"""Diagnose why the 16 priority programs do not reproduce their snRNA-seq PD decrease in the
subtype-resolved spatial (C-SIDE) test. Uses GSE157783 (independent snRNA cohort), the spatial
C-SIDE gene effects, and the Ma author reference profiles.

Outputs (stdout + CSVs in OUT):
  1. per-program gene coverage in the spatial estimates, and snRNA effects of covered vs uncovered genes
  2. per-program PD effect of the same genes in other cell types (neurons, oligodendrocytes, microglia, astro)
  3. per-program lineage attribution of genes from the Ma reference profiles
  4. gene-level agreement between snRNA (GSE157783 astro / MG-AP) and spatial C-SIDE effects
"""
import json, sys, warnings
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"
R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/optimization_20261003"
OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)

programs = json.load(open(R + "/joint_analysis/programs.json"))
pr = pd.read_csv(R + "/mr_parallel_gate_20261003/priority_all_programs.csv")
six = ["raw::CD294__classical_mono__x__astro__up", "raw::CD294__pDC__x__astro__up", "raw::CD30__classical_mono__x__astro__up",
       "raw::CD80__nonclassical_mono__x__astro__up", "raw::CD80__pDC__x__astro__up", "relative::CD274__pDC__x__astro__up"]
ge = pd.read_csv(R + "/spatial_subtype_author_20261003/spatial_cell_state_gene_effects.csv")
spatial_t = {st: g.set_index("gene").t for st, g in ge.groupby("state")}

# ---------------------------------------------------------------- GSE157783 pseudobulk per donor x state
a = ad.read_h5ad(A + "/processed/gse157783/v0.1/gse157783_qc.h5ad")
sym = a.var["gene_symbol"].astype(str).values
uniq = sorted(set(sym)); gidx = {g: i for i, g in enumerate(uniq)}
collapse = sparse.csr_matrix((np.ones(len(sym)), (np.arange(len(sym)), [gidx[g] for g in sym])), shape=(len(sym), len(uniq)))
X = a.X @ collapse  # cells x unique genes, raw counts
obs = a.obs.copy(); obs["donor"] = obs["donor"].astype(str); obs["condition"] = obs["condition"].astype(str); obs["cell_type"] = obs["cell_type"].astype(str)
fine = pd.Series("other", index=obs.index)
for fam in ["brain_astro", "brain_microglia"]:
    m = pd.read_csv(A + f"/interim/v0.1/purification/c157_{fam}.csv", index_col=0)
    if fam == "brain_astro": print("purification map columns:", list(m.columns), "| kept:", int(m.kept.sum()), "| states:", m.loc[m.kept, 'state_pure'].value_counts().to_dict())
    m = m[m.kept]; common = fine.index.intersection(m.index); fine.loc[common] = m.loc[common, "state_pure"].astype(str)
obs["fine"] = fine.values
print("fine states:", obs.fine.value_counts().to_dict())
groups = {"astro(fine)": obs.fine == "astro", "MG-AP(fine)": obs.fine == "microglia_mhc2", "MG-H(fine)": obs.fine == "microglia_homeostatic",
          "Astrocytes(all)": obs.cell_type == "Astrocytes", "Microglia(all)": obs.cell_type == "Microglia",
          "Neurons": obs.cell_type.isin(["DaNs", "Excitatory", "GABA", "Inhibitory", "CADPS2+ neurons"]), "DaNs": obs.cell_type == "DaNs",
          "Oligodendrocytes": obs.cell_type == "Oligodendrocytes", "OPCs": obs.cell_type == "OPCs", "Endothelial": obs.cell_type == "Endothelial cells"}
def welch_t(group_mask):
    donors = sorted(set(obs.donor[group_mask])); rows = []; cond = []
    for d in donors:
        m = (group_mask & (obs.donor == d)).values
        if m.sum() < 25: continue
        s = np.asarray(X[m].sum(0)).ravel(); rows.append(np.log1p(s / max(s.sum(), 1) * 1e4)); cond.append(obs.condition[m].iloc[0] == "PD")
    if len(rows) == 0: return None, (0, 0)
    M = np.stack(rows); c = np.array(cond)
    if c.sum() < 2 or (~c).sum() < 2: return None, (int(c.sum()), int((~c).sum()))
    a_, b_ = M[c], M[~c]
    t = (a_.mean(0) - b_.mean(0)) / np.sqrt(a_.var(0, ddof=1) / len(a_) + b_.var(0, ddof=1) / len(b_) + 1e-12)
    expressed = (M > 0).mean(0) > 0  # gene detected in at least one donor pseudobulk
    t[~expressed] = np.nan
    return pd.Series(t, index=uniq), (int(c.sum()), int((~c).sum()))
T = {}
for name, mask in groups.items():
    t, n = welch_t(mask); T[name] = t; print(f"{name:18s} donors PD/Ctrl = {n} | genes with t: {0 if t is None else int(t.notna().sum())}")

# ---------------------------------------------------------------- Ma reference lineage attribution
rp = pd.read_csv(A + "/results/optimization_20261003/spatial_subtype_author_20261003/inputs/reference_profiles.csv.gz", index_col=0)
lineage = {"astro_fibrous": "astro", "astro_protoplasmic": "astro", "microglia_homeostatic": "microglia", "microglia_mhc2": "microglia",
           "neuron_GABAergic": "neuron", "neuron_MEG3": "neuron", "neuron_dopaminergic": "neuron", "Oligo": "oligo", "OPC": "opc", "Endothelial": "endo", "BAM": "bam", "Tcells": "tcell"}
rpl = rp.T.groupby(lambda c: lineage[c]).sum().T  # genes x lineages (sum of state profiles)
share = rpl.div(rpl.sum(1) + 1e-12, axis=0)
dom = share.idxmax(1)

# ---------------------------------------------------------------- per-program diagnostics
rows = []
for combo in list(pr.combo) + six:
    genes = programs[combo]; bs = combo.split("__x__")[1].split("__")[0]
    st = "astro_fibrous" if bs == "astro" else "microglia_mhc2"; sn_state = "astro(fine)" if bs == "astro" else "MG-AP(fine)"
    cov = [g for g in genes if g in spatial_t[st].index]; unc = [g for g in genes if g not in spatial_t[st].index]
    tsn = T[sn_state]
    def mean_t(gs, series): v = series.reindex(gs).dropna(); return (round(float(v.mean()), 2), len(v)) if len(v) else (np.nan, 0)
    r = dict(program=combo.replace("__x__", " x ").replace("raw::", "").replace("relative::", "rel::"), priority=combo in set(pr.combo), n_genes=len(genes), spatial_covered=len(cov), spatial_uncovered=len(unc))
    r["snRNA_t_covered"], _ = mean_t(cov, tsn); r["snRNA_t_uncovered"], _ = mean_t(unc, tsn)
    sp = spatial_t[st].reindex(cov).dropna(); r["spatial_t_covered"] = round(float(sp.mean()), 2) if len(sp) else np.nan
    both = tsn.reindex(cov).dropna().index.intersection(sp.index)
    r["gene_corr_snRNA_vs_spatial"] = round(float(np.corrcoef(tsn[both], sp[both])[0, 1]), 2) if len(both) > 4 else np.nan
    r["sign_agree_frac"] = round(float(np.mean(np.sign(tsn[both]) == np.sign(sp[both]))), 2) if len(both) > 4 else np.nan
    for name in ["Astrocytes(all)", "Microglia(all)", "Neurons", "DaNs", "Oligodendrocytes", "OPCs", "Endothelial"]:
        r["t_" + name] = mean_t(genes, T[name])[0] if T[name] is not None else np.nan
    d = dom.reindex(genes).dropna(); r["lineage_neuron_frac"] = round(float((d == "neuron").mean()), 2); r["lineage_astro_frac"] = round(float((d == "astro").mean()), 2); r["lineage_microglia_frac"] = round(float((d == "microglia").mean()), 2); r["lineage_oligo_frac"] = round(float((d == "oligo").mean()), 2)
    # mean expression share in the program's own lineage
    own = "astro" if bs == "astro" else "microglia"
    r["own_lineage_share_mean"] = round(float(share.reindex(genes)[own].dropna().mean()), 2)
    rows.append(r)
df = pd.DataFrame(rows)
df.to_csv(OUT + "/diag_programs.csv", index=False)
print("\n=== per-program diagnostics (16 priority + 6 spatially supported) ===")
print(df.to_string())
# gene-level table for the priority programs
gl = []
for combo in pr.combo:
    bs = combo.split("__x__")[1].split("__")[0]; st = "astro_fibrous" if bs == "astro" else "microglia_mhc2"; sn_state = "astro(fine)" if bs == "astro" else "MG-AP(fine)"
    for g in programs[combo]:
        gl.append(dict(program=combo, gene=g, snRNA_t=T[sn_state].get(g, np.nan), spatial_t=spatial_t[st].get(g, np.nan), neuron_t=(T["Neurons"].get(g, np.nan) if T["Neurons"] is not None else np.nan), oligo_t=(T["Oligodendrocytes"].get(g, np.nan) if T["Oligodendrocytes"] is not None else np.nan), dominant_lineage=dom.get(g, "NA"), own_share=share.reindex([g]).iloc[0]["astro" if bs == "astro" else "microglia"] if g in share.index else np.nan))
pd.DataFrame(gl).to_csv(OUT + "/diag_genes.csv", index=False)
print("\nwrote", OUT + "/diag_programs.csv", OUT + "/diag_genes.csv")
