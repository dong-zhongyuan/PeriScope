"""Test the ambient-neuronal-RNA explanation for the matched-brain PD decrease of the priority programs
(GSE157783, 5 PD / 6 control donors), and evaluate lineage-core subsets in snRNA and spatial C-SIDE.
usage: python diag_ambient.py OUTDIR
"""
import json, sys, warnings
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse, stats
warnings.filterwarnings("ignore")
A = "/public/home/mengxl/dzy/pd_product_assets"
R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/optimization_20261003"
OUT = sys.argv[1]; rng = np.random.RandomState(42)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
programs = json.load(open(R + "/joint_analysis/programs.json"))
pr = pd.read_csv(R + "/mr_parallel_gate_20261003/priority_all_programs.csv")
six = ["raw::CD294__classical_mono__x__astro__up", "raw::CD294__pDC__x__astro__up", "raw::CD30__classical_mono__x__astro__up",
       "raw::CD80__nonclassical_mono__x__astro__up", "raw::CD80__pDC__x__astro__up", "relative::CD274__pDC__x__astro__up"]
ge = pd.read_csv(R + "/spatial_subtype_author_20261003/spatial_cell_state_gene_effects.csv")
spatial = {st: g.set_index("gene") for st, g in ge.groupby("state")}

a = ad.read_h5ad(A + "/processed/gse157783/v0.1/gse157783_qc.h5ad")
sym = a.var["gene_symbol"].astype(str).values; uniq = sorted(set(sym)); gidx = {g: i for i, g in enumerate(uniq)}
collapse = sparse.csr_matrix((np.ones(len(sym)), (np.arange(len(sym)), [gidx[g] for g in sym])), shape=(len(sym), len(uniq)))
X = a.X @ collapse
obs = a.obs.copy()
for c in ["donor", "condition", "cell_type"]: obs[c] = obs[c].astype(str)
fine = pd.Series("other", index=obs.index)
for fam in ["brain_astro", "brain_microglia"]:
    m = pd.read_csv(A + f"/interim/v0.1/purification/c157_{fam}.csv", index_col=0); m = m[m.kept]
    common = fine.index.intersection(m.index); fine.loc[common] = m.loc[common, "state_pure"].astype(str)
obs["fine"] = fine.values
neuron_types = ["DaNs", "Excitatory", "GABA", "Inhibitory", "CADPS2+ neurons"]
donors = sorted(set(obs.donor)); is_pd = np.array([obs.condition[obs.donor == d].iloc[0] == "PD" for d in donors])
comp = pd.DataFrame({"neuron_frac": [(obs.cell_type[obs.donor == d].isin(neuron_types)).mean() for d in donors],
                     "DaN_frac": [(obs.cell_type[obs.donor == d] == "DaNs").mean() for d in donors],
                     "PD": is_pd}, index=donors)
print("=== donor composition ===\n", comp.round(3).to_string())
def pseudobulk(mask):
    M = []; keep = []
    for d in donors:
        mm = (mask & (obs.donor == d)).values
        if mm.sum() < 25: M.append(np.full(len(uniq), np.nan)); keep.append(False); continue
        s = np.asarray(X[mm].sum(0)).ravel(); M.append(np.log1p(s / max(s.sum(), 1) * 1e4)); keep.append(True)
    return np.stack(M), np.array(keep)
PB = {"astro": pseudobulk(obs.fine == "astro"), "microglia_mhc2": pseudobulk(obs.fine == "microglia_mhc2"), "neurons": pseudobulk(obs.cell_type.isin(neuron_types)),
      "oligo": pseudobulk(obs.cell_type == "Oligodendrocytes")}
def welch(M, keep):
    a_, b_ = M[keep & is_pd], M[keep & ~is_pd]
    return (a_.mean(0) - b_.mean(0)) / np.sqrt(a_.var(0, ddof=1) / len(a_) + b_.var(0, ddof=1) / len(b_) + 1e-12)
T = {k: pd.Series(welch(*v), index=uniq) for k, v in PB.items()}
mean_expr = {k: pd.Series(np.nanmean(v[0][v[1]], 0), index=uniq) for k, v in PB.items()}
# reference lineage shares
rp = pd.read_csv(A + "/results/optimization_20261003/spatial_subtype_author_20261003/inputs/reference_profiles.csv.gz", index_col=0)
lineage = {"astro_fibrous": "astro", "astro_protoplasmic": "astro", "microglia_homeostatic": "microglia", "microglia_mhc2": "microglia",
           "neuron_GABAergic": "neuron", "neuron_MEG3": "neuron", "neuron_dopaminergic": "neuron", "Oligo": "oligo", "OPC": "opc", "Endothelial": "endo", "BAM": "bam", "Tcells": "tcell"}
share = rp.T.groupby(lambda c: lineage[c]).sum().T; share = share.div(share.sum(1) + 1e-12, axis=0)
# expression-matched null helper (snRNA universe within a state)
def matched_p(state, genes, tseries, expr, B=3000):
    gi = [g for g in genes if g in tseries.index and np.isfinite(tseries[g])]
    if len(gi) < 4: return np.nan, np.nan, len(gi)
    uni = tseries.dropna(); e = expr.reindex(uni.index); bins = pd.qcut(e.rank(method="first"), 20, labels=False)
    obs_stat = uni[gi].mean(); pools = {b: np.where(bins.values == b)[0] for b in np.unique(bins.values)}; gb = bins[gi].values
    null = np.array([np.mean([uni.values[rng.choice(pools[b])] for b in gb]) for _ in range(B)])
    return round(float(obs_stat), 3), round((1 + np.sum(null <= obs_stat)) / (B + 1), 4), len(gi)
rows = []
for combo in list(pr.combo) + six:
    genes = programs[combo]; bs = combo.split("__x__")[1].split("__")[0]; own = "astro" if bs == "astro" else "microglia"; st = "astro_fibrous" if bs == "astro" else "microglia_mhc2"
    sh = share.reindex(genes).dropna()
    tg = T[bs].reindex(sh.index)
    # (a) gene level: does the glial PD decrease scale with the gene's neuronal share?
    ok = tg.notna()
    rho_a = stats.spearmanr(sh.loc[ok, "neuron"], tg[ok]).correlation if ok.sum() > 5 else np.nan
    # (b) donor level: program score in own state vs donor neuron fraction (within-condition partial)
    M, keep = PB[bs]; gi = [gidx[g] for g in genes if g in gidx]
    Z = (M[:, gi] - np.nanmean(M[keep][:, gi], 0)) / (np.nanstd(M[keep][:, gi], 0) + 1e-9); score = np.nanmean(Z, 1)
    d = pd.DataFrame({"score": score, "nf": comp.neuron_frac.values, "pd": is_pd.astype(float)})[keep]
    rho_b = stats.spearmanr(d.score, d.nf).correlation
    resid_s = d.score - d.groupby("pd").score.transform("mean"); resid_n = d.nf - d.groupby("pd").nf.transform("mean")
    rho_b_within = stats.spearmanr(resid_s, resid_n).correlation
    # (c) lineage core: genes the recipient lineage actually expresses
    core = [g for g in sh.index if sh.loc[g, own] >= 0.25 and sh.loc[g, "neuron"] <= 0.5]
    noncore = [g for g in sh.index if g not in set(core)]
    sn_core = matched_p(bs, core, T[bs], mean_expr[bs]); sn_non = matched_p(bs, noncore, T[bs], mean_expr[bs])
    sp = spatial[st]; sp_core = matched_p(st, core, sp.t, sp.mean_log_expression); sp_all = matched_p(st, genes, sp.t, sp.mean_log_expression)
    rows.append(dict(program=combo.replace("__x__", " x ").replace("raw::", "").replace("relative::", "rel::"), priority=combo in set(pr.combo), n=len(genes),
                     rho_glialT_vs_neuronShare=round(float(rho_a), 2), rho_score_vs_neuronFrac=round(float(rho_b), 2), rho_within_condition=round(float(rho_b_within), 2),
                     n_core=len(core), snRNA_core_t=sn_core[0], snRNA_core_p_down=sn_core[1], snRNA_noncore_t=sn_non[0], snRNA_noncore_p_down=sn_non[1],
                     spatial_core_t=sp_core[0], spatial_core_p_down=sp_core[1], spatial_core_n=sp_core[2], spatial_all_p_down=sp_all[1]))
df = pd.DataFrame(rows); df.to_csv(OUT + "/diag_ambient.csv", index=False)
print("\n=== ambient-RNA diagnostics and lineage-core re-tests ===")
print(df.to_string())
print("\nneuron fraction PD vs control (Mann-Whitney):", stats.mannwhitneyu(comp.neuron_frac[is_pd], comp.neuron_frac[~is_pd]))
# global check: across ALL genes, is the astro PD t related to neuron share?
allg = share.index.intersection(T["astro"].dropna().index)
print("all genes: Spearman(astro PD t, neuron share) =", round(stats.spearmanr(share.loc[allg, "neuron"], T["astro"][allg]).correlation, 3), "| n =", len(allg))
print("all genes: Spearman(MG-AP PD t, neuron share) =", round(stats.spearmanr(share.loc[allg, "neuron"], T["microglia_mhc2"][allg]).correlation, 3))
print("all genes: Spearman(neuron PD t, neuron share) =", round(stats.spearmanr(share.loc[allg, "neuron"], T["neurons"][allg]).correlation, 3))
