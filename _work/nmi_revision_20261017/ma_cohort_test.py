"""Third-cohort, subtype-resolved snRNA test of the programs in the Ma et al. 2025 substantia nigra atlas
(author-annotated Fibrous/Protoplasmic astrocytes and Activated/Baseline microglia; PD vs control donors).
Input: donor x sublineage pseudobulk counts exported from the Seurat object. usage: python ma_cohort_test.py DIAGDIR
"""
import json, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
R = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/optimization_20261003"
D = sys.argv[1]; rng = np.random.RandomState(42)
pd.set_option("display.width", 260); pd.set_option("display.max_columns", 40)
programs = json.load(open(R + "/joint_analysis/programs.json")); pr = pd.read_csv(R + "/mr_parallel_gate_20261003/priority_all_programs.csv")
six = ["raw::CD294__classical_mono__x__astro__up", "raw::CD294__pDC__x__astro__up", "raw::CD30__classical_mono__x__astro__up",
       "raw::CD80__nonclassical_mono__x__astro__up", "raw::CD80__pDC__x__astro__up", "relative::CD274__pDC__x__astro__up"]
pb = pd.read_csv(D + "/ma/ma_pseudobulk_counts.csv", index_col=0)  # genes x "Donor|lineage|sublineage"
cells = pd.read_csv(D + "/ma/ma_pseudobulk_cells.csv").set_index("group").n_cells
donors = pd.read_csv(D + "/ma/ma_donors.csv"); cond = donors.set_index("Donor").Condition.astype(str).to_dict()
print("donors:", len(cond), "| PD:", sum(v == "PD" for v in cond.values()), "| Control:", sum(v == "Control" for v in cond.values()))
print("groups:", pb.shape[1], "| genes:", pb.shape[0])
cols = pd.DataFrame([c.split("|") for c in pb.columns], columns=["donor", "lineage", "sub"], index=pb.columns)
def state_matrix(lineage, subs, min_cells=25):
    sel = cols[(cols.lineage == lineage) & (cols["sub"].isin(subs))]
    M = {};
    for d, g in sel.groupby("donor"):
        if cells[g.index].sum() < min_cells: continue
        s = pb[g.index].sum(1).values; M[d] = np.log1p(s / max(s.sum(), 1) * 1e4)
    M = pd.DataFrame(M, index=pb.index).T
    is_pd = np.array([cond.get(d) == "PD" for d in M.index]); known = np.array([d in cond for d in M.index])
    M = M[known]; is_pd = is_pd[known]
    a_, b_ = M.values[is_pd], M.values[~is_pd]
    t = (a_.mean(0) - b_.mean(0)) / np.sqrt(a_.var(0, ddof=1) / len(a_) + b_.var(0, ddof=1) / len(b_) + 1e-12)
    detected = (M.values > 0).sum(0) >= 3; t[~detected] = np.nan
    return dict(t=pd.Series(t, index=pb.index), expr=pd.Series(M.values.mean(0), index=pb.index), nPD=int(is_pd.sum()), nC=int((~is_pd).sum()))
S = {"astro_all": state_matrix("Astrocyte", ["Fibrous", "Protoplasmic"]), "astro_fibrous": state_matrix("Astrocyte", ["Fibrous"]), "astro_protoplasmic": state_matrix("Astrocyte", ["Protoplasmic"]),
     "mg_activated": state_matrix("Myeloid", ["Activated"]), "mg_baseline": state_matrix("Myeloid", ["Baseline"]), "neuron_all": state_matrix("Neuron", ["Dopaminergic", "GABEergic", "MEG3_Interneurons"]), "oligo": state_matrix("Oligo", list(cols[cols.lineage == "Oligo"]["sub"].unique()))}
for k, v in S.items(): print(f"{k:20s} donors PD/Ctrl = {v['nPD']}/{v['nC']} | genes with t: {int(v['t'].notna().sum())}")
def matched_p(d, gset, B=5000):
    t = d["t"].dropna(); e = d["expr"].reindex(t.index); bins = pd.qcut(e.rank(method="first"), 20, labels=False)
    gi = [g for g in gset if g in t.index]
    if len(gi) < 4 or len(gi) / len(gset) < 0.5: return dict(stat=np.nan, p_down=np.nan, n=len(gi))
    stat = t[gi].mean(); pools = {b: np.where(bins.values == b)[0] for b in np.unique(bins.values)}; gb = bins[gi].values; tv = t.values
    null = np.array([np.mean([tv[rng.choice(pools[b])] for b in gb]) for _ in range(B)])
    return dict(stat=round(float(stat), 3), p_down=round((1 + np.sum(null <= stat)) / (B + 1), 4), n=len(gi))
def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p); q = np.empty(n); run = 1.0
    for r in range(n, 0, -1):
        i = o[r - 1]; run = min(run, p[i] * n / r); q[i] = run
    return q
rows = []
for combo in list(pr.combo) + six:
    gs = programs[combo]; bs = combo.split("__x__")[1].split("__")[0]
    states = ["astro_all", "astro_fibrous", "astro_protoplasmic"] if bs == "astro" else ["mg_activated", "mg_baseline"]
    r = dict(program=combo.replace("__x__", " x ").replace("raw::", "").replace("relative::", "rel::"), priority=combo in set(pr.combo), n=len(gs))
    for st in states:
        res = matched_p(S[st], gs); r[st + "_t"] = res["stat"]; r[st + "_p_down"] = res["p_down"]
    r["neuron_t"] = matched_p(S["neuron_all"], gs)["stat"]; r["oligo_t"] = matched_p(S["oligo"], gs)["stat"]
    rows.append(r)
df = pd.DataFrame(rows)
for st in ["astro_all", "astro_fibrous", "astro_protoplasmic", "mg_activated", "mg_baseline"]:
    col = st + "_p_down"
    if col in df:
        m = df.priority.values & df[col].notna().values; q = np.full(len(df), np.nan); q[m] = bh(df.loc[m, col].values); df[st + "_q"] = np.round(q, 4)
df.to_csv(D + "/ma_cohort_test.csv", index=False)
print("\n=== Ma 2025 snRNA cohort, subtype-resolved donor-level test (p_down = PD decrease; q = BH over priority programs within state) ===")
print(df.to_string())
