"""Are the snRNA (Ma, subtype-resolved) and spatial (spot-level) replications complementary, and does the
recipient-lineage expression share of each program explain which modality replicates it?
usage: python ma_spatial_complement.py OUTDIR
"""
import sys, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
W = "/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017"; OUT = sys.argv[1]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
ma = pd.read_csv(W + "/diag/ma_cohort_test.csv"); sp = pd.read_csv(W + "/diag/spot_level_test.csv"); dg = pd.read_csv(W + "/diag/diag_programs.csv")
ma = ma[ma.priority]; sp = sp[sp.priority]; dg = dg[dg.priority]
def norm(x): return x.replace("rel::", "rel::").strip()
rows = []
for _, r in ma.iterrows():
    p = r.program; mg = "microglia_mhc2" in p
    ma_q = r.mg_activated_q if mg else r.astro_fibrous_q
    ma_t = r.mg_activated_t if mg else r.astro_fibrous_t
    alt_t = r.mg_baseline_t if mg else r.astro_protoplasmic_t
    s = sp[sp.program == p]; d = dg[dg.program == p]
    rows.append(dict(program=p, recipient="MG-AP" if mg else "Astro", n_genes=int(r.n),
                     ma_t=ma_t, ma_q=ma_q, ma_alt_state_t=alt_t,
                     spatial_stat=float(s.primary_stat.iloc[0]) if len(s) else np.nan, spatial_q=float(s.primary_q16.iloc[0]) if len(s) else np.nan,
                     own_lineage_share=float(d.own_lineage_share_mean.iloc[0]) if len(d) else np.nan))
t = pd.DataFrame(rows)
t["ma_pass"] = t.ma_q <= 0.05; t["spatial_pass"] = t.spatial_q <= 0.05
t["ma_evaluable"] = t.ma_q.notna()
print(t.round(4).sort_values(["recipient", "program"]).to_string(index=False))
print(f"\nMa pass {int(t.ma_pass.sum())}/{int(t.ma_evaluable.sum())} evaluable | spatial pass {int(t.spatial_pass.sum())}/16")
print(f"replicated in both: {int((t.ma_pass & t.spatial_pass).sum())} | in at least one: {int((t.ma_pass | t.spatial_pass).sum())}/16 | in neither: {int((~t.ma_pass & ~t.spatial_pass).sum())}")
print("\nprograms failing the spatial test and their Ma result:")
print(t[~t.spatial_pass][["program", "ma_t", "ma_q", "spatial_stat", "spatial_q", "own_lineage_share"]].round(4).to_string(index=False))
print("\nprograms failing or not evaluable in Ma and their spatial result:")
print(t[~t.ma_pass][["program", "ma_t", "ma_q", "spatial_stat", "spatial_q", "own_lineage_share"]].round(4).to_string(index=False))
print("\n=== lineage-share mechanism ===")
for col in ["own_lineage_share"]:
    a = t[t.spatial_pass][col].dropna(); b = t[~t.spatial_pass][col].dropna()
    print(f"  recipient-lineage share: spatial-pass median {a.median():.3f} (n={len(a)}) vs spatial-fail median {b.median():.3f} (n={len(b)}) | one-sided MWU P = {stats.mannwhitneyu(a, b, alternative='greater').pvalue:.4f}")
    for rc in ["Astro", "MG-AP"]:
        g = t[t.recipient == rc]
        a = g[g.spatial_pass][col].dropna(); b = g[~g.spatial_pass][col].dropna()
        if len(a) and len(b): print(f"    within {rc}: pass median {a.median():.3f} (n={len(a)}) vs fail median {b.median():.3f} (n={len(b)})")
print("\n=== direction consistency in the matching snRNA state ===")
ev = t[t.ma_evaluable]
print(f"  programs with a negative statistic in their matching Ma state: {int((ev.ma_t < 0).sum())}/{len(ev)}")
for rc in ["Astro", "MG-AP"]:
    g = ev[ev.recipient == rc]; print(f"    {rc}: matching-state t {g.ma_t.min():.3f} to {g.ma_t.max():.3f} (all negative: {bool((g.ma_t<0).all())}) | alternative-state t {g.ma_alt_state_t.min():.3f} to {g.ma_alt_state_t.max():.3f}")
t.to_csv(OUT + "/ma_spatial_complement.csv", index=False); print("\nwrote", OUT + "/ma_spatial_complement.csv")
