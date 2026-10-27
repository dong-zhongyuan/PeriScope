"""MR regional check: for every supported protein-imaging pair, is the imaging association in the surrounding
region centred on the protein instrument, as a shared causal variant implies? Harmonised GWAS Catalog outcome
files are read remotely with tabix; no new data are downloaded.
usage: python mr_region.py OUTDIR
"""
import sys, subprocess, io, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
M = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/blood_brain_MR_20261003"
OUT = sys.argv[1]; WIN = 250000
pd.set_option("display.width", 250)
mr = pd.read_csv(M + "/UKB_primary_MR.csv")
mr = mr[mr.candidate_antigens.isin(["CD22", "CD40", "CD71"]) & mr.selected]
sig = mr[mr.q_BH < 0.05].copy()
print(f"supported pairs (BH q<0.05): {len(sig)}")
COLS = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "variant_id", "ci_upper", "ci_lower", "neg_log10_p_value", "hm_coordinate_conversion", "odds_ratio", "hm_code", "rsid"]
rows = []
for _, r in sig.iterrows():
    url = r.source_url_outcome
    if "harmonised" not in str(url): print(f"  {r.gene}/{r.outcome_trait}: not harmonised, skipped"); continue
    reg = f"{int(r.chromosome)}:{max(1,int(r.position)-WIN)}-{int(r.position)+WIN}"
    out = subprocess.run(["tabix", url, reg], capture_output=True, text=True, timeout=900)
    if out.returncode != 0 or not out.stdout.strip():
        print(f"  {r.gene}/{r.outcome_trait}: tabix empty ({out.stderr.strip()[:80]})"); continue
    d = pd.read_csv(io.StringIO(out.stdout), sep="\t", names=COLS, low_memory=False)
    d = d[pd.to_numeric(d.p_value, errors="coerce").notna()]; d["p_value"] = d.p_value.astype(float)
    d = d[d.p_value > 0]
    lead = d.loc[d.p_value.idxmin()]
    inst = d[d.rsid == r.rsid]
    if inst.empty: inst = d[d.base_pair_location == int(r.position)]
    pi = float(inst.p_value.iloc[0]) if len(inst) else np.nan
    rank = int((d.p_value < pi).sum()) + 1 if np.isfinite(pi) else -1
    dist = abs(int(lead.base_pair_location) - int(r.position))
    print(f"  {r.gene:5s} / {r.outcome_trait:34s} n={len(d):6d} in +/-250 kb | instrument {r.rsid} P={pi:.2e} rank {rank} (top {100*rank/len(d):.2f}%) | regional lead {lead.rsid} P={float(lead.p_value):.2e} at {dist/1000:.1f} kb")
    rows.append(dict(gene=r.gene, antigen=r.candidate_antigens, outcome=r.outcome_trait, accession=r.outcome_accession,
                     instrument=r.rsid, instrument_p=pi, instrument_rank=rank, n_variants=len(d), rank_pct=100 * rank / len(d),
                     regional_lead=lead.rsid, regional_lead_p=float(lead.p_value), distance_kb=dist / 1000.0))
pd.DataFrame(rows).to_csv(OUT + "/mr_regional_check.csv", index=False)
print("\nwrote", OUT + "/mr_regional_check.csv")
