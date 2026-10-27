"""MR directionality: does each instrument explain more variance in the circulating protein than in the imaging
phenotype (Steiger)? Uses only the released exposure and outcome estimates, no new data.
usage: python mr_steiger.py OUTDIR
"""
import sys, warnings
import numpy as np, pandas as pd
from scipy import stats
warnings.filterwarnings("ignore")
M = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/blood_brain_MR_20261003"
OUT = sys.argv[1]; pd.set_option("display.width", 250)
mr = pd.read_csv(M + "/UKB_primary_MR.csv"); sel = mr[mr.selected].copy()
print("columns with a sample size:", [c for c in sel.columns if "n" == c.lower() or "_n" in c.lower() or c.lower().endswith("n")])
print("outcome_n range:", sel.outcome_n.min(), "-", sel.outcome_n.max())
inst = pd.read_csv(M + "/combined_exposure_instruments.csv")
print("exposure instrument columns:", [c for c in inst.columns if "n" in c.lower()][:10])
# UKB-PPP discovery sample size for the Olink assays used here
N_EXP = 34557
print(f"\nexposure N taken as the UKB-PPP European discovery sample ({N_EXP}); exposure r2 from the reported F")
tgt = sel[sel.candidate_antigens.isin(["CD22", "CD40", "CD71"])].copy()
tgt["r2_exposure"] = tgt.F_exposure / (N_EXP - 2 + tgt.F_exposure)
z_out = tgt.beta_outcome / tgt.se_outcome
tgt["r2_outcome"] = z_out ** 2 / (z_out ** 2 + tgt.outcome_n - 2)
tgt["steiger_ratio"] = tgt.r2_exposure / tgt.r2_outcome
# Steiger z test for correlated r (same variant, different traits, independent samples assumed)
def steiger_p(r2a, r2b, na, nb):
    ra, rb = np.sqrt(r2a), np.sqrt(r2b)
    za, zb = np.arctanh(np.clip(ra, 0, .999)), np.arctanh(np.clip(rb, 0, .999))
    se = np.sqrt(1 / (na - 3) + 1 / (nb - 3))
    return stats.norm.sf((za - zb) / se)
tgt["steiger_p"] = steiger_p(tgt.r2_exposure.values, tgt.r2_outcome.values, N_EXP, tgt.outcome_n.values)
sig = tgt[tgt.q_BH < 0.05]
print("\n=== the five supported pairs ===")
print(sig[["gene", "candidate_antigens", "outcome_trait", "F_exposure", "r2_exposure", "r2_outcome", "steiger_ratio", "steiger_p", "q_BH"]].round(6).to_string(index=False))
print(f"\n  all five explain more variance in protein than in imaging: {bool((sig.r2_exposure > sig.r2_outcome).all())}")
print(f"  ratio range {sig.steiger_ratio.min():.0f}x to {sig.steiger_ratio.max():.0f}x | Steiger P max {sig.steiger_p.max():.2e}")
print("\n=== all 48 target-outcome pairs ===")
print(f"  correct direction in {int((tgt.r2_exposure > tgt.r2_outcome).sum())}/{len(tgt)} pairs; smallest ratio {tgt.steiger_ratio.min():.0f}x")
tgt[["gene", "candidate_antigens", "outcome_trait", "F_exposure", "r2_exposure", "r2_outcome", "steiger_ratio", "steiger_p", "q_BH"]].to_csv(OUT + "/mr_steiger.csv", index=False)
print("\nwrote", OUT + "/mr_steiger.csv")
