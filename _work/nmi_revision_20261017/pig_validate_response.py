"""Validate that my blood-side response reproduces the released brain-side convention exactly."""
import warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
P = "/public/home/mengxl/dzy/pd_product/publication/analysis/pd_product/results/pig_external_validation_20261003"
obs = pd.read_csv(P + "/observed_brain_log2_response.csv", index_col=0)
animals = list(obs.columns); lps = np.array([c.endswith("LPS") for c in animals])
brt = pd.read_csv(P + "/inputs/brain_one2one_sum_provided_transcript_TPM.csv", index_col=0)[animals].reindex(obs.index)
def response(T, scale):
    v = T.values.astype(float)
    if scale: v = v / v.sum(0, keepdims=True) * scale
    ctrl = v[:, ~lps].mean(1, keepdims=True)
    return np.log2((v + .1) / (ctrl + .1))
print("released brain response vs my reconstruction, per scaling choice:")
for name, sc in [("scale to 1e4", 1e4), ("scale to 1e6", 1e6), ("no rescaling", None)]:
    R = response(brt, sc)
    d = np.abs(R - obs.values)
    print(f"  {name:16s} max abs diff {d.max():.6f} | mean {d.mean():.6f} | identical={d.max()<1e-6}")
