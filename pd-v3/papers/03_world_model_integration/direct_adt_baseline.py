"""How well can a plain supervised RNA -> ADT model do on the same metric?

The world model reaches median Spearman 0.090 on cell-level ADT while the
measurement ceiling is 0.567. That gap is either architectural (the latent is
optimised for the coupling objective, not for ADT) or intrinsic (per-cell RNA
simply does not determine per-cell ADT). A direct supervised baseline separates
the two: it is an empirical predictability bound for this input/target pair.

Trained only on CITE-seq donors that are not held out; scored with exactly the
original metric on P6/P7/P8. Read-only on the archive; writes into pd-v3.
"""
import json, time
from pathlib import Path
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from scipy.stats import spearmanr

I = Path("/public/home/mengxl/dzy/pd_product_assets/interim/rescreen_20261002")
O = Path("/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003")
OUT = Path("/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/direct_adt_baseline_20261118")
OUT.mkdir(parents=True, exist_ok=True)
dev = "cuda" if torch.cuda.is_available() else "cpu"
torch.set_num_threads(4)

cg = np.load(I / "citeseq.npz")
X, Y = cg["X"], cg["Y"]
donor, state, sub = cg["donor"], cg["state"], cg["subtype"]
is_test, is_val = cg["is_test"], cg["is_val"]
vocab = cg["state_vocab"].astype(str).tolist()
names = cg["donor_names"].astype(str)
registry = [r for r in json.loads((O / "candidate_registry.json").read_text())]

tr = np.where(~is_test & ~is_val)[0]
va = np.where(is_val & ~is_test)[0]
te = np.where(is_test)[0]
print("cells train/val/test: %d / %d / %d" % (len(tr), len(va), len(te)))
print("train donors:", sorted(set(names[donor[tr]])), " test donors:", sorted(set(names[donor[te]])))

Xt = torch.from_numpy(X).float()
Yt = torch.from_numpy(Y).float()


def run(hidden, epochs, lr, seed):
    torch.manual_seed(seed)
    net = nn.Sequential(nn.Linear(X.shape[1], hidden), nn.GELU(), nn.Dropout(0.1),
                        nn.Linear(hidden, hidden), nn.GELU(), nn.Dropout(0.1),
                        nn.Linear(hidden, Y.shape[1])).to(dev)
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=1e-4)
    best, best_state = float("inf"), None
    idx = tr.copy()
    for ep in range(epochs):
        net.train()
        np.random.default_rng(seed * 100 + ep).shuffle(idx)
        for s in range(0, len(idx), 512):
            b = idx[s:s + 512]
            opt.zero_grad()
            loss = nn.functional.mse_loss(net(Xt[b].to(dev)), Yt[b].to(dev))
            loss.backward()
            opt.step()
        net.eval()
        with torch.no_grad():
            vl = float(np.mean([nn.functional.mse_loss(net(Xt[va[s:s + 2048]].to(dev)),
                                                       Yt[va[s:s + 2048]].to(dev)).item()
                                for s in range(0, len(va), 2048)]))
        if vl < best:
            best, best_state = vl, {k: v.detach().clone() for k, v in net.state_dict().items()}
        print("    epoch %d val_mse %.5f%s" % (ep, vl, "  *" if vl == best else ""), flush=True)
    net.load_state_dict(best_state)
    net.eval()
    with torch.no_grad():
        pred = np.concatenate([net(Xt[te[s:s + 2048]].to(dev)).cpu().numpy()
                               for s in range(0, len(te), 2048)])
    return pred, best


def score(pred):
    rows = []
    Yte = Y[te]
    for st in ["cDC", "classical_mono", "nonclassical_mono", "pDC"]:
        si = vocab.index(st)
        for d in np.unique(donor[te]):
            m = (state[te] == si) & (donor[te] == d)
            if m.sum() < 25:
                continue
            for r in registry:
                if r.get("is_control"):
                    continue
                for subtype in np.unique(sub[te][m]):
                    f = m & (sub[te] == subtype)
                    if f.sum() < 25:
                        continue
                    y = Yte[f][:, r["channels"]].mean(1)
                    p = pred[f][:, r["channels"]].mean(1)
                    if np.std(y) < 1e-8 or np.std(p) < 1e-8:
                        continue
                    rho = spearmanr(y, p).statistic
                    if not np.isnan(rho):
                        rows.append((r["target"], st, str(subtype), names[d], int(f.sum()), float(rho)))
    return pd.DataFrame(rows, columns=["target", "blood_state", "citeseq_subtype", "donor", "n", "rho"])


started = time.time()
results = {}
for hidden, epochs, lr in ((512, 6, 1e-3), (1024, 6, 1e-3)):
    tag = "mlp%d" % hidden
    print("  training %s ..." % tag, flush=True)
    pred, vl = run(hidden, epochs, lr, 42)
    t = score(pred)
    t["model"] = tag
    results[tag] = t
    print("  %s: groups=%d  median rho=%.4f  val_mse=%.5f" % (tag, len(t), t.rho.median(), vl), flush=True)

allr = pd.concat(results.values(), ignore_index=True)
allr.to_csv(OUT / "direct_baseline_protein_accuracy.csv", index=False)
summary = dict(
    created_at_utc=pd.Timestamp.utcnow().isoformat(),
    device=dev,
    elapsed_seconds=round(time.time() - started, 1),
    train_cells=int(len(tr)), val_cells=int(len(va)), test_cells=int(len(te)),
    test_donors=sorted(set(names[donor[te]].tolist())),
    metric="identical to evaluate_models.py: Spearman across cells within donor x blood_state x citeseq_subtype, per target",
    median_rho={k: round(float(v.rho.median()), 4) for k, v in results.items()},
    reference_world_model_median_rho=0.0900,
    reference_coupling_off_median_rho=0.0954,
    measurement_ceiling_median=0.5668,
    note="supervised RNA->ADT only; no coupling, no brain side. An empirical "
         "predictability bound for this input/target pair under the same metric.",
)
(OUT / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
print()
print(json.dumps(summary["median_rho"], indent=2))
print("world model 0.0900 | coupling_off 0.0954 | measurement ceiling 0.5668")
print("written:", OUT)
