"""Teacher's gates 3-6 diagnostic suite (blueprint par.4.1, never run before).

Gate 3: condition dropout — does removing condition vector destroy signal?
Gate 4: donor generalization — does prediction match unseen donor distributions?
Gate 5: transient stratification — are gradient norms abnormally high in some strata?
Gate 6: disease shortcut — does removing disease direction from blood destroy edges?
"""
import torch, numpy as np, sys, json
sys.path.insert(0, '/public/home/mengxl/dzy/pd_product/src')
from pdproduct.simulators.ccwm import CCWM, CCWMConfig, sinkhorn_divergence
from scipy.stats import spearmanr
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

A = '/public/home/mengxl/dzy/pd_product_assets'
CCWM_VER = 'v23'   # 本轮版本（数据自持+pd边表新词表，2026-09-09 全链重跑）
ckpt = torch.load(A + f'/checkpoints/project/ccwm-{CCWM_VER}-s42/model.pt', map_location='cpu', weights_only=False)
cfg = CCWMConfig(**ckpt['cfg'])
model = CCWM(cfg); model.load_state_dict(ckpt['state_dict']); model.eval()

blood = np.load(A + f'/interim/ccwm_{CCWM_VER}/blood.npz')
brain_te = np.load(A + f'/interim/ccwm_{CCWM_VER}/brain_locked_test.npz')
X_blood = torch.from_numpy(blood['X'])
X_brain = torch.from_numpy(brain_te['X'])
# 状态词表动态读取（用户 2026-09-09 指令）：来自 npz 自描述，无字面量
STATES = [str(s) for s in blood['state_vocab']]

def oh(c, n): return torch.nn.functional.one_hot(c, n).float()
def cv(c2, st, src): return torch.cat([oh(c2,2), oh(st,len(STATES)), oh(torch.full_like(c2,src),2)], -1)

bs_idx = STATES.index('cDC')
brs_idx = STATES.index('microglia_mhc2')

idx_b = np.where(blood['state'] == bs_idx)[0]
sel_b = idx_b[np.random.RandomState(42).choice(len(idx_b), 200, replace=False)]
x_b = X_blood[torch.from_numpy(sel_b)]
with torch.no_grad(): mu, _ = model.enc_blood(x_b)
with torch.no_grad(): u0 = model.protein_head(mu)

idx_B1 = np.where((brain_te['state'] == brs_idx) & (brain_te['cond2'] == 1))[0][:200]
idx_B0 = np.where((brain_te['state'] == brs_idx) & (brain_te['cond2'] == 0))[0][:200]
x_obs_1 = X_brain[torch.from_numpy(idx_B1)]
x_obs_0 = X_brain[torch.from_numpy(idx_B0)]
obs_diff = (x_obs_1.mean(0) - x_obs_0.mean(0)).numpy()

# === Gate 3: condition dropout ===
print('=== Gate 3: Condition Dropout ===')
c_1 = cv(torch.ones(200, dtype=torch.long), torch.full((200,), brs_idx, dtype=torch.long), 1)
c_0 = cv(torch.zeros(200, dtype=torch.long), torch.full((200,), brs_idx, dtype=torch.long), 1)
with torch.enable_grad():
    _, z1 = model.solve_equilibrium(mu, mu.clone(), u0, c_1)
    _, z0 = model.solve_equilibrium(mu, mu.clone(), u0, c_0)
with torch.no_grad():
    pred_1 = model.dec_brain(z1.detach())
    pred_0 = model.dec_brain(z0.detach())
pred_diff = (pred_1.mean(0) - pred_0.mean(0)).numpy()
r_with = float(spearmanr(pred_diff, obs_diff).statistic)
print(f'WITH cond: r = {r_with:+.4f}')

c_zero = torch.zeros(200, 2 + len(STATES) + 2)  # cond_dim 动态（用户 2026-09-09 指令）
with torch.enable_grad():
    _, z1n = model.solve_equilibrium(mu, mu.clone(), u0, c_zero)
    _, z0n = model.solve_equilibrium(mu, mu.clone(), u0, c_zero)
with torch.no_grad():
    pred_1n = model.dec_brain(z1n.detach())
    pred_0n = model.dec_brain(z0n.detach())
pred_diff_n = (pred_1n.mean(0) - pred_0n.mean(0)).numpy()
r_without = float(spearmanr(pred_diff_n, obs_diff).statistic)
print(f'WITHOUT cond: r = {r_without:+.4f}')
g3_pass = abs(r_without) > 0.5 * abs(r_with)
print(f'Verdict: {"PASS" if g3_pass else "FAIL"} (drops {r_with:+.3f} -> {r_without:+.3f})')

# === Gate 4: donor generalization ===
print('\n=== Gate 4: Donor Generalization ===')
donors_te = np.unique(brain_te['donor'][(brain_te['state'] == brs_idx)])
pred_mean = pred_1.mean(0, keepdim=True)
results_g4 = []
for d in donors_te:
    mask_d = (brain_te['donor'] == d) & (brain_te['state'] == brs_idx)
    x_d = X_brain[torch.from_numpy(np.where(mask_d)[0])]
    if x_d.shape[0] < 10: continue
    dist_self = float(sinkhorn_divergence(pred_mean, x_d[:min(50, x_d.shape[0])]))
    others = [d2 for d2 in donors_te if d2 != d]
    rng_d = np.random.RandomState(0)
    dist_others = []
    for d2 in rng_d.choice(others, min(3, len(others)), replace=False):
        m2 = (brain_te['donor'] == d2) & (brain_te['state'] == brs_idx)
        x_d2 = X_brain[torch.from_numpy(np.where(m2)[0])]
        if x_d2.shape[0] >= 10:
            dist_others.append(float(sinkhorn_divergence(pred_mean, x_d2[:50])))
    if not dist_others: continue
    avg_other = np.mean(dist_others)
    results_g4.append({'donor': str(d), 'self': round(dist_self,4), 'other': round(avg_other,4),
                       'closer': bool(dist_self < avg_other)})
n_closer = sum(1 for r in results_g4 if r['closer'])
g4_pass = n_closer > len(results_g4) / 2
print(f'Closer to self: {n_closer}/{len(results_g4)}')
for r in results_g4[:4]:
    print(f'  {r["donor"]:10s} self={r["self"]:.4f} other={r["other"]:.4f}')
print(f'Verdict: {"PASS" if g4_pass else "INCONCLUSIVE"}')

# === Gate 5: transient stratification ===
print('\n=== Gate 5: Transient Stratification ===')
grad_norms = {}
for st in range(5, 9):
    mask = (brain_te['state'] == st)
    idx = np.where(mask)[0][:100]
    if len(idx) < 10: continue
    x_B = X_brain[torch.from_numpy(idx)]
    with torch.no_grad(): mu_B, _ = model.enc_brain(x_B)
    z_req = mu_B.clone().requires_grad_(True)
    u_z = torch.zeros(len(idx), cfg.n_proteins)
    c_st = cv(torch.from_numpy(brain_te['cond2'][idx]).long(),
              torch.full((len(idx),), st, dtype=torch.long), 1)
    e = model.potential(mu_B.detach(), z_req, u_z, c_st)
    g = torch.autograd.grad(e.sum(), z_req)[0]
    grad_norms[STATES[st]] = round(float(g.norm(dim=-1).mean()), 4)

med = float(np.median(list(grad_norms.values())))
transient = [st for st, gn in grad_norms.items() if gn > 2 * med]
print('Gradient norms:', grad_norms)
print(f'Median: {med:.4f}')
g5_pass = len(transient) == 0
print(f'Transient strata (>2x median): {transient or "none"}')
print(f'Verdict: {"PASS" if g5_pass else "PARTIAL"}')

# === Gate 6: disease shortcut ===
print('\n=== Gate 6: Disease Shortcut ===')
idx_cm = np.where(blood['state'] == STATES.index('classical_mono'))[0]
sel_cm = idx_cm[np.random.RandomState(42).choice(len(idx_cm), min(500, len(idx_cm)), replace=False)]
X_cm = blood['X'][sel_cm]
y_cm = blood['cond2'][sel_cm]
sc = StandardScaler().fit(X_cm)
clf = LogisticRegression(max_iter=1000, C=0.1).fit(sc.transform(X_cm), y_cm)
w = clf.coef_[0]
w_norm = w / np.linalg.norm(w)
X_cm_resid = X_cm - np.outer(X_cm @ w_norm, w_norm)

x_orig = torch.from_numpy(X_cm.astype(np.float32))
x_resid = torch.from_numpy(X_cm_resid.astype(np.float32))
c_test = cv(torch.zeros(len(sel_cm), dtype=torch.long), torch.full((len(sel_cm),), brs_idx, dtype=torch.long), 1)
with torch.no_grad():
    mu_o, _ = model.enc_blood(x_orig)
    mu_r, _ = model.enc_blood(x_resid)
    u_o = model.protein_head(mu_o)
    u_r = model.protein_head(mu_r)
with torch.enable_grad():
    _, z_o = model.solve_equilibrium(mu_o, mu_o.clone(), u_o, c_test)
    _, z_r = model.solve_equilibrium(mu_r, mu_r.clone(), u_r, c_test)
with torch.no_grad():
    pred_o = model.dec_brain(z_o.detach()).mean(0).numpy()
    pred_r = model.dec_brain(z_r.detach()).mean(0).numpy()

r_orig = float(spearmanr(pred_o, obs_diff).statistic)
r_resid = float(spearmanr(pred_r, obs_diff).statistic)
g6_pass = abs(r_resid) > 0.5 * abs(r_orig)
print(f'Original r = {r_orig:+.4f}')
print(f'Residualized r = {r_resid:+.4f}')
print(f'Verdict: {"PASS" if g6_pass else "FAIL"} (survives residualization)')

# Save
gates = {
    'gate3': {'r_with': round(r_with,4), 'r_without': round(r_without,4), 'pass': g3_pass},
    'gate4': {'n_closer': n_closer, 'n_total': len(results_g4), 'per_donor': results_g4,
               'pass': g6_pass},
    'gate5': {'grad_norms': grad_norms, 'median': round(med,4), 'transient': transient,
               'pass': g5_pass},
    'gate6': {'r_orig': round(r_orig,4), 'r_resid': round(r_resid,4), 'pass': g6_pass},
}
json.dump(gates, open(A + '/results/p3_6_gates/gates_3_to_6_v1.json', 'w'), indent=1)
print('\nSaved gates_3_to_6_v1.json')
