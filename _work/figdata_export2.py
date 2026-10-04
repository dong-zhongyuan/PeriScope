# figdata_export2.py — 服务器端补导出 #2(本地写好 scp 上来运行)
# 4e: Kamath 逐基因 Welch t(PD vs Ctrl,基因=全部 GSEA 显著程序集并集)
# 1f: 血端 BlockA X 按 state 子采样 → PCA → UMAP 坐标
# 运行:/public/home/mengxl/dzy/envs/pd_bbm/bin/python figdata_export2.py
import json
import os

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.stats import ttest_ind

A = '/public/home/mengxl/dzy/pd_product_assets'
OUT = '/public/home/mengxl/dzy/pd_product/_work/figdata'
os.makedirs(OUT, exist_ok=True)

# ---------- 4e ----------
gsea = json.load(open(f'{A}/results/p5_gsea/gsea_stable_v1.json'))
prog_genes = {k: sorted({g for t in v for g in t['genes']}) for k, v in gsea.items()}
allg = sorted({g for v in prog_genes.values() for g in v})

ad = sc.read_h5ad(f'{A}/processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad',
                  backed='r')
gidx = {g: i for i, g in enumerate(ad.var['gene_symbol'].astype(str).values)}
present = [g for g in allg if g in gidx]
cond = ad.obs['condition'].astype(str).values
print('condition values:', pd.unique(cond))
pd_mask = np.array([('pd' in v.lower() or 'disease' in v.lower()) for v in cond])

sub = ad[:, [gidx[g] for g in present]]
X = sub.X
Xa = np.asarray(X.todense(), dtype=np.float32) if sp.issparse(X) else np.asarray(X, dtype=np.float32)
rows = []
for j, g in enumerate(present):
    x1, x0 = Xa[pd_mask, j], Xa[~pd_mask, j]
    t, p = ttest_ind(x1, x0, equal_var=False)
    rows.append(dict(gene=g, t=float(t), p=float(p),
                     mean_pd=float(x1.mean()), mean_ctrl=float(x0.mean()),
                     n_pd=int(pd_mask.sum()), n_ctrl=int((~pd_mask).sum())))
dfe = pd.DataFrame(rows)
dfe['programs'] = dfe['gene'].map(
    {g: ';'.join(k for k, v in prog_genes.items() if g in v) for g in present})
dfe.to_csv(f'{OUT}/fig4e_spatial_gene.csv', index=False)
print('4e done:', len(dfe), 'genes,', int(pd_mask.sum()), 'PD /', int((~pd_mask).sum()), 'ctrl cells')

# ---------- 1f ----------
z = np.load(f'{A}/interim/ccwm_v23/blood.npz', allow_pickle=True)
Xb, st, coh = z['X'], z['state'].astype(str), z['cohort'].astype(str)
rng = np.random.default_rng(0)
sel = np.array(sorted(np.concatenate(
    [rng.choice(np.where(st == s)[0], min(400, (st == s).sum()), replace=False)
     for s in np.unique(st)])))
ad2 = sc.AnnData(Xb[sel].astype('float32'))
ad2.obs['state'] = st[sel]
ad2.obs['cohort'] = coh[sel]
sc.pp.pca(ad2, n_comps=30, svd_solver='randomized')
sc.pp.neighbors(ad2, n_neighbors=15)
sc.tl.umap(ad2)
dfu = pd.DataFrame(dict(umap1=ad2.obsm['X_umap'][:, 0],
                        umap2=ad2.obsm['X_umap'][:, 1],
                        state=ad2.obs['state'].values,
                        cohort=ad2.obs['cohort'].values))
dfu.to_csv(f'{OUT}/fig1f_umap.csv', index=False)
print('1f done:', len(dfu), 'cells')

# ---------- 1d/1h: donor x state composition + state marker genes ----------
zb = np.load(f'{A}/interim/ccwm_v23/blood.npz', allow_pickle=True)
Xb, stb, cohb, donb = zb['X'], zb['state'].astype(str), zb['cohort'].astype(str), zb['donor'].astype(str)
vocab = [str(s) for s in zb['state_vocab']]
rows = []
for d in pd.unique(donb):
    m = donb == d
    for s in pd.unique(stb[m]):
        mm = m & (stb == s)
        rows.append(dict(donor=d, cohort=cohb[mm][0], state=vocab[int(s)] if s.isdigit() else s,
                         state_id=int(s), n=int(mm.sum())))
pd.DataFrame(rows).to_csv(f'{OUT}/fig1_donor_state.csv', index=False)
print('1d donors:', len(pd.unique(donb)))

Xf = np.asarray(Xb, dtype=np.float32)
rows = []
for sid in range(len(vocab)):
    m = stb == str(sid)
    if m.sum() < 50:
        continue
    mu_s = Xf[m].mean(0)
    mu_o = Xf[~m].mean(0)
    diff = mu_s - mu_o
    top = np.argsort(-np.abs(diff))[:6]
    for j in top:
        rows.append(dict(state=vocab[sid], state_id=sid, gene_idx=int(j),
                         mean_state=float(mu_s[j]), mean_other=float(mu_o[j]),
                         diff=float(diff[j]), pct_pos_state=float((Xf[m, j] > 0).mean())))
mk = pd.DataFrame(rows)
# gene symbol: 1469 词表序 = sorted(brain∩blood∩cite)，与 fig3_curves 基因序一致
cur = pd.read_csv(f'{OUT}/fig3_curves.csv')
genes = sorted(cur.gene.unique())
mk['gene'] = mk.gene_idx.map(dict(enumerate(genes)))
mk.to_csv(f'{OUT}/fig1_markers.csv', index=False)
print('1h markers:', len(mk), 'genes vocab match:', len(genes))
