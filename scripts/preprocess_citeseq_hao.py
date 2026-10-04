#!/usr/bin/env python3
"""preprocess_citeseq_hao: CITE-seq 配对数据 → CCWM Block B 桥训练集 (公开数据, 一键重建)。

edit 自 ad_omics_v2 同名脚本（基底已删，git 历史可溯）（用户 2026-09-09
指令：项目自持、无外部流入；差异仅路径 + 保护基因边表来源）。

数据源 (公开可获取):
  Hao et al. 2021 Cell CITE-seq PBMC 图谱, GEO GSE164378
  (RNA + 228 种表面蛋白 ADT, 同细胞配对, 3P 批 ~16 万细胞)。
  本脚本只依赖已下载的原始文件 (pd_product_assets/raw/citeseq_hao/GSE164378_RAW.tar +
  sc.meta.data_3P.csv.gz); 缺失时给出 GEO 下载地址, 不自动联网。

处理 (全部确定性, 无随机):
  1. tar 解包 → 读 10x mtx (RNA_3P / ADT_3P), 细胞 barcode 取交集
  2. 细胞类型: 元数据 celltype.l1 → 项目血端标签
     (CD4 T→CD4T, CD8 T→CD8T, B→B, Mono→Mono, NK→NK; 其余→other)
  3. RNA: counts→CPM1e4→log1p; 特征 = 2000 HVG + 保护基因
     (pd_product 自己的边表基因全集: pd_edges_forward.csv 边靶 +
      pd_edges_reverse.csv 边源, 由 p4 剂量-效应与 randko 读出表生成) → 逐基因 z-score
  4. ADT: 逐细胞 CLR (log1p 几何均值化)
  5. 附带保存匹配用表达参考: 逐 (细胞类型 × 基因) log1pCPM 均值,
     逐 (细胞类型 × 蛋白) CLR 均值

输出 (pd_product_assets/processed/citeseq_hao/):
  bridge_data.npz   X(细胞×基因 z-score), Y(细胞×蛋白 CLR), ct(整数码),
                    genes, proteins, ct_labels, expr_mean_ct, adt_mean_ct
裸启动: PYTHONPATH=src python scripts/preprocess_citeseq_hao.py
"""
import os
import gzip
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix

PROJECT = Path("/public/home/mengxl/dzy/pd_product_assets")
DATA = PROJECT / "raw/citeseq_hao"
RAW_TAR = DATA / "GSE164378_RAW.tar"
META = DATA / "GSE164378_sc.meta.data_3P.csv.gz"
RAW_DIR = DATA / "raw_3p"
OUT = PROJECT / "processed/citeseq_hao/bridge_data.npz"
GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/GSE164378/suppl/"

N_HVG = 2000          # 与 GenKI 口径一致的特征数
MIN_COUNTS = 500      # 轻 QC (图谱本身已质控)

CT_MAP = {"CD4 T": "CD4T", "CD8 T": "CD8T", "B": "B", "Mono": "Mono", "NK": "NK"}
CT_LABELS = ["CD4T", "CD8T", "B", "Mono", "NK", "other"]


def log(s):
    print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)


def extract_raw():
    if RAW_DIR.exists() and any(RAW_DIR.iterdir()):
        return
    if not RAW_TAR.exists() or not META.exists():
        raise SystemExit(f"缺原始文件: 请从 GEO 下载到 {DATA}/\n  {GEO_URL}GSE164378_RAW.tar\n"
                         f"  {GEO_URL}GSE164378_sc.meta.data_3P.csv.gz")
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    log(f"解包 {RAW_TAR.name} (1.5GB) ...")
    with tarfile.open(RAW_TAR) as tf:
        for m in tf.getmembers():
            if ("_3P-" in m.name) or m.name.startswith("GSE164378_sc.meta"):
                tf.extract(m, RAW_DIR)


def read_10x(prefix):
    """读 10x 三件套 → (csr cells×features, features 名, barcodes)。"""
    mtx = sorted(RAW_DIR.glob(f"{prefix}-matrix.mtx.gz"))
    feat = sorted(RAW_DIR.glob(f"{prefix}-features.tsv.gz"))
    bc = sorted(RAW_DIR.glob(f"{prefix}-barcodes.tsv.gz"))
    assert mtx and feat and bc, f"{prefix} 三件套缺失: {list(RAW_DIR.iterdir())}"
    log(f"  读 {prefix} ({mtx[0].stat().st_size/1e6:.0f} MB) ...")
    with gzip.open(mtx[0], "rb") as fh:
        M = mmread(fh)
    M = csr_matrix(M).T.tocsr()  # 10x = features×barcodes → 转置
    feats = pd.read_csv(feat[0], sep="\t", header=None, compression="gzip")
    bcs = pd.read_csv(bc[0], sep="\t", header=None, compression="gzip")[0].astype(str).values
    return M, feats, bcs


def main():
    extract_raw()

    # ---- 1. 读 RNA + ADT ----
    log("[1/5] 读 RNA/ADT 矩阵")
    X_rna, rna_feat, rna_bc = read_10x("GSM5008737_RNA_3P")
    Y_adt, adt_feat, adt_bc = read_10x("GSM5008738_ADT_3P")
    rna_sym = rna_feat[1].astype(str).values if rna_feat.shape[1] > 1 else rna_feat[0].astype(str).values
    adt_names = adt_feat[1].astype(str).values if adt_feat.shape[1] > 1 else adt_feat[0].astype(str).values
    log(f"  RNA: {X_rna.shape[0]} 细胞 × {X_rna.shape[1]} 基因; ADT: {Y_adt.shape[0]} 细胞 × {Y_adt.shape[1]} 蛋白")

    # 细胞交集 (RNA∩ADT∩元数据)
    log("[2/5] 细胞交集 + 细胞类型映射")
    meta = pd.read_csv(META, compression="gzip", index_col=0)
    meta.index = meta.index.astype(str)
    ct_l1 = meta["celltype.l1"].astype(str)
    ct_map = ct_l1.map(lambda s: CT_MAP.get(s, "other"))
    # 供体/批次列 (供交叉核对泛化; Hao 图谱为 donor, 缺失则回退 lane)
    donor_col = next((c for c in meta.columns if "donor" in c.lower()), None) or \
        next((c for c in meta.columns if "lane" in c.lower() or "batch" in c.lower()), None)
    donor = meta[donor_col].astype(str) if donor_col else pd.Series("all", index=meta.index)
    log(f"  供体列: {donor_col}, 水平数 {donor.nunique()}")
    common = sorted(set(rna_bc) & set(adt_bc) & set(ct_map.index))
    log(f"  三方交集细胞: {len(common)}")
    rna_pos = pd.Series(range(len(rna_bc)), index=rna_bc)
    adt_pos = pd.Series(range(len(adt_bc)), index=adt_bc)
    X_rna = X_rna[rna_pos[common].values]
    Y_adt = Y_adt[adt_pos[common].values]
    ct = ct_map[common].values
    dn = donor[common].values

    # QC: 总计数过滤
    lib = np.asarray(X_rna.sum(axis=1)).ravel()
    keep = lib >= MIN_COUNTS
    X_rna, Y_adt, ct, dn = X_rna[keep], Y_adt[keep], ct[keep], dn[keep]
    log(f"  QC(>={MIN_COUNTS} counts) 后: {X_rna.shape[0]} 细胞; 类型分布: "
        f"{pd.Series(ct).value_counts().to_dict()}")

    # ---- 2. RNA 归一化 + 特征选择 ----
    log("[3/5] RNA log1p(CPM) + 2000 HVG + 保护基因")
    lib = np.asarray(X_rna.sum(axis=1)).ravel()
    lib[lib == 0] = 1
    Xn = X_rna.multiply(1e4 / lib[:, None]).tocsr()  # multiply 产 COO, 转 CSR 才能列切
    Xn.data = np.log1p(Xn.data)
    mean = np.asarray(Xn.mean(axis=0)).ravel()
    var = np.asarray(Xn.power(2).mean(axis=0)).ravel() - mean ** 2
    disp = var / np.maximum(mean, 1e-8)  # 离散度排序选 HVG (seurat 思路, 确定性)
    hvg = set(np.argsort(-disp)[:N_HVG].tolist())

    protected = set()
    for d in ("forward", "reverse"):
        f = DATA / f"pd_edges_{d}.csv"
        e = pd.read_csv(f)
        col = "target" if d == "forward" else "source"
        protected |= set(e[col].astype(str))
    sym2idx = {}
    for i, s in enumerate(rna_sym):
        sym2idx.setdefault(s, i)  # 重复 symbol 取第一个
    prot_idx = sorted(sym2idx[s] for s in protected if s in sym2idx)
    sel = sorted(hvg | set(prot_idx))
    log(f"  特征: HVG {len(hvg)} + 保护命中 {len(prot_idx)} → 合计 {len(sel)} "
        f"(保护基因 {len(protected)} 个, 命中 {len(prot_idx)})")
    Xn = Xn[:, sel].toarray().astype(np.float32)
    genes = rna_sym[sel]

    # z-score (逐基因, 存均值方差供反演)
    g_mu = Xn.mean(axis=0)
    g_sd = Xn.std(axis=0)
    g_sd[g_sd < 1e-8] = 1.0
    X = ((Xn - g_mu) / g_sd).astype(np.float32)

    # ---- 3. ADT CLR ----
    log("[4/5] ADT CLR 标准化")
    Yl = np.log1p(Y_adt.toarray().astype(np.float32))
    Y = (Yl - Yl.mean(axis=1, keepdims=True)).astype(np.float32)

    # ---- 4. 匹配用参考 (逐细胞类型均值; VK 表达/丰度匹配用) ----
    ct_codes = np.array([CT_LABELS.index(c) for c in ct], dtype=np.int64)
    expr_mean_ct = np.zeros((len(CT_LABELS), len(genes)), dtype=np.float32)
    adt_mean_ct = np.zeros((len(CT_LABELS), Y.shape[1]), dtype=np.float32)
    for k in range(len(CT_LABELS)):
        m = ct_codes == k
        if m.sum():
            expr_mean_ct[k] = Xn[m].mean(axis=0)   # log1pCPM 均值 (未 z-score)
            adt_mean_ct[k] = Y[m].mean(axis=0)     # CLR 均值

    # ---- 5. 保存 ----
    np.savez_compressed(
        OUT, X=X, Y=Y, ct=ct_codes,
        genes=np.asarray(genes, dtype=str), proteins=np.asarray(adt_names, dtype=str),
        ct_labels=np.array(CT_LABELS, dtype=str), donor=np.asarray(dn, dtype=str),
        expr_mean_ct=expr_mean_ct, adt_mean_ct=adt_mean_ct,
        gene_mu=g_mu.astype(np.float32), gene_sd=g_sd.astype(np.float32))
    log(f"[5/5] 已保存 {OUT}: X{X.shape}, Y{Y.shape}, 类型 {len(CT_LABELS)}")


if __name__ == "__main__":
    main()
