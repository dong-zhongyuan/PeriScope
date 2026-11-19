"""Shared definitions for the corrected benchmark. Everything is recomputed from pb18 (donor x state pseudobulks) and stored per-seed predictions."""
import sys, json, pickle, itertools
sys.dont_write_bytecode = True
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import rankdata

OUT = Path('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/benchmark_20261119'); SRC = OUT / 'src'
RES = Path('/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/research_20261118')
BENCH = Path('/public/home/mengxl/dzy/pd_product_assets/results/benchmark_20261002')
ST = ['astro', 'microglia_mhc2']; BLOOD = ['cDC', 'classical_mono', 'nonclassical_mono', 'pDC']; AX = [f'{b}__{s}' for b in BLOOD for s in ST]
TP = ['SN-2142', 'SN-4568', 'SN-4775']; TC = ['SN-3298', 'SN-3322']; TEST = TP + TC
VALP = ['SN-2569', 'SN-4956']; CALP = ['SN-2561', 'SN-3345']
META = pd.read_csv(RES / 'donor_cv/pb18_meta.csv'); _Z = np.load(RES / 'donor_cv/pb18.npz'); MU = _Z['mean'].astype(float); GENES = _Z['genes'].astype(str); NG = MU.shape[1]
KEY = {(r.donor, r.state): i for i, r in META.iterrows()}; COND = dict(zip(META.donor, META.cond)); NCELL = {(r.donor, r.state): int(r.n) for r in META.itertuples()}
SPLIT = META.drop_duplicates('donor').set_index('donor').split.to_dict()
FIT9 = sorted(d for d, s in SPLIT.items() if s == 'train'); D13 = sorted(d for d in SPLIT if d not in TEST); ALL18 = sorted(SPLIT)
assert len(FIT9) == 9 and len(D13) == 13 and set(D13) == set(FIT9) | set(VALP) | set(CALP)

def rz(M):
    """rank-normalise along the last axis: dot products of two such vectors are Spearman correlations"""
    R = rankdata(M, axis=-1); R = R - R.mean(-1, keepdims=True); return R / np.linalg.norm(R, axis=-1, keepdims=True)
def Z(v):
    v = np.asarray(v, float); return (v - v.mean()) / v.std()
def A(ds, s): return np.stack([MU[KEY[(d, s)]] for d in ds])
