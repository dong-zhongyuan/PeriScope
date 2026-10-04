"""Brain competitive validation using discovery_union_v1 programs; 2000 expression-matched random sets."""
import hashlib
from datetime import date
import json
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/public/home/mengxl/dzy/pd_product/src")
from pdproduct.analysis.spatial_overlap import colocalization_test, roi_lognorm

A = "/public/home/mengxl/dzy/pd_product_assets"
H5AD = "/public/home/mengxl/dzy/pd_product_assets/processed/gse253975/v0.1/GSE253975_geomx.h5ad"
MOESM11 = "/public/home/mengxl/dzy/pd_product_assets/raw/ma2025_zenodo/Supplementary Data/41467_2025_62478_MOESM11_ESM.xlsx"
MOESM11_SHEET_A_XML = "xl/worksheets/sheet1.xml"  # SN PD-vs-Control DE, 193 rows (probed 2026-09-04)
KAMATH_H5AD = "/public/home/mengxl/dzy/pd_product_assets/processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad"
V2_JSON = Path(A) / "results" / "p6_spatial" / "spatial_validation_v2.json"

REGISTRY = "/public/home/mengxl/dzy/pd_product_assets/results/figure_program_registry/programs.json"
N_PERM = 2000
SEED = 42
ASTRO_MARKERS = ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"]   # Ma-leg astro-context scoring
MA_ASTRO_TOPFRAC = 0.5   # keep top-50% astro-score ROIs per donor (context matching to Kamath cell-typing)
EXPR_BINS = 20
MIN_MARKER_GENES = 3
CAND_MIN_PRESENT = 4
CROSS_EDGE_MIN_COMBOS = 3

MARKERS = {
    "Microglia_MHC2": ["HLA-DRA", "HLA-DRB1", "CD74", "HLA-DPA1", "HLA-DPB1", "HLA-DQB1",
                       "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DQB2"],
    "Microglia_homeostatic": ["P2RY12", "CX3CR1", "TMEM119", "SALL1", "P2RY13", "GPR34", "SIGLECH"],
    "Astrocyte": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1"],
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MOBP"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "FLT1"],
    "Pericyte": ["PDGFRB", "RGS5", "CSPG4", "NOTCH3"],
}
PRIMARY_BY_BRAIN_STATE = {
    "microglia_mhc2": "Microglia_MHC2",
    "microglia_homeostatic": "Microglia_homeostatic",
    "astro": "Astrocyte",
}
CROSS_EDGE_PRIMARY = {
    "GO:0019886": "Microglia_MHC2",
    "GO:0002495": "Microglia_MHC2",
    "GO:0002478": "Microglia_MHC2",
    "GO:0016064": "Microglia_MHC2",
    "GO:0001774": "Microglia_MHC2",
    "GO:0010594": "Endothelial",
    "GO:0042063": "Astrocyte",
}
CONTEXT_PANELS = ["Oligodendrocyte", "Endothelial", "Pericyte"]


def term_key(term: str) -> str:
    m = re.search(r"\((GO:\d+)\)", term)
    return m.group(1) if m else term


# ---------- shared vocab + dose-curve deltas (same reconstruction as p4) ----------
def load_vocab_and_deltas():
    import scanpy as sc
    brain = sc.read_h5ad("/public/home/mengxl/dzy/pd_product_assets/processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad", backed="r")
    g_brain = brain.var["gene_symbol"].astype(str).values
    brain.file.close()
    bmb = sc.read_h5ad(A + "/processed/moquin_beaudry2025/v0.1/MB2025_pbmc_annotated.h5ad", backed="r")
    g_blood = bmb.var["gene_symbol"].astype(str).values
    bmb.file.close()
    hao = np.load("/public/home/mengxl/dzy/pd_product_assets/processed/citeseq_hao/bridge_data.npz", allow_pickle=True)
    g_cite = hao["genes"].astype(str)
    genes = sorted(set(g_brain) & set(g_blood) & set(g_cite))
    gidx = {g: i for i, g in enumerate(genes)}
    dose = json.load(open(A + "/results/p4_dose_response/dose_response_stable_v1.json"))
    delta = {}
    for edge, rec in dose.items():
        for p, cv in rec["curves"].items():
            xm = np.asarray(cv["x_mean"])
            delta[(p, edge)] = xm[-1] - xm[0]
    return genes, gidx, delta


def load_candidate_sets(panel, gidx=None, delta=None):
    programs = json.load(open(REGISTRY))
    return {'gsea::'+c: dict(source='discovery_union_v1', kind='per_combo',
            protein=c.split('__',1)[0], edge=c.split('__',1)[1],
            brain_state=c.split('__x__')[1], genes=genes)
            for c,genes in programs.items()}


def read_moesm11_sheet(sheet_xml):
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    M = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(MOESM11)
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    ss = ["".join(t.text or "" for t in si.iter(M + "t")) for si in root.findall("m:si", ns)]
    root = ET.fromstring(z.read(sheet_xml))
    rows = []
    for row in root.iter(M + "row"):
        vals = {}
        for c in row.findall("m:c", ns):
            col = "".join(ch for ch in c.get("r") if ch.isalpha())
            v = c.find("m:v", ns)
            if v is None:
                continue
            vals[col] = ss[int(v.text)] if c.get("t") == "s" else v.text
        rows.append(vals)
    hdr = rows[0]
    out = {}
    for r in rows[1:]:
        rec = {hdr.get(k, k): v for k, v in r.items()}
        if "Gene" in rec and rec.get("log2FC") is not None:
            out[rec["Gene"]] = {"log2FC": float(rec["log2FC"]),
                                "p_val": float(rec["p_val"]),
                                "p_val_adj": float(rec["p_val_adj"])}
    return out


def bh_qvalues(pvals):
    items = [(k, p) for k, p in pvals.items() if p is not None]
    m = len(items)
    order = sorted(items, key=lambda x: x[1])
    q, prev = {}, 1.0
    for rank, (k, p) in reversed(list(enumerate(order, start=1))):
        prev = min(prev, p * m / rank)
        q[k] = prev
    return q


KAMATH_SPLIT = "/public/home/mengxl/dzy/pd_product_assets/interim/v0.1/splits/GSE178265_sn_donor_split_v1.json"
MIN_CELLS_PER_DONOR = 25
MIN_DONORS_PER_ARM = 3
BSTATE_TO_CELLTYPE = {"microglia_mhc2": "Microglia",
                      "microglia_homeostatic": "Microglia", "astro": "Astrocyte"}
PANEL_TO_CELLTYPE = {"Microglia_MHC2": "Microglia",
                     "Microglia_homeostatic": "Microglia",
                     "Astrocyte": "Astrocyte", "Endothelial": "Endothelial",
                     "Pericyte": "Pericyte"}


def crosscohort_celltype(s):
    """Frozen Amendment-4 §4 rule: annotated cell type scored for each set."""
    if s["kind"] == "per_combo":
        return BSTATE_TO_CELLTYPE[s["brain_state"]]
    panel = CROSS_EDGE_PRIMARY.get(s["term_key"])
    if panel is not None:
        return PANEL_TO_CELLTYPE[panel]
    from collections import Counter
    bs = Counter(c.split("__x__")[1] for c in s["combos"])
    return BSTATE_TO_CELLTYPE[bs.most_common(1)[0][0]]




def _competitive_test(M, is_pd, set_cols, n_perm=N_PERM, seed=SEED, n_bins=20):
    """Mean Welch-t of the set vs expression-matched random-set null (N_PERM draws)."""
    t = _welch_t(M, is_pd)
    ok = np.isfinite(t)
    expr = M.mean(axis=0)
    bg = np.where(ok)[0]
    qs = np.quantile(expr[bg], np.linspace(0, 1, n_bins + 1))
    bin_idx = np.digitize(expr[bg], qs[1:-1])
    bins = [bg[bin_idx == b] for b in range(n_bins)]
    cols = np.array([c for c in set_cols if ok[c]])
    if len(cols) < 4:
        return {"status": "insufficient_genes", "n_used": len(cols)}
    obs = float(t[cols].mean())
    comp = {}
    for c in cols:
        b = int(np.digitize(expr[c], qs[1:-1]))
        comp[b] = comp.get(b, 0) + 1
    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        draw = [rng.choice(bins[b], size=min(k, len(bins[b])), replace=False)
                for b, k in comp.items() if len(bins[b])]
        null[i] = t[np.concatenate(draw)].mean()
    p_up = float((1 + (null >= obs).sum()) / (1 + n_perm))
    p_down = float((1 + (null <= obs).sum()) / (1 + n_perm))
    return {"status": "ok", "stat": obs, "p_up": p_up, "p_down": p_down,
            "n_genes_used": len(cols),
            "null_mean": float(null.mean()), "null_sd": float(null.std())}


def _welch_t(mat, m_pd):
    """Per-gene Welch t (PD vs Ctrl) for (obs x genes) matrix."""
    a, b = mat[m_pd], mat[~m_pd]
    se = np.sqrt(a.var(axis=0, ddof=1) / a.shape[0]
                 + b.var(axis=0, ddof=1) / b.shape[0])
    return (a.mean(axis=0) - b.mean(axis=0)) / se




def main_testb3():
    """Amendment 9 Part 2: competitive Test B on the 19 per-combo programs."""
    import scanpy as sc

    de = read_moesm11_sheet(MOESM11_SHEET_A_XML)
    made = {"MaDE_up": (sorted(g for g, r in de.items()
                               if r["log2FC"] > 0 and r["p_val_adj"] < 0.05), "p_up"),
            "MaDE_down": (sorted(g for g, r in de.items()
                                 if r["log2FC"] < 0 and r["p_val_adj"] < 0.05), "p_down")}
    out = {"metadata": {
        "date": str(date.today()), "registry": REGISTRY, "registry_sha256": hashlib.sha256(Path(REGISTRY).read_bytes()).hexdigest(), "n_permutations": N_PERM, "n_programs": len(json.load(open(REGISTRY))),
        "version": "discovery_union_v1",
        "statistic": "donor-level log1p-CPM pseudobulk Welch t; set stat = mean t",
        "null": "2000 expression-matched (20-bin) random gene sets, seed 42",
        "ma_leg": "astro-context restricted: top-50% astro-score ROIs per donor (GFAP/AQP4/SLC1A3/ALDH1L1 z-mean), 2026-09-18 optimization, disclosed",
        "wording": "competitive claim only: stronger than expression-matched background; never 'non-zero'",
        "claim_boundary": "C3 brain-side compatibility only; association-level",
    }}

    # ---- Ma leg (donor means of ROI log1p-CPM, full panel) ----
    adata = ad.read_h5ad(H5AD)
    vn = [str(g) for g in adata.var_names]
    gene_index = {g: i for i, g in enumerate(vn)}
    gidx, delta = None, None
    sets_ma = {k: v for k, v in load_candidate_sets(set(vn), gidx, delta).items()
               if v["kind"] == "per_combo"}
    L = roi_lognorm(adata)
    donor = adata.obs["donor"].values.astype(str)
    cond = adata.obs["condition"].values.astype(str)
    donor_cond_ma = {d: cond[donor == d][0] for d in set(donor.tolist())}
    ds_ma = sorted(donor_cond_ma)
    # astro-context restriction (2026-09-18): top-50% astro-score ROIs per donor
    Zs = (L - L.mean(axis=0)) / (L.std(axis=0) + 1e-9)
    ai = [gene_index[g] for g in ASTRO_MARKERS if g in gene_index]
    astro_sc = Zs[:, ai].mean(axis=1)
    del Zs
    M_ma = np.stack([L[(donor == d) & (astro_sc >= np.quantile(astro_sc[donor == d], 1 - MA_ASTRO_TOPFRAC))].mean(axis=0)
                     for d in ds_ma])
    is_pd_ma = np.array([donor_cond_ma[d] == "PD" for d in ds_ma])
    del L
    ma_res = {}
    for sid, s in sets_ma.items():
        r = _competitive_test(M_ma, is_pd_ma,
                              [gene_index[g] for g in s["genes"] if g in gene_index])
        r["n_present"] = len([g for g in s["genes"] if g in gene_index])
        r.update(n_donors=len(ds_ma), n_PD=int(is_pd_ma.sum()), n_Control=int((~is_pd_ma).sum()))
        ma_res[sid] = r
    gate = {name: _competitive_test(M_ma, is_pd_ma,
                                    [gene_index[g] for g in genes if g in gene_index])
            for name, (genes, _) in made.items()}
    gate_pass = all(gate[n].get(made[n][1]) is not None and gate[n][made[n][1]] < 0.05
                    for n in made)
    print(f"Ma leg done; gate: {gate_pass} "
          f"(up {gate['MaDE_up'].get('p_up')}, down {gate['MaDE_down'].get('p_down')})")
    out["metadata"]["positive_control_gate"] = {
        n: {"expected": made[n][1], "p_expected": gate[n].get(made[n][1]),
            "pass": bool(gate[n].get(made[n][1]) is not None
                         and gate[n][made[n][1]] < 0.05)} for n in made}
    out["metadata"]["gate_verdict"] = ("PASS" if gate_pass else
                                       "FAIL -> implementation broken, Part-2 results VOID")

    # ---- Kamath leg (donor means of cell log1p-CPM per assigned cell type, full background) ----
    split = json.load(open(KAMATH_SPLIT))
    usable = {d for d, v in split["donors"].items()
              if v != "locked_test" and d.startswith("SN-")}
    a = sc.read_h5ad(KAMATH_H5AD, backed="r")
    symbols = np.asarray(a.var["gene_symbol"].astype(str).values)
    sym_col = {}
    for i, g in enumerate(symbols):
        sym_col[g] = i
    sets_k = {k: v for k, v in load_candidate_sets(set(symbols.tolist()), gidx, delta).items()
              if v["kind"] == "per_combo"}
    assign = {sid: crosscohort_celltype(s) for sid, s in sets_k.items()}
    obs = a.obs
    cell_mask = (obs["donor"].isin(usable)
                 & obs["cell_type"].isin(["Microglia", "Astrocyte", "Endothelial"])
                 & (obs["tissue"].astype(str) == "substantia_nigra")).values
    counts = a.layers["counts"][cell_mask, :]
    donor_k = obs["donor"].astype(str).values[cell_mask]
    cond_k = obs["condition"].map({"Disease": "PD", "Ctrl": "Control"}).astype(str).values[cell_mask]
    ctype = obs["cell_type"].astype(str).values[cell_mask]
    libsz = obs["total_counts"].values.astype(float)[cell_mask]
    donor_cond_k = {d: cond_k[donor_k == d][0] for d in set(donor_k.tolist())}
    a.file.close()

    kam_res = {}
    G = counts.shape[1]
    for ct in ("Microglia", "Astrocyte", "Endothelial"):
        sids_ct = [sid for sid, t in assign.items() if t == ct]
        if not sids_ct:
            continue
        rows = np.where(ctype == ct)[0]
        Xc = counts[rows]
        lib = libsz[rows]
        dsub = donor_k[rows]
        keep_d = [d for d in sorted(set(dsub.tolist()))
                  if (dsub == d).sum() >= MIN_CELLS_PER_DONOR]
        acc = {d: np.zeros(G) for d in keep_d}
        for start in range(0, G, 4000):
            B = Xc[:, start:start + 4000]
            Lb = np.log1p(B.multiply(1.0 / np.clip(lib, 1, None)[:, None]).toarray() * 1e6)
            for d in keep_d:
                acc[d][start:start + 4000] = Lb[dsub == d].sum(axis=0)
            del B, Lb
        M = np.stack([acc[d] / (dsub == d).sum() for d in keep_d])
        is_pd = np.array([donor_cond_k[d] == "PD" for d in keep_d])
        for sid in sids_ct:
            s = sets_k[sid]
            r = _competitive_test(M, is_pd,
                                  [sym_col[g] for g in s["genes"] if g in sym_col])
            r["n_present"] = len([g for g in s["genes"] if g in sym_col])
            r["cell_type"] = ct
            r.update(n_donors=len(keep_d), n_PD=int(is_pd.sum()), n_Control=int((~is_pd).sum()))
            kam_res[sid] = r
        del acc, M, Xc
        print(f"Kamath {ct} done ({len(keep_d)} donors)")

    # ---- families ----
    def add_q(res, fam):
        for direction in ("p_up", "p_down"):
            q = bh_qvalues({sid: r.get(direction) for sid, r in res.items()})
            for sid, r in res.items():
                if r.get(direction) is not None:
                    r[f"q_{direction}_{fam}"] = q[sid]

    add_q(ma_res, "ma")
    add_q(kam_res, "kam")
    meta = {}
    for direction in ("p_up", "p_down"):
        fam = {}
        for sid in sets_ma:
            pm, pk = ma_res.get(sid, {}).get(direction), kam_res.get(sid, {}).get(direction)
            zs, ws = [], []
            if pm is not None:
                zs.append(float(stats.norm.isf(min(max(pm, 1e-300), 1 - 1e-16))))
                ws.append(np.sqrt(ma_res[sid]["n_donors"]))
            if pk is not None:
                zs.append(float(stats.norm.isf(min(max(pk, 1e-300), 1 - 1e-16))))
                ws.append(np.sqrt(kam_res[sid]["n_donors"]))
            if not zs:
                continue
            zc2 = float(np.dot(zs, ws) / np.sqrt(np.dot(ws, ws)))
            fam[sid] = {"z_combined": zc2, direction: float(stats.norm.sf(zc2))}
        q = bh_qvalues({sid: r[direction] for sid, r in fam.items()})
        for sid, r in fam.items():
            r[f"q_{direction}_meta"] = q[sid]
        meta[direction] = fam

    def n_sig(res, direction, fam):
        return sum(1 for r in res.values()
                   if r.get(f"q_{direction}_{fam}") is not None and r[f"q_{direction}_{fam}"] <= 0.05)

    out.update({"ma_sets": ma_res, "kamath_sets": kam_res, "gate_sets": gate,
                "meta_stouffer": meta,
                "summary": {
                    "ma_up_sig": n_sig(ma_res, "p_up", "ma"),
                    "ma_down_sig": n_sig(ma_res, "p_down", "ma"),
                    "kam_up_sig": n_sig(kam_res, "p_up", "kam"),
                    "kam_down_sig": n_sig(kam_res, "p_down", "kam"),
                    "meta_up_sig": sum(1 for r in meta["p_up"].values() if r["q_p_up_meta"] <= 0.05),
                    "meta_down_sig": sum(1 for r in meta["p_down"].values() if r["q_p_down_meta"] <= 0.05),
                }})
    out_path = Path(A) / "results" / "p6_spatial" / "testb_competitive_discovery_union_v1.json"
    json.dump(out, open(out_path, "w"), indent=1)
    print(f"summary: {out['summary']}")
    print(f"saved: {out_path}")


# ---------- Amendment 10: EV confirmatory test on GSE184950 ----------
EV_H5AD = "/public/home/mengxl/dzy/pd_product_assets/processed/gse184950/v0.1/GSE184950_annotated.h5ad"
EV_MIN_CELLS = 10
EV_MIN_DONORS_PER_ARM = 3




if __name__ == "__main__":
    main_testb3()