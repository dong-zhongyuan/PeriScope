#!/usr/bin/env python
"""Read-only all-target aggregator with a selection-corrected null next to coloc H4.

Stages (default --stage all runs them in order and skips work whose cache is valid):
  panel      build the decoy z panel from the complete BIG40 GWAS files already on disk
  null       per-region selection null for coloc H4 under every null in NULL_SETS, plus the
             decoy calibration pairs for this simulator (cached per region)
  validate   compare with the earlier 122-region results before anything is scaled up
  aggregate  write the delivery tables from the inputs plus the cached null
--refresh    run only `aggregate` (use after the repro or whole-file audit files change).

Every path this process can create, change or remove is confined to OUT by an audit hook.
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"
sys.dont_write_bytecode = True
from pathlib import Path

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
MR = ROOT / "papers/01_mr_imaging_scan"
OUT = MR / "results/all_target_delivery_20261119"


def _inside_out(p):
    try:
        s = os.path.realpath(os.fsdecode(p))
    except Exception:
        return False
    o = os.path.realpath(str(OUT))
    return s == o or s.startswith(o + os.sep) or s == "/dev/null"


def _guard(event, args):
    if event == "open":
        path, _mode, flags = args
        if isinstance(path, int) or path is None or not isinstance(flags, int):
            return
        if flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND) and not _inside_out(path):
            raise PermissionError(f"refusing to write outside {OUT}: {path}")
    elif event in ("os.remove", "os.rmdir", "os.mkdir", "os.truncate", "os.chmod", "os.chown"):
        if not isinstance(args[0], int) and not _inside_out(args[0]):
            raise PermissionError(f"refusing {event} outside {OUT}: {args[0]}")
    elif event in ("os.rename", "os.link", "os.symlink"):
        for a in args[:2]:
            if not _inside_out(a):
                raise PermissionError(f"refusing {event} outside {OUT}: {a}")
    elif event == "subprocess.Popen":
        exe = os.path.basename(os.fsdecode(args[0]))
        if exe != "zcat":
            raise PermissionError(f"refusing to start {exe}")


sys.addaudithook(_guard)

import argparse, gzip, hashlib, io, json, subprocess, time
from multiprocessing import Pool
import numpy as np
import pandas as pd

sys.path.insert(0, str(ROOT))
from pdv3.mr import ivw, wald  # the estimators production used for reverse MR

COLOC = MR / "results/all_target_coloc_20261118"
REV = MR / "results/all_target_reverse_20261118"
COVERAGE = MR / "results/all_target_coverage_20261118"
REPRO = MR / "results/big40_repro_20261119/region_repro.csv.gz"
FU = ROOT / "data/all_target_followup_20261118"
REG = FU / "regional_all"
IMG = REG / "imaging"
ACC = FU / "acceptance_20261118"
AUD = FU / "fast_range_20261118"
RES = ROOT / "papers/03_world_model_integration/results/research_20261118"
INTG = RES / "integration"
MRI = RES / "mr_imaging"
PF6 = RES / "verify-partial_fetch-PF6-a"
INT1A = RES / "verify-integration-INT-1-a"
INT1R = RES / "verify-integration-INT-1-r"
CACHE = OUT / "cache"
NULLDIR = CACHE / "null_v2"
PARTS = CACHE / "panel_parts"
PANEL = CACHE / "panel_z.f32"
PANEL_KEYS = CACHE / "panel_keys.npy"
PANEL_MANIFEST = CACHE / "panel_manifest.json"

RECIPE = "selnull_v2"
GROUP_VERSION = "groups_v2"
SCALES = (1.0, 1.025, 1.03, 1.05, 1.10)
VALID_SCALES = (1.0, 1.025, 1.05, 1.10)  # the scales the earlier 122-region tables were built with
PRIMARY_SCALE = 1.025
WINDOWS = (1000000, 500000)
PRIORS = ((1e-5, "1e5"), (1e-6, "1e6"))
P12_LABEL = {"1e5": "p12_1e-5", "1e6": "p12_1e-6"}
B = 2000
CLEAN_Z = 5.45
MIN_DECOYS = 30
W2 = 0.15 ** 2
BIG40_ROWS = 17103079
MAX_WORKERS = 12
THRESHOLDS = (0.5, 0.8, 0.9)
WINDOW_TAGS = ("w1000k", "w500k")
UNITS = ("region", "exposure", "locus")
LOCUS_LINK_BP = 1000000
CIS_BP = 1000000
# decoy calibration of this simulator: real decoy GWAS against the simulator's own prediction
CAL_ZMIN = 2.5
CAL_DRAWS = 200
CAL_BOOT = 200
CAL_EDGES = np.array([2.5, 3, 3.25, 3.5, 3.75, 4, 4.5, 5.45, np.inf])
CAL_BINS = ("2.5-3", "3-3.25", "3.25-3.5", "3.5-3.75", "3.75-4", "4-4.5", "4.5-5.45", "5.45-99")


def stag(s):
    return "s" + f"{s:.3f}".replace(".", "p")


# Reference range = the nulls INT-1 derived from the data, plus the two unit-variance scales that pass its
# max-|z| check. The other three stay in every table as a light bound (1.0) and heavy bounds (1.05, 1.10).
DRAW_FILES = (("lambda_median_family", INT1A / "null_lambda_median_family.npy", None),
              ("lambda_median_each", INT1A / "null_lambda_median_each.npy", None),
              ("lambda_meanchi2_each", INT1A / "null_lambda_meanchi2_each.npy", None),
              ("realvar_signflip", INT1R / "own_maxnull_draws.npz", "V1_realvar_signflip"),
              ("realvar_winsor4p5", INT1R / "own_maxnull_draws.npz", "V2_realvar_winsor4.5"))
NULL_NAMES = tuple(stag(s) for s in SCALES) + tuple(d[0] for d in DRAW_FILES)
NULL_SCALE = {stag(s): s for s in SCALES}
REFERENCE_NULLS = (stag(1.025), stag(1.03), "lambda_median_family", "lambda_median_each", "lambda_meanchi2_each", "realvar_signflip", "realvar_winsor4p5")
NULL_ROLE = {n: "reference_range" for n in REFERENCE_NULLS}
NULL_ROLE.update({stag(1.0): "light_bound", stag(1.05): "heavy_bound", stag(1.10): "heavy_bound"})
NULL_LABEL = {stag(1.0): "单位方差 ×1.000", stag(1.025): "单位方差 ×1.025", stag(1.03): "单位方差 ×1.030", stag(1.05): "单位方差 ×1.050", stag(1.10): "单位方差 ×1.100",
              "lambda_median_family": "λ(中位数, 按家族)", "lambda_median_each": "λ(中位数, 逐性状)", "lambda_meanchi2_each": "λ(均值 χ², 逐性状)",
              "realvar_signflip": "真实方差(符号翻转)", "realvar_winsor4p5": "真实方差(缩尾 4.5)"}


def null_sets():
    mx = np.sort(np.load(INTG / "atlas/maxnull_draws.npy").astype(np.float64))
    out = [(stag(s), mx * s) for s in SCALES]
    for name, path, key in DRAW_FILES:
        a = np.load(path)[key] if key else np.load(path)
        out.append((name, np.sort(np.asarray(a).astype(np.float64))))
    assert tuple(n for n, _ in out) == NULL_NAMES
    return out


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def atomic_bytes(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".pending{os.getpid()}")
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, path)


def atomic_json(path, obj):
    atomic_bytes(path, (json.dumps(obj, indent=1, sort_keys=True, default=_js) + "\n").encode())


def _js(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def write_csv(df, path):
    path = Path(path)
    if path.suffix == ".gz":
        import io
        s = io.StringIO()
        df.to_csv(s, index=False)
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0, filename="") as g:
            g.write(s.getvalue().encode())
        atomic_bytes(path, buf.getvalue())
    else:
        atomic_bytes(path, df.to_csv(index=False).encode())


def plan():
    return json.loads((REG / "PLAN.json").read_text())


def imaging_plan():
    return json.loads((IMG / "PLAN.json").read_text())


def chrom_label(c):
    c = str(c).lstrip("0")
    return "X" if c in ("23", "X") else c


def allele_key(chrom, pos, a, b):
    return f"{chrom_label(chrom)}:{int(pos)}:" + "_".join(sorted([str(a), str(b)]))


def full_big40_ids():
    return sorted(p.parent.name for p in IMG.glob("BIG40/*/*.txt.gz") if (p.parent / "COMPLETE.json").exists())


def region_tasks():
    p = plan()["regions"]
    ins = pd.read_csv(MR / "results/full_atlas_20261116/full/selected_cis_instruments.csv", low_memory=False).set_index("exposure_id")
    Z = np.load(INTG / "atlas/z_infer.npy", mmap_mode="r")
    ex = {e: i for i, e in enumerate(pd.read_csv(INTG / "atlas/exposures.csv").exposure_id)}
    ph = {e: i for i, e in enumerate(pd.read_csv(INTG / "atlas/phenotypes.csv").phenotype)}
    out = []
    for rid, r in sorted(p.items()):
        e = r["exposure_ids"][0]
        assert len(r["exposure_ids"]) == 1
        i = ins.loc[e]
        assert int(i.pos37) == int(r["center_grch37"])
        za = np.nan
        if e in ex and r["phenotype"] in ph:
            za = float(Z[ex[e], ph[r["phenotype"]]])
        idp = r["phenotype"].split(":")[-1] if r["phenotype"].startswith("BIG40:") else None
        out.append(dict(region_id=rid, phenotype=r["phenotype"], idp=idp, chromosome=r["chromosome"],
                        center=int(r["center_grch37"]), exposure=e, instrument=str(i.variant),
                        inst_key=allele_key(r["chromosome"], i.pos37, i.ea_x, i.oa_x),
                        z_atlas=abs(za) if np.isfinite(za) else np.nan, model_linked=bool(r["model_linked"])))
    return out


# ---------------------------------------------------------------- panel stage
def _keys_one(rid):
    v = pd.read_csv(COLOC / rid / "window_1000000_snp_posteriors.csv.gz", usecols=["variant_key"], keep_default_na=False)
    return v.variant_key.values.astype(str)


def _bigread(idp, cols, names, dtype):
    p = subprocess.Popen(["zcat", str(IMG / "BIG40" / idp / f"{idp}.txt.gz")], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    d = pd.read_csv(p.stdout, sep=" ", header=0, usecols=cols, names=names, dtype=dtype, engine="c", na_filter=False)
    rc = p.wait()
    if rc != 0:
        raise RuntimeError(f"zcat failed for {idp}: {rc}")
    return d


ALLC = ["chr", "rsid", "pos", "a1", "a2", "beta", "se", "lp"]


def _panel_part(a):
    idp, rows_path, pos_path = a
    f = PARTS / f"{idp}.f32"
    src = IMG / "BIG40" / idp / f"{idp}.txt.gz"
    want = json.loads((src.parent / "COMPLETE.json").read_text())["source_sha256"]
    got = sha256(src)
    if want != got:
        return idp, "sha_mismatch", got
    if f.exists():
        return idp, "cached", got
    rows = np.load(rows_path)
    pos = np.load(pos_path)
    d = _bigread(idp, [2, 5, 6], ALLC, {"pos": np.int32, "beta": np.float64, "se": np.float64})
    if len(d) != BIG40_ROWS:
        return idp, f"rows_{len(d)}", got
    if not (d.pos.values[rows] == pos).all():
        return idp, "order_mismatch", got
    z = (d.beta.values[rows] / d.se.values[rows]).astype(np.float32)
    atomic_bytes(f, z.tobytes())
    return idp, "ok", got


def stage_panel(workers):
    if PANEL.exists() and PANEL_KEYS.exists() and PANEL_MANIFEST.exists():
        m = json.loads(PANEL_MANIFEST.read_text())
        if m.get("ids") == full_big40_ids():
            print("panel cached", m["n_keys"], "keys", len(m["ids"]), "phenotypes", flush=True)
            return
    t = time.time()
    CACHE.mkdir(parents=True, exist_ok=True)
    PARTS.mkdir(parents=True, exist_ok=True)
    tasks = region_tasks()
    with Pool(workers) as po:
        ks = po.map(_keys_one, [x["region_id"] for x in tasks], chunksize=8)
    K = np.unique(np.concatenate(ks + [np.array([x["inst_key"] for x in tasks])]))
    print("union keys", len(K), round(time.time() - t), "s", flush=True)
    ids = full_big40_ids()
    ref = ids[0]
    d = _bigread(ref, [0, 2, 3, 4], ALLC, {"chr": "category", "pos": np.int32, "a1": "category", "a2": "category"})
    assert len(d) == BIG40_ROWS
    c = d.chr.astype(str).str.lstrip("0").values
    want = pd.Series(K).str.split(":", n=2, expand=True)
    wset = set(zip(want[0], want[1].astype(int)))
    posv = d.pos.values
    mask = np.fromiter(((a, int(b)) in wset for a, b in zip(c, posv)), dtype=bool, count=len(d))
    sub = np.where(mask)[0]
    a1 = d.a1.astype(str).values[sub]
    a2 = d.a2.astype(str).values[sub]
    lo = np.where(a1 < a2, a1, a2)
    hi = np.where(a1 < a2, a2, a1)
    keys = pd.Series(c[sub]) + ":" + pd.Series(posv[sub]).astype(str) + ":" + pd.Series(lo) + "_" + pd.Series(hi)
    s = pd.Series(sub, index=keys.values)
    dup = int(s.index.duplicated().sum())
    s = s[~s.index.duplicated()]
    hit = s.reindex(K)
    keep = hit.notna().values
    Kk = K[keep]
    rows = hit.values[keep].astype(np.int64)
    POS = posv[rows]
    print("keys found", int(keep.sum()), "missing", int((~keep).sum()), "dup BIG40 keys", dup, round(time.time() - t), "s", flush=True)
    del d, c, posv, mask, a1, a2, lo, hi, keys, s, hit
    rp, pp = CACHE / "panel_rows.npy", CACHE / "panel_pos.npy"
    with open(rp, "wb") as f:
        np.save(f, rows)
    with open(pp, "wb") as f:
        np.save(f, POS)
    with Pool(workers) as po:
        res = []
        for i, r in enumerate(po.imap_unordered(_panel_part, [(x, str(rp), str(pp)) for x in ids])):
            res.append(r)
            if (i + 1) % 25 == 0:
                print("panel files", i + 1, len(ids), round(time.time() - t), "s", flush=True)
    bad = [r for r in res if r[1] not in ("ok", "cached")]
    if bad:
        raise RuntimeError(f"panel source problems: {bad[:5]}")
    M = np.empty((len(Kk), len(ids)), np.float32)
    for j, i in enumerate(ids):
        M[:, j] = np.fromfile(PARTS / f"{i}.f32", dtype=np.float32)
    tmp = CACHE / "panel_z.f32.pending"
    with open(tmp, "wb") as f:
        f.write(M.data)
    os.replace(tmp, PANEL)
    with open(CACHE / "panel_keys.pending.npy", "wb") as f:
        np.save(f, Kk)
    os.replace(CACHE / "panel_keys.pending.npy", PANEL_KEYS)
    atomic_json(PANEL_MANIFEST, dict(ids=ids, n_keys=int(len(Kk)), n_union_keys=int(len(K)), n_missing_in_big40=int((~keep).sum()),
                                     duplicate_big40_keys_dropped=dup, reference_order_file=ref,
                                     source_sha256={r[0]: r[2] for r in res}, shape=[int(len(Kk)), len(ids)], dtype="float32"))
    for i in ids:
        os.remove(PARTS / f"{i}.f32")
    os.rmdir(PARTS)
    print("panel done", M.shape, round(time.time() - t), "s", flush=True)


# ---------------------------------------------------------------- null stage
def lse(a, axis=-1):
    m = a.max(axis=axis, keepdims=True)
    return np.squeeze(m, axis=axis) + np.log(np.exp(a - m).sum(axis=axis))


def h4_matrix(l1, l2, lp1, lp12):
    s1 = lse(l1)
    s2 = lse(l2, axis=1)
    s12 = lse(l1 + l2, axis=1)
    ab = s1 + s2
    d = ab + np.log(-np.expm1(np.minimum(s12 - ab, -1e-12)))
    L = np.stack([np.zeros_like(s2), lp1 + s1 + 0 * s2, lp1 + s2, 2 * lp1 + d, lp12 + s12], 1)
    return np.exp(L[:, 4] - lse(L, axis=1))


def l2_of(z, rV):
    return 0.5 * (np.log1p(-rV) + rV * z * z)


def priors_for(ns):
    p1 = 1 / (ns + 1) if ns * 1e-4 >= 1 else 1e-4
    out = {}
    for p, tag in PRIORS:
        e = 1 / (ns + 1) if ns * p >= 1 else p
        out[tag] = np.log(e)
    return np.log(p1), out


def stats(h, obs):
    return dict(mean=float(h.mean()), q05=float(np.quantile(h, .05)), q50=float(np.quantile(h, .5)), q95=float(np.quantile(h, .95)),
                f50=float((h >= .5).mean()), f80=float((h >= .8).mean()), f90=float((h >= .9).mean()),
                p=float(((h >= obs).sum() + 1) / (len(h) + 1)), n=int(len(h)))


def seed_of(text):
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:8], "little")


def _init_null():
    global KEYS, ZP, IDS, MX, NS
    KEYS = np.load(PANEL_KEYS, mmap_mode="r")
    IDS = np.array(json.loads(PANEL_MANIFEST.read_text())["ids"])
    ZP = np.memmap(PANEL, dtype=np.float32, mode="r", shape=(len(KEYS), len(IDS)))
    MX = np.sort(np.load(INTG / "atlas/maxnull_draws.npy").astype(np.float64))
    NS = null_sets()


def panel_lookup(keys):
    keys = np.asarray(keys).astype(KEYS.dtype)
    i = np.searchsorted(KEYS, keys)
    i[i >= len(KEYS)] = len(KEYS) - 1
    return i, np.asarray(KEYS[i]) == keys


def null_digest():
    m = sha256(PANEL_MANIFEST)
    x = sha256(INTG / "atlas/maxnull_draws.npy")
    d = "|".join(sha256(p) for p in sorted({f[1] for f in DRAW_FILES}))
    return hashlib.sha256(f"{RECIPE}|{B}|{NULL_NAMES}|{CAL_ZMIN}|{CAL_DRAWS}|{m}|{x}|{d}".encode()).hexdigest()


def region_input_digest(rid):
    return hashlib.sha256("|".join(sha256(COLOC / rid / f"window_{w}_snp_posteriors.csv.gz") for w in WINDOWS).encode()).hexdigest()


USE = ["variant_key", "pos37", "beta_x", "standard_error_x", "beta_y", "standard_error_y", "log_abf_protein", "log_abf_imaging"]
CAL_FIELDS = ["decoy", "z_inst"] + [f"H4_{t}" for _, t in PRIORS] + [f"P{int(th * 100)}_{t}" for _, t in PRIORS for th in THRESHOLDS]


def decoy_pairs(rid, names, Cc, Csc, b, l1, rV, lp1, lpp):
    """Every clean decoy with |z| >= CAL_ZMIN at the instrument: its real H4, and what this simulator
    predicts for a null field forced to the same |z| (background drawn from the other clean decoys)."""
    zi = np.abs(Csc)
    hit = np.where(zi >= CAL_ZMIN)[0]
    if not len(hit):
        return []
    rng = np.random.default_rng(seed_of("decoycal|" + rid))
    L = l2_of(Cc[hit], rV)
    real = {tag: h4_matrix(l1, L, lp1, lp) for tag, lp in lpp.items()}
    rows = []
    for q, d in enumerate(hit):
        j = rng.integers(0, len(Cc) - 1, CAL_DRAWS)
        j[j >= d] += 1
        L = l2_of(Cc[j] + np.outer(zi[d] - Csc[j], b), rV)
        r = [str(names[d]), float(zi[d])] + [float(real[t][q]) for _, t in PRIORS]
        for _, t in PRIORS:
            h = h4_matrix(l1, L, lp1, lpp[t])
            r += [float((h >= th).mean()) for th in THRESHOLDS]
        rows.append(r)
    return rows


def null_one(a):
    t, nd = a
    rid = t["region_id"]
    fj, fn = NULLDIR / f"{rid}.json", NULLDIR / f"{rid}.npy"
    rd = region_input_digest(rid)
    if fj.exists() and fn.exists():
        try:
            o = json.loads(fj.read_text())
            if o.get("null_digest") == nd and o.get("region_digest") == rd:
                return rid, "cached"
        except Exception:
            pass
    keep_ids = IDS != t["idp"] if t["idp"] else np.ones(len(IDS), bool)
    out = dict(region_id=rid, null_digest=nd, region_digest=rd, exposure=t["exposure"], instrument=t["instrument"], inst_key=t["inst_key"],
               n_panel_phenotypes=int(keep_ids.sum()), own_phenotype_in_panel=bool(t["idp"] and (IDS == t["idp"]).any()), windows={}, results={})
    sims = np.full((len(WINDOWS), len(PRIORS), len(NS), B), np.nan, dtype=np.float32)
    rng = np.random.default_rng(seed_of("selnull|" + t["exposure"]))
    u = rng.random(B)
    uj = rng.random(B)
    zdraw = lambda S: S[(u * len(S)).astype(int)]  # same u for every null, so nulls differ only through their z* law
    i_inst, ok_inst = panel_lookup([t["inst_key"]])
    status, reason = "computed", ""
    pr = pd.read_csv(COLOC / rid / "coloc_prior_sensitivity.csv")
    for wi, w in enumerate(WINDOWS):
        wtag = f"w{w // 1000}k"
        q = pr[pr.window_bp == w]
        if q.status.iloc[0] != "defined":
            out["windows"][wtag] = dict(status="window_undefined", reason=str(q.undefined_reason.iloc[0]))
            continue
        v = pd.read_csv(COLOC / rid / f"window_{w}_snp_posteriors.csv.gz", usecols=USE, keep_default_na=False,
                        dtype={"variant_key": str})
        ns = len(v)
        idx, ok = panel_lookup(v.variant_key.values)
        l1all = v.log_abf_protein.values.astype(float)
        zy = (v.beta_y / v.standard_error_y).values
        Vy = v.standard_error_y.values ** 2
        lp1, lpp = priors_for(ns)
        info = dict(status="computed", n_snps=ns, n_snps_without_panel_z=int((~ok).sum()))
        rVa = W2 / (W2 + Vy)
        info["max_abs_diff_l2_vs_production"] = float(np.abs(l2_of(zy, rVa) - v.log_abf_imaging.values.astype(float)).max())
        obs_prod = {}
        for p, tag in PRIORS:
            obs_prod[tag] = float(q.iloc[(np.log(q.p12_requested) - np.log(p)).abs().values.argmin()]["PP.H4.abf"])
        for tag, lp in lpp.items():
            h = float(h4_matrix(l1all, l2_of(zy, rVa)[None, :], lp1, lp)[0])
            info[f"recomputed_minus_production_H4_{tag}"] = h - obs_prod[tag]
        if not ok_inst[0]:
            info.update(status="not_computable", reason="instrument_not_in_BIG40_panel")
            out["windows"][wtag] = info
            continue
        Cf = np.nan_to_num(np.asarray(ZP[idx[ok]])[:, keep_ids].T.astype(np.float64))
        nanfrac = float(np.isnan(np.asarray(ZP[idx[ok][:50]])[:, keep_ids]).mean())
        l1 = l1all[ok]
        rV = rVa[ok]
        keys_k = v.variant_key.values[ok]
        hit = np.where(keys_k == t["inst_key"])[0]
        s = int(hit[0]) if len(hit) else -1
        if s >= 0:
            Cs = Cf[:, s]
        else:
            Cs = np.nan_to_num(np.asarray(ZP[i_inst[0]])[keep_ids].astype(np.float64))
        info["instrument_in_coloc_snps"] = bool(np.any(v.variant_key.values == t["inst_key"]))
        info["instrument_in_simulated_snps"] = bool(s >= 0)
        clean = (np.abs(Cf).max(1) < CLEAN_Z) & (np.abs(Cs) < CLEAN_Z)
        info["n_decoys"] = int(len(Cf))
        info["n_clean_decoys"] = int(clean.sum())
        info["panel_nan_frac_sample"] = nanfrac
        if t["idp"] and (IDS == t["idp"]).any():
            own = np.asarray(ZP[idx[ok]])[:, IDS == t["idp"]][:, 0].astype(float)
            info["own_phenotype_max_abs_z_diff"] = float(np.abs(np.abs(own) - np.abs(zy[ok])).max())
        if clean.sum() < MIN_DECOYS:
            info.update(status="not_computable", reason=f"fewer_than_{MIN_DECOYS}_clean_decoys")
            out["windows"][wtag] = info
            continue
        Cc, Csc = Cf[clean], Cs[clean]
        den = float(Csc @ Csc)
        b = (Cc.T @ Csc) / den
        if s >= 0:
            b[s] = 1.0
        jj = (uj * len(Cc)).astype(int)
        zobs_idx = np.where(v.variant_key.values == t["inst_key"])[0]
        z_abs = float(abs(zy[zobs_idx[0]])) if len(zobs_idx) else t["z_atlas"]
        info["instrument_abs_z_observed"] = z_abs
        info["instrument_abs_z_source"] = "coloc_shared_snp" if len(zobs_idx) else ("MR_atlas_infer" if np.isfinite(t["z_atlas"]) else "unavailable")
        l2o = l2_of(zy[ok], rV)
        obs = {tag: float(h4_matrix(l1, l2o[None, :], lp1, lp)[0]) for tag, lp in lpp.items()}
        info["H4_obs_on_simulated_snps"] = obs
        info["H4_production"] = obs_prod
        info["protein_abf_mass_on_dropped_snps"] = float(1 - np.exp(lse(l1) - lse(l1all))) if (~ok).any() else 0.0

        def run(D_fn, n):
            acc = {tag: [] for _, tag in PRIORS}
            for k in range(0, n, 250):
                L2 = l2_of(D_fn(np.arange(k, min(k + 250, n))), rV)
                for tag, lp in lpp.items():
                    acc[tag].append(h4_matrix(l1, L2, lp1, lp))
            return {tag: np.concatenate(x) for tag, x in acc.items()}

        res = {}
        for si, (nm, S) in enumerate(NS):
            zb = zdraw(S)
            H = run(lambda r: Cc[jj[r]] + np.outer(zb[r] - Csc[jj[r]], b), B)
            for pi, (_, tag) in enumerate(PRIORS):
                res[f"{tag}|{nm}"] = stats(H[tag], obs[tag])
                sims[wi, pi, si] = H[tag].astype(np.float32)
        if np.isfinite(z_abs):
            H = run(lambda r: Cc[r] + np.outer(z_abs - Csc[r], b), len(Cc))
            for _, tag in PRIORS:
                res[f"{tag}|cond_obs_z"] = stats(H[tag], obs[tag])
        H = run(lambda r: Cc[r], len(Cc))
        for _, tag in PRIORS:
            res[f"{tag}|decoy_no_selection"] = stats(H[tag], obs[tag])
        if w == 1000000:
            bb = (Cf.T @ Cs) / float(Cs @ Cs)
            shr = np.median(bb * bb) / 0.4549
            bb = np.sign(bb) * np.sqrt(np.maximum(bb * bb - shr, 0))
            if s >= 0:
                bb[s] = 1.0
            ja = (uj * len(Cf)).astype(int)
            zb = zdraw(dict(NS)[stag(PRIMARY_SCALE)])
            H = run(lambda r: Cf[ja[r]] + np.outer(zb[r] - Cs[ja[r]], bb), B)
            for _, tag in PRIORS:
                res[f"{tag}|shrunkLD_alldecoys_{stag(PRIMARY_SCALE)}"] = stats(H[tag], obs[tag])
            out["decoy_pairs"] = decoy_pairs(rid, IDS[keep_ids][clean], Cc, Csc, b, l1, rV, lp1, lpp)
        out["windows"][wtag] = info
        out["results"][wtag] = res
    w1 = out["windows"].get("w1000k", {})
    if w1.get("status") != "computed":
        status, reason = "not_computable", w1.get("reason", "window_1000k_undefined")
    out["status"] = status
    out["reason"] = reason
    NULLDIR.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO()
    np.save(buf, sims)
    atomic_bytes(fn, buf.getvalue())
    atomic_json(fj, out)
    return rid, status


def stage_null(workers, only=None):
    t0 = time.time()
    tasks = region_tasks()
    if only:
        tasks = [x for x in tasks if x["region_id"] in only]
    _init_null()
    nd = null_digest()
    args = [(x, nd) for x in tasks]
    n = {}
    with Pool(workers) as po:
        for i, (rid, st) in enumerate(po.imap_unordered(null_one, args, chunksize=2)):
            n[st] = n.get(st, 0) + 1
            if (i + 1) % 100 == 0 or i + 1 == len(args):
                print("null", i + 1, len(args), n, round(time.time() - t0), "s", flush=True)
    return n


# ---------------------------------------------------------------- validation against the earlier 122-region tables
def _replay_one(a):
    rec, ids = a
    rid = rec["region_id"]
    panel_ix = np.array([int(np.where(IDS == i)[0][0]) for i in ids])
    pid = rec["phenotype"].split(":")[-1]
    keep = np.array([i != pid for i in ids])
    v = pd.read_csv(COLOC / rid / "window_1000000_snp_posteriors.csv.gz", usecols=USE + ["variant_id_x", "variant_id_y"], keep_default_na=False)
    idx, ok = panel_lookup(v.variant_key.values)
    if not ok.all():
        return dict(region_id=rid, replay_status="snps_missing_from_panel", n_missing=int((~ok).sum()))
    l1 = v.log_abf_protein.values.astype(float)
    zx = (v.beta_x / v.standard_error_x).values
    zy = (v.beta_y / v.standard_error_y).values
    Vy = v.standard_error_y.values ** 2
    ns = len(v)
    lp1, lpp = priors_for(ns)
    m = np.where((v.variant_id_y.values == rec["mr_variant_infer"]) | (v.variant_id_x.values == rec["mr_variant_infer"]))[0]
    src = "MR_instrument" if len(m) else "protein_lead_proxy"
    s = int(m[0]) if len(m) else int(np.argmax(np.abs(zx)))
    zs = abs(zy[s])
    C = np.asarray(ZP[idx])[:, panel_ix[keep]].T.astype(np.float64)
    C = C[np.abs(C).max(1) < CLEAN_Z]
    b = (C.T @ C[:, s]) / (C[:, s] @ C[:, s])
    rV = W2 / (W2 + Vy)
    ZC = C + np.outer(zs - C[:, s], b)
    o = dict(region_id=rid, replay_status="ok", replay_instrument_source=src, replay_panel_draws=len(C))
    L2 = l2_of(ZC, rV)
    for tag, lp in lpp.items():
        H = h4_matrix(l1, L2, lp1, lp)
        ob = float(h4_matrix(l1, l2_of(zy, rV)[None, :], lp1, lp)[0])
        o[f"replay_H4_{tag}"] = ob
        o[f"replay_H4null_mean_{tag}"] = float(H.mean())
        o[f"replay_H4null_q95_{tag}"] = float(np.quantile(H, .95))
        o[f"replay_H4null_fracGE08_{tag}"] = float((H >= .8).mean())
        o[f"replay_H4null_P_ge_obs_{tag}"] = float(((H >= ob).sum() + 1) / (len(H) + 1))
    rng = np.random.default_rng(seed_of("replay|" + rid))
    kk = rng.integers(0, len(C), B)
    zb = MX[rng.integers(0, len(MX), B)]
    sims = {}
    for sc in VALID_SCALES:
        acc = []
        for k in range(0, B, 250):
            j = kk[k:k + 250]
            D = C[j] + np.outer(zb[k:k + 250] * sc - C[j, s], b)
            acc.append(h4_matrix(l1, l2_of(D, rV), lp1, lpp["1e5"]))
        sims[stag(sc)] = np.concatenate(acc).astype(np.float32)
    o["_sims"] = sims
    return o


def stage_validate(workers):
    t0 = time.time()
    _init_null()
    Q = pd.read_csv(INTG / "coloc_selection_null_122.csv")
    T = pd.read_csv(INTG / "table_122_regions_raw.csv")[["region_id", "phenotype", "mr_variant_infer", "mr_exposure_id_infer"]]
    ids = [i for i in json.loads((INTG / "stream/ids.json").read_text()) if (INTG / f"stream/z/{i}.npy").exists()]
    with Pool(workers) as po:
        R = po.map(_replay_one, [(r, ids) for r in T.to_dict("records")], chunksize=2)
    sims = {r["region_id"]: r.pop("_sims") for r in R if "_sims" in r}
    D = pd.DataFrame(R).merge(Q, on="region_id", how="left")
    rows = []

    def add(item, earlier, ours, note=""):
        rows.append(dict(item=item, earlier=earlier, this_run=ours, difference=(ours - earlier) if earlier == earlier and ours == ours else np.nan, note=note))

    ok = D[D.replay_status == "ok"]
    add("regions replayed with the earlier recipe (same 254 panel phenotypes, same filters)", 122, len(ok))
    add("regions with the MR instrument among coloc SNPs", int((Q.instrument_source == "MR_instrument").sum()), int((ok.replay_instrument_source == "MR_instrument").sum()))
    for tag in ("1e5", "1e6"):
        for mine, theirs in [("replay_H4", f"H4_{tag}_1000k"), ("replay_H4null_mean", f"H4null_mean_{tag}"), ("replay_H4null_q95", f"H4null_q95_{tag}"),
                             ("replay_H4null_fracGE08", f"H4null_fracGE08_{tag}"), ("replay_H4null_P_ge_obs", f"H4null_P_ge_obs_{tag}")]:
            d = (ok[f"{mine}_{tag}"] - ok[theirs]).abs()
            add(f"max abs difference, {theirs} (deterministic part of the earlier recipe)", 0.0, float(d.max()), f"n={int(d.notna().sum())}")
    S0 = pd.read_csv(INTG / "coloc_null_scale_sensitivity.csv")
    inst = ok[ok.replay_instrument_source == "MR_instrument"].region_id.tolist()
    obs = ok.set_index("region_id")["replay_H4_1e5"].reindex(inst).values
    for sc in VALID_SCALES:
        X = np.stack([sims[r][stag(sc)] for r in inst])
        for th in (.5, .8, .9):
            c = (X >= th).sum(0)
            e = S0[(np.isclose(S0.max_z_scale, sc)) & (np.isclose(S0.threshold, th))]
            add(f"expected regions with H4>={th} among {len(inst)} under selection null, scale {sc}, earlier recipe, fresh Monte Carlo",
                float(e.null_mean.iloc[0]) if len(e) else np.nan, float(c.mean()), f"observed {int((obs >= th).sum())}; Monte Carlo sd of this mean about {c.std() / np.sqrt(B):.2f}")
    per = D[["region_id", "gene", "phenotype", "instrument_source", "replay_status", "replay_instrument_source", "replay_panel_draws", "panel_null_draws",
             "H4_1e5_1000k", "replay_H4_1e5", "H4null_mean_1e5", "replay_H4null_mean_1e5", "H4null_P_ge_obs_1e5", "replay_H4null_P_ge_obs_1e5",
             "H4_1e6_1000k", "replay_H4_1e6", "H4null_mean_1e6", "replay_H4null_mean_1e6"]].copy()
    # same 122 regions under the recipe used for all regions
    prod = []
    for rid in T.region_id:
        f = NULLDIR / f"{rid}.json"
        if not f.exists():
            continue
        o = json.loads(f.read_text())
        wi = o["windows"].get("w1000k", {})
        rs = o["results"].get("w1000k", {})
        rec = dict(region_id=rid, all_region_status=o["status"], all_region_instrument_in_coloc_snps=wi.get("instrument_in_coloc_snps"), all_region_clean_decoys=wi.get("n_clean_decoys"))
        for k in ("cond_obs_z", stag(1.0), stag(PRIMARY_SCALE)):
            x = rs.get(f"1e5|{k}")
            if x:
                rec[f"all_region_mean_{k}"] = x["mean"]
                rec[f"all_region_p_{k}"] = x["p"]
                rec[f"all_region_f80_{k}"] = x["f80"]
        prod.append(rec)
    if prod:
        P = pd.DataFrame(prod)
        per = per.merge(P, on="region_id", how="left")
        J = per[per.instrument_source == "MR_instrument"]
        if "all_region_mean_cond_obs_z" in J:
            add("Pearson r, conditional-null mean H4, earlier 254-panel vs all-region recipe (301 panel), 101 regions", 1.0, float(J.H4null_mean_1e5.corr(J.all_region_mean_cond_obs_z)))
            add("mean of conditional-null mean H4, 101 regions", float(J.H4null_mean_1e5.mean()), float(J.all_region_mean_cond_obs_z.mean()), "earlier 254-panel vs all-region recipe")
            add("max abs difference of conditional-null mean H4, 101 regions", 0.0, float((J.H4null_mean_1e5 - J.all_region_mean_cond_obs_z).abs().max()), "panel grows from 254 to 301 phenotypes, so this is not expected to be zero")
        for sc in (1.0, PRIMARY_SCALE):
            c = f"all_region_f80_{stag(sc)}"
            if c in J:
                e = S0[(np.isclose(S0.max_z_scale, sc)) & (np.isclose(S0.threshold, .8))]
                add(f"expected H4>=0.8 count among 101 under selection null, scale {sc}, all-region recipe", float(e.null_mean.iloc[0]), float(J[c].sum()), "sum over regions of P(H4>=0.8 | selection null)")
        Pr = per[per.instrument_source != "MR_instrument"]
        add("regions the earlier run handled with a protein-lead proxy, now conditioned on the real instrument through the panel or marked not computable", len(Pr), int(Pr.all_region_status.notna().sum()), "see per-region file")
    L = pd.read_csv(MRI / "null_with_ld_by_region.csv")
    rr = []
    for rid in L.region_id:
        f = NULLDIR / f"{rid}.json"
        if f.exists():
            o = json.loads(f.read_text())
            x = o["results"].get("w500k", {})
            if f"1e5|{stag(1.0)}" in x:
                rr.append(dict(region_id=rid, p10=x[f"1e5|{stag(1.0)}"]["p"], p1025=x[f"1e5|{stag(PRIMARY_SCALE)}"]["p"], m10=x[f"1e5|{stag(1.0)}"]["mean"], m1025=x[f"1e5|{stag(PRIMARY_SCALE)}"]["mean"], f1025=x[f"1e5|{stag(PRIMARY_SCALE)}"]["f80"]))
    if rr:
        J = L.merge(pd.DataFrame(rr), on="region_id")
        add("regions shared with the earlier genotype-LD null (mr_imaging s16, 500 kb, p12=1e-5)", len(L), len(J))
        add("Spearman, selection p: s16 vs this run (scale 1.025)", 1.0, float(J.emp_p.corr(J.p1025, method="spearman")), "s16 uses EUR503 genotype LD and its own max-|z| null")
        add("mean expected H4 under selection", float(J.null_mean_H4.mean()), float(J.m1025.mean()), "s16 vs this run at scale 1.025")
        add("expected count H4>=0.8", float(J.null_p80.sum()), float(J.f1025.sum()), "s16 vs this run at scale 1.025")
    write_csv(per, OUT / "validation_122_per_region.csv")
    write_csv(pd.DataFrame(rows), OUT / "validation_122_summary.csv")
    print(pd.DataFrame(rows).round(5).to_string(), flush=True)
    print("validate done", round(time.time() - t0), "s", flush=True)


# ---------------------------------------------------------------- aggregate stage
AUDIT_LINES = {}


def _jsonl(p):
    """Read complete lines only; the audit is still appending while this runs."""
    out = {}
    n = 0
    if Path(p).exists():
        raw = Path(p).read_bytes()
        end = raw.rfind(b"\n")
        for line in raw[:end + 1].decode().splitlines() if end >= 0 else []:
            if line.strip():
                r = json.loads(line)
                out[str(r.get("unit"))] = r
                n += 1
    AUDIT_LINES[str(p)] = n
    return out


def load_audits():
    A = {"protein": {}, "imaging": {}}
    for kind in A:
        for src, p in (("earlier_stratified_sample", AUD / f"audit/results_{kind}.jsonl"), ("targeted_H4_20261119", AUD / f"audit_targeted_H4_20261119/results_{kind}.jsonl")):
            for u, r in _jsonl(p).items():
                A[kind].setdefault(u, []).append((src, r))
    m = json.loads((AUD / "audit_targeted_H4_20261119/MANIFEST.json").read_text())
    return A, {k: set(v) for k, v in m["units"].items()}


def audit_status(kind, unit, route, audits, targeted):
    if route == "full_file":
        return "whole_file_read_at_acquisition", ""
    if route != "fast_range_v1_20261118":
        return "route_" + str(route), ""
    recs = audits[kind].get(str(unit), [])
    if not recs:
        return ("range_not_audited_in_targeted_queue" if str(unit) in targeted[kind] else "range_not_audited"), ""
    states = []
    for src, r in recs:
        if "outputs_match" not in r or "rows_match" not in r:
            states.append("incomplete")
            continue
        if kind == "protein":
            good = r.get("outputs_match") is True and r.get("rows_match") is True and r.get("published_full_md5_verified") is True
        else:
            good = r.get("outputs_match") is True and r.get("rows_match") is True and r.get("gzip_crc_eof_verified") is True and r.get("full_numeric_validation") is True
        states.append("match" if good else "mismatch")
    src = "+".join(s for s, _ in recs)
    if "mismatch" in states:
        return "range_audited_mismatch", src
    if "incomplete" in states:
        return "range_audit_record_incomplete", src
    return "range_audited_match", src


REV_NUM = ["source_peaks", "LD_selected_before_outcome_join", "usable_instruments", "n_instruments", "beta", "se", "p", "se_random", "p_random", "ci_lower", "ci_upper", "F", "F_min", "Q", "Q_p"]
NOCIS_NUM = ["beta", "se", "p", "se_random", "p_random", "ci_lower", "ci_upper", "F", "F_min", "Q", "Q_p"]
CLUMP_ELIGIBLE = ("retained_by_exposure_only_LD_pruning", "LD_pruned", "LD_not_estimable")
REGION_CENTER = {}


def fmt_p12(x):
    return f"{x:.3e}"


PENDING_X = set()


def _chr(s):
    return s.astype(str).str.lstrip("0").replace({"23": "X"})


def reclump(au, pairs, exclude):
    """Production's exposure-only greedy clump (r2 <= 0.001 within 10 Mb), optionally without some peaks."""
    r = {frozenset((a, b)): c for a, b, c in zip(pairs.variant_a, pairs.variant_b, pairs.r)}
    kept = []
    x = au[au.ld_selection_status.isin(CLUMP_ELIGIBLE) & ~au.variant_key.isin(exclude)].sort_values(["discovery_neglog10_p", "variant_key"], ascending=[False, True])
    for k, c, p in zip(x.variant_key, x.chromosome, x.position_grch37):
        good = True
        for k2, c2, p2 in kept:
            if c == c2 and abs(p - p2) <= 10000000:
                v = r.get(frozenset((k, k2)), np.nan)
                if not np.isfinite(v) or v * v > 0.001:
                    good = False
                    break
        if good:
            kept.append((k, c, p))
    return {k for k, _, _ in kept}


def mr_estimate(v):
    if len(v) == 1:
        r = v.iloc[0]
        return wald(r.beta_imaging, r.se_imaging, r.protein_beta_aligned, r.protein_se)
    return ivw(v.beta_imaging.values, v.se_imaging.values, v.protein_beta_aligned.values, v.protein_se.values, independent=True)


def reverse_cis_columns(rid, e):
    """Which instruments of each reverse estimate lie inside the protein's own 1 Mb window, and the estimate without them.
    Cis peaks are removed before the exposure-only clump, so a peak they had pruned can come back."""
    chrom, center = REGION_CENTER[rid]
    rd = REV / rid
    ra = pd.read_csv(rd / "all_peak_reverse_ratios_and_coverage.csv", low_memory=False)
    au = pd.read_csv(rd / "all_peaks_exposure_only_LD_audit.csv", low_memory=False)
    o = {}
    for x in (ra, au):
        x["dist"] = (x.position_grch37 - center).abs() if len(x) else []
        x["is_cis"] = (_chr(x.chromosome) == chrom) & (x.dist <= CIS_BP) if len(x) else []
    cis_keys = set(au.variant_key[au.is_cis]) if len(au) else set()
    o["rev_n_published_peaks_in_cis_window"] = len(cis_keys)
    eligible_cis = bool(len(au) and au[au.is_cis].ld_selection_status.isin(CLUMP_ELIGIBLE).any())
    prod_kept = set(au.variant_key[au.ld_selection_status.eq("retained_by_exposure_only_LD_pruning")]) if len(au) else set()
    nocis_kept = prod_kept
    o["rev_reclump_reproduces_production"] = True
    if eligible_cis:
        pairs = pd.read_csv(rd / "all_peak_pair_reference_LD.csv", usecols=["variant_a", "variant_b", "r"])
        o["rev_reclump_reproduces_production"] = reclump(au, pairs, set()) == prod_kept
        nocis_kept = reclump(au, pairs, cis_keys)
    dmax_b, dmax_p, n_with = 0.0, 0.0, 0
    for r in e.to_dict("records"):
        t = f"rev_{r['palindrome_policy']}_{r['imaging_effect_set']}"
        g = ra[ra.palindrome_policy.eq(r["palindrome_policy"]) & ra.imaging_effect_set.eq(r["imaging_effect_set"]) & ra.status.eq("estimated")] if len(ra) else ra
        v = g[g.selected_after_LD.astype(str).eq("True")] if len(g) else g
        est = r["status"] == "estimated"
        n_cis = int(v.is_cis.sum()) if len(v) else 0
        o[t + "_n_cis_instruments"] = n_cis if est else 0
        o[t + "_n_cis_instruments_500kb"] = int((v.is_cis & (v.dist <= 500000)).sum()) if est and len(v) else 0
        o[t + "_evidence_tier"] = ("B" if n_cis else "A") if est else ""
        rec = dict(status=r["status"], method="none", n=0, added=0)
        vals = {}
        if est:
            m = mr_estimate(v)
            dmax_b = max(dmax_b, abs(m["beta"] - r["beta"]))
            dmax_p = max(dmax_p, abs(m["p"] - r["p"]))
            n_with += int(n_cis > 0)
            w = g[g.variant_key.isin(nocis_kept)]
            rec.update(n=len(w), added=int((~w.variant_key.isin(set(v.variant_key))).sum()))
            if len(w):
                vals = mr_estimate(w)
                rec.update(status="estimated", method=vals["method"])
            else:
                rec.update(status="not_estimated_no_instrument_outside_cis_window")
        o[t + "_nocis_status"] = rec["status"]
        o[t + "_nocis_method"] = rec["method"]
        o[t + "_nocis_n_instruments"] = rec["n"]
        o[t + "_nocis_n_added_by_reclump"] = rec["added"]
        for c in NOCIS_NUM:
            o[f"{t}_nocis_{c}"] = vals.get(c, np.nan)
        o[t + "_nocis_evidence_tier"] = ("A" if (n_cis == 0 and rec["added"] == 0) else "A*") if rec["status"] == "estimated" and est else ""
    o["rev_n_estimates_with_cis_instrument"] = n_with
    o["rev_max_abs_recomputed_minus_production_beta"] = dmax_b
    o["rev_max_abs_recomputed_minus_production_p"] = dmax_p
    return o


def read_region(rid):
    o = dict(region_id=rid)
    d = COLOC / rid
    q = pd.read_csv(d / "coloc_prior_sensitivity.csv", keep_default_na=False, na_values=["NA"])
    grid = sorted(q.p12_requested.unique())
    for w, g in q.groupby("window_bp"):
        wt = f"w{int(w) // 1000}k"
        g = g.sort_values("p12_requested")
        o[f"coloc_{wt}_status"] = g.status.iloc[0]
        o[f"coloc_{wt}_undefined_reason"] = g.undefined_reason.iloc[0] if str(g.undefined_reason.iloc[0]) != "nan" else ""
        o[f"coloc_{wt}_nsnps"] = int(g.nsnps.iloc[0])
        o[f"coloc_{wt}_p1"] = float(g.p1.iloc[0])
        o[f"coloc_{wt}_sdY_protein"] = float(g.sdY_protein_estimated.iloc[0]) if pd.notna(g.sdY_protein_estimated.iloc[0]) else np.nan
        o[f"coloc_{wt}_n_priors_with_changed_p12"] = int((~np.isclose(g.p12, g.p12_requested)).sum())
        for r in g.to_dict("records"):
            p = fmt_p12(r["p12_requested"])
            for h in range(5):
                o[f"coloc_{wt}_p12_{p}_PPH{h}"] = float(r[f"PP.H{h}.abf"])
    v = pd.read_csv(d / "coloc_scale_diagnostics.csv")
    for r in v.itertuples(index=False):
        wt = f"w{int(r.window_bp) // 1000}k"
        o[f"coloc_{wt}_min_p_protein"] = float(r.min_p_protein)
        o[f"coloc_{wt}_min_p_imaging"] = float(r.min_p_imaging)
    rd = REV / rid
    e = pd.read_csv(rd / "reverse_direction_estimates.csv")
    assert len(e) == 4
    vc = e.status.value_counts().to_dict()
    o["rev_status_counts"] = json.dumps({k: int(vc[k]) for k in sorted(vc)}, sort_keys=True)
    o["rev_n_estimated"] = int(vc.get("estimated", 0))
    o["rev_n_no_usable_instrument"] = int(vc.get("not_estimated_no_usable_instrument", 0))
    o["rev_n_no_published_peak"] = int(vc.get("not_estimated_no_published_peak", 0))
    o["rev_n_other_status"] = 4 - o["rev_n_estimated"] - o["rev_n_no_usable_instrument"] - o["rev_n_no_published_peak"]
    for r in e.to_dict("records"):
        t = f"rev_{r['palindrome_policy']}_{r['imaging_effect_set']}"
        o[t + "_status"] = r["status"]
        o[t + "_method"] = r["method"]
        for c in REV_NUM:
            o[f"{t}_{c}"] = r.get(c, np.nan)
    rc = pd.read_csv(rd / "reference_coverage.csv")
    x = rc.variant_key.astype(str).str.startswith("23:")
    o["rev_n_published_peaks"] = int(len(rc))
    o["rev_n_chrX_peaks"] = int(x.sum())
    o["rev_n_chrX_peaks_reference_matched"] = int((x & rc.reference_status.astype(str).eq("matched")).sum())
    o["rev_n_chrX_peaks_listed_pending_in_request_table"] = int(rc.variant_key[x].astype(str).isin(PENDING_X).sum())
    ra = pd.read_csv(rd / "all_peak_reverse_ratios_and_coverage.csv", usecols=["variant_key", "selected_after_LD", "status"])
    u = ra[ra.variant_key.astype(str).str.startswith("23:") & ra.selected_after_LD.astype(str).eq("True") & ra.status.eq("estimated")]
    o["rev_n_chrX_instruments_in_estimates"] = int(u.variant_key.nunique())
    o["rev_chrX_status"] = "no_chrX_peaks" if not x.any() else ("chrX_instruments_in_estimates_X_dosage_convention_not_author_reviewed" if len(u) else "chrX_peaks_present_none_in_estimates")
    o["rev_reference_status_counts"] = json.dumps({k: int(n) for k, n in sorted(rc.reference_status.astype(str).value_counts().items())}, sort_keys=True)
    o.update(reverse_cis_columns(rid, e))
    return o, grid


def _hash_region(rid):
    d, r = COLOC / rid, REV / rid
    return dict(region_id=rid, coloc_prior_sha256=sha256(d / "coloc_prior_sensitivity.csv"), coloc_complete_sha256=sha256(d / "COMPLETE.json"),
                window_1000k_sha256=sha256(d / "window_1000000_snp_posteriors.csv.gz"), window_500k_sha256=sha256(d / "window_500000_snp_posteriors.csv.gz"),
                reverse_estimates_sha256=sha256(r / "reverse_direction_estimates.csv"), reverse_reference_sha256=sha256(r / "reference_coverage.csv"), reverse_complete_sha256=sha256(r / "COMPLETE.json"),
                reverse_ratios_sha256=sha256(r / "all_peak_reverse_ratios_and_coverage.csv"), reverse_ld_audit_sha256=sha256(r / "all_peaks_exposure_only_LD_audit.csv"), reverse_ld_pairs_sha256=sha256(r / "all_peak_pair_reference_LD.csv"))


def _read_null(rid):
    f = NULLDIR / f"{rid}.json"
    if not f.exists():
        return None
    return json.loads(f.read_text())


def null_columns(o):
    r = dict(region_id=o["region_id"], selnull_status=o["status"], selnull_reason=o["reason"], selnull_n_panel_phenotypes=o["n_panel_phenotypes"], selnull_own_phenotype_in_panel=o["own_phenotype_in_panel"])
    for wt in ("w1000k", "w500k"):
        i = o["windows"].get(wt, {})
        for k in ("n_snps", "n_snps_without_panel_z", "n_decoys", "n_clean_decoys", "instrument_in_coloc_snps", "instrument_in_simulated_snps",
                  "instrument_abs_z_observed", "instrument_abs_z_source", "protein_abf_mass_on_dropped_snps", "own_phenotype_max_abs_z_diff"):
            r[f"selnull_{wt}_{k}"] = i.get(k, np.nan)
        for tag in ("1e5", "1e6"):
            r[f"selnull_{wt}_{tag}_H4_obs_used"] = i.get("H4_obs_on_simulated_snps", {}).get(tag, np.nan)
            r[f"selnull_{wt}_{tag}_recomputed_minus_production_H4"] = i.get(f"recomputed_minus_production_H4_{tag}", np.nan)
        for key, x in o["results"].get(wt, {}).items():
            tag, name = key.split("|")
            for f, lab in (("p", "p"), ("mean", "expected_H4"), ("q95", "H4_q95"), ("f50", "P_H4ge050"), ("f80", "P_H4ge080"), ("f90", "P_H4ge090")):
                r[f"selnull_{wt}_{tag}_{name}_{lab}"] = x[f]
    return r


def load_sims(rids):
    return np.stack([np.load(NULLDIR / f"{r}.npy") for r in rids])


def assign_loci(T):
    """Single-linkage loci: protein instruments on one chromosome whose positions chain within LOCUS_LINK_BP."""
    u = T[["chromosome_label", "center_grch37"]].drop_duplicates().copy()
    u["o"] = [23 if c == "X" else int(c) for c in u.chromosome_label]
    u = u.sort_values(["o", "center_grch37"]).reset_index(drop=True)
    new = (u.o != u.o.shift()) | (u.center_grch37 - u.center_grch37.shift() > LOCUS_LINK_BP)
    u["k"] = new.cumsum()
    g = u.groupby("k").center_grch37.agg(["min", "max"])
    u["locus_id"] = [f"chr{c}:{g.loc[k, 'min']}-{g.loc[k, 'max']}" for c, k in zip(u.chromosome_label, u.k)]
    m = {(c, p): l for c, p, l in zip(u.chromosome_label, u.center_grch37, u.locus_id)}
    return [m[(c, p)] for c, p in zip(T.chromosome_label, T.center_grch37)]


GX = {}


def group_masks(T):
    m = {"all": np.ones(len(T), bool)}
    for v in (True, False):
        m[f"model_linked={v}"] = (T.model_linked == v).values
    for c, lab in (("protein_route", "protein_route"), ("imaging_route", "imaging_route")):
        for v in sorted(T[c].dropna().unique()):
            m[f"{lab}={v}"] = (T[c] == v).values
    for v in (True, False):
        m[f"instrument_in_coloc_snps={v}"] = (T.selnull_w1000k_instrument_in_coloc_snps == v).values
    for a in (True, False):
        for pr in sorted(T.protein_route.dropna().unique()):
            for ir in sorted(T.imaging_route.dropna().unique()):
                m[f"model_linked={a}&protein_route={pr}&imaging_route={ir}"] = ((T.model_linked == a) & (T.protein_route == pr) & (T.imaging_route == ir)).values
    return m


def setup_groups(T):
    n = len(T)
    GX.clear()
    GX["all"] = load_sims(T.region_id.tolist())
    GX["masks"] = group_masks(T)
    GX["exposure"] = pd.factorize(T.exposure_id)[0]
    GX["locus"] = pd.factorize(T.locus_id)[0]
    GX["obs"] = np.full((2, 2, n), np.nan)
    GX["p"] = np.full((2, 2, len(NULL_NAMES), n), np.nan)
    GX["cond"] = {}
    for wi, wt in enumerate(WINDOW_TAGS):
        for pi, (_, tag) in enumerate(PRIORS):
            GX["obs"][wi, pi] = pd.to_numeric(T[f"selnull_{wt}_{tag}_H4_obs_used"], errors="coerce").values
            for ki, nm in enumerate(NULL_NAMES):
                GX["p"][wi, pi, ki] = pd.to_numeric(T[f"selnull_{wt}_{tag}_{nm}_p"], errors="coerce").values
            for th in THRESHOLDS:
                c = f"selnull_{wt}_{tag}_cond_obs_z_P_H4ge{int(th * 100):03d}"
                GX["cond"][(wi, pi, th)] = pd.to_numeric(T[c], errors="coerce").values if c in T else np.full(n, np.nan)
    # z* of each exposure and draw, and one thinning uniform per exposure and draw for the calibrated counts
    ex = list(pd.unique(T.exposure_id))
    U, V = np.empty((len(ex), B)), np.empty((len(ex), B), np.float32)
    for i, e in enumerate(ex):
        U[i] = np.random.default_rng(seed_of("selnull|" + e)).random(B)
        V[i] = np.random.default_rng(seed_of("thin|" + e)).random(B)
    GX["V"] = V
    GX["zbin"] = np.empty((len(NULL_NAMES), len(ex), B), np.int8)
    for ki, (_, S) in enumerate(null_sets()):
        GX["zbin"][ki] = np.clip(np.searchsorted(CAL_EDGES, S[(U * len(S)).astype(int)], side="right") - 1, 0, len(CAL_BINS) - 1)


def _reduce(codes):
    order = np.argsort(codes, kind="stable")
    c = codes[order]
    return order, np.flatnonzero(np.r_[True, c[1:] != c[:-1]])


def _thr(th):
    return f"{int(th * 100):03d}"


def group_task(a):
    gname, wi, pi, ki = a
    nm, tag, wt = NULL_NAMES[ki], PRIORS[pi][1], WINDOW_TAGS[wi]
    m = GX["masks"][gname]
    idx = np.flatnonzero(m & np.isfinite(GX["obs"][wi, pi]))
    base = dict(group=gname, window=wt, p12=P12_LABEL[tag], null_set=nm, null_scale=NULL_SCALE.get(nm, np.nan), null_role=NULL_ROLE[nm], n_regions=int(m.sum()), n_selection_null_computed=len(idx))
    if not len(idx):
        return [dict(base, unit=u, n_units=0) for u in UNITS]
    X = GX["all"][idx, wi, pi, ki]
    obs = GX["obs"][wi, pi][idx].astype(np.float32)
    rp = GX["p"][wi, pi, ki][idx]
    rng = np.random.default_rng(seed_of(f"dep|{gname}|{wt}|{tag}|{nm}"))
    lu, linv = np.unique(GX["locus"][idx], return_inverse=True)
    Ud = rng.random((len(lu), B), dtype=np.float32)[linv]
    Pr = {th: (X >= th).mean(1).astype(np.float32) for th in THRESHOLDS}
    rows = []
    for unit in UNITS:
        if unit == "region":
            Xu, ou, pu, anyp, order, st = X, obs, rp, rp, None, None
        else:
            order, st = _reduce(GX[unit][idx])
            Xu = np.maximum.reduceat(X[order], st, axis=0)
            ou = np.maximum.reduceat(obs[order], st)
            pu = ((Xu >= ou[:, None]).sum(1) + 1) / (B + 1)
            anyp = np.minimum.reduceat(rp[order], st)
        rec = dict(base, unit=unit, n_units=len(ou))
        for th in THRESHOLDS:
            t_ = _thr(th)
            c = (Xu >= th).sum(0)
            o = int((ou >= th).sum())
            D = Ud < Pr[th][:, None]
            if unit != "region":
                D = np.logical_or.reduceat(D[order], st, axis=0)
            cd = D.sum(0)
            rec.update({f"observed_H4ge{t_}": o, f"expected_H4ge{t_}_under_selection": float(c.mean()), f"mc_se_expected_H4ge{t_}": float(c.std(ddof=1) / np.sqrt(B)),
                        f"null_q025_H4ge{t_}": float(np.quantile(c, .025)), f"null_q975_H4ge{t_}": float(np.quantile(c, .975)), f"P_null_count_ge_observed_H4ge{t_}": float(((c >= o).sum() + 1) / (B + 1)),
                        f"expected_H4ge{t_}_within_locus_dependent": float(cd.mean()), f"null_q025_H4ge{t_}_within_locus_dependent": float(np.quantile(cd, .025)),
                        f"null_q975_H4ge{t_}_within_locus_dependent": float(np.quantile(cd, .975)), f"P_null_count_ge_observed_H4ge{t_}_within_locus_dependent": float(((cd >= o).sum() + 1) / (B + 1))})
            if unit == "region":
                rec[f"expected_H4ge{t_}_given_observed_instrument_z"] = float(np.nansum(GX["cond"][(wi, pi, th)][idx]))
        hi = ou >= .8
        rec.update(selection_p_median=float(np.median(pu)), n_selection_p_le_0p05=int((pu <= .05).sum()), n_selection_p_le_0p01=int((pu <= .01).sum()),
                   selection_p_median_among_H4ge080=float(np.median(pu[hi])) if hi.any() else np.nan, n_H4ge080_with_selection_p_le_0p05=int((hi & (pu <= .05)).sum()),
                   n_units_any_region_selection_p_le_0p05=int((anyp <= .05).sum()), n_H4ge080_units_any_region_selection_p_le_0p05=int((hi & (anyp <= .05)).sum()))
        rows.append(rec)
    return rows


# ---------------------------------------------------------------- decoy calibration
def apply_weight_rule(w):
    """Same convention as INT-1-r: the three lowest |z| bins are left at 1, they carry almost no simulated hits and the
    ratios there are unstable. Clean decoys cannot reach the last bin, so it takes the ratio of the bin below."""
    w = np.array(w, float)
    w[..., :3] = 1.0
    w[..., 7] = w[..., 6]
    w[~np.isfinite(w)] = 1.0
    return w


def own_calibration(T, nl):
    rows = []
    for o in nl:
        for r in o.get("decoy_pairs", []):
            rows.append([o["region_id"]] + r)
    D = pd.DataFrame(rows, columns=["region_id"] + CAL_FIELDS)
    loc = dict(zip(T.region_id, GX["locus"]))
    D["locus"] = D.region_id.map(loc)
    D["bin"] = np.clip(np.searchsorted(CAL_EDGES, D.z_inst.values, side="right") - 1, 0, len(CAL_BINS) - 1)
    write_csv(D.drop(columns=["locus"]), OUT / "decoy_pairs_this_simulator.csv.gz")
    L = len(loc) and int(max(loc.values())) + 1
    out, tab = {}, []
    for pi, (_, tag) in enumerate(PRIORS):
        for th in THRESHOLDS:
            k = int(th * 100)
            hl, ml = np.zeros((L, 8)), np.zeros((L, 8))
            np.add.at(hl, (D.locus.values, D.bin.values), (D[f"H4_{tag}"].values >= th).astype(float))
            np.add.at(ml, (D.locus.values, D.bin.values), D[f"P{k}_{tag}"].values)
            h, mm = hl.sum(0), ml.sum(0)
            rng = np.random.default_rng(seed_of(f"calboot|{tag}|{th}"))
            wt = rng.multinomial(L, np.ones(L) / L, size=CAL_BOOT)
            with np.errstate(divide="ignore", invalid="ignore"):
                rep = (wt @ hl) / (wt @ ml)
                raw = h / mm
            out[(pi, th)] = dict(point=apply_weight_rule(raw), reps=apply_weight_rule(rep))
            for b, lab in enumerate(CAL_BINS):
                s = D[D.bin == b]
                tab.append(dict(calibration_source="this_simulator", p12=P12_LABEL[tag], threshold=th, z_bin=lab, n_pairs=len(s), n_regions=s.region_id.nunique(), n_loci=s.locus.nunique(), real_hits=float(h[b]), model_expected=float(mm[b]),
                                ratio_raw=float(raw[b]) if np.isfinite(raw[b]) else np.nan, ratio_boot_q025=float(np.nanquantile(rep[:, b], .025)) if np.isfinite(rep[:, b]).any() else np.nan,
                                ratio_boot_q975=float(np.nanquantile(rep[:, b], .975)) if np.isfinite(rep[:, b]).any() else np.nan, weight_used=float(out[(pi, th)]["point"][b])))
    return out, tab


def int1r_calibration():
    out, tab = {}, []
    for th in THRESHOLDS:
        C = pd.read_csv(INT1R / f"decoy_calibration_H4ge{th}.csv")
        assert tuple(C.z_bin) == CAL_BINS
        h, m = C.real_hits.values.astype(float), C.model_expected.values.astype(float)
        rng = np.random.default_rng(seed_of(f"int1r|{th}"))
        rep = rng.gamma(h + .5, size=(CAL_BOOT, len(h))) / m
        pt = h / m
        pw, rw = pt.copy(), rep.copy()
        pw[:3] = 1
        rw[:, :3] = 1
        out[(0, th)] = dict(point=pw, reps=rw)
        for b, lab in enumerate(CAL_BINS):
            tab.append(dict(calibration_source="INT1r_bins", p12="p12_1e-5", threshold=th, z_bin=lab, n_pairs=int(C.pairs[b]), n_regions=int(C.regions[b]), n_loci=np.nan, real_hits=float(h[b]), model_expected=float(m[b]),
                            ratio_raw=float(pt[b]), ratio_boot_q025=float(np.quantile(rep[:, b], .025)), ratio_boot_q975=float(np.quantile(rep[:, b], .975)), weight_used=float(pw[b])))
    return out, tab


CAL_GROUPS = ("all", "model_linked=True", "model_linked=False")


def calib_task(a):
    gname, pi, ki, th, src = a
    idx = np.flatnonzero(GX["masks"][gname] & np.isfinite(GX["obs"][0, pi]))
    W = GX["cal"][src][(pi, th)]
    X = GX["all"][idx, 0, pi, ki] >= th
    obs = GX["obs"][0, pi][idx] >= th
    r, b = np.nonzero(X)
    ex = GX["exposure"][idx][r]
    # exposure factor codes follow first appearance, the same order as U and V
    bn = GX["zbin"][ki][ex, b]
    v = GX["V"][ex, b]
    rng = np.random.default_rng(seed_of(f"calcount|{gname}|{pi}|{ki}|{th}|{src}"))
    rows = []
    for unit in UNITS:
        codes = np.arange(len(idx)) if unit == "region" else GX[unit][idx]
        uq, uinv = np.unique(codes, return_inverse=True)
        nu = len(uq)
        oo = np.zeros(nu, bool)
        np.logical_or.at(oo, uinv, obs)
        o = int(oo.sum())
        ur = uinv[r]

        def count(w):
            wv = w[bn]
            keep = v < np.minimum(wv, 1)
            A = np.zeros((nu, B), bool)
            A[ur[keep], b[keep]] = True
            c = A.sum(0)
            up = wv > 1
            if up.any():
                E = np.zeros((nu, B), np.float32)
                np.maximum.at(E, (ur[up], b[up]), (wv[up] - 1).astype(np.float32))
                c = c + rng.poisson(E.sum(0))
            return c

        c0 = count(np.ones(len(CAL_BINS)))
        cp = count(W["point"])
        pool = np.stack([count(w) for w in W["reps"]])
        Em = pool.mean(1)
        rows.append(dict(group=gname, unit=unit, n_units=nu, window="w1000k", p12=P12_LABEL[PRIORS[pi][1]], threshold=th, null_set=NULL_NAMES[ki], null_role=NULL_ROLE[NULL_NAMES[ki]], calibration_source=src, observed=o,
                         expected_uncalibrated=float(c0.mean()), P_uncalibrated=float(((c0 >= o).sum() + 1) / (B + 1)), expected_calibrated=float(cp.mean()),
                         expected_calibrated_q025_ratio_uncertainty=float(np.quantile(Em, .025)), expected_calibrated_q975_ratio_uncertainty=float(np.quantile(Em, .975)),
                         null_count_q025_with_ratio_uncertainty=float(np.quantile(pool, .025)), null_count_q975_with_ratio_uncertainty=float(np.quantile(pool, .975)),
                         P_calibrated_with_ratio_uncertainty=float(((pool >= o).mean(1).mean() * B + 1) / (B + 1))))
    return rows


def unit_selection_columns(T):
    """Exposure-level and locus-level selection p at 1 Mb, p12=1e-5: the unit's best H4 against the same maximum over its regions in each draw."""
    obs = GX["obs"][0, 0]
    idx = np.flatnonzero(np.isfinite(obs))
    out = {}
    for unit in ("exposure", "locus"):
        codes = GX[unit]
        order, st = _reduce(codes[idx])
        ou = np.maximum.reduceat(obs[idx][order].astype(np.float32), st)
        ucode = codes[idx][order][st]
        res = {}
        for ki, nm in enumerate(NULL_NAMES):
            if nm not in REFERENCE_NULLS:
                continue
            Xu = np.maximum.reduceat(GX["all"][idx, 0, 0, ki][order], st, axis=0)
            res[nm] = (((Xu >= ou[:, None]).sum(1) + 1) / (B + 1), (Xu >= .8).mean(1))
        pm = np.stack([res[n][0] for n in REFERENCE_NULLS])
        fm = np.stack([res[n][1] for n in REFERENCE_NULLS])
        pmap = {}
        for lab, arr in ((f"{unit}_max_H4_w1000k_1e5", ou.astype(float)), (f"{unit}_selection_p_H4", res[stag(PRIMARY_SCALE)][0]), (f"{unit}_selection_p_H4_ref_min", pm.min(0)), (f"{unit}_selection_p_H4_ref_max", pm.max(0)),
                         (f"{unit}_selection_P_any_H4ge080_ref_min", fm.min(0)), (f"{unit}_selection_P_any_H4ge080_ref_max", fm.max(0))):
            full = np.full(codes.max() + 1, np.nan)
            full[ucode] = arr
            out[lab] = full[codes]
    return pd.DataFrame(out, index=T.index)


def null_admissibility():
    """INT-1-a's max-|z| check: real effects only push an assay's best |z| up, so the observed distribution function should not
    lie above a null's in its lower half. Recomputed here for every null, with INT-1-a's quantile grid."""
    best = pd.read_csv(INTG / "atlas/assay_maxT.csv").best_abs_z.values
    xg = np.quantile(best, [.02, .05, .1, .2, .3, .4, .5])
    Fo = np.array([(best <= x).mean() for x in xg])
    se = np.sqrt(Fo * (1 - Fo) / len(best))
    prev = pd.read_csv(INT1A / "atlas_null_admissibility.csv").set_index("null").z_of_excess.to_dict()
    alias = {stag(1.0): "analyst_scale_1.0", stag(1.025): "analyst_scale_1.025", stag(1.05): "analyst_scale_1.05", stag(1.10): "analyst_scale_1.1"}
    rows = []
    for nm, S in null_sets():
        F0 = np.searchsorted(S, xg) / len(S)
        rows.append(dict(null_set=nm, label=NULL_LABEL[nm], null_role=NULL_ROLE[nm], n_draws=len(S), null_q05=float(np.quantile(S, .05)), null_median=float(np.median(S)), null_q95=float(np.quantile(S, .95)),
                         max_excess_Fobs_minus_Fnull=float((Fo - F0).max()), z_of_excess=float(((Fo - F0) / se).max()), z_of_excess_in_INT1a_table=prev.get(alias.get(nm, nm), np.nan)))
    return pd.DataFrame(rows)


def reverse_cis_summary(T):
    rows = []
    long = []
    for pol in ("drop", "infer"):
        for eff in ("discovery", "replication"):
            t = f"rev_{pol}_{eff}"
            d = pd.DataFrame(dict(region_id=T.region_id, exposure_id=T.exposure_id, model_linked=T.model_linked, policy=pol, effect_set=eff, status=T[t + "_status"], p=T[t + "_p"], p_random=T[t + "_p_random"],
                                  n_cis=T[t + "_n_cis_instruments"], ns=T[t + "_nocis_status"], pn=T[t + "_nocis_p"], pnr=T[t + "_nocis_p_random"]))
            long.append(d[d.status == "estimated"])
    E = pd.concat(long, ignore_index=True)
    E["p_re"] = E.p_random.fillna(E.p)
    E["pn_re"] = E.pnr.fillna(E.pn)

    def rate(s):
        return float((s < .05).mean()) if len(s) else np.nan

    for scope, m in (("all", np.ones(len(E), bool)), ("model_linked", E.model_linked.values), ("not_model_linked", ~E.model_linked.values)):
        for (pol, eff), g in list(E[m].groupby(["policy", "effect_set"])) + [(("all_four", "all_four"), E[m])]:
            a, c = g[g.n_cis == 0], g[g.n_cis > 0]
            kept = c[c.ns == "estimated"]
            rows.append(dict(scope=scope, palindrome_policy=pol, imaging_effect_set=eff, n_estimates=len(g), n_without_cis_instrument=len(a), n_with_cis_instrument=len(c), n_regions_with_cis_instrument=c.region_id.nunique(), n_exposures_with_cis_instrument=c.exposure_id.nunique(),
                             frac_p_lt_0p05_without_cis_fixed=rate(a.p), frac_p_lt_0p05_without_cis_random=rate(a.p_re), frac_p_lt_0p05_with_cis_fixed=rate(c.p), frac_p_lt_0p05_with_cis_random=rate(c.p_re),
                             n_with_cis_still_estimable_outside_cis=len(kept), n_with_cis_no_instrument_outside_cis=len(c) - len(kept), frac_p_lt_0p05_same_estimates_before_removal_fixed=rate(kept.p), frac_p_lt_0p05_same_estimates_before_removal_random=rate(kept.p_re),
                             frac_p_lt_0p05_after_removal_fixed=rate(kept.pn), frac_p_lt_0p05_after_removal_random=rate(kept.pn_re)))
    return pd.DataFrame(rows)


def range_summary(G, C):
    rows = []
    for gname in CAL_GROUPS:
        for unit in UNITS:
            for wt in WINDOW_TAGS:
                for _, tag in PRIORS:
                    lab = P12_LABEL[tag]
                    g = G[(G.group == gname) & (G.unit == unit) & (G.window == wt) & (G.p12 == lab)]
                    r = g[g.null_role == "reference_range"]
                    if not len(r):
                        continue
                    for th in THRESHOLDS:
                        t_ = _thr(th)
                        e, pp, pd_ = r[f"expected_H4ge{t_}_under_selection"], r[f"P_null_count_ge_observed_H4ge{t_}"], r[f"P_null_count_ge_observed_H4ge{t_}_within_locus_dependent"]
                        rec = dict(group=gname, unit=unit, n_units=int(r.n_units.iloc[0]), window=wt, p12=lab, threshold=th, observed=int(r[f"observed_H4ge{t_}"].iloc[0]),
                                   expected_ref_min=float(e.min()), expected_ref_max=float(e.max()), P_ref_min=float(pp.min()), P_ref_max=float(pp.max()), P_within_locus_dependent_ref_min=float(pd_.min()), P_within_locus_dependent_ref_max=float(pd_.max()))
                        for nm in (stag(1.0), stag(1.05), stag(1.10)):
                            x = g[g.null_set == nm]
                            rec[f"expected_{nm}"] = float(x[f"expected_H4ge{t_}_under_selection"].iloc[0])
                            rec[f"P_{nm}"] = float(x[f"P_null_count_ge_observed_H4ge{t_}"].iloc[0])
                        c = C[(C.group == gname) & (C.unit == unit) & (C.window == wt) & (C.p12 == lab) & np.isclose(C.threshold, th) & (C.calibration_source == "this_simulator") & (C.null_role == "reference_range")] if wt == "w1000k" else C.iloc[0:0]
                        if len(c):
                            rec.update(expected_calibrated_ref_min=float(c.expected_calibrated.min()), expected_calibrated_ref_max=float(c.expected_calibrated.max()),
                                       expected_calibrated_ratio_uncertainty_low=float(c.expected_calibrated_q025_ratio_uncertainty.min()), expected_calibrated_ratio_uncertainty_high=float(c.expected_calibrated_q975_ratio_uncertainty.max()),
                                       P_calibrated_ref_min=float(c.P_calibrated_with_ratio_uncertainty.min()), P_calibrated_ref_max=float(c.P_calibrated_with_ratio_uncertainty.max()))
                        rows.append(rec)
    return pd.DataFrame(rows)


def stage_aggregate(workers):
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    pl = plan()
    regs = pl["regions"]
    iph = imaging_plan()["phenotypes"]
    rids = sorted(regs)
    assert len(rids) == 1375
    cov = pd.read_csv(ACC / "per_region_coverage.csv", keep_default_na=False)
    assert len(cov) == 1375 and cov.region_id.is_unique
    nulls_present = [r for r in rids if (NULLDIR / f"{r}.json").exists() and (NULLDIR / f"{r}.npy").exists()]
    if len(nulls_present) != len(rids):
        raise SystemExit(f"null cache incomplete: {len(nulls_present)}/{len(rids)}; run --stage null first")
    global PENDING_X
    pk = pd.read_csv(REG / "reverse_LD/all_unique_peak_coverage.csv", usecols=["variant_key", "reference_request_status"], low_memory=False)
    PENDING_X = set(pk.variant_key[pk.reference_request_status.eq("pending_X_ploidy_method_review")].astype(str))
    REGION_CENTER.clear()
    REGION_CENTER.update({r: (chrom_label(regs[r]["chromosome"]), int(regs[r]["center_grch37"])) for r in rids})
    with Pool(workers) as po:
        rr = po.map(read_region, rids, chunksize=8)
        nl = po.map(_read_null, rids, chunksize=16)
        hashes = po.map(_hash_region, rids, chunksize=8)
    grid = rr[0][1]
    C = pd.DataFrame([r[0] for r in rr])
    N = pd.DataFrame([null_columns(o) for o in nl])
    ins = pd.read_csv(MR / "results/full_atlas_20261116/full/selected_cis_instruments.csv", low_memory=False).set_index("exposure_id")
    maxT = pd.read_csv(INTG / "atlas/assay_maxT.csv").set_index("exposure_id")
    cal = pd.read_csv(MRI / "assay_level_calibration.csv.gz").set_index("exposure_id")
    units = pd.read_csv(PF6 / "route_ledger_units.csv", low_memory=False, dtype={"unit": str})
    uu = {(r.kind, str(r.unit)): r for r in units.itertuples(index=False)}
    audits, targeted = load_audits()
    rows = []
    for rid in rids:
        r = regs[rid]
        e = r["exposure_ids"][0]
        ph = iph[r["phenotype"]]
        i = ins.loc[e]
        u_p = uu[("protein", r["accession"])]
        o = dict(region_id=rid, accession=r["accession"], gene=r["gene"], chromosome=r["chromosome"], chromosome_label=chrom_label(r["chromosome"]), center_grch37=int(r["center_grch37"]),
                 phenotype=r["phenotype"], phenotype_description=r["phenotype_description"], family=ph["family"], idp=ph["idp"], model_linked=bool(r["model_linked"]),
                 chrX=chrom_label(r["chromosome"]) == "X", exposure_id=e, selection_policies="|".join(r["policies"]), instrument_variant=i.variant, instrument_pos_grch37=int(i.pos37),
                 instrument_effect_allele=i.ea_x, instrument_other_allele=i.oa_x, instrument_class=i.variant_class, instrument_is_palindromic=bool(i.is_palindromic_snv), instrument_MHC=i.MHC,
                 reverse_peak_source_status=r["reverse_peak_source_status"])
        o["protein_route"] = u_p.route
        o["protein_whole_file_read"] = bool(u_p.whole_file_read)
        o["protein_md5_verified_at_acquisition"] = u_p.md5_verified if pd.notna(u_p.md5_verified) else np.nan
        st, src = audit_status("protein", r["accession"], u_p.route, audits, targeted)
        o["protein_audit_status"], o["protein_audit_source"] = st, src
        if ph["family"] == "BIG40":
            u_i = uu[("imaging", ph["idp"])]
            o["imaging_route"] = u_i.route
            o["imaging_whole_file_read"] = bool(u_i.whole_file_read)
            st, src = audit_status("imaging", ph["idp"], u_i.route, audits, targeted)
        else:
            c = json.loads((IMG / ph["family"] / ph["idp"] / "COMPLETE.json").read_text())
            o["imaging_route"] = "qsm_archive_member_full_read"
            o["imaging_whole_file_read"] = True
            st, src = ("whole_archive_member_read_gzip_crc_verified" if c.get("gzip_crc_eof_verified") else "qsm_member_read_crc_not_recorded"), ""
        o["imaging_audit_status"], o["imaging_audit_source"] = st, src
        o["both_units_whole_file_or_audit_match"] = all(s in ("whole_file_read_at_acquisition", "range_audited_match", "whole_archive_member_read_gzip_crc_verified") for s in (o["protein_audit_status"], o["imaging_audit_status"]))
        if e in cal.index:
            o.update(assay_p=float(cal.loc[e, "assay_p"]), assay_q=float(cal.loc[e, "assay_q"]), assay_max_abs_z=float(cal.loc[e, "max_abs_z"]), assay_q_status="available")
        else:
            o.update(assay_p=np.nan, assay_q=np.nan, assay_max_abs_z=np.nan, assay_q_status="exposure_absent_from_drop_policy_calibration")
        o.update(p_maxT=float(maxT.loc[e, "p_maxT"]), q_maxT_all_assays=float(maxT.loc[e, "q_maxT_all_assays"]), maxT_best_abs_z=float(maxT.loc[e, "best_abs_z"]),
                 maxT_best_phenotype=maxT.loc[e, "best_phenotype"], region_phenotype_is_maxT_best=bool(maxT.loc[e, "best_phenotype"] == r["phenotype"]))
        rows.append(o)
    I = pd.DataFrame(rows)
    cc = cov[["region_id", "protein_state", "protein_reason", "protein_error", "imaging_state", "imaging_window_present", "coloc_state", "coloc_reason", "coloc_window_status", "reverse_state", "reverse_reason", "reverse_no_published_peak", "reverse_status_counts", "blocked_by"]].rename(columns={"reverse_status_counts": "acceptance_reverse_status_counts"})
    T = I.merge(cc, on="region_id", validate="one_to_one").merge(C, on="region_id", validate="one_to_one").merge(N, on="region_id", validate="one_to_one")
    T["region_complete"] = (T.protein_state == "complete") & (T.imaging_state == "complete") & (T.coloc_state == "complete") & (T.reverse_state == "complete")
    # counting units
    T["n_regions_for_exposure"] = T.groupby("exposure_id").region_id.transform("size")
    T["n_accessions_for_exposure"] = T.groupby("exposure_id").accession.transform("nunique")
    T["n_exposures_for_accession"] = T.groupby("accession").exposure_id.transform("nunique")
    T["locus_id"] = assign_loci(T)
    T["n_regions_in_locus"] = T.groupby("locus_id").region_id.transform("size")
    T["n_exposures_in_locus"] = T.groupby("locus_id").exposure_id.transform("nunique")
    p5 = fmt_p12(1e-5)
    p6 = fmt_p12(1e-6)
    p3 = fmt_p12(3e-5)
    for wt in WINDOW_TAGS:
        for lab, p in (("1e5", p5), ("1e6", p6), ("3e5", p3)):
            for h in range(5):
                T[f"H{h}_{wt}_{lab}"] = T[f"coloc_{wt}_p12_{p}_PPH{h}"]
    ps = stag(PRIMARY_SCALE)
    T["selection_p_H4"] = T[f"selnull_w1000k_1e5_{ps}_p"]
    T["expected_H4_under_selection"] = T[f"selnull_w1000k_1e5_{ps}_expected_H4"]
    T["selection_P_H4ge080"] = T[f"selnull_w1000k_1e5_{ps}_P_H4ge080"]
    for lab, suf in (("selection_p_H4", "p"), ("expected_H4_under_selection", "expected_H4"), ("selection_P_H4ge080", "P_H4ge080")):
        m = T[[f"selnull_w1000k_1e5_{n}_{suf}" for n in REFERENCE_NULLS]]
        T[lab + "_ref_min"] = m.min(axis=1)
        T[lab + "_ref_max"] = m.max(axis=1)
    T["expected_H4_given_observed_instrument_z"] = T.get("selnull_w1000k_1e5_cond_obs_z_expected_H4", np.nan)
    T["selection_p_H4_given_observed_instrument_z"] = T.get("selnull_w1000k_1e5_cond_obs_z_p", np.nan)
    T["selection_null_computable"] = T.selnull_status == "computed"
    T["selection_null_evidence_tier"] = np.where(T.selection_null_computable, np.where(T.model_linked, "B", "A*"), "")
    T["instrument_placement"] = np.where(T.selnull_w1000k_instrument_in_coloc_snps == True, "instrument_in_coloc_snps", np.where(T.selection_null_computable, "instrument_outside_coloc_snps_conditioned_through_panel", "not_computable"))
    # repro columns
    T["repro_available"] = False
    repro_cols = []
    repro_note = "region_repro.csv.gz not present when this table was written"
    if REPRO.exists():
        Rp = pd.read_csv(REPRO, low_memory=False)
        if "region_id" in Rp and Rp.region_id.is_unique:
            ren = {c: (c if c.startswith("repro_") else "repro_" + c) for c in Rp.columns if c != "region_id"}
            Rp = Rp.rename(columns=ren)
            repro_cols = list(ren.values())
            T = T.merge(Rp, on="region_id", how="left", validate="one_to_one")
            T["repro_available"] = T.region_id.isin(Rp.region_id)
            repro_note = f"merged from region_repro.csv.gz, {int(T.repro_available.sum())} of 1375 regions"
        else:
            repro_note = "region_repro.csv.gz present but has no unique region_id column; not merged"
    T = T.sort_values("region_id").reset_index(drop=True)
    assert len(T) == 1375 and T.region_id.is_unique
    # ---- groups, calibration, unit-level selection columns
    setup_groups(T)
    h = hashlib.sha256(f"{GROUP_VERSION}|{null_digest()}|{B}|{THRESHOLDS}|{CAL_BOOT}|{NULL_NAMES}".encode())
    h.update(GX["obs"].tobytes())
    h.update(GX["locus"].tobytes())
    h.update(GX["exposure"].tobytes())
    for k in sorted(GX["masks"]):
        h.update(k.encode())
        h.update(GX["masks"][k].tobytes())
    for th in THRESHOLDS:
        h.update(sha256(INT1R / f"decoy_calibration_H4ge{th}.csv").encode())
    key = h.hexdigest()[:20]
    fG, fC, fB = CACHE / f"groups_{key}.csv.gz", CACHE / f"calibrated_{key}.csv.gz", CACHE / f"calbins_{key}.csv.gz"
    if fG.exists() and fC.exists() and fB.exists():
        G, Cal, Bins = pd.read_csv(fG), pd.read_csv(fC), pd.read_csv(fB)
        print("group tables reused from cache", key, flush=True)
    else:
        for old in list(CACHE.glob("groups_*.csv.gz")) + list(CACHE.glob("calibrated_*.csv.gz")) + list(CACHE.glob("calbins_*.csv.gz")):
            os.remove(old)
        own, tab1 = own_calibration(T, nl)
        ir, tab2 = int1r_calibration()
        GX["cal"] = {"this_simulator": own, "INT1r_bins": ir}
        tasks = [(g, wi, pi, ki) for g in GX["masks"] for wi in range(2) for pi in range(2) for ki in range(len(NULL_NAMES))]
        ct = [(g, pi, ki, th, "this_simulator") for g in CAL_GROUPS for pi in range(2) for ki in range(len(NULL_NAMES)) for th in THRESHOLDS]
        ct += [(g, 0, ki, th, "INT1r_bins") for g in CAL_GROUPS for ki in range(len(NULL_NAMES)) for th in THRESHOLDS]
        with Pool(workers) as po:
            G = pd.DataFrame([r for rs in po.imap(group_task, tasks, chunksize=4) for r in rs])
            print("group rows", len(G), round(time.time() - t0), "s", flush=True)
            Cal = pd.DataFrame([r for rs in po.imap(calib_task, ct, chunksize=2) for r in rs])
            print("calibrated rows", len(Cal), round(time.time() - t0), "s", flush=True)
        Bins = pd.DataFrame(tab1 + tab2)
        write_csv(G, fG)
        write_csv(Cal, fC)
        write_csv(Bins, fB)
    U = unit_selection_columns(T)
    T = pd.concat([T, U], axis=1)
    first = ["region_id", "accession", "gene", "chromosome", "chromosome_label", "chrX", "center_grch37", "phenotype", "phenotype_description", "family", "idp", "model_linked", "exposure_id", "n_regions_for_exposure", "n_accessions_for_exposure", "locus_id", "n_regions_in_locus", "n_exposures_in_locus", "selection_policies",
             "region_complete", "protein_state", "imaging_state", "coloc_state", "reverse_state", "protein_route", "imaging_route", "protein_audit_status", "imaging_audit_status", "H4_w1000k_1e5", "H4_w500k_1e5", "H4_w1000k_1e6", "H4_w500k_1e6",
             "assay_q", "p_maxT", "selection_p_H4", "selection_p_H4_ref_min", "selection_p_H4_ref_max", "expected_H4_under_selection", "expected_H4_under_selection_ref_min", "expected_H4_under_selection_ref_max", "selection_P_H4ge080", "selection_P_H4ge080_ref_min", "selection_P_H4ge080_ref_max",
             "exposure_max_H4_w1000k_1e5", "exposure_selection_p_H4", "exposure_selection_p_H4_ref_min", "exposure_selection_p_H4_ref_max", "locus_max_H4_w1000k_1e5", "locus_selection_p_H4", "locus_selection_p_H4_ref_min", "locus_selection_p_H4_ref_max",
             "expected_H4_given_observed_instrument_z", "selection_p_H4_given_observed_instrument_z", "selection_null_computable", "selection_null_evidence_tier", "selnull_status", "selnull_reason", "instrument_placement", "rev_n_estimated", "rev_n_no_usable_instrument", "rev_n_no_published_peak", "rev_n_estimates_with_cis_instrument", "rev_chrX_status", "repro_available"]
    T = T[first + [c for c in T.columns if c not in first]]
    write_csv(T, OUT / "all_target_region_table.csv.gz")
    write_csv(G, OUT / "h4_vs_selection_by_model_route.csv")
    write_csv(Cal, OUT / "h4_vs_selection_decoy_calibrated.csv")
    write_csv(Bins, OUT / "decoy_calibration_bins.csv")
    RS = range_summary(G, Cal)
    write_csv(RS, OUT / "h4_reference_range_summary.csv")
    AD = null_admissibility()
    write_csv(AD, OUT / "null_admissibility.csv")
    # exposure and locus tables
    H4 = "H4_w1000k_1e5"
    ex_rows, lo_rows = [], []
    for e, g in T.groupby("exposure_id", sort=True):
        g = g.sort_values(H4, ascending=False)
        ex_rows.append(dict(exposure_id=e, gene="|".join(sorted(set(g.gene))), phenotype=g.phenotype.iloc[0], family=g.family.iloc[0], idp=g.idp.iloc[0], model_linked=bool(g.model_linked.any()), locus_id=g.locus_id.iloc[0], instrument_variant=g.instrument_variant.iloc[0], n_regions=len(g),
                            accessions="|".join(g.accession), H4_by_accession_w1000k_1e5="|".join(f"{x:.4g}" for x in g[H4]), H4_max=float(g[H4].max()), H4_min=float(g[H4].min()), accessions_on_both_sides_of_H4_0p8=bool(g[H4].max() >= .8 and g[H4].min() < .8), p_maxT=float(g.p_maxT.iloc[0]), assay_q=g.assay_q.iloc[0],
                            exposure_selection_p_H4=g.exposure_selection_p_H4.iloc[0], exposure_selection_p_H4_ref_min=g.exposure_selection_p_H4_ref_min.iloc[0], exposure_selection_p_H4_ref_max=g.exposure_selection_p_H4_ref_max.iloc[0],
                            exposure_selection_P_any_H4ge080_ref_min=g.exposure_selection_P_any_H4ge080_ref_min.iloc[0], exposure_selection_P_any_H4ge080_ref_max=g.exposure_selection_P_any_H4ge080_ref_max.iloc[0]))
    for l, g in T.groupby("locus_id", sort=False):
        b = g.sort_values(H4, ascending=False).iloc[0]
        lo_rows.append(dict(locus_id=l, chromosome_label=g.chromosome_label.iloc[0], first_instrument_pos=int(g.center_grch37.min()), last_instrument_pos=int(g.center_grch37.max()), n_regions=len(g), n_exposures=g.exposure_id.nunique(), n_instruments=g.instrument_variant.nunique(), model_linked_any=bool(g.model_linked.any()), genes="|".join(sorted(set(g.gene))),
                            H4_max=float(g[H4].max()), region_with_H4_max=b.region_id, gene_with_H4_max=b.gene, locus_selection_p_H4=g.locus_selection_p_H4.iloc[0], locus_selection_p_H4_ref_min=g.locus_selection_p_H4_ref_min.iloc[0], locus_selection_p_H4_ref_max=g.locus_selection_p_H4_ref_max.iloc[0]))
    EX, LO = pd.DataFrame(ex_rows), pd.DataFrame(lo_rows).sort_values("locus_id")
    write_csv(EX, OUT / "all_target_exposure_table.csv.gz")
    write_csv(LO, OUT / "all_target_locus_table.csv.gz")
    long_rows = []
    for o in nl:
        for wt, res in o["results"].items():
            for key_, x in res.items():
                tag, name = key_.split("|")
                long_rows.append(dict(region_id=o["region_id"], window=wt, p12=P12_LABEL[tag], null_set=name, **x))
    write_csv(pd.DataFrame(long_rows), OUT / "selection_null_by_region_long.csv.gz")
    # reverse tables
    rv = []
    for pol in ("drop", "infer"):
        for eff in ("discovery", "replication"):
            c = f"rev_{pol}_{eff}_status"
            for scope, m in (("all", T.index == T.index), ("model_linked", T.model_linked), ("chrX_protein_region", T.chrX)):
                vc = T.loc[m, c].value_counts()
                for k, n in vc.items():
                    rv.append(dict(palindrome_policy=pol, imaging_effect_set=eff, scope=scope, status=k, n_regions=int(n)))
    for scope, m in (("all", T.index == T.index), ("model_linked", T.model_linked), ("chrX_protein_region", T.chrX)):
        for k, n in T.loc[m, "rev_chrX_status"].value_counts().items():
            rv.append(dict(palindrome_policy="any", imaging_effect_set="chrX_peak_status", scope=scope, status=k, n_regions=int(n)))
    write_csv(pd.DataFrame(rv), OUT / "reverse_status_table.csv")
    RC = reverse_cis_summary(T)
    write_csv(RC, OUT / "reverse_cis_summary.csv")
    # coverage summary
    cs = []

    def add(section, item, value):
        cs.append(dict(section=section, item=item, value=value))

    add("regions", "planned", len(T))
    add("regions", "unique_region_id_rows", int(T.region_id.nunique()))
    add("regions", "complete_all_four_stages", int(T.region_complete.sum()))
    add("regions", "model_linked", int(T.model_linked.sum()))
    add("regions", "chrX", int(T.chrX.sum()))
    add("units", "regions", len(T))
    add("units", "AGES_protein_accessions", int(T.accession.nunique()))
    add("units", "UKB_exposures", int(T.exposure_id.nunique()))
    add("units", "genes", int(T.gene.nunique()))
    add("units", "instruments", int(T.instrument_variant.nunique()))
    add("units", "loci_instruments_chained_within_1Mb", int(T.locus_id.nunique()))
    for k, n in T.groupby("exposure_id").size().value_counts().sort_index().items():
        add("units", f"exposures_with_{k}_regions", int(n))
    add("units", "exposures_model_linked", int(T[T.model_linked].exposure_id.nunique()))
    add("units", "exposures_not_model_linked", int(T[~T.model_linked].exposure_id.nunique()))
    add("units", "accessions_with_more_than_one_exposure", int(T[T.n_exposures_for_accession > 1].accession.nunique()))
    for c in ("family", "protein_state", "imaging_state", "coloc_state", "reverse_state", "protein_route", "imaging_route", "protein_audit_status", "imaging_audit_status", "coloc_w1000k_status", "coloc_w500k_status", "assay_q_status", "selnull_status", "selnull_reason", "instrument_placement", "rev_chrX_status"):
        for k, n in T[c].fillna("").astype(str).value_counts().items():
            add(c, k if k else "(blank)", int(n))
    add("assay", "regions_with_assay_q", int(T.assay_q.notna().sum()))
    add("assay", "regions_with_p_maxT", int(T.p_maxT.notna().sum()))
    add("assay", "regions_where_phenotype_is_best_in_maxT_table", int(T.region_phenotype_is_maxT_best.sum()))
    add("reverse_chrX", "regions_with_chrX_published_peaks", int((T.rev_n_chrX_peaks > 0).sum()))
    add("reverse_chrX", "chrX_peak_links_total", int(T.rev_n_chrX_peaks.sum()))
    add("reverse_chrX", "chrX_peak_links_reference_matched", int(T.rev_n_chrX_peaks_reference_matched.sum()))
    add("reverse_chrX", "chrX_peak_links_listed_pending_in_request_table", int(T.rev_n_chrX_peaks_listed_pending_in_request_table.sum()))
    add("reverse_chrX", "regions_with_chrX_instrument_in_estimates", int((T.rev_n_chrX_instruments_in_estimates > 0).sum()))
    add("reverse_cis", "regions_with_published_peak_in_cis_window", int((T.rev_n_published_peaks_in_cis_window > 0).sum()))
    add("reverse_cis", "regions_with_estimate_using_cis_instrument", int((T.rev_n_estimates_with_cis_instrument > 0).sum()))
    add("reverse_cis", "estimates_using_cis_instrument", int(T.rev_n_estimates_with_cis_instrument.sum()))
    add("reverse_cis", "max_abs_recomputed_minus_production_beta", float(T.rev_max_abs_recomputed_minus_production_beta.max()))
    add("reverse_cis", "max_abs_recomputed_minus_production_p", float(T.rev_max_abs_recomputed_minus_production_p.max()))
    add("reverse_cis", "regions_where_reclump_differs_from_production", int((~T.rev_reclump_reproduces_production.astype(bool)).sum()))
    add("repro", "regions_with_repro_row", int(T.repro_available.sum()))
    add("repro", "note", repro_note)
    add("units", "protein_units_in_targeted_audit_queue", len(targeted["protein"]))
    add("units", "imaging_units_in_targeted_audit_queue", len(targeted["imaging"]))
    add("units", "protein_targeted_audit_results_present", sum(1 for v in audits["protein"].values() if any(s.startswith("targeted") for s, _ in v)))
    add("units", "imaging_targeted_audit_results_present", sum(1 for v in audits["imaging"].values() if any(s.startswith("targeted") for s, _ in v)))
    for th in THRESHOLDS:
        add("coloc_w1000k_p12_1e-5", f"regions_H4_ge_{th}", int((T[H4] >= th).sum()))
        add("coloc_w1000k_p12_1e-5", f"exposures_H4_ge_{th}_in_any_region", int(T[T[H4] >= th].exposure_id.nunique()))
        add("coloc_w1000k_p12_1e-5", f"loci_H4_ge_{th}_in_any_region", int(T[T[H4] >= th].locus_id.nunique()))
    hi = T[T[H4] >= .8]
    add("coloc_w1000k_p12_1e-5", "regions_H4_ge_0.8_both_units_whole_file_or_audit_match", int(hi.both_units_whole_file_or_audit_match.sum()))
    add("coloc_w1000k_p12_1e-5", "regions_H4_ge_0.8_at_least_one_range_unit_without_whole_file_audit", int((~hi.both_units_whole_file_or_audit_match).sum()))
    add("coloc_w1000k_p12_1e-5", "exposures_with_accessions_on_both_sides_of_H4_0.8", int(EX.accessions_on_both_sides_of_H4_0p8.sum()))
    add("selection_null", "snps_in_coloc_windows_without_panel_z_1000k", int(T.selnull_w1000k_n_snps_without_panel_z.sum()))
    add("selection_null", "min_clean_decoys_1000k", int(T.selnull_w1000k_n_clean_decoys.min()))
    add("selection_null", "max_abs_recomputed_minus_production_H4_1000k_p12_1e-5", float(T.selnull_w1000k_1e5_recomputed_minus_production_H4.abs().max()))
    add("selection_null", "max_abs_own_phenotype_z_difference_panel_vs_production", float(T.selnull_w1000k_own_phenotype_max_abs_z_diff.max()))
    add("selection_null", "decoy_pairs_for_calibration", int(sum(len(o.get("decoy_pairs", [])) for o in nl)))
    write_csv(pd.DataFrame(cs), OUT / "coverage_summary.csv")
    # identities
    H = pd.DataFrame(hashes).sort_values("region_id")
    write_csv(H, OUT / "per_region_input_hashes.csv.gz")
    named = [REG / "PLAN.json", IMG / "PLAN.json", ACC / "ACCEPTANCE.json", ACC / "per_region_coverage.csv", REG / "all_worklist_execution_links.csv", COVERAGE / "MANIFEST.json",
             PF6 / "route_ledger_units.csv", PF6 / "route_ledger_regions.csv", AUD / "audit/results_protein.jsonl", AUD / "audit/results_imaging.jsonl", AUD / "audit_targeted_H4_20261119/results_protein.jsonl",
             AUD / "audit_targeted_H4_20261119/results_imaging.jsonl", AUD / "audit_targeted_H4_20261119/MANIFEST.json", INTG / "atlas/maxnull_draws.npy", INTG / "atlas/assay_maxT.csv", INTG / "atlas/z_infer.npy", INTG / "atlas/exposures.csv", INTG / "atlas/phenotypes.csv",
             MRI / "assay_level_calibration.csv.gz", MRI / "region_table_with_assay_q.csv.gz", MRI / "null_with_ld_by_region.csv", REG / "reverse_LD/all_unique_peak_coverage.csv", REG / "reverse_supplement/PLAN.json", REG / "reverse_supplement/reference_LD/chr23/MANIFEST.json",
             INT1A / "atlas_null_admissibility.csv", INT1A / "null_lambda_median_family.npy", INT1A / "null_lambda_median_each.npy", INT1A / "null_lambda_meanchi2_each.npy", INT1R / "own_maxnull_draws.npz", INT1R / "own_null_final_counts.csv",
             INT1R / "decoy_calibration_H4ge0.5.csv", INT1R / "decoy_calibration_H4ge0.8.csv", INT1R / "decoy_calibration_H4ge0.9.csv", MRI / "s12_null_expectation.py", MRI / "s14_ld_prepare.py", MRI / "s16_null_with_ld.py",
             INTG / "coloc_selection_null_122.csv", INTG / "coloc_expected_under_selection_null.csv", INTG / "coloc_null_scale_sensitivity.csv", INTG / "DELIVERY_table_122_model_linked_regions.csv", INTG / "coloc_null.py", INTG / "null_scale.py",
             MR / "results/full_atlas_20261116/full/selected_cis_instruments.csv", MR / "run_one_all_target_coloc.py", MR / "all_target_coloc.R", MR / "run_one_all_target_reverse.py", MR / "reverse_followup_core.py", ROOT / "pdv3/mr.py", Path(__file__).resolve()]
    if REPRO.exists():
        named.append(REPRO)
    ident = dict(script=str(Path(__file__).resolve()), recipe=RECIPE, group_tables=GROUP_VERSION, group_cache_key=key, null_digest=null_digest(), null_sets=list(NULL_NAMES), reference_nulls=list(REFERENCE_NULLS), files={str(p): sha256(p) for p in named}, per_region_hash_table="per_region_input_hashes.csv.gz",
                 per_region_hash_table_sha256=sha256(OUT / "per_region_input_hashes.csv.gz"), audit_result_lines=dict(AUDIT_LINES),
                 panel=json.loads(PANEL_MANIFEST.read_text()), panel_z_sha256=sha256(PANEL), repro_status=repro_note,
                 pinned_not_edited=dict(acquire_all_target_protein_regions=sha256(MR / "acquire_all_target_protein_regions.py")))
    atomic_json(OUT / "SOURCE_IDENTITIES.json", ident)
    pass  # the internal README writer is not part of the public copy
    print("aggregate done", T.shape, round(time.time() - t0), "s", flush=True)
    return T, G


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["all", "panel", "null", "validate", "aggregate"], default="all")
    ap.add_argument("--refresh", action="store_true", help="only rewrite outputs from cached null results and current inputs")
    ap.add_argument("--workers", type=int, default=MAX_WORKERS)
    ap.add_argument("--only", default="", help="comma-separated region ids (null stage only)")
    ap.add_argument("--guard-selftest", action="store_true", help="confirm that writes outside the output directory are refused, then exit")
    a = ap.parse_args()
    if a.guard_selftest:
        t = MR / "results/_aggregate_guard_selftest_should_not_exist"
        try:
            open(t, "w").close()
        except PermissionError as e:
            print("refused as intended:", e)
            print("exists after attempt:", t.exists())
            return
        raise SystemExit("guard failed")
    w = min(a.workers, MAX_WORKERS)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.refresh:
        stage_aggregate(w)
        return
    if a.stage in ("all", "panel"):
        stage_panel(w)
    if a.stage in ("all", "null"):
        stage_null(w, set(a.only.split(",")) if a.only else None)
    if a.stage in ("all", "validate"):
        stage_validate(w)
    if a.stage in ("all", "aggregate"):
        stage_aggregate(w)


if __name__ == "__main__":
    main()
