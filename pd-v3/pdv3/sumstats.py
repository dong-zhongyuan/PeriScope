"""Canonical summary-statistic loading with full-row validation audit.

Columns: variant_id, effect_allele, other_allele, beta, se; optional eaf, p,
n, build, chromosome, position, protein, outcome. Rename source-specific fields
explicitly via read_sumstats(column_map={source: canonical}). No inferred beta
signs, genome builds, units, phenotype IDs, or silent duplicate selection.
"""
from pathlib import Path
import numpy as np
import pandas as pd

REQUIRED = ("variant_id", "effect_allele", "other_allele", "beta", "se")


def validate_sumstats(df, key_cols=("variant_id",)):
    """Return a copy with valid/validation_errors, preserving every row.

    Schema errors raise; row errors accumulate. Duplicate keys mark ALL copies
    invalid, even identical copies. Include protein/outcome in key_cols for long
    tables. Optional missing EAF/P/N are retained, with no imputation.
    """
    missing = set(REQUIRED).union(key_cols) - set(df.columns)
    if missing:
        raise ValueError(f"missing required columns: {sorted(missing)}")
    if df.columns.duplicated().any():
        raise ValueError("duplicate column names")
    out = df.copy()
    issues = [[] for _ in range(len(out))]

    def mark(mask, reason):
        for i in np.flatnonzero(np.asarray(mask, dtype=bool)):
            issues[i].append(reason)

    for col in key_cols:
        mark(out[col].isna() | out[col].astype(str).str.strip().eq(""), f"missing_{col}")
    mark(out.duplicated(list(key_cols), keep=False), "duplicate_key")
    for col in ("effect_allele", "other_allele"):
        out[col] = out[col].astype("string").str.upper().str.strip()
        mark(~out[col].isin(list("ACGT")), f"invalid_{col}")
    mark(out.effect_allele.eq(out.other_allele).fillna(False), "identical_alleles")
    for col in ("beta", "se", "eaf", "p", "n", "position"):
        if col not in out:
            continue
        raw = out[col]
        numeric = pd.to_numeric(raw, errors="coerce").to_numpy(dtype=float, na_value=np.nan)
        mark(raw.notna() & ~np.isfinite(numeric), f"invalid_{col}")
        if col in {"beta", "se"}:
            mark(~np.isfinite(numeric), f"missing_or_nonfinite_{col}")
        if col == "se":
            mark(numeric <= 0, "nonpositive_se")
        if col in {"p", "eaf"}:
            mark((numeric < 0) | (numeric > 1), f"out_of_range_{col}")
        if col in {"n", "position"}:
            mark((numeric <= 0) | (np.isfinite(numeric) & (numeric != np.floor(numeric))), f"invalid_{col}_integer")
        out[col] = numeric
    out["validation_errors"] = [";".join(dict.fromkeys(x)) for x in issues]
    out["valid"] = [not x for x in issues]
    return out


def read_sumstats(path, *, sep=None, column_map=None, key_cols=("variant_id",), **kwargs):
    """Read CSV/TSV (including compression) and attach a full validation audit.

    `sep=None` chooses TSV for .tsv[.gz] or .txt[.gz], CSV otherwise. Explicitly
    pass sep for other sources. Invalid records remain in the returned table.
    """
    path = Path(path)
    if sep is None:
        sep = "\t" if any(s in path.suffixes for s in (".tsv", ".txt")) else ","
    out = pd.read_csv(path, sep=sep, **kwargs)
    if column_map:
        missing = set(column_map) - set(out.columns)
        if missing:
            raise ValueError(f"mapped source columns missing: {sorted(missing)}")
        out = out.rename(columns=column_map)
    return validate_sumstats(out, key_cols=key_cols)
