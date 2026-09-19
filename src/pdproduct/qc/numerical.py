"""Numerical QC gates for step outputs (G4/G6)."""
from __future__ import annotations

import numpy as np


def check_numeric_array(arr, name: str = "array", *, finite: bool = True) -> dict:
    a = np.asarray(arr, dtype=np.float64)
    report = {
        "name": name,
        "shape": list(a.shape),
        "n_nan": int(np.isnan(a).sum()),
        "n_inf": int(np.isinf(a).sum()),
        "min": float(np.nanmin(a)) if a.size else None,
        "max": float(np.nanmax(a)) if a.size else None,
    }
    ok = report["n_nan"] == 0 and (not finite or report["n_inf"] == 0)
    report["ok"] = bool(ok)
    return report


def check_correlation_value(rho: float, name: str = "rho") -> dict:
    ok = np.isfinite(rho) and -1.0 <= rho <= 1.0
    return {"name": name, "value": rho, "ok": bool(ok)}


def qc_outputs(reports: list[dict]) -> dict:
    return {"ok": all(r["ok"] for r in reports), "reports": reports}
