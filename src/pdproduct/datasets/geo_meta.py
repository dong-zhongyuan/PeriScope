"""GEO series/sample metadata parsing (one-source-of-truth for donor facts).

Parses NCBI GEO SOFT text (acc.cgi?targ=gsm&form=text) into structured
sample records. Dataset preprocessing must use these records — not README
claims — for donor/condition attribution (HANDOFF_SERVER_V1 §10.4).
"""
from __future__ import annotations

import json
import re
import urllib.request
from pathlib import Path

GEO_GSM_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?targ=gsm&form=text&view=brief&acc={accession}"
)


def fetch_gsm_text(accession: str, dest_dir: str | Path) -> Path:
    dest = Path(dest_dir) / f"{accession}_gsm.txt"
    if dest.exists():
        return dest
    url = GEO_GSM_URL.format(accession=accession)
    req = urllib.request.Request(url, headers={"User-Agent": "pdproduct/0.1"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
        f.write(r.read())
    return dest


def parse_gsm_text(path: str | Path) -> list[dict]:
    """Parse GEO SOFT text into per-sample records."""
    samples: list[dict] = []
    current: dict | None = None
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    samples.append(current)
                current = {"accession": line.split("=", 1)[1].strip(), "characteristics": {}}
            elif current is not None and line.startswith("!Sample_"):
                key, _, value = line.partition(" = ")
                key = key.removeprefix("!Sample_")
                if key == "characteristics_ch1":
                    ck, _, cv = value.partition(":")
                    current["characteristics"][ck.strip()] = cv.strip()
                elif key in current and key != "characteristics":
                    current[key] += " | " + value  # repeated single-valued fields
                else:
                    current[key] = value
    if current is not None:
        samples.append(current)
    return samples


def human_samples(samples: list[dict]) -> list[dict]:
    out = []
    for s in samples:
        source = str(s.get("source_name_ch1", ""))
        title = str(s.get("title", ""))
        organism = str(s.get("organism_ch1", ""))
        if "Human" in source or "Homo" in source or "human" in title or organism == "Homo sapiens":
            out.append(s)
    return out


def _condition_of(ch: dict) -> str:
    for key, value in ch.items():
        kl = key.lower()
        if "disease" in kl or "status" in kl or "diagnosis" in kl or "stage" in kl:
            return f"{value}"
    return "?"


def _donor_of(ch: dict) -> str:
    for key, value in ch.items():
        if "donor" in key.lower():
            return value
    return "?"


def summarize(accession: str, cache_dir: str | Path) -> dict:
    """Fetch (or reuse cached) metadata, parse, save a structured JSON."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    raw = fetch_gsm_text(accession, cache_dir)
    samples = parse_gsm_text(raw)
    human = human_samples(samples)
    result = {
        "accession": accession,
        "n_samples_total": len(samples),
        "n_samples_human": len(human),
        "human_samples": [
            {
                "accession": s["accession"],
                "title": s.get("title", ""),
                "source": s.get("source_name_ch1", ""),
                "characteristics": s["characteristics"],
            }
            for s in human
        ],
    }
    out = cache_dir / f"{accession}_human_meta.json"
    with open(out, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2, ensure_ascii=False)
    result["meta_json"] = str(out)
    return result


# --- GSE178265-specific donor extraction (patterns verified against data) ---

_PAT_CN = re.compile(r"pPDCN(\d+)")            # caudate: pPDCN4340DAPIA030419
_PAT_SN = re.compile(r"pPDsHSrSNxi(\d+)d")     # SN: pPDsHSrSNxi3345d200429DAPIA


def gse178265_donor_of(title_or_barcode_prefix: str) -> str | None:
    """Map a sample title / barcode prefix to a donor id.

    Returns e.g. 'CN-4340' (caudate) or 'SN-3345' (substantia nigra);
    None for non-human or unmatched patterns.
    """
    m = _PAT_SN.search(title_or_barcode_prefix)
    if m:
        return f"SN-{m.group(1)}"
    m = _PAT_CN.search(title_or_barcode_prefix)
    if m:
        return f"CN-{m.group(1)}"
    return None
