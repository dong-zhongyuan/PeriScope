"""Acquire an explicitly selected BIG40 IDP; verify before publishing a download.

Uses only Python's standard library. Resumable downloads keep an HTTP identity
sidecar and .part suffix until size, SHA256, gzip CRC and row validation finish.
The full source gzip is retained; canonical first/last-row previews are NOT the
full outcome table. No MR, filtering, or biological interpretation is performed.
"""
import argparse
import collections
import csv
import fcntl
import gzip
import hashlib
from html.parser import HTMLParser
from http.client import IncompleteRead
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys
import time
import urllib.request
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parents[2]
BASE = "https://open.oxcin.ox.ac.uk/ukbiobank/big40/"
INDEX_URL = BASE + "BIG40-IDPs_v4/IDPs.html"
SOURCE_HEADER = "chr rsid pos a1 a2 beta se pval(-log10)".split()
CANONICAL = ["variant_id", "effect_allele", "other_allele", "beta", "se", "p",
             "n", "build", "chromosome", "position", "outcome", "neg_log10_p"]
COHORTS = {"discovery": ("stats", "disc"), "replication": ("repro", "rep"),
           "combined": ("stats33k", "all")}


def utc():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path, obj):
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(obj, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    tmp.replace(path)


def request(url, *, method="GET", headers=None):
    return urllib.request.urlopen(urllib.request.Request(
        url, method=method, headers={"Accept-Encoding": "identity",
        "User-Agent": "pd-v3-BIG40-acquisition/1.0", **(headers or {})}), timeout=30)


def acquire(url, target, validator=None, attempts=5):
    """Validate a complete staged file before atomically naming it as final."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with request(url, method="HEAD") as response:
        identity = {"url": url, "resolved_url": response.url,
                    "http_size_bytes": int(response.headers["Content-Length"]),
                    "etag": response.headers.get("ETag"),
                    "last_modified": response.headers.get("Last-Modified"),
                    "content_type": response.headers.get("Content-Type")}
    part = target.with_name(target.name + ".part")
    resume = target.with_name(target.name + ".part.json")
    receipt_path = target.with_name(target.name + ".download.json")
    size = identity["http_size_bytes"]
    if target.exists():
        if not receipt_path.exists():
            raise RuntimeError(f"Existing final file has no receipt: {target}")
        old = json.loads(receipt_path.read_text())
        if old["http"] != identity or target.stat().st_size != size or sha256(target) != old["sha256"]:
            raise RuntimeError(f"Existing final file or upstream identity changed: {target}")
        print(f"Verified cached file: {target.name}", flush=True)
        return old
    if part.exists():
        if not resume.exists() or json.loads(resume.read_text())["http"] != identity:
            raise RuntimeError(f"Staged file lacks matching upstream identity: {part}")
        if part.stat().st_size > size:
            raise RuntimeError(f"Staged file exceeds upstream length: {part}")
        if not identity["etag"] and not identity["last_modified"]:
            raise RuntimeError("Cannot safely resume without an HTTP validator")
    else:
        atomic_json(resume, {"http": identity, "started_at_utc": utc()})
    remaining = size - (part.stat().st_size if part.exists() else 0)
    if shutil.disk_usage(target.parent).free < remaining + (100 << 20):
        raise RuntimeError("Insufficient free disk space for the selected file")
    started = time.monotonic()
    for attempt in range(attempts):
        offset = part.stat().st_size if part.exists() else 0
        if offset == size:
            break
        headers = {}
        if offset:
            headers = {"Range": f"bytes={offset}-",
                       "If-Range": identity["etag"] or identity["last_modified"]}
        try:
            with request(url, headers=headers) as response:
                if offset:
                    expected = f"bytes {offset}-{size-1}/{size}"
                    if response.status != 206 or response.headers.get("Content-Range") != expected:
                        raise RuntimeError("Server did not honor the exact resume range; no bytes appended")
                elif response.status != 200:
                    raise RuntimeError(f"Unexpected initial HTTP status {response.status}")
                for key, header in (("etag", "ETag"), ("last_modified", "Last-Modified")):
                    if identity[key] and response.headers.get(header) != identity[key]:
                        raise RuntimeError("HTTP source changed between HEAD and GET")
                if int(response.headers.get("Content-Length", size-offset)) != size-offset:
                    raise RuntimeError("Unexpected response length")
                last_update = time.monotonic()
                with part.open("ab" if offset else "wb") as handle:
                    while True:
                        chunk = response.read(1 << 20)
                        if not chunk:
                            break
                        if handle.tell() + len(chunk) > size:
                            raise RuntimeError("Received more bytes than HTTP size")
                        handle.write(chunk)
                        if time.monotonic() - last_update >= 15:
                            print(f"{target.name}: {handle.tell():,}/{size:,} bytes; "
                                  f"elapsed {time.monotonic()-started:.0f}s", flush=True)
                            last_update = time.monotonic()
                    handle.flush()
                    os.fsync(handle.fileno())
            if part.stat().st_size != size:
                raise OSError(f"Incomplete transfer: {part.stat().st_size}/{size}")
            break
        except (OSError, EOFError, IncompleteRead) as error:
            print(f"Transfer attempt {attempt+1}/{attempts}: {type(error).__name__}: {error}", flush=True)
            if attempt + 1 == attempts:
                raise
            time.sleep(min(2 ** attempt, 10))
    if not part.exists() or part.stat().st_size != size:
        raise RuntimeError("Transfer is incomplete; retained as .part")
    print(f"Validating {target.name} ({size:,} bytes)", flush=True)
    validation = validator(part) if validator else {"size_verified": True}
    receipt = {"http": identity, "download_started_at_utc": json.loads(resume.read_text())["started_at_utc"],
               "completed_at_utc": utc(), "size_bytes": size, "sha256": sha256(part),
               "validation": validation}
    # The source remains staged if either validation or hashing raises.
    atomic_json(receipt_path, receipt)
    part.replace(target)
    resume.unlink()
    print(f"Complete: {target.name}; SHA256={receipt['sha256']}", flush=True)
    return receipt


class HTMLTable(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows, self.row, self.value, self.active = [], [], "", False
        self.spans = {}
        self.colspan = self.rowspan = 1

    def fill_spans(self):
        while len(self.row) in self.spans:
            column = len(self.row)
            value, remaining = self.spans.pop(column)
            self.row.append(value)
            if remaining > 1:
                self.spans[column] = (value, remaining - 1)

    def handle_starttag(self, tag, attrs):
        if tag == "tr":
            self.row = []
        if tag in ("td", "th"):
            self.fill_spans()
            attributes = dict(attrs)
            self.colspan = int(attributes.get("colspan", 1))
            self.rowspan = int(attributes.get("rowspan", 1))
            self.active, self.value = True, ""

    def handle_data(self, data):
        if self.active:
            self.value += data

    def handle_endtag(self, tag):
        if tag in ("td", "th"):
            self.active = False
            value = self.value.strip()
            for _ in range(self.colspan):
                column = len(self.row)
                self.row.append(value)
                if self.rowspan > 1:
                    self.spans[column] = (value, self.rowspan - 1)
        if tag == "tr":
            self.fill_spans()
            self.rows.append(self.row)


def parse_index(path):
    table = HTMLTable()
    table.feed(path.read_text(encoding="utf-8"))
    header = next(row for row in table.rows if "Pheno" in row and "N(disc)" in row)
    keys = ["idp", "ukb_id", "short_name", "units", "type", "category_id", "category",
            "description", "n_discovery", "n_x_par_discovery", "n_x_nonpar_discovery",
            "n_replication", "n_x_par_replication", "n_x_nonpar_replication",
            "n_combined", "n_x_par_combined", "n_x_nonpar_combined", "heritability", "heritability_se"]
    labels = ["Pheno", "UKB ID", "IDP short name", "Units", "Type", "Cat.", "Category name",
              "IDP description", "N(disc)", "Npar(disc)", "Nnonpar(disc)",
              "N(rep)", "Npar(rep)", "Nnonpar(rep)", "N(all)", "Npar(all)",
              "Nnonpar(all)", "Heritability", "Heritability(SE)"]
    indices = [header.index(label) for label in labels]
    records = []
    for row in table.rows:
        if len(row) > indices[0] and re.fullmatch(r"\d{4}", row[indices[0]]):
            if len(row) != len(header):
                raise ValueError("IDP row width differs from header")
            records.append(dict(zip(keys, (row[i] for i in indices))))
    ids = [x["idp"] for x in records]
    if len(ids) != len(set(ids)) or not ids:
        raise ValueError("Empty or duplicated IDP index")
    return sorted(records, key=lambda x: x["idp"])


def chromosome_number(value):
    return 23 if value.upper() in {"X", "0X"} else int(value)


def canonical_row(fields, meta, cohort):
    chrom, variant, pos, oa, ea, beta, se, lp = fields
    chromosome = chromosome_number(chrom)
    return dict(zip(CANONICAL, [variant, ea, oa, float(beta), float(se),
        10 ** (-float(lp)), int(meta[f"n_{cohort}"]) if 1 <= chromosome <= 22 else None,
        "GRCh37", str(chromosome), int(pos), f"BIG40:release2:{cohort}:{meta['idp']}", float(lp)]))


def validate_sumstats(path, meta, cohort):
    counts, errors = collections.Counter(), collections.Counter()
    first = last = None
    minimum, maximum, previous = {}, {}, {}
    rows = indels = p_underflow = position_decreases = 0
    with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
        header = handle.readline().split()
        if header != SOURCE_HEADER:
            raise ValueError(f"Unexpected summary-statistic header: {header}")
        for line in handle:
            rows += 1
            fields = line.split()
            if len(fields) != len(header):
                raise ValueError(f"Wrong column count on data row {rows}")
            chrom, variant, pos, oa, ea, beta, se, lp = fields
            c, p = chromosome_number(chrom), int(pos)
            b, s, l = float(beta), float(se), float(lp)
            reasons = []
            if not (1 <= c <= 23 and p > 0): reasons.append("invalid_coordinate")
            if not math.isfinite(b): reasons.append("invalid_beta")
            if not math.isfinite(s) or s <= 0: reasons.append("invalid_se")
            if not math.isfinite(l) or l < 0: reasons.append("invalid_neg_log10_p")
            if not oa or not ea or oa == ea: reasons.append("invalid_alleles")
            errors.update(reasons)
            counts[str(c)] += 1
            minimum[str(c)] = min(minimum.get(str(c), p), p)
            maximum[str(c)] = max(maximum.get(str(c), p), p)
            position_decreases += p < previous.get(str(c), p)
            previous[str(c)] = p
            indels += len(oa) != 1 or len(ea) != 1
            p_underflow += math.isfinite(l) and l > 300 and 10 ** (-l) == 0
            if first is None: first = fields
            last = fields
            if rows % 2_000_000 == 0:
                print(f"Validated {rows:,} source rows", flush=True)
    if rows == 0 or errors:
        raise ValueError(f"Source validation failed: rows={rows}, errors={dict(errors)}")
    return {"gzip_crc_and_eof_verified": True, "source_header": header,
            "canonical_header": CANONICAL, "data_rows": rows,
            "chromosome_rows": dict(sorted(counts.items(), key=lambda x: int(x[0]))),
            "min_position_by_chromosome": minimum, "max_position_by_chromosome": maximum,
            "indel_or_multibase_rows": indels, "p_underflow_rows_if_float64": p_underflow,
            "position_decreases_within_chromosome": position_decreases,
            "row_errors": dict(errors), "duplicate_variant_ids": "not_checked_full_file",
            "first_source_row": first, "last_source_row": last,
            "first_canonical_row": canonical_row(first, meta, cohort),
            "last_canonical_row": canonical_row(last, meta, cohort),
            "x_chromosome_n": "Not assigned by position; use catalog PAR/nonPAR counts with an explicit boundary definition"}


def write_csv(path, rows, fields):
    tmp = path.with_name(path.name + ".tmp")
    with tmp.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--idp", required=True, help="Exactly one BIG40 phenotype number, e.g. 0342")
    parser.add_argument("--cohort", choices=COHORTS, default="discovery")
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data/big40")
    args = parser.parse_args()
    if not re.fullmatch(r"\d{1,4}", args.idp): parser.error("--idp must be a 1–4 digit number")
    idp = args.idp.zfill(4)
    dest = args.data_dir.resolve()
    download, extracted = dest / "download", dest / "extracted"
    download.mkdir(parents=True, exist_ok=True)
    extracted.mkdir(parents=True, exist_ok=True)
    with (dest / ".acquire.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        print(f"Free disk bytes: {shutil.disk_usage(dest).free:,}", flush=True)
        index = download / "IDPs_v4.html"
        index_receipt = acquire(INDEX_URL, index,
            lambda p: {"parsed_idp_rows": len(parse_index(p)), "utf8_verified": True})
        catalog = parse_index(index)
        meta = next((row for row in catalog if row["idp"] == idp), None)
        if meta is None: raise ValueError(f"IDP {idp} absent from complete source directory")
        catalog_path = extracted / "IDPs_v4.csv"
        write_csv(catalog_path, catalog, list(catalog[0]))
        # Preserve the official page that defines build, A2, AF and effect units.
        methods = download / "BIG40_release2_index.html"
        methods_receipt = acquire(BASE, methods,
            lambda p: {"utf8_verified": bool(p.read_text(encoding="utf-8"))})
        subdir = COHORTS[args.cohort][0]
        filename = f"{idp}.txt.gz"
        target = download / subdir / filename
        receipt = acquire(BASE + f"release2/{subdir}/{filename}", target,
                          lambda p: validate_sumstats(p, meta, args.cohort))
        v = receipt["validation"]
        prefix = f"{idp}_{args.cohort}"
        validation_path = extracted / (prefix + ".validation.json")
        atomic_json(validation_path, v)
        preview_path = extracted / (prefix + ".first_last_rows.csv")
        write_csv(preview_path, [v["first_canonical_row"], v["last_canonical_row"]], CANONICAL)
        manifest_path = dest / "MANIFEST.json"
        manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {"datasets": {}}
        manifest["written_at_utc"] = utc()
        manifest["index"] = {"path": str(index.relative_to(dest)), **index_receipt,
                             "parsed_csv": str(catalog_path.relative_to(dest)),
                             "parsed_csv_sha256": sha256(catalog_path)}
        manifest["methods_page"] = {"path": str(methods.relative_to(dest)), **methods_receipt}
        manifest["datasets"][prefix] = {
            "path": str(target.relative_to(dest)), **receipt, "catalog_record": meta,
            "build": "GRCh37/hg19", "effect_allele_source_column": "a2",
            "other_allele_source_column": "a1", "effect_scale": "SD of deconfounded phenotype per a2 allele",
            "n_autosomal": int(meta[f"n_{args.cohort}"]), "n_source": INDEX_URL,
            "build_allele_scale_source": BASE, "license": "Public summary download; explicit data-file license UNVERIFIED",
            "doi": "10.1038/s41593-021-00826-4",
            "canonical_preview": str(preview_path.relative_to(dest)),
            "canonical_preview_sha256": sha256(preview_path),
            "validation_path": str(validation_path.relative_to(dest)),
            "validation_sha256": sha256(validation_path),
            "script": str(Path(__file__).relative_to(ROOT)), "script_sha256": sha256(__file__),
            "command": sys.argv, "runtime_python": sys.version,
            "notes": ["Preview has only two rows, not the full outcome table.",
                      "Negative log10 P is retained to avoid information loss under float64 underflow.",
                      "Source indels are retained; downstream SNP eligibility is a separate audit."]}
        atomic_json(manifest_path, manifest)
        print(json.dumps({"complete": True, "idp": idp, "cohort": args.cohort,
                          "n_autosomal": int(meta[f"n_{args.cohort}"]),
                          "source_rows": v["data_rows"], "manifest": str(manifest_path)}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
