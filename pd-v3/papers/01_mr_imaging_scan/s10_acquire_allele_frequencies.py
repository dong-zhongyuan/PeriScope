#!/usr/bin/env python3
"""Download BIG40 variant frequencies and retain all rows at selected cis positions.

One network connection at a time. The complete source gzip is retained and read
through EOF (including gzip CRC validation); no allele or AF filtering is used.
"""
import collections
import csv
import fcntl
import gzip
import json
import math
import os
from pathlib import Path
import platform
import time

from s07_acquire_big40 import acquire, atomic_json, chromosome_number, request, sha256, utc


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/big40_variants"
EXPOSURES = ROOT / "data/ukb_ppp/extracted/cis_primary_selected.csv"
SOURCE_URL = "https://open.oxcin.ox.ac.uk/ukbiobank/big40/release2/variants.txt.gz"
METHODS_URL = "https://open.oxcin.ox.ac.uk/ukbiobank/big40/"
SOURCE = DATA / "download/variants.txt.gz"
HEADER = "chr rsid pos a1 a2 af info".split()
FIELDS = ["chr", "pos37", "oa_y", "ea_y", "eaf_y", "info", "variant_y", "source_row",
          "chr_raw", "af_raw", "af_status", "info_status", "build", "source_y"]
MISSING = {"", ".", "NA", "NAN", "NULL", "NONE", "-", "N/A"}


def value_status(raw, frequency=False):
    if raw.upper() in MISSING:
        return "missing"
    try:
        value = float(raw)
    except ValueError:
        return "not_numeric"
    if not math.isfinite(value):
        return "missing_nonfinite"
    if frequency and not 0 <= value <= 1:
        return "outside_0_1"
    return "available"


def scan_full_source(path, positions):
    counts, af_statuses, info_statuses = collections.Counter(), collections.Counter(), collections.Counter()
    matched = []
    rows = 0
    first = last = None
    last_update = time.monotonic()
    with gzip.open(path, "rt", encoding="ascii", errors="strict") as handle:
        header = handle.readline().split()
        if header != HEADER:
            raise ValueError(f"Unexpected source header: {header}")
        for source_row, line in enumerate(handle, 2):
            fields = line.split()
            if len(fields) != len(HEADER):
                raise ValueError(f"Unexpected source column count on physical line {source_row}")
            chrom_raw, rsid, pos_raw, a1, a2, af, info = fields
            chrom, pos = chromosome_number(chrom_raw), int(pos_raw)
            if not 1 <= chrom <= 23 or pos <= 0:
                raise ValueError(f"Invalid coordinate on physical line {source_row}")
            af_status, info_status = value_status(af, frequency=True), value_status(info)
            af_statuses[af_status] += 1
            info_statuses[info_status] += 1
            counts[str(chrom)] += 1
            rows += 1
            first = fields if first is None else first
            last = fields
            if (chrom, pos) in positions:
                matched.append(dict(zip(FIELDS, [
                    chrom, pos, a1, a2, af, info, rsid, source_row,
                    chrom_raw, af, af_status, info_status, "GRCh37", SOURCE_URL,
                ])))
            if time.monotonic() - last_update >= 15:
                print(f"Read {rows:,} source rows; retained {len(matched):,} coordinate matches", flush=True)
                last_update = time.monotonic()
    if not rows:
        raise ValueError("Source table is empty")
    # Reading to EOF above checks every gzip member's trailer CRC and length.
    validation = {
        "gzip_crc_verified": True, "complete_stream_read": True,
        "source_header": header, "source_data_rows": rows,
        "source_row_definition": "physical line number including header as line 1",
        "chromosome_counts": dict(counts), "af_status_counts": dict(af_statuses),
        "info_status_counts": dict(info_statuses), "first_row": first, "last_row": last,
        "retained_coordinate_rows": len(matched), "positions_requested": len(positions),
        "filters": {"coordinate_membership_only": True, "af": None, "info": None,
                    "alleles": None, "rsid": None},
    }
    return validation, matched


def atomic_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.replace(path)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    SOURCE.parent.mkdir(parents=True, exist_ok=True)
    with (DATA / ".acquire.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        input_digest = sha256(EXPOSURES)
        with EXPOSURES.open(newline="") as handle:
            exposures = list(csv.DictReader(handle))
        if len({r["exposure_id"] for r in exposures}) != len(exposures):
            raise ValueError("Selected exposure IDs are not unique")
        positions = {(chromosome_number(r["chrom"]), int(r["pos37"])) for r in exposures}
        if not positions:
            raise ValueError("Selected instrument table is empty")

        # Preserve the small official page describing AF and allele orientation.
        methods_path = DATA / "download/BIG40_methods.html"
        methods_receipt = DATA / "download/BIG40_methods.json"
        if not methods_path.exists():
            with request(METHODS_URL) as response:
                content = response.read()
                methods_http = {"url": METHODS_URL, "resolved_url": response.url,
                                "retrieved_at_utc": utc(), "http_status": response.status}
            methods_path.write_bytes(content)
            atomic_json(methods_receipt, {**methods_http, "sha256": sha256(methods_path)})
        else:
            recorded = json.loads(methods_receipt.read_text())
            if recorded["sha256"] != sha256(methods_path):
                raise ValueError("Cached official methods page SHA256 mismatch")

        collected = {}

        def validator(path):
            validation, matches = scan_full_source(path, positions)
            collected.update(validation=validation, matches=matches)
            return validation

        # acquire() downloads sequentially, validates the complete staged gzip,
        # hashes it, and only then publishes it without the .part suffix.
        receipt = acquire(SOURCE_URL, SOURCE, validator=validator)
        if not collected:
            validator(SOURCE)  # A cached gzip still needs extraction for this plan.
        if input_digest != sha256(EXPOSURES):
            raise RuntimeError("Exposure plan changed during extraction; rerun against a stable input")
        matches = collected["matches"]
        by_position = collections.defaultdict(list)
        for row in matches:
            by_position[row["chr"], row["pos37"]].append(row)
        coverage = []
        for row in exposures:
            chrom, pos = chromosome_number(row["chrom"]), int(row["pos37"])
            hits = by_position.get((chrom, pos), [])
            coverage.append({
                "exposure_id": row["exposure_id"], "gene": row["gene"], "assay": row["assay"],
                "chr": chrom, "pos37": pos, "variant_x": row["variant"],
                "matched_source_rows": len(hits),
                "available_frequency_rows": sum(h["af_status"] == "available" for h in hits),
                "coordinate_status": "found" if hits else "coordinate_not_in_source",
                "source_rows": ";".join(str(h["source_row"]) for h in hits),
            })
        extracted = DATA / "extracted/cis_instrument_frequencies.csv"
        coverage_path = DATA / "extracted/instrument_frequency_coverage.csv"
        atomic_csv(extracted, matches, FIELDS)
        atomic_csv(coverage_path, coverage, list(coverage[0]))
        summary = {
            **collected["validation"], "selected_assays": len(exposures),
            "selected_unique_positions": len(positions),
            "positions_found": len(by_position), "positions_missing": len(positions) - len(by_position),
            "assays_with_coordinate_match": sum(r["matched_source_rows"] > 0 for r in coverage),
            "assays_without_coordinate_match": sum(r["matched_source_rows"] == 0 for r in coverage),
            "positions_with_multiple_source_rows": sum(len(v) > 1 for v in by_position.values()),
            "retained_af_statuses": dict(collections.Counter(r["af_status"] for r in matches)),
            "retained_chrX_rows": sum(r["chr"] == 23 for r in matches),
            "allele_rule": "oa_y=a1; ea_y=a2; eaf_y=af, the a2 frequency, as documented by BIG40",
            "missing_frequency_rule": "retain original token; no imputation or filtering",
            "mapping_rule": "coordinate match only; all same-position source rows retained, no allele matching yet",
            "frequency_population": "BIG40 release2 variants table; exact cohort denominator is not given on the official page",
            "intended_use": "parallel allele-frequency sensitivity; no change to default palindromic-SNP exclusion",
        }
        summary_path = DATA / "extracted/extraction_summary.json"
        atomic_json(summary_path, summary)
        outputs = [extracted, coverage_path, summary_path]
        atomic_json(DATA / "MANIFEST.json", {
            "generated_at_utc": utc(), "python": platform.python_version(),
            "source": receipt, "source_path": str(SOURCE.relative_to(ROOT)),
            "official_methods": {"url": METHODS_URL, "file": str(methods_path.relative_to(ROOT)),
                                 "sha256": sha256(methods_path)},
            "input": {"path": str(EXPOSURES.relative_to(ROOT)), "sha256": input_digest},
            "code": {str(p.relative_to(ROOT)): sha256(p) for p in
                     [Path(__file__).resolve(), Path(__file__).with_name("s07_acquire_big40.py").resolve()]},
            "outputs": {str(p.relative_to(ROOT)): {"sha256": sha256(p), "bytes": p.stat().st_size}
                        for p in outputs},
            "network_connections_at_once": 1,
        })
        print(json.dumps({k: summary[k] for k in ["source_data_rows", "retained_coordinate_rows",
              "selected_assays", "selected_unique_positions", "positions_found", "positions_missing",
              "positions_with_multiple_source_rows", "retained_af_statuses", "retained_chrX_rows"]}, indent=2))


if __name__ == "__main__":
    main()
