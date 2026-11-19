#!/usr/bin/env python3
"""Exact-key annotation of official peaks from the already migrated BIG40 variant table."""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/big40_peak_catalog"
VARIANTS = ROOT / "data/big40_variants/download/variants.txt.gz"
PEAKS = OUT / "all_10889_source_peaks.csv"

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        while block := f.read(8*1024*1024):
            h.update(block)
    return h.hexdigest()

def main():
    assert sha(VARIANTS) == "152517ada2fb5c572b26959d184f099369e521eeac681c01b3a3ca94b1e7299a"
    verification = json.loads((OUT / "VERIFICATION.json").read_text())
    assert sha(PEAKS) == verification["output_sha256"][PEAKS.name]
    with PEAKS.open() as f:
        peaks = list(csv.DictReader(f))
    def key(row):
        return (row["source_chr"].lstrip("0"),) + tuple(row[c] for c in ["source_variant_id", "position_grch37", "other_allele", "effect_allele"])
    wanted = {key(row) for row in peaks}
    found = {}
    # Exhaust the gzip stream to validate its CRC; no early exit after last match.
    with gzip.open(VARIANTS, "rt") as f:
        assert next(f).split() == ["chr", "rsid", "pos", "a1", "a2", "af", "info"]
        count = 0
        for count, line in enumerate(f, 1):
            fields = line.split()
            assert len(fields) == 7
            variant_key = (fields[0].lstrip("0"),) + tuple(fields[1:5])
            if variant_key in wanted:
                assert variant_key not in found
                assert 0 <= float(fields[5]) <= 1 and math.isfinite(float(fields[6])) and float(fields[6]) >= 0
                found[variant_key] = (count+1, fields[5], fields[6])
        assert count == 17103079
    rows = []
    for row in peaks:
        match = found.get(key(row))
        rows.append(dict(row, frequency_match="exact_id_coordinate_and_alleles" if match else "not_found",
                         variant_source_line=match[0] if match else "",
                         effect_allele_frequency=match[1] if match else "",
                         imputation_info=match[2] if match else ""))
    output = OUT / "all_10889_source_peaks_with_frequencies.csv"
    with output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    with output.open() as f:
        stored = list(csv.DictReader(f))
    assert len(stored) == len(peaks) == 10889
    for original, annotated in zip(peaks, stored, strict=True):
        assert all(annotated[k] == v for k,v in original.items())
        if key(original) in found:
            assert (int(annotated["variant_source_line"]), annotated["effect_allele_frequency"], annotated["imputation_info"]) == found[key(original)]
        else:
            assert annotated["frequency_match"] == "not_found"
    result = {
        "complete": True, "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "variant_rows_scanned": count, "gzip_crc_checked_by_complete_read": True,
        "peak_rows_preserved": len(peaks), "unique_peak_variant_keys": len(wanted),
        "matched_unique_variant_keys": len(found), "unmatched_unique_variant_keys": len(wanted-found.keys()),
        "matched_peak_rows": sum(row["frequency_match"] != "not_found" for row in rows),
        "key": "chromosome with leading zeros removed (0X=X), source variant id, GRCh37 position, a1, a2; no allele swaps or complements",
        "frequency_definition": "BIG40 official variant-table a2 frequency; not per-IDP/SNP sample-specific AF",
        "reverse_mr_computed": False, "ld_independence_verified": False,
        "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in [VARIANTS, PEAKS, Path(__file__)]},
        "output_sha256": {output.name: sha(output)},
    }
    (OUT / "FREQUENCY_VERIFICATION.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")}, indent=2))

if __name__ == "__main__":
    main()
