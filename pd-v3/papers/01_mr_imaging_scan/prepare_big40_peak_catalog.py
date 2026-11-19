#!/usr/bin/env python3
"""Index the complete official BIG40 local-peak table without selecting MR instruments."""
from pathlib import Path
import csv
import hashlib
import json
import math
from collections import Counter
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/big40_peak_catalog"
SOURCE = OUT / "Table1.txt"
META = ROOT / "data/big40_atlas/phenotypes_3935.csv"
EXPECTED_SOURCE_SHA = "36812cc2fe42df69c5b1e4e4e270fc3e4e19681ba7a658ffec0568794fa34f2e"
EXPECTED_HEADER = "rsid pheno chr pos a1 a2 p-value beta se p-value(repro) beta(repro) se(repro) nominal".split()
FIELD_MAP = {
    "rsid": "source_variant_id", "pheno": "source_pheno", "chr": "source_chr",
    "pos": "position_grch37", "a1": "other_allele", "a2": "effect_allele",
    "p-value": "discovery_neglog10_p", "beta": "discovery_beta",
    "se": "discovery_se", "p-value(repro)": "replication_neglog10_p",
    "beta(repro)": "replication_beta", "se(repro)": "replication_se",
    "nominal": "source_nominal_replication_flag",
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def main():
    assert sha(SOURCE) == EXPECTED_SOURCE_SHA
    receipt = json.loads((OUT / "DOWNLOAD.json").read_text())
    assert receipt["sha256"] == sha(SOURCE) and receipt["bytes"] == SOURCE.stat().st_size
    lines = SOURCE.read_text().splitlines()
    assert lines[0].split() == EXPECTED_HEADER
    raw = [dict(zip(EXPECTED_HEADER, line.split(), strict=True)) for line in lines[1:]]
    assert len(raw) == 10889
    with META.open() as f:
        metadata = list(csv.DictReader(f))
    assert len(metadata) == 3935
    by_id = {row["idp"]: row for row in metadata}
    assert len(by_id) == 3935
    assert sum(row["is_qc"] == "true" for row in metadata) == 16
    counts, auto_counts, x_counts = Counter(), Counter(), Counter()
    records, keys = [], set()
    for source_line, row in enumerate(raw, 2):
        assert row["pheno"].startswith("V")
        idp = row["pheno"][1:]
        assert len(idp) == 4 and idp in by_id
        m = by_id[idp]
        chrom = row["chr"].lstrip("0")
        if chrom == "23":
            chrom = "X"
        assert chrom in [str(i) for i in range(1, 23)] + ["X"]
        assert int(row["pos"]) > 0
        for key in ["p-value", "beta", "se", "p-value(repro)", "beta(repro)", "se(repro)"]:
            assert math.isfinite(float(row[key]))
        assert float(row["p-value"]) >= 7.5
        assert float(row["p-value(repro)"]) >= 0
        assert float(row["se"]) > 0 and float(row["se(repro)"]) > 0
        assert row["nominal"] in {"0", "1"}
        key = (idp, chrom, row["pos"], row["a1"], row["a2"])
        assert key not in keys
        keys.add(key)
        counts[idp] += 1
        (x_counts if chrom == "X" else auto_counts)[idp] += 1
        record = {"source_line": source_line, "idp": idp, "phenotype": m["phenotype"],
                  "chromosome": chrom, "family": m["family"], "is_qc": m["is_qc"]}
        record.update({new: row[old] for old, new in FIELD_MAP.items()})
        record.update({"ld_independence_verified": "false", "selected_for_reverse_mr": "false"})
        records.append(record)
    coverage = []
    for m in metadata:
        idp = m["idp"]
        row = {key: m[key] for key in ["idp", "phenotype", "ukb_id", "short_name", "family", "is_qc",
               "n_discovery", "n_x_par_discovery", "n_x_nonpar_discovery", "n_replication",
               "n_x_par_replication", "n_x_nonpar_replication"]}
        row.update({"source_peak_rows": counts[idp], "autosomal_peak_rows": auto_counts[idp],
                    "x_peak_rows": x_counts[idp],
                    "source_peak_status": "reported_local_peak" if counts[idp] else "no_reported_peak_at_source_filters",
                    "ld_independence_verified": "false", "reverse_mr_computed": "false"})
        coverage.append(row)
    write_csv(OUT / "all_10889_source_peaks.csv", records)
    write_csv(OUT / "all_3935_phenotype_coverage.csv", coverage)
    # Independent line/token comparison of every exported source field, including X and replication.
    with (OUT / "all_10889_source_peaks.csv").open() as f:
        exported = list(csv.DictReader(f))
    cells = 0
    for line_number, (line, row) in enumerate(zip(lines[1:], exported, strict=True), 2):
        assert int(row["source_line"]) == line_number
        for source_index, source_key in enumerate(EXPECTED_HEADER):
            assert line.split()[source_index] == row[FIELD_MAP[source_key]]
            cells += 1
    with (OUT / "all_3935_phenotype_coverage.csv").open() as f:
        exported_coverage = list(csv.DictReader(f))
    assert {r["idp"] for r in exported_coverage} == set(by_id)
    assert sum(int(r["source_peak_rows"]) for r in exported_coverage) == len(exported)
    check_counts = Counter(r["idp"] for r in exported)
    assert all(int(r["source_peak_rows"]) == check_counts[r["idp"]] for r in exported_coverage)
    assert all(int(r["autosomal_peak_rows"]) + int(r["x_peak_rows"]) == int(r["source_peak_rows"]) for r in exported_coverage)
    result = {
        "complete": True, "kind": "public_source_catalog_only", "verified_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_rows": len(records), "source_cells_exactly_preserved": cells,
        "phenotypes_in_original_denominator": len(metadata),
        "phenotypes_with_source_peaks": len(counts), "phenotypes_without_source_peaks": 3935-len(counts),
        "qc_phenotypes_retained": 16,
        "qc_source_peak_rows": sum(counts[m["idp"]] for m in metadata if m["is_qc"] == "true"),
        "autosomal_source_peak_rows": sum(auto_counts.values()), "x_source_peak_rows": sum(x_counts.values()),
        "source_nominal_replication_rows": sum(int(r["nominal"]) for r in raw),
        "unique_source_variant_ids": len({r["rsid"] for r in raw}),
        "reverse_mr_computed": False, "ld_independence_verified": False,
        "QSM2022_included": False,
        "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in [SOURCE, META, OUT / "DOWNLOAD.json", Path(__file__)]},
        "output_sha256": {p.name: sha(p) for p in [OUT / "all_10889_source_peaks.csv", OUT / "all_3935_phenotype_coverage.csv"]},
    }
    (OUT / "VERIFICATION.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")}, indent=2))

if __name__ == "__main__":
    main()
