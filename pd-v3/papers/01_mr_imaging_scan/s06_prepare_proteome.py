#!/usr/bin/env python3
"""Extract every Sun 2023 ST9 cis primary pQTL and the complete ST3 catalogue.

Uses the complete publisher XLSX, with no outcome inputs and no allele filtering.
Run with the project hh/env Python. --download fetches the source if absent.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import posixpath
import re
import sys
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/ukb_ppp"
URL = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-023-06592-6/MediaObjects/41586_2023_6592_MOESM3_ESM.xlsx"
XLSX = DATA / "download/41586_2023_6592_MOESM3_ESM.xlsx"
META = DATA / "download/sun2023_download_metadata.json"
VERIFIED = ROOT / "reference/sun2023_three_leads_verified_20261115.json"
OLD = ROOT / "reference/periscope_mr_20261115/results/combined_exposure_instruments.csv"
EXPECTED_SHA256 = "e174964b81e0d16f84e485593d4b349e9c24790f60983724d4f9f37681938cc7"
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PID = "UKBPPP ProteinID"
VID = "Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)"
PAV_TERMS = {
    "missense_variant", "stop_gained", "stop_lost", "start_lost",
    "frameshift_variant", "inframe_insertion", "inframe_deletion",
    "protein_altering_variant",
}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def download():
    """Commit a download only after length, ZIP CRC and fixed source hash checks."""
    XLSX.parent.mkdir(parents=True, exist_ok=True)
    partial = XLSX.with_suffix(".xlsx.partial")
    started = now()
    req = urllib.request.Request(URL, headers={"User-Agent": "pd-v3-source-verification/1.0"})
    with urllib.request.urlopen(req, timeout=90) as response:
        if response.status != 200:
            raise RuntimeError(f"Full download expected HTTP 200, got {response.status}")
        headers = dict(response.headers.items())
        with partial.open("wb") as handle:
            while chunk := response.read(1024 * 1024):
                handle.write(chunk)
    if partial.stat().st_size != 34411394:
        raise ValueError("Unexpected source byte count; do not treat a partial or new version as verified")
    with zipfile.ZipFile(partial) as archive:
        if archive.testzip() is not None:
            raise ValueError("Source ZIP CRC check failed")
    digest = sha256(partial)
    if digest != EXPECTED_SHA256:
        raise ValueError("Publisher XLSX changed; review new source before replacing the pinned input")
    os.replace(partial, XLSX)
    write_json(META, {
        "url": URL, "download_started_at_utc": started,
        "download_finished_at_utc": now(), "filename": XLSX.name,
        "size_bytes": XLSX.stat().st_size, "sha256": digest, "http_status": 200,
        "http_headers": headers, "full_file": True, "zip_crc_verified": True,
    })


def sheet_members(archive):
    rels = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    targets = {r.get("Id"): r.get("Target") for r in rels}
    workbook = ET.fromstring(archive.read("xl/workbook.xml"))
    return {
        s.get("name"): posixpath.normpath("xl/" + targets[s.get("{" + REL_NS + "}id")])
        for s in workbook.findall("m:sheets/m:sheet", NS)
    }


def read_table(xml, shared, header_row):
    """Read full sheetData, not the workbook's stale autofilter/named ranges."""
    rows = []
    headers = None
    for row in ET.fromstring(xml).findall("m:sheetData/m:row", NS):
        number = int(row.get("r"))
        cells = {}
        for cell in row:
            col = re.sub(r"[0-9]+$", "", cell.get("r"))
            value = cell.find("m:v", NS)
            if cell.get("t") == "s":
                cells[col] = shared[int(value.text)] if value is not None else ""
            elif cell.get("t") == "inlineStr":
                cells[col] = "".join(t.text or "" for t in cell.findall(".//m:t", NS))
            else:
                cells[col] = value.text if value is not None else ""
        if number == header_row:
            headers = cells
            if len(set(headers.values())) != len(headers):
                raise ValueError("Duplicate source headers")
        elif number > header_row and any(cells.values()):
            if headers is None or not cells.get("A"):
                raise ValueError(f"Unexpected source table row {number}")
            extra = set(cells) - set(headers)
            if any(cells[c] for c in extra):
                raise ValueError(f"Unmapped source columns in row {number}: {extra}")
            rows.append({"original_row": number, **{h: cells.get(c, "") for c, h in headers.items()}})
    if not rows:
        raise ValueError("Empty source table")
    return rows


def write_csv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def verify_with_openpyxl(raw_tables):
    """Independent reader check of every extracted source value, without edits."""
    import openpyxl
    workbook = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    checked = {}
    try:
        for sheet, rows in raw_tables.items():
            headers = list(rows[0])[1:]
            extracted = {r["original_row"]: r for r in rows}
            found, cell_count = set(), 0
            header_row = 3 if sheet == "ST3" else 5
            for number, values in enumerate(workbook[sheet].iter_rows(values_only=True), 1):
                if number == header_row and list(values[:len(headers)]) != headers:
                    raise ValueError(f"Independent reader header mismatch: {sheet}")
                if number <= header_row or not any(v is not None for v in values):
                    continue
                if number not in extracted:
                    raise ValueError(f"Independent reader found an unextracted row: {sheet} {number}")
                found.add(number)
                for header, value in zip(headers, values):
                    raw = extracted[number][header]
                    if isinstance(value, (int, float)):
                        equal = math.isclose(float(raw), value, rel_tol=1e-15, abs_tol=0)
                    else:
                        equal = raw == (value if value is not None else "")
                    if not equal:
                        raise ValueError(f"Independent reader cell mismatch: {sheet} {number} {header}")
                    cell_count += 1
            if found != set(extracted):
                raise ValueError(f"Independent reader row count mismatch: {sheet}")
            checked[sheet] = {"rows": len(found), "cells": cell_count}
    finally:
        workbook.close()
    return {"reader": "openpyxl", "version": openpyxl.__version__, "checks": checked,
            "numeric_comparison": "relative tolerance 1e-15 against independently decoded numeric cells; strings exact",
            "source_sha256": sha256(XLSX), "passed": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="Download the pinned XLSX if it is absent")
    parser.add_argument("--verify-openpyxl", action="store_true", help="Cross-check every source cell with the independent read-only reader")
    args = parser.parse_args()
    if not XLSX.exists() and args.download:
        download()
    if not XLSX.exists() or not META.exists():
        raise FileNotFoundError("Full XLSX and download metadata required; run --download if absent")
    metadata = json.loads(META.read_text())
    digest = sha256(XLSX)
    if digest != EXPECTED_SHA256 or digest != metadata["sha256"]:
        raise ValueError("Source SHA256 mismatch")
    if XLSX.stat().st_size != 34411394 or not metadata["full_file"]:
        raise ValueError("Full source identity not verified")
    excerpt = json.loads(VERIFIED.read_text())
    with zipfile.ZipFile(XLSX) as archive:
        if archive.testzip() is not None:
            raise ValueError("Source ZIP CRC check failed")
        members = sheet_members(archive)
        shared = ["".join(t.text or "" for t in si.findall(".//m:t", NS))
                  for si in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        xmls = {sheet: archive.read(members[sheet]) for sheet in ("ST3", "ST9")}
        member_hashes = {members[s]: hashlib.sha256(b).hexdigest() for s, b in xmls.items()}
        for member, expected in excerpt["source"]["zip_member_hashes"].items():
            if member_hashes[member] != expected:
                raise ValueError("Previously verified ST9 member does not match full workbook")
        st3 = read_table(xmls["ST3"], shared, 3)
        st9 = read_table(xmls["ST9"], shared, 5)
    if len(st3) != 2941 or len(st9) != 14287:
        raise ValueError("Complete table row counts differ from verified source version")
    lookup = {r[PID]: r for r in st3}
    if len(lookup) != len(st3) or any(r[PID] not in lookup for r in st9):
        raise ValueError("ST3 assay identity duplicate or missing ST9 mapping")
    by_row = {r["original_row"]: r for r in st9}
    three_checks = []
    for rec in excerpt["records"]:
        source = by_row[rec["excel_row"]]
        if any(source[k] != v for k, v in rec["raw_fields"].items()):
            raise ValueError(f"Previously verified raw row changed: {rec['protein']}")
        three_checks.append({"protein": rec["protein"], "original_row": rec["excel_row"],
                             "all_previously_recorded_raw_fields_exactly_equal": True})

    mapping = [{
        "exposure_id": r[PID], "gene": r["Gene symbol"], "assay": r["Olink ID"],
        "assay_target": r["Assay Target"], "uniprot": r["UniProt"],
        "protein_panel": r["Protein panel"], "gene_chrom38": r["Gene CHROM"],
        "gene_start38": r["Gene start"], "gene_end38": r["Gene end"],
        "dilution_factor": r["Dilution factor"],
        "percent_samples_below_lod_raw": r["% of samples below LOD"],
        "coefficient_variation_median_percent_raw": r["Coefficient of variation (median %)"],
        "percent_data_failing_qc_raw": r["% of data failing QC"],
        "icc_sample_age_raw": r["ICC (interclass correlation) of sample age "],
        "st3_original_row": r["original_row"], "source_sheet": "ST3", "source_url": URL,
    } for r in st3]

    instruments = []
    for r in st9:
        if r["cis/trans"] != "cis":
            continue
        chrom, pos37, a0, a1, *suffix = r[VID].split(":")
        chrom_numeric = {"X": "23", "Y": "24", "MT": "25"}.get(chrom, chrom)
        if chrom_numeric != r["CHROM"] or suffix != ["imp", "v1"]:
            raise ValueError(f"Unexpected variant ID row {r['original_row']}")
        beta, se = Decimal(r["BETA (discovery, wrt. A1)"]), Decimal(r["SE (discovery)"])
        eaf = Decimal(r["A1FREQ (discovery)"])
        if not beta.is_finite() or not se.is_finite() or se <= 0 or not 0 <= eaf <= 1:
            raise ValueError(f"Unexpected discovery statistics row {r['original_row']}")
        assay = lookup[r[PID]]
        rsid = r["rsID"] if re.fullmatch(r"rs\d+", r["rsID"]) else "NA"
        consequence = r["Annotated gene consequence"]
        is_pav = bool(set(re.split(r"[,;&|]", consequence)) & PAV_TERMS)
        snv = len(a0) == len(a1) == 1 and set(a0 + a1) <= set("ACGT")
        instruments.append({
            "exposure_id": r[PID], "gene": assay["Gene symbol"], "assay": assay["Olink ID"],
            "variant": rsid if rsid != "NA" else r[VID], "chrom": int(chrom_numeric),
            "pos37": int(pos37), "pos38": int(r["GENPOS (hg38)"]),
            "ea_x": a1, "oa_x": a0, "beta_x": r["BETA (discovery, wrt. A1)"],
            "se_x": r["SE (discovery)"], "eaf_x": r["A1FREQ (discovery)"], "n_x": "NA",
            "source_x": "Sun2023_ST9_discovery", "units_x": "reported protein analysis unit",
            "original_row": r["original_row"], "source_sheet": "ST9", "source_url": URL,
            "variant_id37": r[VID], "variant_chrom_raw": chrom, "rsid": rsid, "rsid_raw": r["rsID"],
            "assay_target": assay["Assay Target"], "uniprot": assay["UniProt"],
            "protein_panel": assay["Protein panel"], "st3_original_row": assay["original_row"],
            "neglog10_p_x": r["log10(p) (discovery)"], "F_x": str(float((beta / se) ** 2)),
            "abs_z_x": str(float(abs(beta / se))),
            "beta_x_replication": r["BETA (replication)"], "se_x_replication": r["SE (replication)"],
            "eaf_x_replication": r["A1FREQ (replication)"], "neglog10_p_x_replication": r["log10(p) (replication)"],
            "cis_trans_published": r["cis/trans"], "cis_gene_published": r["cis gene"],
            "region_id": r["Region ID"], "region_start": r["Region Start"], "region_end": r["Region End"],
            "MHC": r["MHC"], "pav_consequence": consequence,
            "pav_coding_consequence": is_pav, "pav_annotation_scope": "sentinel_consequence_only_not_LD_linkage",
            "annotated_gene": r["Bioinfomatic annotated gene"], "ensembl_gene_id": r["Ensembl gene ID"],
            "biotype": r["Biotype"], "distance_to_gene": r["Distance to gene"],
            "CADD_phred": r["CADD_phred"], "SIFT": r["SIFT"], "PolyPhen": r["PolyPhen"],
            "PHAST_Phylop_score": r["PHAST Phylop_score"], "FitCons_score": r["FitCons_score"], "IMPACT": r["IMPACT"],
            "is_palindromic_snv": snv and {a0, a1} in ({"A", "T"}, {"C", "G"}),
            "variant_class": "snv" if snv else "indel" if len(a0) != len(a1) else "other",
            "n_x_status": "not_reported_per_assay_or_variant_in_ST9",
            "units_x_status": "not_relabelled_as_SD_from_ST9_alone",
        })
    grouped = defaultdict(list)
    for row in instruments:
        grouped[row["exposure_id"]].append(row)
    for rows in grouped.values():
        ranked = sorted(rows, key=lambda r: (-abs(Decimal(r["beta_x"]) / Decimal(r["se_x"])), r["original_row"]))
        for rank, row in enumerate(ranked, 1):
            row.update(selection_rank=rank, selected_for_single_instrument=rank == 1,
                       selection_rule="descending_abs_discovery_beta_over_se_then_ascending_original_row")
    selected = [r for r in instruments if r["selected_for_single_instrument"]]
    if len(instruments) != 1955 or len(selected) != 1954:
        raise ValueError("Cis primary coverage mismatch")
    selected_lookup = {r["exposure_id"]: r for r in selected}
    coverage = []
    for row in mapping:
        tool = selected_lookup.get(row["exposure_id"])
        coverage.append({**row,
            "n_cis_primary": len(grouped[row["exposure_id"]]), "has_cis_primary": tool is not None,
            "instrument_status": "published_cis_primary_available" if tool else "no_published_cis_primary",
            "exclusion_reason": "" if tool else "no_published_cis_primary",
            "selected_variant": tool["variant"] if tool else "NA",
            "selected_original_row": tool["original_row"] if tool else "NA",
        })

    with OLD.open(newline="") as handle:
        old_rows = list(csv.DictReader(handle))
    for check, rec in zip(three_checks, excerpt["records"]):
        new = selected_lookup[rec["raw_fields"][PID]]
        old = [r for r in old_rows if r["assay"] == new["exposure_id"] and r["rsid"] == new["rsid"]]
        if len(old) != 1:
            raise ValueError(f"Legacy comparison row not unique for {new['exposure_id']}")
        old = old[0]
        for target, legacy in [("beta_x", "beta_exposure"), ("se_x", "se_exposure"), ("eaf_x", "eaf_exposure")]:
            if not math.isclose(float(new[target]), float(old[legacy]), rel_tol=1e-12, abs_tol=1e-14):
                raise ValueError(f"Legacy numeric mismatch {new['gene']} {target}")
        for target, legacy in [("ea_x", "effect_allele"), ("oa_x", "other_allele"), ("pos37", "pos37"), ("pos38", "pos")]:
            if str(new[target]) != old[legacy]:
                raise ValueError(f"Legacy allele/coordinate mismatch {new['gene']} {target}")
        check.update(legacy_beta_se_eaf_and_alleles_positions_match=True,
                     selected_variant=new["variant"], effect_allele=new["ea_x"], other_allele=new["oa_x"])

    out = DATA / "extracted"
    tables = {"ST3_raw.csv": st3, "ST9_raw.csv": st9, "assay_mapping.csv": mapping,
              "assay_coverage.csv": coverage, "cis_primary_all.csv": instruments, "cis_primary_selected.csv": selected}
    for name, rows in tables.items():
        write_csv(out / name, rows)
    independent = verify_with_openpyxl({"ST3": st3, "ST9": st9}) if args.verify_openpyxl else None
    if independent:
        write_json(out / "independent_reader_verification.json", independent)
    multi = {k: [{"variant": r["variant"], "original_row": r["original_row"],
                  "abs_z_x": r["abs_z_x"], "selected": r["selected_for_single_instrument"]} for r in v]
             for k, v in grouped.items() if len(v) > 1}
    summary = {
        "st3_assays": len(st3), "st9_all_primary_rows": len(st9),
        "st9_cis_primary_rows": len(instruments), "assays_with_cis_primary": len(selected),
        "assays_without_cis_primary": len(st3) - len(selected),
        "multi_cis_primary_assays": multi,
        "cis_primary_variant_classes": dict(Counter(r["variant_class"] for r in instruments)),
        "cis_primary_palindromic_snv_rows_preserved": sum(r["is_palindromic_snv"] for r in instruments),
        "cis_primary_without_rsid_rows_preserved": sum(r["rsid"] == "NA" for r in instruments),
        "cis_primary_pav_coding_consequence_rows": sum(r["pav_coding_consequence"] for r in instruments),
        "three_lead_checks": three_checks, "selection_uses_outcomes": False,
        "additional_p_value_or_F_filter_applied": False, "allele_filters_applied": False,
        "per_row_sample_sizes_imputed": False,
        "source_member_sha256": member_hashes,
    }
    write_json(out / "extraction_summary.json", summary)
    schema = {"canonical_columns": list(instruments[0]), "coverage_columns": list(coverage[0]),
              "identity": "exposure_id is full original UKBPPP ProteinID; assay is Olink ID; composite and repeated gene targets are kept",
              "coordinate_rule": "pos37 from ST9 variant ID; pos38 from GENPOS (hg38); gene bounds from ST3 GRCh38; chrom uses ST9 numeric CHROM (23=X), with variant_chrom_raw preserved",
              "allele_rule": "ea_x=A1, oa_x=A0 from ST9 variant ID; BETA (discovery, wrt. A1), SE and A1FREQ retained",
              "raw_values": "Raw CSV stores exact XML scalar text (shared strings resolved); raw blanks stay blank",
              "unknown_values": "Canonical n_x is NA for every row; ST9 does not supply per-row n",
              "p_values": "positive published log10(p) retained as neglog10_p_x; no exponentiation underflow",
              "pav_rule": "pav_coding_consequence flags a protein-altering sentinel consequence only; no LD proxy annotation inferred",
              "pav_terms": sorted(PAV_TERMS),
              "selection": "For each full assay ID, maximum absolute discovery beta/SE; exact ties by original Excel row",
              "units": "Reported protein analysis unit retained; no claim that ST9 alone defines an SD scale",
              "full_plan": "All ST3 assays remain in coverage including absence of a published cis primary; this absence is not a null association",
              "scope": "Published primary pQTLs, not complete genome-wide or regional summary statistics; not LD-independent multiple-IV input"}
    write_json(out / "schema.json", schema)
    outputs = [out / name for name in tables] + [out / "extraction_summary.json", out / "schema.json"]
    if independent:
        outputs.append(out / "independent_reader_verification.json")
    write_json(DATA / "MANIFEST.json", {
        "generated_at_utc": now(), "command": sys.argv, "python": platform.python_version(),
        "source": metadata, "source_member_sha256": member_hashes,
        "code": {str(Path(__file__).relative_to(ROOT)): sha256(__file__)},
        "validation_inputs": {str(p.relative_to(ROOT)): sha256(p) for p in (VERIFIED, OLD)},
        "outputs": {str(p.relative_to(ROOT)): {"sha256": sha256(p), "size_bytes": p.stat().st_size} for p in outputs},
        "table_row_counts": {name: len(rows) for name, rows in tables.items()},
        "portal_used": False, "independent_reader_verification_performed": independent is not None,
    })
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
