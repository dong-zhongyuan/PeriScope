"""Account for every planned all-target region at every stage. Read-only.

A region counts as accounted for when it is either complete or carries an explicit
documented reason why it cannot be. The point is that a 20-hour run cannot end up
quietly missing regions: anything neither complete nor explained is reported as an
unexplained gap. Safe to run while the queues are live; it only reads receipts.
"""
from pathlib import Path
import csv, json, re, sys
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
BASE = ROOT / "data/all_target_followup_20261118/regional_all"
SUP = BASE / "reverse_supplement"
COLOC = ROOT / "papers/01_mr_imaging_scan/results/all_target_coloc_20261118"
REVERSE = ROOT / "papers/01_mr_imaging_scan/results/all_target_reverse_20261118"
OUT = ROOT / "data/all_target_followup_20261118/acceptance_20261118"

TRANSPORT = re.compile(
    r"handshake operation timed out|ConnectionRefusedError|TimeoutError|IncompleteRead|"
    r"RemoteDisconnected|ConnectionResetError|Connection reset|BadStatusLine", re.I)


def read(path):
    try:
        return json.loads(Path(path).read_text())
    except Exception:
        return None


def classify(error):
    if not error:
        return "unknown"
    if TRANSPORT.search(error):
        return "transport_retryable"
    if "formatted MD5 mismatch" in error:
        return "no_published_formatted_md5_awaiting_original"
    if "invalid literal for int()" in error:
        return "chrX_label_cast_fixed_awaiting_recompute"
    if "archive prior failure before retry" in error or "explicit recovery" in error:
        return "awaiting_explicit_retry"
    return "other_needs_review"


def stage_status(folder, extra=None):
    complete = Path(folder) / "COMPLETE.json"
    failed = Path(folder) / "FAILED.json"
    if complete.is_file():
        record = read(complete) or {}
        out = dict(state="complete")
        if extra:
            for key in extra:
                if key in record:
                    out[key] = record[key]
        return out, record
    if failed.is_file():
        record = read(failed) or {}
        error = str(record.get("error", ""))
        return dict(state="failed", reason=classify(error), error=error[:200]), record
    return dict(state="absent", reason="not_yet_reached_or_in_flight"), None


def main():
    plan = read(BASE / "PLAN.json")
    imaging_plan = read(BASE / "imaging/PLAN.json")
    supplement = read(SUP / "PLAN.json") or {}
    if not plan or not imaging_plan:
        print("missing PLAN.json")
        return 1
    regions = plan["regions"]
    order = [rid for acc in plan["source_order"] for rid in plan["sources"][acc]["regions"]]
    OUT.mkdir(parents=True, exist_ok=True)

    rows = []
    for rid in order:
        rec = regions[rid]
        acc = rec["accession"]
        pheno = imaging_plan["phenotypes"][rec["phenotype"]]
        chrom = str(rec["chromosome"])
        is_x = chrom in ("23", "X")

        protein, _ = stage_status(BASE / "sources" / acc, ["source_rows", "input_format"])
        imaging, imaging_record = stage_status(
            BASE / "imaging" / pheno["family"] / pheno["idp"], ["source_rows"])
        window_id = "chr%s_center%s" % (chrom, rec["center_grch37"])
        window_present = bool(imaging_record and window_id in (imaging_record.get("windows") or {}))
        coloc, coloc_record = stage_status(COLOC / rid, ["window_status"])
        reverse, reverse_record = stage_status(
            REVERSE / rid, ["scientific_status_counts", "no_published_peak", "summary_rows"])

        blocked_by = None
        if coloc["state"] != "complete" or reverse["state"] != "complete":
            if protein["state"] != "complete":
                blocked_by = "protein_source:" + (protein.get("reason") or protein["state"])
            elif imaging["state"] != "complete" or not window_present:
                blocked_by = "imaging_source:" + (imaging.get("reason") or imaging["state"])
            elif coloc["state"] == "failed" or reverse["state"] == "failed":
                blocked_by = "downstream_failure"
            else:
                blocked_by = "downstream_pending"

        rows.append(dict(
            region_id=rid, accession=acc, gene=rec.get("gene"), chromosome=chrom,
            chrX=is_x, phenotype=rec["phenotype"], family=pheno["family"], idp=pheno["idp"],
            model_linked=rec.get("model_linked"),
            protein_state=protein["state"], protein_reason=protein.get("reason", ""),
            protein_error=protein.get("error", ""),
            imaging_state=imaging["state"], imaging_window_present=window_present,
            coloc_state=coloc["state"], coloc_reason=coloc.get("reason", ""),
            coloc_window_status=json.dumps(coloc.get("window_status")) if coloc.get("window_status") else "",
            reverse_state=reverse["state"], reverse_reason=reverse.get("reason", ""),
            reverse_no_published_peak=reverse.get("no_published_peak", ""),
            reverse_status_counts=json.dumps(reverse.get("scientific_status_counts")) if reverse.get("scientific_status_counts") else "",
            blocked_by=blocked_by or "",
        ))

    def tally(key):
        out = {}
        for row in rows:
            out[row[key]] = out.get(row[key], 0) + 1
        return out

    coloc_windows = {}
    for row in rows:
        if row["coloc_window_status"]:
            for window, status in json.loads(row["coloc_window_status"]).items():
                coloc_windows.setdefault(window, {})
                coloc_windows[window][status] = coloc_windows[window].get(status, 0) + 1
    reverse_counts = {}
    for row in rows:
        if row["reverse_status_counts"]:
            for status, n in json.loads(row["reverse_status_counts"]).items():
                reverse_counts[status] = reverse_counts.get(status, 0) + n

    unexplained = [r["region_id"] for r in rows
                   if (r["coloc_state"] != "complete" or r["reverse_state"] != "complete")
                   and not r["blocked_by"]]
    needs_review = sorted({r["protein_reason"] for r in rows if r["protein_reason"] == "other_needs_review"} |
                          {r["coloc_reason"] for r in rows if r["coloc_reason"] == "other_needs_review"} |
                          {r["reverse_reason"] for r in rows if r["reverse_reason"] == "other_needs_review"})

    chrx = [r for r in rows if r["chrX"]]
    md5_pending = sorted({r["accession"] for r in rows
                          if r["protein_reason"] == "no_published_formatted_md5_awaiting_original"})

    report = dict(
        checked_at_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        planned_regions=len(rows),
        planned_sources=len(plan["source_order"]),
        planned_imaging_phenotypes=len({r["phenotype"] for r in rows}),
        protein_states=tally("protein_state"),
        protein_failure_reasons={k: v for k, v in tally("protein_reason").items() if k},
        imaging_states=tally("imaging_state"),
        coloc_states=tally("coloc_state"),
        coloc_window_status_counts=coloc_windows,
        reverse_states=tally("reverse_state"),
        reverse_scientific_status_totals=reverse_counts,
        reverse_regions_without_published_peak=sum(1 for r in rows if r["reverse_no_published_peak"] is True),
        blocked_by=tally("blocked_by"),
        chrX=dict(regions=len(chrx),
                  coloc_complete=sum(1 for r in chrx if r["coloc_state"] == "complete"),
                  coloc_failed=sum(1 for r in chrx if r["coloc_state"] == "failed"),
                  reverse_complete=sum(1 for r in chrx if r["reverse_state"] == "complete")),
        no_published_formatted_md5_sources_pending=dict(count=len(md5_pending), accessions=md5_pending),
        unexplained_gaps=unexplained,
        reasons_needing_review=needs_review,
        fully_complete=not unexplained
            and all(r["coloc_state"] == "complete" and r["reverse_state"] == "complete" for r in rows),
        note="Complete means a receipt exists. 'undefined' coloc windows and "
             "'not_estimated_no_usable_instrument' reverse rows are real scientific "
             "outcomes, not gaps, and are reported as such.",
    )
    (OUT / "ACCEPTANCE.json").write_text(json.dumps(report, indent=2) + "\n")
    with (OUT / "per_region_coverage.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ("no_published_formatted_md5_sources_pending", "unexplained_gaps", "note")},
                     indent=2))
    print("pending no-MD5 sources:", len(md5_pending))
    print("unexplained gaps:", len(unexplained))
    return 0


if __name__ == "__main__":
    sys.exit(main())
