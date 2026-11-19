"""Check that the chrX fix left every autosomal colocalisation result unchanged.

For each region that already had a COMPLETE.json before the fix, compare the
recomputed outputs_sha256 against the archived manifest. INPUT_AUDIT.json is
expected to differ because it embeds the dependency sha256 list, which now carries
the edited worker and the added normaliser dependency; every scientific output file
must be byte-identical. chrX regions had no pre-fix receipt, so they are reported
separately as newly computed.
"""
from pathlib import Path
import hashlib, json, os, sys, time
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
OUT = ROOT / "papers/01_mr_imaging_scan/results/all_target_coloc_20261118"
ARCHIVE = OUT / "pre_chrX_fix_20261118"
REPORT = ARCHIVE / "REGRESSION.json"
POLL = 120
DEADLINE_HOURS = 36
EXPECTED_DIFFERENT = {"INPUT_AUDIT.json"}


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def dump(path, payload):
    tmp = Path(path).with_name(Path(path).name + ".pending")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    tmp.replace(path)


def main():
    manifest = json.loads((ARCHIVE / "MANIFEST.json").read_text())["regions"]
    plan = json.loads((ROOT / "data/all_target_followup_20261118/regional_all/PLAN.json").read_text())
    chrx = sorted(
        rid for rid, rec in plan["regions"].items()
        if str(rec["chromosome"]) in ("23", "X")
    )
    state = dict(
        pid=os.getpid(),
        started_at_utc=utc(),
        purpose="autosomal byte-identity after the chrX chromosome-label fix",
        pre_fix_regions=len(manifest),
        chrX_regions_planned=len(chrx),
        expected_to_differ=sorted(EXPECTED_DIFFERENT),
        identical=[],
        differing={},
        failed_after_fix={},
        chrX_computed={},
        pending=sorted(manifest),
    )
    dump(REPORT, state)
    deadline = time.time() + DEADLINE_HOURS * 3600
    while time.time() < deadline:
        for rid in list(state["pending"]):
            folder = OUT / rid
            complete = folder / "COMPLETE.json"
            failed = folder / "FAILED.json"
            if complete.is_file():
                new = json.loads(complete.read_text())["outputs_sha256"]
                old = manifest[rid]["outputs_sha256"] or {}
                changed = {
                    name: dict(before=value, after=new.get(name))
                    for name, value in old.items()
                    if new.get(name) != value and name not in EXPECTED_DIFFERENT
                }
                missing = [name for name in old if name not in new]
                if changed or missing:
                    state["differing"][rid] = dict(changed=changed, missing=missing,
                                                   chromosome=manifest[rid]["chromosome"])
                else:
                    state["identical"].append(rid)
                state["pending"].remove(rid)
            elif failed.is_file():
                record = json.loads(failed.read_text())
                state["failed_after_fix"][rid] = str(record.get("error"))[:300]
                state["pending"].remove(rid)
        for rid in chrx:
            if rid in state["chrX_computed"]:
                continue
            complete = OUT / rid / "COMPLETE.json"
            failed = OUT / rid / "FAILED.json"
            if complete.is_file():
                record = json.loads(complete.read_text())
                audit = json.loads((OUT / rid / "INPUT_AUDIT.json").read_text())
                state["chrX_computed"][rid] = dict(
                    window_status=record.get("window_status"),
                    shared_snps=audit.get("shared_snps"),
                    match_counts=audit.get("match_counts"),
                    chromosome_scope_recorded=bool(audit.get("chromosome_scope")),
                )
            elif failed.is_file():
                state["chrX_computed"][rid] = dict(
                    failed=str(json.loads(failed.read_text()).get("error"))[:300]
                )
        state.update(
            updated_at_utc=utc(),
            identical_count=len(state["identical"]),
            differing_count=len(state["differing"]),
            pending_count=len(state["pending"]),
            chrX_done=len(state["chrX_computed"]),
        )
        dump(REPORT, state)
        if not state["pending"] and len(state["chrX_computed"]) == len(chrx):
            state.update(
                finished_at_utc=utc(),
                complete=True,
                autosomal_identity_holds=not state["differing"] and not state["failed_after_fix"],
            )
            dump(REPORT, state)
            return 0
        time.sleep(POLL)
    state.update(finished_at_utc=utc(), complete=False, note="deadline reached with regions still pending")
    dump(REPORT, state)
    return 1


if __name__ == "__main__":
    sys.exit(main())
