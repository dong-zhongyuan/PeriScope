"""Include chrX in the all-target colocalisation, per the 2026-11-18 decision.

The frozen worker asserted `chromosome.astype(int)` and built variant keys with
`int(chromosome)`, so the 26 chrX regions could only raise
ValueError("invalid literal for int() with base 10: 'X'"). The two sides also label
the same chromosome differently: the AGES protein extracts carry "X", the BIG40
imaging windows carry "23". A naive string comparison would therefore have produced
zero shared SNPs and a silent "undefined", so the fix normalises both sides with
`normchr` from acquire_all_target_protein_regions.py, the same function the frozen
selection step already used to decide which variants belong to each region.

Autosomal keys are delegated to the untouched frozen `key()` so they are identical
by construction. prepare_new_scan_coloc.py, all_target_coloc.R, claudia.R, the prior
grid, the windows and every effect prior are unchanged.

Every completed region's COMPLETE.json records this worker's sha256 among its inputs,
so all of them are archived here and recomputed; verify_chrx_coloc_regression.py then
checks the recomputed autosomal outputs against the archived sha256 values.
"""
from pathlib import Path
import hashlib, json, shutil, sys
import datetime

ROOT = Path("/public/home/mengxl/dzy/pd-v3")
SCRIPTS = ROOT / "papers/01_mr_imaging_scan"
BASE = ROOT / "data/all_target_followup_20261118/regional_all"
OUT = ROOT / "papers/01_mr_imaging_scan/results/all_target_coloc_20261118"
WORKER = SCRIPTS / "run_one_all_target_coloc.py"
ARCHIVE = OUT / "pre_chrX_fix_20261118"
TAG = "chrX_int_cast_20261118"


def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(8 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dump(path, payload):
    tmp = Path(path).with_name(Path(path).name + ".pending")
    tmp.write_text(json.dumps(payload, indent=2) + "\n")
    tmp.replace(path)


def running(pid):
    try:
        raw = Path("/proc/%d/cmdline" % int(pid)).read_bytes()
    except Exception:
        return False
    return bool(raw) and b"run_all_target_coloc.py" in raw


EDITS = [
    (
        "from prepare_new_scan_coloc import key,write\n",
        "from prepare_new_scan_coloc import key as autosomal_key,write\n"
        "from acquire_all_target_protein_regions import normchr\n"
        "def chrom_label(v):\n"
        " if isinstance(v,float) and float(v).is_integer():v=int(v)\n"
        " return normchr(v)\n"
        "def region_key(row):\n"
        ' """Frozen key for numeric chromosomes; normchr label for chrX."""\n'
        " label=chrom_label(row.chromosome)\n"
        " if label.isdigit():return autosomal_key(row)\n"
        " a,b=row.effect_allele,row.other_allele\n"
        " if len(a)!=1 or len(b)!=1 or a not in 'ACGT' or b not in 'ACGT' or a==b:return None\n"
        " return f'{label}:{int(row.base_pair_location)}:'+'_'.join(sorted([a,b]))\n",
    ),
    (
        "  assert frame.chromosome.astype(int).eq(int(rec[\"chromosome\"])).all()\n",
        "  assert frame.chromosome.map(chrom_label).eq(chrom_label(rec[\"chromosome\"])).all()\n",
    ),
    (
        "  frame[\"variant_key\"]=frame.apply(key,axis=1) if len(frame) else pd.Series(dtype=object)\n",
        "  frame[\"variant_key\"]=frame.apply(region_key,axis=1) if len(frame) else pd.Series(dtype=object)\n",
    ),
    (
        " assert shared.chromosome_x.eq(shared.chromosome_y).all()\n",
        " assert shared.chromosome_x.map(chrom_label).eq(shared.chromosome_y.map(chrom_label)).all()\n",
    ),
    (
        "  decision_threshold=None)\n return audit\n",
        "  decision_threshold=None)\n"
        " if not chrom_label(rec[\"chromosome\"]).isdigit():\n"
        "  audit[\"chromosome_scope\"]=(\"chrX included on the author-published per-variant statistics as given: \"\n"
        "   \"no ploidy recoding, no sex-stratified reanalysis, and sdY.est uses the reported assay N and protein \"\n"
        "   \"frequency exactly as for autosomes; the reverse-MR chrX ploidy question remains open\")\n"
        " return audit\n",
    ),
    (
        "  dependencies=[Path(__file__),Path(__file__).with_name(\"all_target_coloc.R\"),Path(__file__).with_name(\"prepare_new_scan_coloc.py\"),ROOT/\"reference/coloc/claudia.R\",\n",
        "  dependencies=[Path(__file__),Path(__file__).with_name(\"all_target_coloc.R\"),Path(__file__).with_name(\"prepare_new_scan_coloc.py\"),\n"
        "   Path(__file__).with_name(\"acquire_all_target_protein_regions.py\"),ROOT/\"reference/coloc/claudia.R\",\n",
    ),
]


def main():
    controller = OUT / "CONTROLLER.json"
    data = json.loads(controller.read_text()) if controller.is_file() else None
    if data and running(data.get("pid")):
        print("coloc controller %s is still alive; stop it first" % data.get("pid"))
        return 1
    plan = json.loads((BASE / "PLAN.json").read_text())
    regions = plan["regions"]
    ARCHIVE.mkdir(parents=True, exist_ok=True)

    report = dict(
        applied_at_utc=utc(),
        decision="include chrX in all-target coloc; fix the chromosome cast only",
        worker_sha256_before=sha(WORKER),
        unchanged=[
            "papers/01_mr_imaging_scan/prepare_new_scan_coloc.py",
            "papers/01_mr_imaging_scan/all_target_coloc.R",
            "reference/coloc/claudia.R",
            "data/all_target_followup_20261118/regional_all/PLAN.json",
            "data/all_target_followup_20261118/regional_all/imaging/PLAN.json",
        ],
        controller_state_archived={},
        completed_regions_archived={},
        chrX_failures_archived={},
        edits_applied=[],
    )

    # 1. archive the stopped controller state and the pre-edit worker
    for name in ("CONTROLLER.json", "controller.log"):
        source = OUT / name
        if source.is_file():
            destination = ARCHIVE / name
            shutil.move(str(source), str(destination))
            report["controller_state_archived"][name] = sha(destination)
    worker_copy = ARCHIVE / (WORKER.name + ".pre_chrX_fix")
    shutil.copy2(str(WORKER), str(worker_copy))
    report["controller_state_archived"][worker_copy.name] = sha(worker_copy)

    # 2. archive every completed region receipt, then clear it so it recomputes
    manifest = {}
    for folder in sorted(OUT.iterdir()):
        if not folder.is_dir() or folder.name not in regions:
            continue
        complete = folder / "COMPLETE.json"
        if complete.is_file():
            record = json.loads(complete.read_text())
            manifest[folder.name] = dict(
                chromosome=regions[folder.name]["chromosome"],
                phenotype=regions[folder.name]["phenotype"],
                window_status=record.get("window_status"),
                outputs_sha256=record.get("outputs_sha256"),
                input_sha256=record.get("input_sha256"),
                receipt_sha256=sha(complete),
            )
            kept = folder / "COMPLETE_pre_chrX_fix.json"
            shutil.copy2(str(complete), str(kept))
            complete.unlink()
            report["completed_regions_archived"][folder.name] = sha(kept)
        failed = folder / "FAILED.json"
        if failed.is_file():
            record = json.loads(failed.read_text())
            error = str(record.get("error", ""))
            if "invalid literal for int()" not in error:
                report["chrX_failures_archived"][folder.name] = "left_in_place:" + error[:80]
                continue
            target = folder / ("attempt_%s" % TAG)
            target.mkdir(exist_ok=True)
            moved = {}
            for name in ("FAILED.json", "RUN.json", "worker.log"):
                source = folder / name
                if source.is_file():
                    destination = target / name
                    shutil.move(str(source), str(destination))
                    moved[name] = sha(destination)
            dump(target / "ARCHIVE.json", dict(archived_at_utc=utc(), previous_error=error, files=moved))
            report["chrX_failures_archived"][folder.name] = str(target.relative_to(ROOT))
    dump(ARCHIVE / "MANIFEST.json", dict(
        created_at_utc=utc(),
        purpose="pre-fix per-region receipts; recomputed autosomal outputs must match these sha256",
        regions=manifest,
    ))
    report["completed_regions_count"] = len(manifest)

    # 3. apply the edits
    text = WORKER.read_text()
    for old, new in EDITS:
        if old not in text:
            print("edit anchor not found, nothing written:\n%r" % old[:120])
            return 2
        if text.count(old) != 1:
            print("edit anchor is not unique: %r" % old[:120])
            return 2
        text = text.replace(old, new)
        report["edits_applied"].append(old.strip()[:120])
    WORKER.write_text(text)
    report["worker_sha256_after"] = sha(WORKER)
    dump(ARCHIVE / "CHRX_FIX.json", report)
    print(json.dumps(dict(
        completed_regions_archived=len(manifest),
        chrX_failures_archived=len(report["chrX_failures_archived"]),
        worker_before=report["worker_sha256_before"],
        worker_after=report["worker_sha256_after"],
    ), indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
