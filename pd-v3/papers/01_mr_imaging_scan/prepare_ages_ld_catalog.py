#!/usr/bin/env python3
"""Inventory a published AGES LD snapshot; download only its small SNP indices."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import collections
import csv
import hashlib
import json
import time
import urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "data/ages_ld_catalog"
SOURCE = ROOT / "data/ages_public/LD_RESOURCE.json"

def now():
    return datetime.now(timezone.utc).isoformat()

def digest(raw, method="sha256"):
    return hashlib.new(method, raw).hexdigest()

def write_json(path, value):
    temp = path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(value, indent=2)+"\n")
    temp.replace(path)

def main():
    source = json.loads(SOURCE.read_text())
    assert source["id"] == 5711426
    files = source["files"]
    assert len({f["key"] for f in files}) == len(files) == 119
    indices = [f for f in files if f["key"].endswith(".z")]
    matrices = {f["key"].removesuffix(".bcor"): f for f in files if f["key"].endswith(".bcor")}
    assert len(indices) == 59 and len(matrices) == 60
    assert sum(f["size"] for f in indices) == 5326961
    (OUT / "download").mkdir(parents=True, exist_ok=True)
    state = {"status": "acquiring_SNP_indices", "complete": False, "started_at_utc": now(),
             "published_files": 119, "planned_indices": 59, "planned_index_bytes": 5326961,
             "matrices_downloaded_by_this_script": 0, "source_sha256": digest(SOURCE.read_bytes()),
             "code_sha256": digest(Path(__file__).read_bytes()), "completed": [], "errors": {}}
    write_json(OUT / "RUN.json", state)
    def acquire(f):
        p = OUT / "download" / f["key"]
        expected = f["checksum"].removeprefix("md5:")
        receipt_path = p.with_suffix(p.suffix+".receipt.json")
        if p.exists():
            raw = p.read_bytes()
            assert len(raw) == f["size"] and digest(raw,"md5") == expected
            assert receipt_path.exists(), "Existing file lacks acquisition receipt"
            return f["key"]
        for attempt in range(2):
            try:
                with urllib.request.urlopen(f["links"]["self"], timeout=35) as response:
                    raw = response.read(f["size"]+1)
                    assert len(raw) == f["size"] and digest(raw,"md5") == expected
                    assert raw.splitlines()[0].split() == [b"rsid",b"chromosome",b"position",b"allele1",b"allele2",b"maf"]
                    tmp = p.with_suffix(p.suffix+".part")
                    tmp.write_bytes(raw); tmp.replace(p)
                    write_json(receipt_path, {"source_record":5711426,"file":f,"received_at_utc":now(),
                        "final_url":response.url,"headers":dict(response.headers),"sha256":digest(raw),"source_md5_verified":True})
                return f["key"]
            except Exception:
                if attempt: raise
                time.sleep(2)
    with ThreadPoolExecutor(max_workers=6) as pool:
        futures = {pool.submit(acquire,f):f["key"] for f in indices}
        for future in as_completed(futures):
            name=futures[future]
            try: state["completed"].append(future.result())
            except Exception as e: state["errors"][name]=repr(e)
            state["updated_at_utc"]=now()
            write_json(OUT / "RUN.json",state)
            if len(state["completed"])%10==0 or state["errors"]:
                print(json.dumps({"completed":len(state["completed"]),"errors":state["errors"]}),flush=True)
    if state["errors"]:
        state["status"]="source_acquisition_incomplete";write_json(OUT/"RUN.json",state)
        raise RuntimeError(state["errors"])
    rows=[]
    for f in sorted(indices,key=lambda f:f["key"]):
        raw=(OUT/"download"/f["key"]).read_bytes()
        assert len(raw)==f["size"] and digest(raw,"md5")==f["checksum"].split(":")[1]
        lines=raw.decode().splitlines()
        assert lines[0].split()==["rsid","chromosome","position","allele1","allele2","maf"]
        parsed=[line.split() for line in lines[1:]]
        assert parsed and all(len(line)==6 for line in parsed)
        keys=[tuple(line[:5]) for line in parsed]
        assert len(set(keys))==len(keys)
        chromosomes=sorted({line[1] for line in parsed})
        assert len(chromosomes)==1
        positions=[int(line[2]) for line in parsed]
        mafs=[float(line[5]) for line in parsed]
        assert min(positions)>0 and all(0<=m<=.5 for m in mafs)
        locus=f["key"].removesuffix(".z");matrix=matrices[locus]
        rows.append({"locus":locus,"chromosome":chromosomes[0],"minimum_position":min(positions),"maximum_position":max(positions),
                     "SNP_rows":len(parsed),"unique_SNP_keys":len(set(keys)),"minimum_MAF":min(mafs),"maximum_MAF":max(mafs),
                     "index_file":f["key"],"index_bytes":f["size"],"index_md5":f["checksum"],"index_sha256":digest(raw),
                     "matrix_file":matrix["key"],"matrix_bytes":matrix["size"],"matrix_md5":matrix["checksum"],
                     "matrix_url":matrix["links"]["self"],"matrix_downloaded":False,"signed_LD_verified":False})
    covered={row["locus"] for row in rows}
    for locus,matrix in sorted(matrices.items()):
        if locus not in covered:
            row={k:"" for k in rows[0]};row.update(locus=locus,matrix_file=matrix["key"],matrix_bytes=matrix["size"],
                matrix_md5=matrix["checksum"],matrix_url=matrix["links"]["self"],matrix_downloaded=False,signed_LD_verified=False)
            rows.append(row)
    output=OUT/"published_region_coverage.csv"
    with output.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n");writer.writeheader();writer.writerows(rows)
    with output.open() as f: stored=list(csv.DictReader(f))
    assert len(stored)==60 and {row["locus"] for row in stored}==set(matrices)
    for row in stored:
        if not row["index_file"]:continue
        # Separate DictReader pass verifies regions against each original index.
        with (OUT/"download"/row["index_file"]).open() as f:
            original=list(csv.DictReader((" ".join(line.split()) for line in f),delimiter=" "))
        assert len(original)==int(row["SNP_rows"])
        assert {x["chromosome"] for x in original}=={row["chromosome"]}
        assert min(int(x["position"]) for x in original)==int(row["minimum_position"])
        assert max(int(x["position"]) for x in original)==int(row["maximum_position"])
    verification={"complete":True,"verified_at_utc":now(),"published_record":5711426,"record_revision":source.get("revision"),
        "source_snapshot_sha256":digest(SOURCE.read_bytes()),"downloaded_indices":len(indices),"downloaded_index_bytes":5326961,
        "all_published_index_sizes_and_MD5_verified":True,"index_regions_by_chromosome":dict(collections.Counter(row["chromosome"] for row in rows if row["index_file"])),
        "total_index_SNP_rows":sum(row["SNP_rows"] for row in rows if row["index_file"]),
        "matrix_without_published_index":sorted(set(matrices)-covered),"matrices_downloaded":0,"signed_LD_verified":False,
        "index_has_beta_or_se":False,"index_has_sample_N":False,"all_772_study_loci_covered":False,
        "output_sha256":{output.name:digest(output.read_bytes())},"code_sha256":digest(Path(__file__).read_bytes())}
    write_json(OUT/"VERIFICATION.json",verification)
    state.update(status="complete_source_coverage_catalog",complete=True,finished_at_utc=now())
    write_json(OUT/"RUN.json",state)
    print(json.dumps(verification,indent=2),flush=True)

if __name__=="__main__": main()
