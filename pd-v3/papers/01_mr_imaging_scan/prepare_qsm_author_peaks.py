"""Read-only extraction of the complete published QSM/T2* Table 2, retaining all source rows."""
from pathlib import Path
import csv,hashlib,json,re,zipfile,xml.etree.ElementTree as ET
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[2];BASE=ROOT/"data/all_target_followup_20261118/qsm_published_peaks"
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 path=BASE/"41593_2022_1074_MOESM3_ESM.xlsx";before=sha(path)
 receipt=json.loads((BASE/"SINGLE_FILE_ACQUISITION.json").read_text())
 successful=[v for v in receipt["attempts"] if v["complete"]]
 assert successful and all(v["sha256"]==before for v in successful)
 ns={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main"};sheets={}
 with zipfile.ZipFile(path) as z:
  assert z.testzip() is None
  strings=["".join(t.text or "" for t in x.findall(".//m:t",ns)) for x in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si",ns)]
  rels={x.attrib["Id"]:x.attrib["Target"] for x in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
  for s in ET.fromstring(z.read("xl/workbook.xml")).findall("m:sheets/m:sheet",ns):
   target=rels[s.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]]
   target=target.lstrip("/") if target.startswith("/") else "xl/"+target;rows={}
   for row in ET.fromstring(z.read(target)).findall("m:sheetData/m:row",ns):
    vals={}
    for cell in row:
     assert cell.find("m:f",ns) is None,"Source formula requires separate handling"
     v=cell.find("m:v",ns)
     if v is not None:
      value=v.text
      if cell.attrib.get("t")=="s":value=strings[int(value)]
      vals[re.sub("[0-9]","",cell.attrib["r"])]=value
     elif cell.attrib.get("t")=="inlineStr":vals[re.sub("[0-9]","",cell.attrib["r"]) ]="".join(t.text or "" for t in cell.findall(".//m:t",ns))
    if vals:rows[int(row.attrib["r"])]=vals
   sheets[s.attrib["name"]]=rows
 rows=sheets["Supplementary Table 2"];header=rows[7]
 expected={"A":"Cluster index","B":"pheno","C":"Modality","D":"Region","E":"Side","F":"rsid","G":"chr","H":"pos","I":"a1","J":"a2","K":"af","L":"info","N":"beta","O":"se","Q":"beta.repro","R":"se.repro","S":"nominal"}
 assert all(header[k]==v for k,v in expected.items())
 assert "discovery" in rows[1]["A"] and "replication" in rows[1]["A"] and "a2 allele" in rows[1]["A"]
 records=[]
 for rownum,cell in rows.items():
  if rownum<=7:continue
  assert set(cell)==set("ABCDEFGHIJKLMNOPQRST"),(rownum,list(cell))
  assert re.fullmatch("V[0-9]{4}",cell["B"])
  rec=dict(source_sheet="Supplementary Table 2",source_row=rownum,source_range=f"A{rownum}:T{rownum}",
   cluster_index=int(cell["A"]),idp=cell["B"][1:],source_pheno=cell["B"],modality=cell["C"],region=cell["D"],side=cell["E"],
   source_variant_id=cell["F"],chromosome="X" if cell["G"]=="X" else str(int(cell["G"])),
   position_grch37=int(cell["H"]),other_allele=cell["I"],effect_allele=cell["J"],
   effect_allele_frequency=float(cell["K"]),imputation_info=float(cell["L"]),discovery_neglog10_p=float(cell["M"]),
   discovery_beta=float(cell["N"]),discovery_se=float(cell["O"]),replication_neglog10_p=float(cell["P"]),
   replication_beta=float(cell["Q"]),replication_se=float(cell["R"]),source_nominal_replication_flag=int(cell["S"]),
   author_beta_signs=cell["T"],frequency_status="author_table_af; effect_allele_orientation_pending_reference_check",
   ld_independence_verified=False,selected_for_reverse_mr=False)
  assert rec["position_grch37"]>0 and rec["effect_allele"]!=rec["other_allele"]
  assert 0<=rec["effect_allele_frequency"]<=1 and 0<=rec["imputation_info"]<=1
  assert rec["discovery_se"]>0 and rec["replication_se"]>0 and rec["discovery_neglog10_p"]>7.5
  assert np.isfinite([rec[k] for k in ["effect_allele_frequency","imputation_info","discovery_neglog10_p","discovery_beta","discovery_se","replication_neglog10_p","replication_beta","replication_se"]]).all()
  rec["phenotype"]="OxfordQSM2022:discovery:"+rec["idp"]
  rec["variant_key"]=("23" if rec["chromosome"]=="X" else rec["chromosome"])+":"+str(rec["position_grch37"])+":"+ "_".join(sorted([rec["effect_allele"],rec["other_allele"]]))
  records.append(rec)
 table=pd.DataFrame(records);assert len(table)==517 and not table.duplicated(["phenotype","variant_key"]).any()
 table.to_csv(BASE/"all_517_author_peaks.csv",index=False)
 planpath=ROOT/"data/all_target_followup_20261118/regional_all/imaging/PLAN.json";plan=json.loads(planpath.read_text())
 requested={k:v for k,v in plan["phenotypes"].items() if v["family"]=="OxfordQSM2022"}
 subset=table[table.phenotype.isin(requested)].copy();subset.to_csv(BASE/"all_requested_QSM_author_peaks.csv",index=False)
 coverage=[]
 for phenotype,meta in requested.items():
  n=int(table.phenotype.eq(phenotype).sum())
  for rid in meta["regions"]:coverage.append(dict(region_id=rid,phenotype=phenotype,source_peak_rows=n,
   status="published_Table2_available" if n else "no_entry_in_complete_author_Table2",reverse_MR_computed=False))
 pd.DataFrame(coverage).to_csv(BASE/"all_22_requested_region_coverage.csv",index=False)
 third=sheets["Supplementary Table 3"];thirdrows=sum(k>7 for k in third)
 assert sha(path)==before
 report=dict(complete=True,source_url=successful[0]["url"],workbook_sha256=before,source_sheet="Supplementary Table 2",source_range="A8:T524",
  source_caption=rows[1]["A"],source_rows=len(table),unique_IDPs=table.idp.nunique(),requested_IDPs=len(requested),
  requested_region_rows=len(coverage),requested_author_peak_rows=len(subset),requested_unique_variants=subset.variant_key.nunique(),
  source_peak_threshold="Published Table 2 lead associations with -log10P > 7.5; not exhaustive GWAS IV discovery",
  effect_allele="a2, explicitly defined in source caption",frequency_orientation="Pending full reference check before infer policy",
  replication_beta_and_SE_available=True,replication_archive_download_required_for_published_peaks=False,
  no_new_variant_discovery=True,no_P_or_target_cutoff_applied_to_source_table=True,full_original_workbook_unchanged=True,
  not_added_to_frozen_running_PLAN=True,followup="Separate source supplement; existing acquisition plans and code hashes stay frozen",
  supplementary_Table3_discrepancy=dict(caption_claimed_associations=407,actual_nonempty_data_rows=thirdrows,
   action="Preserved source discrepancy; narrower Table 3 is not used for peak selection"),
  input_sha256={str(planpath.relative_to(ROOT)):sha(planpath)},
  outputs_sha256={p.name:sha(p) for p in [BASE/"all_517_author_peaks.csv",BASE/"all_requested_QSM_author_peaks.csv",BASE/"all_22_requested_region_coverage.csv"]},
  code_sha256=sha(__file__))
 (BASE/"AUTHOR_PEAK_MANIFEST.json").write_text(json.dumps(report,indent=2)+"\n")
 print(json.dumps({k:report[k] for k in ["source_rows","unique_IDPs","requested_IDPs","requested_region_rows","requested_author_peak_rows","requested_unique_variants","supplementary_Table3_discrepancy"]}),flush=True)
if __name__=="__main__":main()
