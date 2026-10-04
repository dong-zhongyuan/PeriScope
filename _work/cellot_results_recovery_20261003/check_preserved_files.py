"""Verify existing non-CellOT artifacts against the pre-recovery backup."""
import hashlib,json
from pathlib import Path
B=Path('/public/home/mengxl/dzy/pd_recovery_backup_20261003_cellot')
D=Path('/public/home/mengxl/dzy/pd_product_assets/results')
allowed={
 'benchmark_20261002':{'comparison_summary.csv','baseline_metrics.csv','seed_summary.csv','axis_summary.csv','benchmark_evaluation_summary.csv','training_summary.csv','manifest.json'},
 'pig_external_validation_20261003':{'primary_summary.csv','response_metrics.csv','output_manifest.json','all_six_axis_metrics.csv','leave_one_animal_out.csv','paired_animal_bootstrap.csv','primary_bootstrap_intervals.csv','prediction_diagnostics.csv','periscope_comparator_bootstrap.csv'},
 'blood_brain_MR_20261003':set(),'optimization_20261003':set()}
def digest(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(4*1024*1024),b''):h.update(c)
 return h.hexdigest()
counts={};changes=[];missing=[]
for group,excluded in allowed.items():
 n=0
 for old in sorted((B/group).rglob('*')):
  if not old.is_file():continue
  rel=old.relative_to(B/group)
  new=D/group/rel
  if str(rel) in excluded:
   if old.suffix=='.csv' and (not new.is_file() or not new.read_bytes().startswith(old.read_bytes())):changes.append(group+'/'+str(rel)+' [original prefix changed]')
   continue
  if not new.is_file():missing.append(group+'/'+str(rel));continue
  if old.stat().st_size!=new.stat().st_size or digest(old)!=digest(new):changes.append(group+'/'+str(rel))
  n+=1
 counts[group]=n
result=dict(identical_existing_files=counts,changed=changes,missing=missing,status='passed' if not changes and not missing else 'failed')
Path(__file__).with_name('preservation_check.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
