from pathlib import Path
import pandas as pd,json,shutil
P=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
O=P/'spatial_subtype_author_20261003/diagnostics';O.mkdir(exist_ok=True)
branches=['spatial_subtype_20261003','spatial_subtype_markerpanel_20261003','spatial_subtype_paired_markers_20261003','spatial_subtype_corrected_20261003','spatial_subtype_author_20261003']
rows=[]
for b in branches:
 d=P/b;q=d/'state_identifiability_QC.csv'
 if not q.exists():continue
 x=pd.read_csv(q);x.insert(0,'reference_version',b);rows.append(x)
 if b!=branches[-1]:
  dest=O/b;dest.mkdir(exist_ok=True)
  for name in ['protocol.json','state_identifiability_QC.csv','input_report.json']:
   if (d/name).exists():shutil.copy2(d/name,dest/name)
if rows:pd.concat(rows,ignore_index=True).to_csv(O/'reference_QC_comparison.csv',index=False)
(O/'selection_note.json').write_text(json.dumps(dict(final_reference='Author-annotated SN subtypes with independent reference and QC donors',selection_based_on='reference metadata, biological annotation, heldout mixture recovery; no candidate p/q or pass counts',failed_reference_use='diagnostics only; never merged as biological candidate evidence'),indent=2))
