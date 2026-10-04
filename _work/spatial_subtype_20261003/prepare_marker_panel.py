"""Reference-only composition markers, reserved from program evaluation."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
P=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
base=P/'spatial_subtype_20261003'
out=P/'spatial_subtype_markerpanel_20261003';inp=out/'inputs';inp.mkdir(parents=True,exist_ok=True)
assert not (out/'samples').exists(),'Do not alter a reference after fitting.'
x=pd.read_csv(base/'inputs/reference_profiles.csv.gz',index_col=0)
eligible=~x.index.str.startswith(('MT-','RPL','RPS'))
chosen=[]
for st in x:
    other=x.drop(columns=st).max(axis=1)
    rank=np.log2((x[st]+1e-6)/(other+1e-6))
    ids=rank[eligible & (x[st]>=1e-4)].sort_values(ascending=False).head(30).index
    chosen.extend(dict(gene=g,state=st,contrast='versus_max_other_type',log2FC=rank[g]) for g in ids)
for st,other in [('microglia_homeostatic','microglia_mhc2'),('microglia_mhc2','microglia_homeostatic')]:
    rank=np.log2((x[st]+1e-6)/(x[other]+1e-6))
    ids=rank[eligible & (x[st]>=5e-5)].sort_values(ascending=False).head(40).index
    chosen.extend(dict(gene=g,state=st,contrast='versus_other_microglia_state',log2FC=rank[g]) for g in ids)
markers=pd.DataFrame(chosen);markers.to_csv(out/'reference_marker_selection.csv',index=False)
reserved=set(markers.gene)
evaluation=set((base/'inputs/evaluation_genes.txt').read_text().splitlines())-reserved
for p in (base/'inputs').iterdir():
    if p.name not in ['deconvolution_genes.txt','evaluation_genes.txt']:(inp/p.name).symlink_to(p)
(inp/'deconvolution_genes.txt').write_text('\n'.join(sorted(reserved))+'\n')
(inp/'evaluation_genes.txt').write_text('\n'.join(sorted(evaluation))+'\n')
protocol=json.loads((base/'protocol.json').read_text())
protocol.update(version='20261003_spatial_subtype_v3_reserved_markers',
  amendment='V2 held-out composition QC showed microglia confusion after excluding the entire model vocabulary. Reserve reference-only markers and omit these genes from program evaluation; candidate outcomes not inspected.',
  deconvolution_genes='30 reference-enriched markers/type against max other type (mean fraction>=1e-4), plus 40 markers/microglia state against the other state (>=5e-5); no MT/RPL/RPS; markers reserved from evaluation',
  adoption='Global C7 adoption requires usable QC and donor coverage for each recipient state represented in the prior A6 set; otherwise report a spatial-supported subset and retain A6 as the global screen')
(out/'protocol.json').write_text(json.dumps(protocol,indent=2))
rows=[]
for tag,rel in [('raw','programs.json'),('relative','target_specific/programs.json')]:
    for key,genes in json.loads((P/rel).read_text()).items():
        n=len(set(genes)&evaluation)
        rows.append(dict(combo=tag+'::'+key,total_genes=len(genes),evaluation_genes=n,coverage=n/len(genes)))
coverage=pd.DataFrame(rows);coverage.to_csv(out/'marker_reservation_program_coverage.csv',index=False)
(out/'input_report.json').write_text(json.dumps(dict(
  source_input_directory=str(base/'inputs'),n_reserved_markers=len(reserved),n_evaluation_genes=len(evaluation),
  program_coverage_eligible=int(((coverage.evaluation_genes>=4)&(coverage.coverage>=.5)).sum()),
  n_programs=len(coverage),selection_uses_candidate_statistics=False,
  protocol_sha256=hashlib.sha256((out/'protocol.json').read_bytes()).hexdigest()),indent=2))
print((out/'input_report.json').read_text())
