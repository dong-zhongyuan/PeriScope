"""Reference-only composition markers, reserved from program evaluation."""
from pathlib import Path
import json,hashlib,sys
import numpy as np,pandas as pd
import h5py
from anndata.io import read_elem
from anndata._core.sparse_dataset import sparse_dataset
from scipy import sparse
sys.path.insert(0,'/public/home/mengxl/dzy/pd_product/_work/optimization_20261003')
from input_utils import canonical_symbols
P=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
base=P/'spatial_subtype_20261003'
out=P/'spatial_subtype_paired_markers_20261003';inp=out/'inputs';inp.mkdir(parents=True,exist_ok=True)
assert not (out/'samples').exists(),'Do not alter a reference after fitting.'
x=pd.read_csv(base/'inputs/reference_profiles.csv.gz',index_col=0)
states=['microglia_homeostatic','microglia_mhc2']
m=pd.read_csv(base/'inputs/reference_selected_cells.csv.gz',index_col=0)
m=m[m.state.isin(states)]
paired=sorted(set(m[m.state.eq(states[0])].donor)&set(m[m.state.eq(states[1])].donor))
m=m[m.donor.isin(paired)]
brain=P.parents[1]/'processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad'
with h5py.File(brain) as f:
    obs=read_elem(f['obs']);var=read_elem(f['var'])
    ii=np.sort(obs.index.get_indexer(m.index));R=sparse_dataset(f['layers/counts'])[ii].astype(float)
    m=m.loc[obs.index[ii]]
genes=canonical_symbols(var.gene_symbol.astype(str));gi={g:i for i,g in enumerate(x.index)}
sel=[i for i,g in enumerate(genes) if g in gi]
C=sparse.csr_matrix((np.ones(len(sel)),(sel,[gi[genes[i]] for i in sel])),shape=(len(genes),len(x)))
R=(R@C).tocsr();means={}
for st in states:
    means[st]=[]
    for donor in paired:
        ids=np.where(m.state.eq(st)&m.donor.eq(donor))[0]
        means[st].append(np.asarray(R[ids].multiply((1/m.nUMI.iloc[ids].to_numpy())[:,None]).mean(0)).ravel())
    means[st]=np.stack(means[st]);x[st]=means[st].mean(0)
x.to_csv(inp/'reference_profiles.csv.gz')
sex={'XIST','TSIX','UTY','KDM5D','DDX3Y','EIF1AY','RPS4Y1','RPS4Y2','USP9Y','ZFY','TMSB4Y','NLGN4Y','LINC00278','TXLNGY','PRKY','AMELY','PCDH11Y','SRY'}
eligible=~x.index.str.startswith(('MT-','RPL','RPS','TTTY','RBMY','TSPY','DAZ','BPY2','HSFY','CDY','VCY'))&~x.index.isin(sex)
chosen=[]
for st in x:
    other=x.drop(columns=st).max(axis=1)
    rank=np.log2((x[st]+1e-6)/(other+1e-6))
    ids=rank[eligible & (x[st]>=1e-4)].sort_values(ascending=False).head(30).index
    chosen.extend(dict(gene=g,state=st,contrast='versus_max_other_type',log2FC=rank[g]) for g in ids)
for st,other in [('microglia_homeostatic','microglia_mhc2'),('microglia_mhc2','microglia_homeostatic')]:
    paired_fc=np.log2((means[st]+1e-6)/(means[other]+1e-6))
    rank=pd.Series(np.median(paired_fc,axis=0),index=x.index)
    consistent=pd.Series((paired_fc>0).sum(0)>=3,index=x.index)
    ids=rank[eligible & consistent & (x[st]>=5e-5)].sort_values(ascending=False).head(40).index
    chosen.extend(dict(gene=g,state=st,contrast='within_donor_versus_other_microglia_state',log2FC=rank[g]) for g in ids)
markers=pd.DataFrame(chosen);markers.to_csv(out/'reference_marker_selection.csv',index=False)
reserved=set(markers.gene)
evaluation=set((base/'inputs/evaluation_genes.txt').read_text().splitlines())-reserved
for p in (base/'inputs').iterdir():
    if p.name not in ['reference_profiles.csv.gz','deconvolution_genes.txt','evaluation_genes.txt']:(inp/p.name).symlink_to(p)
(inp/'deconvolution_genes.txt').write_text('\n'.join(sorted(reserved))+'\n')
(inp/'evaluation_genes.txt').write_text('\n'.join(sorted(evaluation))+'\n')
protocol=json.loads((base/'protocol.json').read_text())
protocol.update(version='20261003_spatial_subtype_v4_paired_reserved_markers',
  amendment='V2 excluded critical microglia genes. V3 marker inventory revealed sex genes because homeostatic reference donors were all female while MHC-II used both sexes. Pair reference donors for both microglia profiles and select within-donor markers; exclude explicit sex markers. No candidate outcomes inspected.',
  paired_microglia_reference_donors=paired,
  excluded_sex_marker_symbols=sorted(sex),excluded_sex_marker_prefixes=['TTTY','RBMY','TSPY','DAZ','BPY2','HSFY','CDY','VCY'],
  deconvolution_genes='30 reference-enriched markers/type against max other type (mean fraction>=1e-4), plus 40 markers/microglia state by median within-donor log2FC with positive sign in >=3/4 paired donors (>=5e-5); no MT/RPL/RPS or explicit sex markers; markers reserved from evaluation',
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
