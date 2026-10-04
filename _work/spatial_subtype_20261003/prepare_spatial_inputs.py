"""Export real counts/reference signatures and explicitly labelled QC mixtures."""
from pathlib import Path
import sys,json,hashlib,gzip
import numpy as np,pandas as pd,h5py
from scipy import sparse,io
from anndata.io import read_elem,sparse_dataset
sys.path.insert(0,'/public/home/mengxl/dzy/pd_product/_work/optimization_20261003')
from input_utils import canonical_symbols

A=Path('/public/home/mengxl/dzy/pd_product_assets');P=A/'results/optimization_20261003'
O=P/'spatial_subtype_20261003';D=O/'inputs';D.mkdir(parents=True,exist_ok=True)
rng=np.random.default_rng(20261003)
protocol={
 'version':'20261003_spatial_subtype_v2','method':'RCTD full mode; C-SIDE intercept per spatial donor; donor-level program tests',
 'reference':'GSE178265 substantia nigra training donors; donor-balanced signatures; microglia restricted to original Microglia lineage, original Mono_macro retained as separate nuisance; Hao P1-P4 T/B/NK nuisance profiles',
 'pre_fit_amendment':'Reference inventory showed negligible T/B coverage and original Mono_macro cells absorbed into purified microglia labels. Add nuisance lineages before any RCTD or candidate result is examined.',
 'states':['astro','microglia_homeostatic','microglia_mhc2'],
 'resolution':'matches existing model recipient states; astro is a cell type, microglia split into two marker-defined states; no individual-cell spatial measurement claimed',
 'deconvolution_genes':'exclude entire 4741-gene model vocabulary, not just selected candidate genes; counts from remaining shared genes',
 'quality_mixtures':'400 RNA mixtures from brain calibration donors plus Hao P5 immune backgrounds, disjoint from reference donors; 1500 UMIs each; synthetic mixtures are QC only, never candidate evidence',
 'quality_rule_per_state':{'fraction_spearman_min':.5,'fraction_RMSE_max':.15,'minimum_QC_mixtures':100},
 'RCTD':'full mode, default gene/fold thresholds, UMI_min=100; no case/control labels used in deconvolution',
 'CSIDE':{'cell_type_threshold':25,'gene_threshold':1e-5,'weight_threshold':.8,'doublet_mode':False,'sigma_gene':True,'normalize_expr':False,
          'expression':'intercept model fit separately for each donor; all sufficiently represented mixture types included'},
 'program_test':{'minimum_donors_each_arm':3,'minimum_genes':4,'minimum_program_gene_fraction':.5,
     'gene_values':'C-SIDE estimated cell-state log-expression; nonconverged/nonfinite estimates missing',
     'competitive_null':'5000 expression-matched random gene sets,20 bins,seed42; whole6849program family and both directions',
     'C7':'state passes QC; estimable donor contrast; competitive q<=0.05; spatial within-state disease direction agrees with pre-existing Kamath same-state direction',
     'missing':'unresolved, not biological negative; no substitution of mixed tissue or different microglia state'},
 'adoption':'quality evaluated before inspecting target outcomes; seven-criterion results use unchanged C1-C6 and independently report C7 pass/not_supported/unresolved',
 'reference_source':'https://github.com/dmcable/spacexr',
}
pp=O/'protocol.json'
if pp.exists() and json.loads(pp.read_text())!=protocol:
    assert not (O/'samples').exists(), 'Do not amend after fitting'
    (O/'protocol_v1_before_reference_inventory.json').write_text(pp.read_text())
pp.write_text(json.dumps(protocol,indent=2))

def collapse_symbols(gs):
    gs=canonical_symbols(gs);unique=sorted(set(gs));ix={g:i for i,g in enumerate(unique)}
    return unique,sparse.csr_matrix((np.ones(len(gs)),(np.arange(len(gs)),[ix[g] for g in gs])),shape=(len(gs),len(unique)))
def write_mtx(name,x):
    with gzip.open(D/(name+'.mtx.gz'),'wb') as f:io.mmwrite(f,x.T.astype(np.int32),field='integer')

brain=A/'processed/gse178265/v0.1/GSE178265_sn_annotated.h5ad'
with h5py.File(brain) as f:obs=read_elem(f['obs']);var=read_elem(f['var'])
labels=obs.cell_type.astype(str).copy()
for fam in ['brain_astro','brain_microglia']:
    m=pd.read_csv(A/f'interim/v0.1/purification/c178_{fam}.csv',index_col=0);m=m[m.kept];ix=labels.index.intersection(m.index)
    labels.loc[ix]=m.loc[ix,'state_pure'].astype(str)
labels.loc[obs.cell_type.astype(str).eq('Mono_macro')]='Mono_macro'
split=json.loads((A/'interim/v0.1/splits/GSE178265_sn_donor_split_v1.json').read_text())['donors']
roles=obs.donor.astype(str).map(split);valid=obs.tissue.astype(str).eq('substantia_nigra')&~labels.isin(['unassigned','Microglia','Astrocyte'])
meta=obs[['donor','condition']].copy();meta['state']=labels;meta['role']=roles
refrows=[];support=[]
for (donor,state),g in meta[valid&roles.eq('train')].groupby(['donor','state'],observed=True):
    rr=obs.index.get_indexer(g.index)
    if len(rr)<25:continue
    chosen=rng.choice(rr,min(250,len(rr)),replace=False);refrows.extend(chosen)
    support.append(dict(donor=donor,state=state,n_available=len(rr),n_selected=len(chosen)))
support=pd.DataFrame(support);good=support.groupby('state').donor.nunique();types=sorted(good[good>=2].index)
refrows=np.sort([i for i in refrows if labels.iloc[i] in types]);refmeta=meta.iloc[refrows].copy()
calrows=np.where(valid&roles.eq('calibration')&labels.isin(types))[0]
with h5py.File(brain) as f:
    raw=sparse_dataset(f['layers/counts']);R=raw[refrows].astype(float);C=raw[calrows].astype(float)
genes,col=collapse_symbols(var.gene_symbol.astype(str));R=(R@col).tocsr();C=(C@col).tocsr()
refmeta['nUMI']=np.asarray(R.sum(1)).ravel();calmeta=meta.iloc[calrows].copy();calmeta['nUMI']=np.asarray(C.sum(1)).ravel()
# Include lymphocyte backgrounds so immune transcripts need not be assigned to
# the microglia components. These profiles do not become new candidate axes.
rawdir=A/'raw/citeseq_hao/raw_3p'
hf=pd.read_csv(next(rawdir.glob('GSM5008737_RNA_3P-features.tsv.gz')),sep='\t',header=None)
hb=pd.read_csv(next(rawdir.glob('GSM5008737_RNA_3P-barcodes.tsv.gz')),sep='\t',header=None)[0].astype(str)
hm=pd.read_csv(A/'raw/citeseq_hao/GSE164378_sc.meta.data_3P.csv.gz',index_col=0).loc[hb]
hmap={'CD4 T':'T_CD4_background','CD8 T':'T_CD8_background','B':'B_background','NK':'NK_background'}
hstate=hm['celltype.l1'].map(hmap);hd=hm['donor'].astype(str)
hi=[];hcal=[]
for (donor,st),g in pd.DataFrame({'donor':hd.to_numpy(),'state':hstate.to_numpy()}).dropna().groupby(['donor','state']):
    if donor not in ['P1','P2','P3','P4','P5']:continue
    chosen=rng.choice(g.index.to_numpy(),min(150,len(g)),replace=False)
    (hcal if donor=='P5' else hi).extend(chosen)
with gzip.open(next(rawdir.glob('GSM5008737_RNA_3P-matrix.mtx.gz')),'rb') as f:H=io.mmread(f).T.tocsr()
hgenes=canonical_symbols(hf.iloc[:,1].astype(str));gindex={g:i for i,g in enumerate(genes)}
select=[i for i,g in enumerate(hgenes) if g in gindex]
hcollapse=sparse.csr_matrix((np.ones(len(select)),(select,[gindex[hgenes[i]] for i in select])),shape=(len(hgenes),len(genes)))
for which,ixs in [('reference',np.sort(hi)),('calibration',np.sort(hcal))]:
    X=(H[ixs]@hcollapse).tocsr();lib=np.asarray(H[ixs].sum(1)).ravel()
    md=pd.DataFrame({'donor':['Hao:'+d for d in hd.iloc[ixs]],'condition':'healthy','state':hstate.iloc[ixs].to_numpy(),
        'role':which,'nUMI':lib},index=['Hao:'+b for b in hb.iloc[ixs]])
    if which=='reference':R=sparse.vstack([R,X],format='csr');refmeta=pd.concat([refmeta,md])
    else:C=sparse.vstack([C,X],format='csr');calmeta=pd.concat([calmeta,md])
types=sorted(refmeta.state.unique());del H
profiles={}
for state in types:
    means=[]
    for donor in sorted(refmeta.loc[refmeta.state.eq(state),'donor'].unique()):
        ix=np.where(refmeta.state.eq(state)&refmeta.donor.eq(donor))[0]
        means.append(np.asarray(R[ix].multiply((1/refmeta.nUMI.iloc[ix].to_numpy())[:,None]).mean(0)).ravel())
    profiles[state]=np.stack(means).mean(0)
profiles=pd.DataFrame(profiles,index=genes)
spatial=A/'processed/gse253975/v0.1/GSE253975_geomx.h5ad'
with h5py.File(spatial) as f:sobs=read_elem(f['obs']);svar=read_elem(f['var']);S=sparse_dataset(f['X'])[:].astype(float)
sg,coll=collapse_symbols(svar.index.astype(str));S=(S@coll).tocsr();sobs['nUMI']=np.asarray(S.sum(1)).ravel()
common=sorted(set(genes)&set(sg));ri=np.array([genes.index(g) for g in common]);si=np.array([sg.index(g) for g in common])
profiles=profiles.loc[common];profiles.index.name='gene';profiles.to_csv(D/'reference_profiles.csv.gz')
model=np.load(A/'interim/rescreen_20261002/brain_train.npz')['genes'].astype(str).tolist();program_genes=set()
for rel in ['programs.json','target_specific/programs.json']:program_genes.update(set().union(*map(set,json.loads((P/rel).read_text()).values())))
assert program_genes<=set(model)
deconv=[g for g in common if g not in set(model) and not g.startswith(('MT-','RPL','RPS'))]
pd.Series(common).to_csv(D/'genes.txt',index=False,header=False)
pd.Series(deconv).to_csv(D/'deconvolution_genes.txt',index=False,header=False)
pd.Series(sorted(set(common)&set(model))).to_csv(D/'evaluation_genes.txt',index=False,header=False)
coords=pd.read_csv('/public/home/mengxl/dzy/pd_product/_work/optimization_20261003/reference/spatial_coordinates.csv').set_index('spot_id')
coords=coords.loc[sobs.index];assert (coords.donor.astype(str).to_numpy()==sobs.donor.astype(str).to_numpy()).all()
sobs['x']=coords.pxl_col_in_fullres.to_numpy();sobs['y']=coords.pxl_row_in_fullres.to_numpy();sobs.index.name='barcode'
sobs.to_csv(D/'spatial_metadata.csv');write_mtx('spatial_counts',S[:,si])
# A small actual-cell Reference object satisfies package bookkeeping; the full
# donor-balanced profile matrix above is supplied to create.RCTD explicitly.
thin=np.concatenate([np.where(refmeta.state.eq(t))[0][:30] for t in types]);rm=refmeta.iloc[thin].copy();rm.index.name='barcode'
rm.to_csv(D/'reference_metadata.csv');write_mtx('reference_counts',R[thin][:,ri])
support.to_csv(D/'reference_donor_support.csv',index=False);refmeta.to_csv(D/'reference_selected_cells.csv.gz')

# Calibration-only pseudo-spots: known RNA mixture proportions, realistic depth.
rows=[];truth=[];mixmeta=[];caltypes=calmeta.state.to_numpy();targets=protocol['states'];lib=calmeta.nUMI.to_numpy()
for k in range(400):
    present=[t for t in types if (caltypes==t).sum()>=5]
    picked=list(rng.choice(present,size=min(5,len(present)),replace=False))
    if k<300:
        target=targets[k%3]
        if target in present and target not in picked:picked[0]=target
    fractions=rng.dirichlet(np.ones(len(picked))*.7);prob=np.zeros(len(genes));true={t:0. for t in types}
    for t,w in zip(picked,fractions):
        choices=rng.choice(np.where(caltypes==t)[0],5,replace=False)
        per=np.asarray(C[choices].multiply((1/lib[choices])[:,None]).mean(0)).ravel()
        per/=per.sum()
        prob+=w*per;true[t]+=float(w)
    counts=rng.multinomial(1500,prob/prob.sum());rows.append(sparse.csr_matrix(counts[ri]))
    bc=f'QC_{k:04d}';truth.append(dict(barcode=bc,**true));mixmeta.append(dict(barcode=bc,x=k%20,y=k//20,nUMI=1500,condition='QC_synthetic_only',donor='calibration_mixture'))
write_mtx('qc_counts',sparse.vstack(rows));pd.DataFrame(truth).to_csv(D/'qc_truth.csv',index=False);pd.DataFrame(mixmeta).to_csv(D/'qc_metadata.csv',index=False)
report=dict(status='complete',n_spots=len(sobs),n_reference_cells=len(refmeta),reference_donors=sorted(refmeta.donor.unique()),
  calibration_donors=sorted(calmeta.donor.unique()),reference_states=types,n_common_genes=len(common),n_deconv_genes=len(deconv),
  n_model_genes=len(model),n_program_genes=len(program_genes),deconv_program_overlap=len(set(deconv)&program_genes),
  reference_state_donors=refmeta.groupby('state').donor.nunique().to_dict(),protocol_sha256=hashlib.sha256(pp.read_bytes()).hexdigest())
(O/'input_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
