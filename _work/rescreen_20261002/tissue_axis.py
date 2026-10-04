"""Observed blood-brain disease concordance, independent of target nomination."""
from pathlib import Path
import json,re,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','2')
import h5py,numpy as np,pandas as pd
from anndata.io import read_elem,sparse_dataset
from scipy import sparse
from concurrent.futures import ProcessPoolExecutor
from input_utils import canonical_symbols
from tissue_axis_helpers import boot_pair

A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/rescreen_20261002';D=O/'tissue_axis';D.mkdir(exist_ok=True)
FILES={'bmb':'moquin_beaudry2025/v0.1/MB2025_pbmc_annotated.h5ad','b223':'gse223138/v0.1/GSE223138_pbmc_annotated.h5ad',
       'c178':'gse178265/v0.1/GSE178265_sn_annotated.h5ad','c157':'gse157783/v0.1/gse157783_qc.h5ad',
       'ma':'gse253975/v0.1/GSE253975_geomx.h5ad','pdd':'gse184950/v0.1/GSE184950_annotated.h5ad'}
PD={'Disease','PD','PDD','Early PD','Late PD'};CONTROL={'Control','Ctrl','CTRL','Normal'}

def metadata(cohort):
    with h5py.File(A/'processed'/FILES[cohort]) as f:return read_elem(f['obs']),read_elem(f['var'])

def prepare(cohort,genes,brain_cov):
    obs,var=metadata(cohort);symbols=canonical_symbols(var.gene_symbol.astype(str) if 'gene_symbol' in var else var.index)
    gi={g:i for i,g in enumerate(genes)};ri=[i for i,g in enumerate(symbols) if g in gi]
    collapse=sparse.csr_matrix((np.ones(len(ri)),(ri,[gi[symbols[i]] for i in ri])),shape=(len(symbols),len(genes)))
    typed=pd.Series('other',index=obs.index)
    if cohort in ['bmb','b223','c178','c157']:
        families=['blood_myeloid'] if cohort.startswith('b') else ['brain_astro','brain_microglia']
        for family in families:
            m=pd.read_csv(A/'interim/v0.1/purification'/f'{cohort}_{family}.csv',index_col=0);m=m[m.kept];ix=typed.index.intersection(m.index);typed.loc[ix]=m.loc[ix,'state_pure'].astype(str)
    elif cohort=='ma':typed[:]='whole_spatial'
    else:typed.loc[obs.cell_type.astype(str)=='Astrocyte']='astro_broad_PDD'
    ok=typed.ne('other')&obs.condition.astype(str).isin(PD|CONTROL)
    if cohort=='c178':ok &= obs.tissue.astype(str).eq('substantia_nigra')
    records=obs[['donor','condition']].copy();records['state']=typed
    counts=records[ok].groupby(['state','donor'],observed=True).size()
    keys=[(str(s),str(d)) for (s,d),n in counts.items() if n>=(1 if cohort=='ma' else 25)]
    kidx={k:i for i,k in enumerate(keys)};row_group=np.array([kidx.get((str(s),str(d)),-1) for s,d in zip(typed,obs.donor)])
    rows=np.where(ok.to_numpy()&(row_group>=0))[0];M=np.zeros((len(keys),len(genes)));N=np.zeros(len(keys),int)
    with h5py.File(A/'processed'/FILES[cohort]) as f:
        x=sparse_dataset(f['layers/counts'] if 'layers/counts' in f else f['X'])
        for start in range(0,len(rows),4096):
            rr=rows[start:start+4096];raw=x[rr].astype(float);lib=np.asarray(raw.sum(1)).ravel().clip(1)
            b=(raw@collapse).multiply((1e4/lib)[:,None]).tocsr();b.data=np.log1p(b.data)
            for group in np.unique(row_group[rr]):
                mask=row_group[rr]==group;M[group]+=np.asarray(b[mask].sum(0)).ravel();N[group]+=int(mask.sum())
    blocks={};summary=[]
    for state in sorted({s for s,d in keys}):
        ii=[i for i,k in enumerate(keys) if k[0]==state];donors=[keys[i][1] for i in ii];conditions=[];cv=[]
        for d in donors:
            row=obs[obs.donor.astype(str)==d].iloc[0];conditions.append(int(str(row.condition) in PD))
            if cohort=='bmb':cv.append([float(row.Age),1. if row.Sex=='M' else 0. if row.Sex=='F' else np.nan])
            elif cohort=='c178':cv.append(brain_cov.get(d.split('-')[-1],[np.nan]*3))
        conditions=np.array(conditions);means=M[ii]/N[ii,None]
        cov=np.array(cv,float) if cv else None
        if cov is not None:
            if not np.isfinite(cov).all():raise ValueError(('missing source covariates',cohort,state))
            cov=cov[:,cov.std(0)>1e-10];cov=(cov-cov.mean(0))/(cov.std(0)+1e-12)
        rec=dict(cohort=cohort,state=state,n_PD=int(sum(conditions==1)),n_control=int(sum(conditions==0)),donors=donors,covariates=['age','sex'] if cohort=='bmb' else ['age','sex','PMI'] if cohort=='c178' else [])
        summary.append(rec)
        if min(sum(conditions==0),sum(conditions==1))<3:continue
        blocks[state]=dict(means=means,conds=conditions,donors=donors,cov=cov)
        np.savez_compressed(D/f'{cohort}__{state}.npz',means=means,conditions=conditions,donors=donors,genes=genes,covariates=np.empty((len(ii),0)) if cov is None else cov)
    (D/f'{cohort}_support.json').write_text(json.dumps(summary,indent=2));print('TISSUE PROFILES',cohort,list(blocks),flush=True)
    return blocks

def compare(task):
    label,a,bs,adjusted=task
    variable=a['means'].var(0)>1e-12
    variable &= np.logical_or.reduce([b['means'].var(0)>1e-12 for b in bs])
    aa=dict(a,means=a['means'][:,variable]);bb=[dict(b,means=b['means'][:,variable]) for b in bs]
    try:
        res=boot_pair([aa],bb,cov_a=a['cov'] if adjusted else None,cov_b_list=[b['cov'] if adjusted else None for b in bs],n_boot=800,seed=8613)
        # These bootstrap tail fractions are descriptive support, not calibrated null p-values.
        for suffix in ['w','u']:
            if 'boot_p_pos_'+suffix in res:res['bootstrap_nonpositive_fraction_'+suffix]=res.pop('boot_p_pos_'+suffix)
    except (ValueError,np.linalg.LinAlgError) as e:res=dict(status='nonidentified_adjusted_contrast',reason=str(e))
    return dict(analysis=label,covariate_adjusted=adjusted,n_genes=int(variable.sum()),blood_donors=len(a['donors']),brain_donors=sum(len(b['donors']) for b in bs),**res)

def main():
    sets={}
    for cohort in FILES:
        obs,var=metadata(cohort);sets[cohort]=set(canonical_symbols(var.gene_symbol.astype(str) if 'gene_symbol' in var else var.index))
    genes=sorted(set.intersection(*sets.values())-{'','nan','NA','None'})
    raw=json.loads((A/'cache/geo_meta/GSE178265_human_meta.json').read_text());bc={}
    for sample in raw['human_samples']:
        ch=sample.get('characteristics',{});m=re.search(r'(\d{4})',sample.get('title',''))
        if m:bc[m.group(1)]=[float(ch.get('age','nan')),1. if ch.get('Sex')=='Male' else 0. if ch.get('Sex')=='Female' else np.nan,float(ch.get('pmi','nan'))]
    blocks={c:prepare(c,genes,bc) for c in FILES};tasks=[]
    for blood,a in blocks['bmb'].items():
        for brain in sorted(set(blocks['c178'])|set(blocks['c157'])):
            bs=[blocks[c][brain] for c in ['c178','c157'] if brain in blocks[c]]
            for adjusted in [False,True]:tasks.append((blood+'__x__brain_meta_'+brain,a,bs,adjusted))
        for c in ['ma','pdd']:
            for name,b in blocks[c].items():
                for adjusted in [False,True]:tasks.append((blood+'__x__'+c+'_'+name,a,[b],adjusted))
        if blood in blocks['b223']:
            for adjusted in [False,True]:tasks.append((blood+'__blood_replication',a,[blocks['b223'][blood]],adjusted))
    with ProcessPoolExecutor(4) as p:results=list(p.map(compare,tasks))
    (O/'tissue_axis_associations.json').write_text(json.dumps(dict(n_common_genes=len(genes),genes=genes,results=results,
        primary_cell_states='existing fine purification maps; PDD astrocytes retain original broad annotation and are a separate phenotype sensitivity',
        raw_input='canonical-symbol count sums followed by log1p(CP10k); donor means; at least three donors per disease arm',
        adjustment='age and sex in MB2025; age, sex and PMI in Kamath; unavailable covariates are not imputed in the other cohorts',
        scope='observed tissue concordance, not a target nomination criterion; includes training donors',
        bootstrap='800 within-condition donor resampling attempts; nuisance adjustment and cohort aggregation repeated; nonidentified resamples excluded and counts retained'),indent=2))
    pd.DataFrame(results).to_csv(O/'tissue_axis_associations.csv',index=False);print('TISSUE AXIS COMPLETE',len(results),flush=True)
if __name__=='__main__':main()
