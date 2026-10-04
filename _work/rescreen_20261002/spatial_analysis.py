"""New-program spatial statistics with the audited barcode-to-image coordinates."""
from pathlib import Path
import json,gzip,re
import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree
from scipy import sparse,stats
from input_utils import canonical_symbols
from randko import bh

A=Path('/public/home/mengxl/dzy/pd_product_assets');O=A/'results/rescreen_20261002';REF=Path(__file__).parent/'reference'

def threshold_sensitivity(Z,S,names,gi,markers,coord,dn,donors):
    from itertools import combinations
    conditions=np.array([coord.loc[dn==d,'condition'].iloc[0]=='PD' for d in donors]);n1=int(conditions.sum());n0=len(donors)-n1
    weights=[]
    for selected in combinations(range(len(donors)),n1):
        mask=np.zeros(len(donors),bool);mask[list(selected)]=True;weights.append(mask/n1-(~mask)/n0)
    weights=np.stack(weights);rows=[];effects=[]
    for state,marker_genes in markers.items():
        pj=[i for i,name in enumerate(names) if name.split('__x__')[1].split('__')[0]==state]
        mi=[gi[g] for g in marker_genes if g in gi]
        if not pj or len(mi)<3:continue
        context=Z[:,mi].mean(1)
        for fraction in [.3,.5,.7,1.]:
            M=[];sizes=[]
            for d in donors:
                ix=np.where(dn==d)[0];threshold=np.quantile(context[ix],1-fraction);keep=ix[context[ix]>=threshold]
                M.append(S[keep][:,pj].mean(0));sizes.append(len(keep))
            M=np.stack(M);effect=M[conditions].mean(0)-M[~conditions].mean(0);null=weights@M
            p=(np.abs(null)>=np.abs(effect)[None,:]-1e-12).mean(0)
            for k,j in enumerate(pj):
                effects.append(dict(combo=names[j],retained_fraction=fraction,effect_PD_minus_control=float(effect[k]),p_exact_two_sided=float(p[k]),n_permutations=len(weights)))
                for di,d in enumerate(donors):rows.append(dict(combo=names[j],retained_fraction=fraction,donor=d,condition='PD' if conditions[di] else 'Control',score=float(M[di,k]),n_spots=sizes[di]))
    df=pd.DataFrame(effects,columns=['combo','retained_fraction','effect_PD_minus_control','p_exact_two_sided','n_permutations']);df['q_all_program_thresholds']=bh(df.p_exact_two_sided)
    df.to_csv(O/'spatial_threshold_sensitivity.csv',index=False)
    pd.DataFrame(rows,columns=['combo','retained_fraction','donor','condition','score','n_spots']).to_csv(O/'spatial_threshold_donor_scores.csv',index=False)
    (O/'spatial_threshold_definition.json').write_text(json.dumps(dict(thresholds=[.3,.5,.7,1.],context='matching brain-state marker score, computed within each donor',scaling='fixed gene mean and SD from all ten full-ROI donor profiles for every threshold',primary='full-ROI competitive validation; no threshold selected by results',score_test='exact two-sided donor-label permutations; BH over all new program-threshold combinations'),indent=2))

def main():
    coord=pd.read_csv(REF/'spatial_coordinates.csv').set_index('spot_id');look=pd.read_csv(REF/'visium_v1_barcode_grid.csv').set_index('barcode')
    checks=[]
    for f in sorted((REF/'alignments').glob('*.gz')):
        gsm=f.name.split('_')[0];j=json.load(gzip.open(f,'rt'));grid=pd.DataFrame(j['oligo']);q=coord[coord.sample_gsm==gsm]
        lookup={(int(r.row),int(r.col)):r for r in grid.itertuples()}
        assert len(q)==int(grid.tissue.sum())
        for sid,r in q.iterrows():
            bc=re.search(r'([ACGT]{16})[.-]1$',sid).group(1);v=look.loc[bc];truth=lookup[int(v['row'])-1,int(v['col'])-1]
            assert truth.tissue and np.isclose(truth.imageX,r.pxl_col_in_fullres) and np.isclose(truth.imageY,r.pxl_row_in_fullres)
        checks.append(dict(sample=gsm,n=len(q),mismatches=0))
    a=ad.read_h5ad(A/'processed/gse253975/v0.1/GSE253975_geomx.h5ad');coord=coord.loc[a.obs_names]
    assert coord.index.is_unique and np.array_equal(coord.donor.astype(str),a.obs.donor.astype(str))
    gs=canonical_symbols(a.var_names);programs=json.loads((O/'programs.json').read_text());names=list(programs)
    markers={'astro':['GFAP','AQP4','SLC1A3','ALDH1L1'],'microglia_mhc2':['HLA-DRA','HLA-DRB1','CD74','HLA-DPA1','HLA-DPB1'],
             'microglia_homeostatic':['P2RY12','CX3CR1','TMEM119','SALL1','P2RY13','GPR34']}
    chosen=sorted((set().union(*(set(v) for v in programs.values()))|set.union(*(set(v) for v in markers.values())))&set(gs))
    gi={g:i for i,g in enumerate(chosen)};cols=[i for i,g in enumerate(gs) if g in gi]
    collapse=sparse.csr_matrix((np.ones(len(cols)),(cols,[gi[gs[i]] for i in cols])),shape=(len(gs),len(chosen)))
    lib=np.asarray(a.X.sum(1)).ravel().clip(1);L=np.log1p((a.X@collapse).multiply((1e6/lib)[:,None]).toarray()).astype(np.float32)
    dn=a.obs.donor.astype(str).to_numpy();donors=sorted(set(dn));M=np.stack([L[dn==d].mean(0) for d in donors]);Z=(L-M.mean(0))/(M.std(0,ddof=1)+1e-8)
    columns=[];valid=[]
    for name in names:
        ix=[gi[g] for g in programs[name] if g in gi]
        if len(ix)>=4:valid.append(name);columns.append(Z[:,ix].mean(1))
    S=np.stack(columns,axis=1) if columns else np.empty((len(a),0));rows=[];rng=np.random.default_rng(8127);nperm=499
    threshold_sensitivity(Z,S,valid,gi,markers,coord,dn,donors)
    for sample in sorted(set(coord.sample_gsm)):
        ix=np.where(coord.sample_gsm.to_numpy()==sample)[0];q=coord.iloc[ix]
        xy=q[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy();dist,neigh=cKDTree(xy).query(xy,k=7)
        spacing=np.median(dist[:,1]);rr=np.repeat(np.arange(len(ix)),6);cc=neigh[:,1:].ravel();keep=dist[:,1:].ravel()<=1.4*spacing
        W=sparse.csr_matrix((np.ones(keep.sum()),(rr[keep],cc[keep])),shape=(len(ix),len(ix)));W=W.maximum(W.T);W.setdiag(0);W.eliminate_zeros()
        z=S[ix]-S[ix].mean(0);den=np.sum(z*z,axis=0);norm=len(ix)/W.sum()
        obs=norm*np.sum(z*(W@z),axis=0)/(den+1e-12);exceed=np.zeros(len(valid),int)
        for _ in range(nperm):
            zz=z[rng.permutation(len(ix))];v=norm*np.sum(zz*(W@zz),axis=0)/(den+1e-12);exceed+=v>=obs
        for j,name in enumerate(valid):
            state=name.split('__x__')[1].split('__')[0];mi=[gi[g] for g in markers[state] if g in gi]
            rho=float(stats.spearmanr(S[ix,j],Z[ix][:,mi].mean(1)).statistic) if len(mi)>=3 else np.nan
            independent_members=[gi[g] for g in programs[name] if g in gi and g not in markers[state]]
            rho_nonoverlap=float(stats.spearmanr(Z[ix][:,independent_members].mean(1),Z[ix][:,mi].mean(1)).statistic) if len(mi)>=3 and len(independent_members)>=4 else np.nan
            rows.append(dict(combo=name,sample=sample,donor=q.donor.iloc[0],condition=q.condition.iloc[0],n_spots=len(ix),moran_I=float(obs[j]),
                             p_spatial_shuffle=float((1+exceed[j])/(1+nperm)),context_marker_rho=rho,context_rho_excluding_shared_genes=rho_nonoverlap,n_nonoverlap_program_genes=len(independent_members),n_context_markers=len(mi)))
        print('spatial',sample,len(ix),'programs',len(valid),flush=True)
    df=pd.DataFrame(rows) if rows else pd.DataFrame(columns=['combo','sample','moran_I','p_spatial_shuffle']);df['q_all_program_sections']=bh(df.p_spatial_shuffle)
    df.to_csv(O/'spatial_localization.csv',index=False)
    np.savez_compressed(O/'spatial_program_scores.npz',scores=S,programs=valid,spot_ids=a.obs_names.to_numpy(dtype=str),
                        xy=coord[['pxl_col_in_fullres','pxl_row_in_fullres']].to_numpy(),donors=dn)
    (O/'spatial_coordinate_checks.json').write_text(json.dumps(dict(checks=checks,n_total=len(coord),status='passed',
        statistic='Moran I on symmetric six-neighbor image graph, reject edges beyond 1.4x median nearest spot distance',
        null='499 within-section spot-label permutations; BH over all program-section tests',
        scores='gene-wise standardized by all-donor mean profiles, average program genes'),indent=2))

if __name__=='__main__':
    from evidence_cache import hit,finish,acquire
    lock=acquire(__file__)
    if not hit(__file__):
        main()
        finish(__file__,['spatial_localization.csv', 'spatial_program_scores.npz', 'spatial_coordinate_checks.json','spatial_threshold_sensitivity.csv','spatial_threshold_donor_scores.csv','spatial_threshold_definition.json'])
    else:print('IDENTICAL NEW-RUN EVIDENCE ALREADY COMPLETE',flush=True)
