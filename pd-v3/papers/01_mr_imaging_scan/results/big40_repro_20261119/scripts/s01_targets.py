"""Step 1: freeze the query list BEFORE any replication value is read. Nothing here touches repro data.
Groups: (1) MR hit pairs under every atlas q definition, labelled; (2) coloc H4>=0.8 regions; (3) all 1375 all-target regions; (4) matched non-hit null plus a random unselected panel."""
import pandas as pd, numpy as np, json, sys, os, time, gzip
from multiprocessing import Pool
ROOT='/public/home/mengxl/dzy/pd-v3/'; MR=ROOT+'papers/01_mr_imaging_scan/'; RES=MR+'results/'; O=RES+'big40_repro_20261119/'; A=RES+'full_atlas_20261116/full/'
RR=ROOT+'papers/03_world_model_integration/results/research_20261118/'; BASE=ROOT+'data/all_target_followup_20261118/regional_all/'; COL=RES+'all_target_coloc_20261118/'
sys.path.insert(0,O+'scripts'); import rfetch as RF
SEED=20261119; KCTRL=6; NRANDOM=1000; t0=time.time()
def log(*a): print(round(time.time()-t0),*a,flush=True)
Q=['q_bh_global','q_by_global','q_storey_lambda_0p25','q_storey_lambda_0p5','q_storey_lambda_0p75','q_bh_cis_family','q_by_cis_family','q_bh_family','q_bh_estimated_m']
PAL={'AT','TA','CG','GC'}
ph=pd.read_csv(A+'phenotypes.csv',low_memory=False); ph=ph.loc[:,~ph.columns.duplicated()][['phenotype','idp','family','is_qc','short_name','n_discovery','n_replication','heritability']].copy(); ph['idp']=ph.idp.astype(str)
ph['has_repro']=ph.phenotype.str.startswith('BIG40:release2:discovery:')&ph.n_replication.notna()
def bh_est(lp):
    o=np.argsort(lp,kind='stable'); m=len(lp); l=lp[o]+np.log(m)-np.log(np.arange(1,m+1)); l=np.minimum.accumulate(l[::-1])[::-1]; q=np.empty(m); q[o]=np.exp(np.minimum(l,0)); return q
U={}
for pol in ['drop','infer']:
    d=pd.read_pickle(O+f'work/est_{pol}.pkl')
    for c in ['exposure_id','phenotype','gene','variant','family','ea_x','oa_x','ea_y','oa_y']: d[c]=d[c].astype(str)
    d['q_bh_estimated_m']=bh_est(d.log_p.values); d['z_disc']=d.beta_y_aligned/d.se_y; d['az']=d.z_disc.abs(); d['key']=d.exposure_id+'|'+d.phenotype
    d=d.merge(ph[['phenotype','is_qc','has_repro','n_discovery','n_replication']],on='phenotype',how='left'); U[pol]=d; log(pol,len(d))
# best non-QC pair per assay and policy
best={}
for pol,d in U.items():
    x=d[d.is_qc!=True]; best[pol]=set(x.loc[x.groupby('exposure_id').log_p.idxmin(),'key'])
log('best pairs',{k:len(v) for k,v in best.items()})
lab={}; qv={}
for pol,d in U.items():
    for c in Q:
        h=d[d[c]<.05]
        for k in h.key.values: lab.setdefault(k,[]).append(f'{pol}:{c}')
ac=pd.read_csv(RR+'mr_imaging/assay_level_calibration.csv.gz'); a2=pd.read_csv(RR+'verify-mr_imaging-F2-r/assay_level_recomputed.csv.gz')
assay_tabs={'assay_q_analyst':set(ac.exposure_id[ac.assay_q<.05]),'assay_q_F2r':set(a2.exposure_id[a2.assay_q<.05])}
d=U['drop']; x=d[(d.is_qc!=True)]; bb=x.loc[x.groupby('exposure_id').log_p.idxmin()]
for name,S in assay_tabs.items():
    for k in bb.key[bb.exposure_id.isin(S)].values: lab.setdefault(k,[]).append(f'drop:{name}_best_pair')
    log(name,'assays q<.05',len(S))
keys=pd.Index(sorted(lab)); log('union hit pairs (all families)',len(keys))
I=U['infer']; H=I[I.key.isin(keys)].copy(); log('in infer rows',len(H)); assert len(H)==len(keys)
H['hit_labels']=H.key.map(lambda k:';'.join(sorted(set(lab[k]))))
for pol,d in U.items():
    m=d[d.key.isin(keys)].set_index('key')
    for c in ['q_bh_global','q_by_global','q_bh_family','q_bh_cis_family','q_by_cis_family','q_bh_estimated_m']: H[f'{c}__{pol}']=H.key.map(m[c])
    H[f'best_pair_of_assay__{pol}']=H.key.isin(best[pol])
H['palindromic']=(H.ea_x+H.oa_x).isin(PAL); H['maf']=np.minimum(H.eaf_y,1-H.eaf_y)
H.drop(columns=[c for c in Q]).to_csv(O+'work/hits_all.csv.gz',index=False)
log('hit pairs by family',H.family.value_counts().to_dict())
# ---- catalog lookup
cat=RF.catalog(); cache={}
def loc(ch,pos,e,o):
    k=(int(ch),int(pos),e,o)
    if k not in cache:
        r,n,a1,a2=RF.find_row(cat,*k)
        if r is None:   # strand-complement representation
            C=str.maketrans('ACGT','TGCA'); r2,n2,b1,b2=RF.find_row(cat,int(ch),int(pos),e.translate(C),o.translate(C))
            cache[k]=(r2,n,b1,b2,'complement_alleles' if r2 is not None else 'not_in_catalog')
        else: cache[k]=(r,n,a1,a2,'exact_alleles')
    return cache[k]
def annotate(t,ech='chrom',ep='pos37',ee='ea',eo='oa'):
    L=[loc(c,p,e,o) for c,p,e,o in zip(t[ech],t[ep],t[ee],t[eo])]
    t['catalog_row']=[x[0] for x in L]; t['n_catalog_rows_at_pos']=[x[1] for x in L]; t['a1']=[x[2] for x in L]; t['a2']=[x[3] for x in L]; t['catalog_match']=[x[4] for x in L]; return t
base=['pair_type','exposure_id','gene','phenotype','family','is_qc','variant','chrom','pos37','ea','oa','eaf','n_disc','n_rep','beta_disc','se_disc','z_disc','palindromic','discovery_source']
def from_atlas(x,typ):
    t=pd.DataFrame(dict(pair_type=typ,exposure_id=x.exposure_id,gene=x.gene,phenotype=x.phenotype,family=x.family,is_qc=x.is_qc,variant=x.variant,chrom=x.chrom.astype(int),pos37=x.pos37.astype(np.int64),ea=x.ea_x,oa=x.oa_x,eaf=x.eaf_y_aligned,n_disc=x.n_y,n_rep=x.n_replication,beta_disc=x.beta_y_aligned,se_disc=x.se_y,z_disc=x.z_disc,palindromic=(x.ea_x+x.oa_x).isin(PAL),discovery_source='atlas_full_family_GCST_harmonised'))
    return t
Th=from_atlas(H,'hit'); Th['hit_labels']=H.hit_labels.values; Th['log_p_mr']=H.log_p.values; Th['F_x']=H.F.values
for c in [c for c in H.columns if '__' in c]: Th[c]=H[c].values
Th=annotate(Th); log('hit catalog',Th.catalog_match.value_counts().to_dict(),'palindromic',int(Th.palindromic.sum()))
# ---- matched null: pool = drop-policy estimated BIG40 pairs that are not a hit under ANY definition/policy. No palindromes in this pool by construction.
P=U['drop']; P=P[P.has_repro&~P.key.isin(keys)].copy(); P['maf']=np.minimum(P.eaf_y,1-P.eaf_y); P=P[P.maf.notna()]; log('null pool',len(P),'max|z|',round(P.az.max(),3))
MB=[0,.005,.01,.05,.1,.2,.35,.5001]; mb=lambda v:np.digitize(v,MB)-1; P['mb']=mb(P.maf.values); P['isx']=(P.chrom>=23)
Tt=H[H.phenotype.isin(ph.phenotype[ph.has_repro])].copy(); Tt['mb']=mb(np.nan_to_num(Tt.maf.values,nan=.25)); Tt['isx']=Tt.chrom>=23
Tt=Tt.sort_values('az',ascending=False); rng=np.random.default_rng(SEED)
Pa=P[P.az>=1.5].sort_values('az').reset_index(drop=True); log('pool |z|>=1.5',len(Pa))
az=Pa.az.values; used=np.zeros(len(Pa),bool); ex=Pa.exposure_id.values; var=Pa.variant.values
idxs={}
for k,g in Pa.groupby(['family','mb','isx']).indices.items(): idxs[('L1',)+k]=g
for k,g in Pa.groupby(['family','isx']).indices.items(): idxs[('L2',)+k]=g
for k,g in Pa.groupby(['isx']).indices.items(): idxs[('L3',k if not isinstance(k,tuple) else k[0])]=g
def pick(arr,t,n,exp,v):
    a=az[arr]; p=np.searchsorted(a,t); lo=p-1; hi=p; out=[]
    while len(out)<n and (lo>=0 or hi<len(arr)):
        dl=abs(a[lo]-t) if lo>=0 else 9e9; dh=abs(a[hi]-t) if hi<len(arr) else 9e9
        if dl<=dh: j=arr[lo]; lo-=1
        else: j=arr[hi]; hi+=1
        if used[j] or ex[j]==exp or var[j]==v: continue
        out.append(j)
    return out
rows=[]; 
for h in Tt.itertuples():
    got=[]; lev=[]
    for L,k in [('L1',('L1',h.family,h.mb,h.isx)),('L2',('L2',h.family,h.isx)),('L3',('L3',h.isx))]:
        if len(got)>=KCTRL or k not in idxs: continue
        g=pick(idxs[k],h.az,KCTRL-len(got),h.exposure_id,h.variant)
        for j in g: used[j]=True
        got+=g; lev+=[L]*len(g)
    for j,l in zip(got,lev): rows.append((j,h.key,l,h.az))
N=pd.DataFrame(rows,columns=['j','matched_hit','match_level','matched_hit_absz']); C=Pa.iloc[N.j.values].reset_index(drop=True)
Tn=from_atlas(C,'null_matched'); Tn['matched_hit']=N.matched_hit.values; Tn['match_level']=N.match_level.values; Tn['matched_hit_absz']=N.matched_hit_absz.values; Tn['z_gap']=C.az.values-N.matched_hit_absz.values
Tn['assay_has_any_hit']=C.exposure_id.isin(set(H.exposure_id)).values; Tn['log_p_mr']=C.log_p.values; Tn['F_x']=C.F.values; Tn=annotate(Tn)
log('null matched',len(Tn),Tn.match_level.value_counts().to_dict(),'|z| quantiles',Tn.z_disc.abs().quantile([0,.5,.9,1]).round(2).to_dict(),'mean z gap',round(Tn.z_gap.mean(),3))
R=P[(P.is_qc!=True)].sample(NRANDOM,random_state=SEED); Tr=from_atlas(R,'null_random_unselected'); Tr['log_p_mr']=R.log_p.values; Tr['F_x']=R.F.values; Tr['assay_has_any_hit']=R.exposure_id.isin(set(H.exposure_id)).values; Tr=annotate(Tr)
# ---- regions
plan=json.load(open(BASE+'PLAN.json')); ip=json.load(open(BASE+'imaging/PLAN.json')); inst=pd.read_csv(A+'selected_cis_instruments.csv',low_memory=False).set_index('exposure_id')
w_of={}
for wid,w in ip['windows'].items():
    for r in w['regions']: w_of.setdefault(r,[]).append(wid)
def region_job(rid):
    r=plan['regions'][rid]; ph_=r['phenotype']; out=dict(region_id=rid,accession=r['accession'],gene=r['gene'],phenotype=ph_,chrom=23 if str(r['chromosome']) in('X','23') else int(r['chromosome']),center=int(r['center_grch37']),model_linked=bool(r['model_linked']),exposure_ids=';'.join(r['exposure_ids']),policies=';'.join(r['policies']),n_exposures=len(r['exposure_ids']))
    fam='BIG40' if ph_.startswith('BIG40') else 'QSM'; out['source_family']=fam
    cs=pd.read_csv(COL+rid+'/coloc_prior_sensitivity.csv'); x=cs[np.isclose(cs.p12_requested,1e-5)]; h={int(a):b for a,b in zip(x.window_bp,x['PP.H4.abf'])}; out['H4_500k_p12_1e5']=h.get(500000,np.nan); out['H4_1M_p12_1e5']=h.get(1000000,np.nan); out['H4_max_p12_1e5']=np.nanmax(list(h.values())) if h else np.nan; out['H4_max_any_prior']=cs['PP.H4.abf'].max(); out['coloc_status']=';'.join(sorted(set(cs.status)))
    W=500000 if (out['H4_500k_p12_1e5']>=out['H4_1M_p12_1e5']) else 1000000; out['lead_window_bp']=W
    for wname,cols in [(W,None)]:
        s=pd.read_csv(COL+rid+f'/window_{W}_snp_posteriors.csv.gz',usecols=['variant_key','pos37','variant_id_y','effect_allele_y','other_allele_y','effect_allele_frequency_y','beta_y','standard_error_y','p_value_y','n_y','beta_x','standard_error_x','effect_allele_x','other_allele_x','palindromic','snp_pp_H4_conditional'])
        if len(s)==0: continue
        a=s.loc[s.snp_pp_H4_conditional.idxmax()]; b=s.loc[s.p_value_y.idxmin()]
        for tag,v in [('coloc_lead',a),('imaging_lead',b)]:
            out.update({f'{tag}_{k}':v[c] for k,c in dict(variant='variant_id_y',pos37='pos37',ea='effect_allele_y',oa='other_allele_y',eaf='effect_allele_frequency_y',beta='beta_y',se='standard_error_y',n='n_y',pal='palindromic',snp_pp='snp_pp_H4_conditional',beta_x='beta_x',se_x='standard_error_x',ea_x='effect_allele_x').items()})
    return out
def center_job(a):
    rid,idp,wid,chrom,center,ea,oa=a
    try:
        f=ROOT+json.load(open(BASE+f'imaging/BIG40/{idp}/COMPLETE.json'))['windows'][wid]['output']
        w=pd.read_csv(f,usecols=['chromosome','base_pair_location','effect_allele','other_allele','beta','standard_error','effect_allele_frequency','n_y','rsid','source_data_row'])
        m=w[(w.base_pair_location==center)&(((w.effect_allele==ea)&(w.other_allele==oa))|((w.effect_allele==oa)&(w.other_allele==ea)))]
        if len(m)==0: return dict(region_id=rid,center_found=False)
        v=m.iloc[0]; s=1 if v.effect_allele==ea else -1
        return dict(region_id=rid,center_found=True,c_beta=s*v.beta,c_se=v.standard_error,c_eaf=v.effect_allele_frequency if s==1 else 1-v.effect_allele_frequency,c_n=v.n_y,c_rsid=v.rsid,c_row=int(v.source_data_row),c_native_a2=v.effect_allele)
    except Exception as e: return dict(region_id=rid,center_found=False,c_err=repr(e)[:120])
if __name__=='__main__':
    with Pool(10) as p: G=pd.DataFrame(p.map(region_job,list(plan['regions']),8))
    log('regions',len(G),G.source_family.value_counts().to_dict(),'H4max>=.8',int((G.H4_max_p12_1e5>=.8).sum()))
    pf=pd.read_csv(RR+'verify-partial_fetch-PF6-a/coloc_hits_H4_ge_0.8_by_route.csv'); s=set(G.region_id[G.H4_max_p12_1e5>=.8]); log('PF6 list',len(pf),'same set',set(pf.region_id)==s)
    G['first_exposure']=G.exposure_ids.str.split(';').str[0]; G['ea']=G.first_exposure.map(inst.ea_x); G['oa']=G.first_exposure.map(inst.oa_x); G['inst_variant']=G.first_exposure.map(inst.variant); G['inst_pos']=G.first_exposure.map(inst.pos37); log('center==instrument pos',int((G.inst_pos==G.center).sum()),'of',len(G))
    G['idp']=G.phenotype.str.split(':').str[-1]; B=G[G.source_family=='BIG40']
    jobs=[(r.region_id,r.idp,[w for w in w_of[r.region_id] if w in ip['phenotypes'][r.phenotype]['windows']][0],r.chrom,r.center,r.ea,r.oa) for r in B.itertuples()]
    with Pool(10) as p: Cc=pd.DataFrame(p.map(center_job,jobs,8))
    G=G.merge(Cc,on='region_id',how='left').merge(ph[['phenotype','family','is_qc','n_discovery','n_replication','has_repro']],on='phenotype',how='left'); log('center found',G.center_found.value_counts(dropna=False).to_dict(),'err',G.get('c_err',pd.Series(dtype=str)).notna().sum())
    G.to_csv(O+'work/regions_table.csv.gz',index=False)
    Tg=[]
    for r in G[(G.source_family=='BIG40')].itertuples():
        pal=lambda a,b:(a+b) in PAL
        def rec(pt,var,pos,ea,oa,eaf,beta,se,n,src):
            return dict(pair_type=pt,region_id=r.region_id,exposure_id=r.first_exposure,gene=r.gene,phenotype=r.phenotype,family=r.family,is_qc=r.is_qc,variant=var,chrom=r.chrom,pos37=int(pos),ea=ea,oa=oa,eaf=eaf,n_disc=n,n_rep=r.n_replication,beta_disc=beta,se_disc=se,z_disc=beta/se if se and se>0 else np.nan,palindromic=pal(ea,oa),discovery_source=src,H4_max_p12_1e5=r.H4_max_p12_1e5,H4_500k_p12_1e5=r.H4_500k_p12_1e5,H4_1M_p12_1e5=r.H4_1M_p12_1e5,model_linked=r.model_linked)
        if r.center_found==True: Tg.append(rec('region_center',r.c_rsid,r.center,r.ea,r.oa,r.c_eaf,r.c_beta,r.c_se,r.c_n,'native_stats_window'))
        else: Tg.append(rec('region_center',r.inst_variant,r.center,r.ea,r.oa,np.nan,np.nan,np.nan,r.n_discovery,'center_not_in_window'))
        for tag in ['coloc_lead','imaging_lead']:
            if pd.notna(getattr(r,f'{tag}_pos37',np.nan)):
                Tg.append(rec(f'region_{tag}',getattr(r,f'{tag}_variant'),getattr(r,f'{tag}_pos37'),getattr(r,f'{tag}_ea'),getattr(r,f'{tag}_oa'),getattr(r,f'{tag}_eaf'),getattr(r,f'{tag}_beta'),getattr(r,f'{tag}_se'),getattr(r,f'{tag}_n'),'coloc_shared_snps_native'))
    Tg=annotate(pd.DataFrame(Tg)); log('region query rows',len(Tg),Tg.pair_type.value_counts().to_dict(),Tg.catalog_match.value_counts().to_dict())
    qs=G[G.source_family!='BIG40'][['region_id','accession','gene','phenotype','H4_max_p12_1e5']]; qs.to_csv(O+'work/regions_not_testable.csv',index=False)
    T=pd.concat([Th,Tn,Tr,Tg],ignore_index=True,sort=False); T['idp']=T.phenotype.str.split(':').str[-1]
    T['testable']=T.phenotype.isin(ph.phenotype[ph.has_repro]); T['pair_id']=np.arange(len(T)); T.to_pickle(O+'work/query_table.pkl',protocol=4)
    log('query table',len(T),T.groupby('pair_type').testable.agg(['size','sum']).to_dict('index'))
    pri={'region_coloc_lead':0,'region_center':0,'region_imaging_lead':0,'hit':1,'null_matched':2,'null_random_unselected':3}
    Q_=T[T.testable&T.catalog_row.notna()].copy(); Q_['catalog_row']=Q_.catalog_row.astype(np.int64); Q_['prio']=Q_.pair_type.map(pri)
    Fx=Q_.groupby(['idp','catalog_row']).agg(prio=('prio','min'),n_rep=('n_rep','first'),chrom=('chrom','first'),pos37=('pos37','first'),a1=('a1','first'),a2=('a2','first')).reset_index(); Fx['release']='repro'
    Dp=Q_[Q_.palindromic&(Q_.pair_type=='hit')]; Fd=Dp.groupby(['idp','catalog_row']).agg(prio=('prio','min'),n_rep=('n_disc','first'),chrom=('chrom','first'),pos37=('pos37','first'),a1=('a1','first'),a2=('a2','first')).reset_index(); Fd['release']='stats'; Fd['prio']=1
    F_=pd.concat([Fx,Fd],ignore_index=True); F_.to_csv(O+'work/fetch_targets.csv.gz',index=False)
    log('fetch targets',len(F_),F_.groupby(['release','prio']).size().to_dict(),'distinct repro files',F_[F_.release=='repro'].idp.nunique())
    log('DONE')
