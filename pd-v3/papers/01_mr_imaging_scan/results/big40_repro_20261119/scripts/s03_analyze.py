"""Step 3: merge the range-read values with the frozen query table and compute every requested statistic. No value is dropped silently; every queried pair stays in repro_pairs.csv.gz."""
import pandas as pd, numpy as np, json, glob, os, sys, time, re
from scipy import stats, optimize
from scipy.interpolate import BSpline
ROOT='/public/home/mengxl/dzy/pd-v3/'; MR=ROOT+'papers/01_mr_imaging_scan/'; RES=MR+'results/'; O=RES+'big40_repro_20261119/'; W=O+'work/'; A=RES+'full_atlas_20261116/full/'
RR=ROOT+'papers/03_world_model_integration/results/research_20261118/'; SUB=sys.argv[1] if len(sys.argv)>1 else 'fetch'; FINAL=len(sys.argv)>2 and sys.argv[2]=='final'
t0=time.time(); rng=np.random.default_rng(20261119); NB=4000; NF=20000
def log(*a): print(round(time.time()-t0),*a,flush=True)
pd.set_option('display.width',250,'display.max_columns',60,'display.max_rows',300)
# ------------------------------------------------------------------ merge
T=pd.read_pickle(W+'query_table.pkl'); T['catalog_row']=T.catalog_row.astype('float'); T['idp']=T.idp.astype(str)
rec=[]
for f in glob.glob(W+SUB+'/part_*.jsonl'):
    for l in open(f):
        try: rec.append(json.loads(l))
        except Exception: pass
Fx=pd.DataFrame(rec); Fx['rank']=(Fx.status!='ok').astype(int)+0.5*(Fx.digit_status=='interval_wide'); Fx=Fx.sort_values('rank',kind='stable').drop_duplicates(['release','idp','row']); log('fetch rows',len(Fx),Fx.groupby('release').status.value_counts().to_dict())
Rp=Fx[Fx.release=='repro'].rename(columns={'row':'catalog_row','status':'fetch_status','beta':'beta_rep_a2','se':'se_rep','mlogp':'mlogp_rep','z':'z_rep_a2','z_lo':'z_rep_a2_lo','z_hi':'z_rep_a2_hi','z_width':'z_rep_width','n_unknown':'n_inferred_chars','n_feasible_values':'n_feasible','tries':'read_tries','cluster_bytes':'fetch_cluster_bytes','cluster_size':'fetch_cluster_rows'})
keep=['idp','catalog_row','fetch_status','beta_rep_a2','se_rep','mlogp_rep','z_rep_a2','z_rep_a2_lo','z_rep_a2_hi','z_rep_width','digit_status','n_inferred_chars','n_feasible','read_tries','fetch_cluster_bytes','fetch_cluster_rows','note']
T=T.merge(Rp[[c for c in keep if c in Rp.columns]],on=['idp','catalog_row'],how='left')
Dn=Fx[(Fx.release=='stats')&(Fx.status=='ok')][['idp','row','beta','se','mlogp']].rename(columns={'row':'catalog_row','beta':'beta_native_disc_a2','se':'se_native_disc','mlogp':'mlogp_native_disc'}); T=T.merge(Dn,on=['idp','catalog_row'],how='left')
T['status']=np.select([~T.testable,T.catalog_row.isna(),T.fetch_status=='ok',T.fetch_status.isna()],['not_testable_no_repro_file','variant_not_in_catalog','ok','not_fetched'],default=T.fetch_status)
# allele alignment, per variant
sg=np.where((T.a2==T.ea)&(T.a1==T.oa),1,np.where((T.a1==T.ea)&(T.a2==T.oa),-1,0)); T['sign_a2_to_ea']=sg
T['allele_check']=np.select([T.catalog_match=='not_in_catalog',T.catalog_match=='complement_alleles',sg==1,sg==-1],['variant_not_in_catalog','complement_alleles_flagged','effect_allele_is_catalog_a2','effect_allele_is_catalog_a1'],default='unresolved')
T['palindromic_flag']=T.palindromic.astype(bool)
T['z_disc_ea']=T.z_disc; T['discovery_orientation']='catalog_alleles_match'
nat=T.beta_native_disc_a2.notna(); T.loc[nat,'z_disc_native_ea']=(T.sign_a2_to_ea*T.beta_native_disc_a2/T.se_native_disc)[nat]
pal=T.palindromic_flag&(T.pair_type=='hit')
T.loc[pal&nat&(np.sign(T.z_disc_native_ea)==np.sign(T.z_disc_ea)),'discovery_orientation']='palindromic_harmonised_equals_native'
T.loc[pal&nat&(np.sign(T.z_disc_native_ea)!=np.sign(T.z_disc_ea)),'discovery_orientation']='palindromic_harmonised_opposite_to_native'
T.loc[pal&~nat&T.testable,'discovery_orientation']='palindromic_native_discovery_missing'; T.loc[pal&~nat&~T.testable,'discovery_orientation']='palindromic_not_testable_no_repro_file'
T.loc[T.palindromic_flag&T.pair_type.str.startswith('region'),'discovery_orientation']='palindromic_native_window_used'
T['z_disc_used']=np.where(pal&nat,T.z_disc_native_ea,T.z_disc_ea); T['beta_disc_used']=np.where(pal&nat,T.sign_a2_to_ea*T.beta_native_disc_a2,T.beta_disc); T['se_disc_used']=np.where(pal&nat,T.se_native_disc,T.se_disc)
ok=(T.status=='ok')&(T.sign_a2_to_ea!=0)&T.z_disc_used.notna()&~((T.pair_type=='hit')&T.palindromic_flag&~nat)&(T.digit_status!='interval_wide')
T['analysable']=ok; T['beta_rep_ea']=T.sign_a2_to_ea*T.beta_rep_a2; T['z_rep_ea']=T.sign_a2_to_ea*T.z_rep_a2; T['absz_disc']=T.z_disc_used.abs(); s=np.sign(T.z_disc_used)
T['z_rep_signed']=T.z_rep_ea*s; T['z_rep_signed_lo']=np.minimum(T.z_rep_a2_lo*T.sign_a2_to_ea*s,T.z_rep_a2_hi*T.sign_a2_to_ea*s); T['z_rep_signed_hi']=np.maximum(T.z_rep_a2_lo*T.sign_a2_to_ea*s,T.z_rep_a2_hi*T.sign_a2_to_ea*s)
T['same_sign']=np.where(T.analysable,(T.z_rep_signed>0),np.nan); T['one_sided_p']=stats.norm.sf(T.z_rep_signed); T['z_rep_signed_gt_1p645']=np.where(T.analysable,T.z_rep_signed>stats.norm.isf(.05),np.nan)
T['is_x']=T.chrom>=23; scale=np.sqrt(T.n_rep/T.n_disc); scale=scale.where(~T.is_x); T['scale_n']=scale; T['scale_se']=T.se_disc_used/T.se_rep; T['scale']=T.scale_n.fillna(T.scale_se)
T['expected_z_n']=T.absz_disc*T.scale_n; T['expected_z_se']=T.beta_disc_used.abs()/T.se_rep; T['expected_z']=T.expected_z_n.fillna(T.expected_z_se)
T['maf']=np.minimum(T.eaf,1-T.eaf)
prev=pd.read_csv(RR+'verify-mr_imaging-F2-r/replication_raw.csv',low_memory=False)[['exposure_id','phs']].rename(columns={'phs':'phenotype'}); p2=pd.concat([pd.read_csv(RR+f'mr_imaging/{f}') for f in ['repro_results.csv','repro_results_retry.csv']])[['exposure_id','phenotype']]
pk=set(zip(prev.exposure_id,prev.phenotype))|set(zip(p2.exposure_id,p2.phenotype)); T['read_in_earlier_F2r_or_s07']=[(a,b) in pk for a,b in zip(T.exposure_id,T.phenotype)]
# locus clusters: variants within 1 Mb on a chromosome are chained into one locus
v=T.drop_duplicates(['chrom','pos37'])[['chrom','pos37']].sort_values(['chrom','pos37']); nl=(v.chrom.diff().ne(0)|v.pos37.diff().gt(1_000_000)).cumsum(); v['locus']=nl.values; T=T.merge(v,on=['chrom','pos37'],how='left')
# ------------------------------------------------------------------ winner's curse: Tweedie empirical Bayes on the full scan's z distribution (discovery only)
Dp=pd.read_pickle(W+'est_drop.pkl')[['exposure_id','phenotype','family','beta_y_aligned','se_y','variant','q_bh_global']]; Dp['phenotype']=Dp.phenotype.astype(str); Dp['family']=Dp.family.astype(str); Dp=Dp[Dp.phenotype.str.startswith('BIG40')&(Dp.family!='QC')]; Dp['z']=(Dp.beta_y_aligned/Dp.se_y).astype(float); log('scan z rows',len(Dp))
knots=np.array([0,.5,1,1.5,2,2.5,3,3.5,4,4.5,5,6,8,10.]); k3=np.r_[[0]*3,knots,[10.]*3]
def bs(x,der=0):
    n=len(k3)-4; B=np.empty((len(x),n))
    for i in range(n):
        c=np.zeros(n); c[i]=1; sp=BSpline(k3,c,3,extrapolate=True); B[:,i]=sp.derivative(der)(x) if der else sp(x)
    return B
def fit_density(az):
    az=np.abs(az); e=np.arange(0,10.0001,.05); h,_=np.histogram(az[az<10],e); m=(e[:-1]+e[1:])/2; X=bs(m); beta=np.zeros(X.shape[1]); beta[:]=np.log(h.mean()+1)
    for it in range(60):
        eta=np.clip(X@beta,-30,30); mu=np.exp(eta); W_=mu; z_=eta+(h-mu)/mu; XtW=X.T*W_; nb=np.linalg.solve(XtW@X+1e-6*np.eye(X.shape[1]),XtW@z_)
        if np.abs(nb-beta).max()<1e-7: beta=nb; break
        beta=nb
    return beta
GRID=np.arange(0,12.0001,.01)
def curve(beta,zs):
    """posterior mean of |true z| given |observed z| on a grid. Beyond the point with only 200 scan values left the fitted density is not trusted: no further shrinkage is added there (slope 1). Forced monotone and within [0,|z|]."""
    az=np.sort(np.abs(zs)); t=min(az[-200] if len(az)>200 else 3.0,9.99); x=np.minimum(GRID,9.99); raw=GRID+bs(x,1)@beta; i=int(np.searchsorted(GRID,t)); m=raw.copy(); m[i:]=raw[i]+(GRID[i:]-GRID[i]); m=np.maximum.accumulate(m); return np.clip(m,0,GRID),t
def tweedie(cv,az):
    az=np.abs(np.asarray(az,float)); return np.where(az>12,az-(12-cv[-1]),np.interp(az,GRID,cv))
fams=Dp.family.value_counts(); bg,tg=curve(fit_density(Dp.z.values),Dp.z.values); bf={}; tails={'global':float(tg)}
for f in fams.index:
    if fams[f]>=100000:
        zz=Dp.z.values[Dp.family.values==f]; bf[f],tails[f]=curve(fit_density(zz),zz); tails[f]=float(tails[f])
log('tweedie fits',len(bf)+1,'tail points',{k:round(v,2) for k,v in tails.items()})
tw=np.full(len(T),np.nan); a=T.absz_disc.values; tw_g=tweedie(bg,np.nan_to_num(a)); tw=tw_g.copy()
for f,b in bf.items():
    m=(T.family.values==f); tw[m]=tweedie(b,np.nan_to_num(a[m]))
tw[np.isnan(a)]=np.nan; tw_g[np.isnan(a)]=np.nan; T['mu_tweedie_family']=tw; T['mu_tweedie_global']=tw_g
T['expected_z_tweedie']=T.mu_tweedie_family*T.scale; T['expected_z_tweedie_global']=T.mu_tweedie_global*T.scale
gr=np.arange(0,10.01,.25); pd.DataFrame(dict(absz=gr,mu_global=tweedie(bg,gr),**{('mu_'+f.replace(' ','_')):tweedie(b,gr) for f,b in bf.items()})).round(4).to_csv(O+'winners_curse_tweedie_curve.csv',index=False)
# conditional MLE: truncation at the definition's own |z| cutoff
H=T[T.pair_type=='hit']; lab_t={}
for k,g in H[H.hit_labels.notna()].assign(l=H.hit_labels.str.split(';')).explode('l').groupby('l'):
    if 'assay_q' in k: continue
    if k.split(':')[1]=='q_bh_family': lab_t[k]=g.groupby('family').absz_disc.min().to_dict()
    else: lab_t[k]=float(g.absz_disc.min())
def cut(lbls,fam):
    v=[]
    for l in str(lbls).split(';'):
        x=lab_t.get(l)
        if x is None: continue
        v.append(x.get(fam,np.nan) if isinstance(x,dict) else x)
    v=[y for y in v if y==y]; return min(v) if v else np.nan
T['wc_cut_any']=[cut(l,f) if t=='hit' else np.nan for l,f,t in zip(T.hit_labels,T.family,T.pair_type)]
def cmle(z,c):
    z=np.asarray(z,float); c=np.asarray(c,float); out=np.full(len(z),np.nan)
    for i,(x,k) in enumerate(zip(z,c)):
        if not(x==x and k==k): continue
        f=lambda m:m+np.exp(stats.norm.logpdf(k-m)-stats.norm.logsf(k-m))-x
        if f(0.0)>=0: out[i]=0.0; continue
        out[i]=optimize.brentq(f,0.0,x+1e-9) if f(x)>0 else x
    return out
# ------------------------------------------------------------------ phenotype correlation for pairs sharing a variant
Zm=Dp[Dp.exposure_id.isin(Dp.groupby('exposure_id').z.apply(lambda s:s.abs().max()).pipe(lambda m:m[m<4.5]).index)]; ex=Zm.exposure_id.astype('category'); pc=Zm.phenotype.astype('category'); Zt=np.zeros((len(ex.cat.categories),len(pc.cat.categories)),'float32'); Zt[ex.cat.codes,pc.cat.codes]=Zm.z.values.astype('float32'); Zt=(Zt-Zt.mean(0))/(Zt.std(0)+1e-9); Rm=(Zt.T@Zt)/Zt.shape[0]; pidx={p:i for i,p in enumerate(pc.cat.categories)}; log('phenotype z-correlation matrix',Rm.shape,'from clean assays',Zt.shape[0]); del Zm
def var_inflation(g):
    tot=float(len(g)); shared=0
    for v_,x in g.groupby('variant',sort=False):
        if len(x)<2: continue
        idx=[pidx.get(p,-1) for p in x.phenotype]; keep=[i for i,j in enumerate(idx) if j>=0]
        if len(keep)<2: continue
        I=np.array([idx[i] for i in keep]); s_=np.sign(x.z_disc_used.values[keep]).astype('float32'); M=Rm[np.ix_(I,I)]; np.fill_diagonal(M,0.0) if False else None; val=float(s_@M@s_-(s_**2*np.diag(M)).sum()); tot+=val; shared+=len(keep)
    return tot,shared
# ------------------------------------------------------------------ group statistics
def poisson_ratio(N,D,cl,B=NB):
    u,inv=np.unique(cl,return_inverse=True); K=len(u); Sn=np.bincount(inv,N,K); Sd=np.bincount(inv,D,K); w=rng.poisson(1,(B,K)).astype('float32'); a=w@Sn; b=w@Sd
    with np.errstate(all='ignore'): r=a/b
    return np.nanquantile(r,[.025,.975])
def flips(sgnA,B=NF):
    return None
def flip_test(x,cl,pos=None,B=NF):
    """sign flip whole clusters: returns p for (mean x), for sign concordance and for one-sided significance."""
    u,inv=np.unique(cl,return_inverse=True); K=len(u); X=np.bincount(inv,x,K); Sg=np.bincount(inv,np.sign(x),K); Pp=np.bincount(inv,x>1.6448536,K); Pn=np.bincount(inv,x<-1.6448536,K)
    obs=np.array([X.sum(),Sg.sum(),Pp.sum()]); c=np.zeros(3); done=0
    while done<B:
        n=min(2000,B-done); f=rng.integers(0,2,(n,K),dtype=np.int8)*2-1; c[0]+=(f@X>=obs[0]-1e-9).sum(); c[1]+=(f@Sg>=obs[1]-1e-9).sum(); c[2]+=(((f>0)@Pp+(f<0)@Pn)>=obs[2]-1e-9).sum(); done+=n
    return (1+c)/(B+1)
def cvar(x,cl):
    x=np.asarray(x,float); n=len(x); m=x.mean(); e=x-m; u,inv=np.unique(cl,return_inverse=True); s=np.bincount(inv,e,len(u)); return float((s**2).sum()/n**2)
def icc(x,cl):
    x=np.asarray(x,float); u,inv=np.unique(cl,return_inverse=True); k=len(u); n=len(x); ni=np.bincount(inv,minlength=k); 
    if k<2 or n<=k: return np.nan
    gm=x.mean(); mi=np.bincount(inv,x,k)/ni; msb=(ni*(mi-gm)**2).sum()/(k-1); msw=((x-mi[inv])**2).sum()/(n-k); n0=(n-(ni**2).sum()/n)/(k-1)
    return float((msb-msw)/(msb+(n0-1)*msw)) if (msb+(n0-1)*msw)>0 else np.nan
def S(g_all,group,scope='overall',tier='A*',desc='',wc=None):
    g=g_all[g_all.analysable]; n=len(g); r=dict(group=group,scope=scope,tier=tier,description=desc,n_queried=len(g_all),n_analysable=n,n_failed_or_missing=int((g_all.status.isin(['failed','error_head','not_fetched','cap_reached','no_remote_file'])).sum()),n_not_testable=int((g_all.status=='not_testable_no_repro_file').sum()),n_palindromic=int(g.palindromic_flag.sum()),n_inferred_interval_rows=int((g.z_rep_width>0).sum()))
    if n<3: return r
    x=g.z_rep_signed.values; e=g.expected_z.values; es=g.expected_z_se.values; A_=g.exposure_id.astype(str).values; V=g.variant.astype(str).values; L=g.locus.values; P=g.phenotype.astype(str).values
    r.update(n_assays=len(set(A_)),n_variants=len(set(V)),n_loci=len(set(L)),n_phenotypes=len(set(P)),median_absz_discovery=float(g.absz_disc.median()),mean_absz_discovery=float(g.absz_disc.mean()))
    c=int((x>0).sum()); o=int((x>stats.norm.isf(.05)).sum()); r.update(same_sign_n=c,same_sign_frac=c/n,same_sign_binom_p_naive=float(stats.binomtest(c,n,.5,alternative='greater').pvalue),one_sided_p05_n=o,one_sided_p05_frac=o/n,one_sided_binom_p_naive=float(stats.binomtest(o,n,.05,alternative='greater').pvalue),expected_chance_one_sided=.05*n)
    fp=flip_test(x,L); r.update(mean_z_flip_p_by_locus=float(fp[0]),same_sign_flip_p_by_locus=float(fp[1]),one_sided_flip_p_by_locus=float(fp[2]))
    r.update(mean_z_rep_signed=float(x.mean()),expected_z_n_scaled=float(np.nanmean(g.expected_z_n)),expected_z=float(np.mean(e)),expected_z_se_based=float(np.mean(es)),ratio=float(x.mean()/np.mean(e)),ratio_se_based=float(x.mean()/np.mean(es)))
    lo,hi=poisson_ratio(x,e,A_); r.update(ratio_ci_lo_assay=float(lo),ratio_ci_hi_assay=float(hi)); lo,hi=poisson_ratio(x,e,L); r.update(ratio_ci_lo_locus=float(lo),ratio_ci_hi_locus=float(hi))
    et=g.expected_z_tweedie.values; lo,hi=poisson_ratio(x,et,A_); r.update(expected_z_after_wc_tweedie=float(np.mean(et)),ratio_after_wc_tweedie=float(x.mean()/np.mean(et)),ratio_after_wc_tweedie_lo=float(lo),ratio_after_wc_tweedie_hi=float(hi),expected_same_sign_frac_after_wc=float(stats.norm.cdf(et).mean()),expected_one_sided_frac_after_wc=float(stats.norm.sf(1.6448536-et).mean()),expected_one_sided_frac_before_wc=float(stats.norm.sf(1.6448536-e).mean()),expected_same_sign_frac_before_wc=float(stats.norm.cdf(e).mean()),ratio_after_wc_tweedie_global_prior=float(x.mean()/g.expected_z_tweedie_global.mean()))
    if wc is not None:
        cc=wc(g); m=cmle(g.absz_disc.values,cc); fin=~np.isnan(m); 
        if fin.sum()>=3:
            ec=m[fin]*g.scale.values[fin]; lo,hi=poisson_ratio(x[fin],ec,A_[fin]); r.update(n_cmle=int(fin.sum()),expected_z_after_wc_cmle=float(ec.mean()),ratio_after_wc_cmle=float(x[fin].mean()/ec.mean()),ratio_after_wc_cmle_lo=float(lo),ratio_after_wc_cmle_hi=float(hi),expected_one_sided_frac_after_cmle=float(stats.norm.sf(1.6448536-ec).mean()))
    bd=g.beta_disc_used.abs().values; br=(g.beta_rep_ea*np.sign(g.z_disc_used)).values; w=1/g.se_rep.values**2; r.update(ivw_slope_rep_on_disc=float((bd*br*w).sum()/(bd*bd*w).sum()),ivw_slope_se=float(1/np.sqrt((bd*bd*w).sum())))
    s2=x.var(ddof=1); iid=s2/n; va=cvar(x,A_); vv=cvar(x,V); vl=cvar(x,L); vp=cvar(x,P); vb=cvar(x,np.char.add(np.char.add(V.astype(str),'|'),P.astype(str))); two=max(vv+vp-vb,1e-12); two_l=max(vl+vp-cvar(x,np.char.add(L.astype(str),np.char.add('|',P.astype(str)))),1e-12)
    r.update(sd_z_rep_signed=float(np.sqrt(s2)),se_mean_iid=float(np.sqrt(iid)),se_mean_cluster_assay=float(np.sqrt(va)),se_mean_cluster_variant=float(np.sqrt(vv)),se_mean_cluster_locus=float(np.sqrt(vl)),se_mean_cluster_phenotype=float(np.sqrt(vp)),se_mean_twoway_locus_phenotype=float(np.sqrt(two_l)),design_effect_assay=va/iid,design_effect_variant=vv/iid,design_effect_locus=vl/iid,design_effect_phenotype=vp/iid,design_effect_twoway_locus_phenotype=two_l/iid,n_eff_twoway=float(s2/two_l),icc_variant=icc(x,V),icc_locus=icc(x,L),icc_phenotype=icc(x,P),icc_assay=icc(x,A_))
    tv,sh=var_inflation(g); r.update(null_variance_inflation_shared_variant_phenotype_corr=tv/n,pairs_sharing_variant=sh,stouffer_z_naive=float(x.sum()/np.sqrt(n)),stouffer_z_adj_two_way=float(x.mean()/np.sqrt(two_l)),stouffer_p_adj_two_way=float(stats.norm.sf(x.mean()/np.sqrt(two_l))))
    return r
Q_LAB=['q_bh_global','q_by_global','q_storey_lambda_0p25','q_storey_lambda_0p5','q_storey_lambda_0p75','q_bh_cis_family','q_by_cis_family','q_bh_family','q_bh_estimated_m']
H=T[T.pair_type=='hit'].copy(); H['l']=H.hit_labels.fillna('').map(lambda s:set(s.split(';')))
NQ=lambda d:d[d.is_qc!=True]
def lab_mask(lb): return H.l.map(lambda s:lb in s)
def wc_for(lb):
    x=lab_t.get(lb)
    if x is None: return None
    return (lambda g:g.family.map(x).values) if isinstance(x,dict) else (lambda g:np.full(len(g),x))
rows=[]
rows.append(S(NQ(H),'hit_ANY_definition',desc='Union of all pair-level q<0.05 definitions and both palindrome policies, non-QC BIG40',wc=lambda g:g.wc_cut_any.values))
rows.append(S(NQ(H[~H.read_in_earlier_F2r_or_s07]),'hit_ANY_definition__not_read_before',tier='A',desc='as above, excluding pairs whose repro value was already read in the earlier checks',wc=lambda g:g.wc_cut_any.values))
rows.append(S(H[H.is_qc==True],'hit_ANY_definition__QC_phenotypes_only',desc='Hit pairs on the 16 QC phenotypes'))
rows.append(S(H[H.palindromic_flag],'hit_ANY_definition__palindromic_only',desc='Palindromic variants (only reachable under the infer policy); discovery sign taken from the native BIG40 row'))
for pol in ['drop','infer']:
    for q in Q_LAB:
        lb=f'{pol}:{q}'; m=lab_mask(lb); rows.append(S(NQ(H[m]),f'hit__{lb}',desc=f'pair-level {q}<0.05, policy {pol}, non-QC BIG40',wc=wc_for(lb)))
for nm in ['drop:assay_q_analyst_best_pair','drop:assay_q_F2r_best_pair']:
    rows.append(S(NQ(H[lab_mask(nm)]),f'hit__{nm}',desc='best non-QC pair of assays with assay-level max|z| q<0.05'))
B=H[H.best_pair_of_assay__drop==True]; hit133=lab_mask('drop:q_bh_global'); rows.append(S(NQ(B[lab_mask('drop:q_bh_global').reindex(B.index,fill_value=False)&(B.chrom<=22)]),'F2r_reproduction_check__best_pair_of_q_bh_global_hit_assays_autosomal',desc='same set as the earlier 126-pair check (its ratio was 0.484)',wc=wc_for('drop:q_bh_global')))
G=T[T.pair_type.str.startswith('region')]; 
for pt,nm in [('region_center','instrument'),('region_coloc_lead','coloc_lead_snp'),('region_imaging_lead','imaging_lead_snp')]:
    g=G[G.pair_type==pt]
    rows.append(S(g,f'region_all__{nm}',desc='All BIG40 all-target regions (best pair per assay and policy)'))
    rows.append(S(g[g.H4_max_p12_1e5>=.8],f'coloc_H4ge0.8__{nm}',desc='max of 500kb and 1Mb PP.H4 at p12=1e-5 >=0.8'))
    rows.append(S(g[g.H4_1M_p12_1e5>=.8],f'coloc_H4ge0.8_1Mb__{nm}',desc='1Mb window H4>=0.8 at p12=1e-5')); rows.append(S(g[g.H4_500k_p12_1e5>=.8],f'coloc_H4ge0.8_500kb__{nm}',desc='500kb window H4>=0.8 at p12=1e-5'))
    rows.append(S(g[g.H4_max_p12_1e5<.8],f'coloc_H4lt0.8__{nm}',desc='regions with max H4<0.8'))
    rows.append(S(g[g.model_linked==True],f'region_model_linked__{nm}',desc='model-linked regions'))
N=T[T.pair_type=='null_matched']; Rn=T[T.pair_type=='null_random_unselected']
rows.append(S(N,'null_matched',desc='non-hit pairs (no definition, no policy), matched to hits on phenotype class and MAF bin, nearest available |z|')); rows.append(S(N[N.assay_has_any_hit==False],'null_matched__assays_without_any_hit',desc='matched null restricted to assays with no hit pair')); rows.append(S(Rn,'null_random_unselected',desc='1000 random non-hit non-QC pairs, no matching on |z|'))
for lo,hi in [(1.5,2.5),(2.5,3),(3,3.5),(3.5,3.85),(3.85,4.0),(4.0,4.2),(4.2,4.4)]:
    rows.append(S(N[(N.absz_disc>=lo)&(N.absz_disc<hi)],f'null_matched__absz_{lo}_{hi}',scope=f'absz_bin',desc='matched null by discovery |z|'))
for lo,hi in [(3.85,4.0),(4.0,4.2),(4.2,4.4),(4.4,4.6),(4.6,5),(5,6),(6,8),(8,100)]:
    rows.append(S(NQ(H[(H.absz_disc>=lo)&(H.absz_disc<hi)]),f'hit_ANY__absz_{lo}_{hi}',scope='absz_bin',desc='hit pairs by discovery |z|',wc=lambda g:g.wc_cut_any.values))
for nm,g0,wc in [('hit_ANY_definition',NQ(H),lambda g:g.wc_cut_any.values),('hit__drop:q_bh_global',NQ(H[lab_mask('drop:q_bh_global')]),wc_for('drop:q_bh_global')),('coloc_H4ge0.8__coloc_lead_snp',G[(G.pair_type=='region_coloc_lead')&(G.H4_max_p12_1e5>=.8)],None),('region_all__instrument',G[G.pair_type=='region_center'],None),('null_matched',N,None)]:
    for f,g in g0.groupby('family'):
        if len(g)>=10: rows.append(S(g,nm,scope=f'phenotype_class={f}',desc='per-phenotype-class breakdown',wc=wc))
topl=NQ(H)[NQ(H).analysable].groupby('locus').size().sort_values(ascending=False); TOP=int(topl.index[0]); TOPV=NQ(H)[(NQ(H).locus==TOP)].variant.value_counts().index[0]; log('top hit locus',TOP,TOPV,int(topl.iloc[0]))
def one_per(d,by): return d[d.analysable].sort_values('absz_disc',ascending=False).drop_duplicates(by)
hq=lab_mask('drop:q_bh_global')
rows.append(S(one_per(NQ(H),'locus'),'hit_ANY_definition__one_pair_per_locus',desc='strongest discovery pair at each 1 Mb-chained locus; near-independent subset',wc=lambda g:g.wc_cut_any.values))
rows.append(S(one_per(NQ(H),'exposure_id'),'hit_ANY_definition__one_pair_per_assay',desc='strongest discovery pair of each assay',wc=lambda g:g.wc_cut_any.values))
rows.append(S(one_per(NQ(H[hq]),'locus'),'hit__drop:q_bh_global__one_pair_per_locus',desc='strongest pair per locus among drop-policy BH-global hits',wc=wc_for('drop:q_bh_global')))
rows.append(S(NQ(H[H.locus!=TOP]),'hit_ANY_definition__excluding_top_locus',desc=f'without the locus holding most hit pairs (lead variant {TOPV})',wc=lambda g:g.wc_cut_any.values)); rows.append(S(NQ(H[hq&(H.locus!=TOP)]),'hit__drop:q_bh_global__excluding_top_locus',desc=f'without the locus holding most hit pairs (lead variant {TOPV})',wc=wc_for('drop:q_bh_global')))
rows.append(S(N[N.locus!=TOP],'null_matched__excluding_top_locus',desc='matched null without the same locus'))
rows.append(S(one_per(N,'locus'),'null_matched__one_pair_per_locus',desc='strongest matched-null pair per locus'))
for nm,g0,wc in [('hit__drop:q_bh_global',NQ(H[hq]),wc_for('drop:q_bh_global')),('coloc_H4ge0.8__coloc_lead_snp',G[(G.pair_type=='region_coloc_lead')&(G.H4_max_p12_1e5>=.8)],None),('region_all__instrument',G[G.pair_type=='region_center'],None),('null_matched',N,None),('null_random_unselected',Rn,None)]:
    rows.append(S(g0[~g0.read_in_earlier_F2r_or_s07],nm+'__not_read_before',tier='A',desc='excluding pairs whose repro value was read in the earlier F2-r / s07 checks',wc=wc))
Y=pd.DataFrame(rows); Y.to_csv(O+'repro_summary.csv',index=False); log('summary rows',len(Y))
dc=['group','scope','n_analysable','n_assays','n_variants','n_loci','n_phenotypes','pairs_sharing_variant','icc_assay','icc_variant','icc_locus','icc_phenotype','design_effect_assay','design_effect_variant','design_effect_locus','design_effect_phenotype','design_effect_twoway_locus_phenotype','n_eff_twoway','null_variance_inflation_shared_variant_phenotype_corr','se_mean_iid','se_mean_cluster_assay','se_mean_cluster_locus','se_mean_cluster_phenotype','se_mean_twoway_locus_phenotype','stouffer_z_naive','stouffer_z_adj_two_way']
Y[Y.scope=='overall'][[c for c in dc if c in Y.columns]].to_csv(O+'dependence_summary.csv',index=False)
sh=[]
for nm,g in [('hit_ANY_definition',NQ(H)),('hit__drop:q_bh_global',NQ(H[hq])),('coloc_H4ge0.8__coloc_lead_snp',G[(G.pair_type=='region_coloc_lead')&(G.H4_max_p12_1e5>=.8)]),('region_all__instrument',G[G.pair_type=='region_center']),('null_matched',N)]:
    g=g[g.analysable]; a=g.groupby('variant').size(); b=g.groupby('phenotype').size(); sh.append(dict(group=nm,n_pairs=len(g),n_variants=len(a),max_pairs_one_variant=int(a.max()),top_variant=str(a.idxmax()),share_of_pairs_in_top5_variants=float(a.sort_values(ascending=False).head(5).sum()/len(g)),pairs_with_variant_shared_by_2plus=int(a[a>1].sum()),n_phenotypes=len(b),max_pairs_one_phenotype=int(b.max()),pairs_with_phenotype_shared_by_2plus=int(b[b>1].sum())))
pd.DataFrame(sh).to_csv(O+'pair_sharing_counts.csv',index=False)
nt=T[~T.testable][['pair_type','exposure_id','gene','phenotype','family','variant','chrom','pos37','hit_labels']].assign(reason='no BIG40 repro file: QSM/T2* phenotypes come from OxfordQSM2022'); rq=pd.read_csv(W+'regions_not_testable.csv').rename(columns={'H4_max_p12_1e5':'H4_max'}).assign(pair_type='region',reason='no BIG40 repro file: QSM/T2* phenotypes come from OxfordQSM2022'); pd.concat([nt,rq],sort=False).to_csv(O+'not_testable_pairs.csv',index=False)
# ---- hit vs null inside the overlapping |z| band, with assay-level bootstrap
bands=[(3.85,4.0),(4.0,4.2),(4.2,4.4)]; R_=[]
h=NQ(H)[NQ(H).analysable]; n0=N[N.analysable]
for nname,n_ in [('all_matched_null',n0),('null_from_assays_without_any_hit',n0[n0.assay_has_any_hit==False])]:
  for lo,hi in bands+[(3.85,4.4)]:
    a_=h[(h.absz_disc>=lo)&(h.absz_disc<hi)]; b_=n_[(n_.absz_disc>=lo)&(n_.absz_disc<hi)]
    if len(a_)<5 or len(b_)<5: continue
    def mn(d):
        u,i=np.unique(d.exposure_id.astype(str),return_inverse=True); return np.bincount(i,d.z_rep_signed.values,len(u)),np.bincount(i,minlength=len(u))
    sa,ca=mn(a_); sb,cb=mn(b_); W1=rng.poisson(1,(NB,len(sa))); W2=rng.poisson(1,(NB,len(sb))); da=(W1@sa)/(W1@ca); db=(W2@sb)/(W2@cb); dd=da-db
    R_.append(dict(null_set=nname,absz_band=f'{lo}-{hi}',n_hit=len(a_),n_hit_assays=len(sa),n_null=len(b_),n_null_assays=len(sb),hit_mean_zrep=a_.z_rep_signed.mean(),null_mean_zrep=b_.z_rep_signed.mean(),difference=a_.z_rep_signed.mean()-b_.z_rep_signed.mean(),diff_ci_lo=float(np.nanquantile(dd,.025)),diff_ci_hi=float(np.nanquantile(dd,.975)),hit_same_sign=(a_.z_rep_signed>0).mean(),null_same_sign=(b_.z_rep_signed>0).mean(),hit_one_sided=(a_.z_rep_signed>1.6449).mean(),null_one_sided=(b_.z_rep_signed>1.6449).mean(),hit_mean_expected=a_.expected_z.mean(),null_mean_expected=b_.expected_z.mean(),hit_ratio=a_.z_rep_signed.mean()/a_.expected_z.mean(),null_ratio=b_.z_rep_signed.mean()/b_.expected_z.mean()))
pd.DataFrame(R_).to_csv(O+'hit_vs_null_same_z_band.csv',index=False)
# ---- BIG40's own genome-wide peaks as an external yardstick (no fetching)
p=pd.read_csv(ROOT+'data/big40_peak_catalog/all_10889_source_peaks.csv',low_memory=False); p['c']=pd.to_numeric(p.chromosome,errors='coerce'); q=p[(p.is_qc==False)&p.replication_beta.notna()&p.replication_se.notna()&(p.c<=22)].copy(); q['zd']=q.discovery_beta/q.discovery_se; q['zs']=q.replication_beta/q.replication_se*np.sign(q.zd); q['e']=q.discovery_beta.abs()/q.replication_se
P_=[]
for lo,hi in [(0,5.53),(5.53,6),(6,7),(7,10),(10,20),(20,1e9)]:
    g=q[(q.zd.abs()>=lo)&(q.zd.abs()<hi)]
    if len(g): P_.append(dict(absz_bin=f'{lo}-{hi}',n=len(g),same_sign=(g.zs>0).mean(),one_sided=(g.zs>1.6449).mean(),mean_zrep=g.zs.mean(),expected=g.e.mean(),ratio=g.zs.mean()/g.e.mean()))
pd.DataFrame(P_).to_csv(O+'big40_peak_yardstick.csv',index=False)
# ---- overlap
Ov=[]
for pol in ['drop','infer']:
    pp=pd.read_csv(RES+f'precision_overlap_20261116/full/pair_precision_{pol}.csv.gz',usecols=['exposure_id','phenotype','status','conditional_se','Fieller95_rho0_lower','Fieller95_rho0_upper','Fieller95_grid_lower_min','Fieller95_grid_upper_max','delta_SE_grid_min','delta_SE_grid_max']); hh=H[(H.is_qc!=True)&H.l.map(lambda s:any(x.startswith(pol+':') for x in s))]; m=hh.merge(pp,on=['exposure_id','phenotype'],how='left'); 
    Ov.append(dict(policy=pol,n_hit_pairs=len(hh),n_matched_precision_rows=int(m.Fieller95_rho0_lower.notna().sum()),frac_rho0_interval_excludes_0=float(((m.Fieller95_rho0_lower>0)|(m.Fieller95_rho0_upper<0)).mean()),frac_all_rho_grid_intervals_exclude_0=float(((m.Fieller95_grid_lower_min>0)|(m.Fieller95_grid_upper_max<0)).mean()),median_deltaSE_ratio_grid_max_to_rho0=float((m.delta_SE_grid_max/m.conditional_se).median()),F_min=float(m.F_x.min()),F_p01=float(m.F_x.quantile(.01)),F_median=float(m.F_x.median()),max_weak_iv_relative_bias_at_full_overlap_1_over_Fmin=float(1/m.F_x.min())))
pd.DataFrame(Ov).to_csv(O+'sample_overlap_assessment.csv',index=False)
rr=T[(T.pair_type=='null_random_unselected')&T.analysable]; r_=np.corrcoef(rr.z_disc_used,rr.z_rep_ea)[0,1] if len(rr)>10 else np.nan
ov2=dict(random_panel_n=len(rr),random_panel_corr_z_disc_z_rep=float(r_) if r_==r_ else None,random_panel_corr_se=float(1/np.sqrt(max(len(rr)-3,1))),random_panel_mean_product=float((rr.z_disc_used*rr.z_rep_ea).mean()) if len(rr) else None)
# ---- region file
g0=G.copy(); plan=pd.read_csv(W+'regions_table.csv.gz',low_memory=False); reg=plan[['region_id','accession','gene','phenotype','idp','family','chrom','center','model_linked','H4_500k_p12_1e5','H4_1M_p12_1e5','H4_max_p12_1e5','source_family']].copy()
def side(pt,pre):
    x=g0[g0.pair_type==pt]; c=dict(variant='variant',pos37='pos37',ea='ea',oa='oa',palindromic_flag='palindromic',z_disc_used='z_disc',beta_rep_ea='repro_beta',se_rep='repro_se',z_rep_ea='repro_z_on_ea',z_rep_signed='repro_z',same_sign='repro_same_sign',expected_z='expected_z',expected_z_se='expected_z_se_based',expected_z_tweedie='expected_z_after_wc_tweedie',one_sided_p='repro_one_sided_p',status='repro_status',digit_status='repro_digit_status',z_rep_width='repro_z_interval_width',discovery_orientation='discovery_orientation',allele_check='allele_check')
    x=x[['region_id']+list(c)].rename(columns=c); 
    return x.rename(columns={k:pre+k for k in x.columns if k!='region_id'}) if pre else x
M=reg.merge(side('region_center',''),on='region_id',how='left').merge(side('region_coloc_lead','coloc_lead_'),on='region_id',how='left').merge(side('region_imaging_lead','imaging_lead_'),on='region_id',how='left')
M.loc[M.source_family!='BIG40','repro_status']='not_testable_no_repro_file'; M['repro_same_sign']=M.repro_same_sign.map({1.0:True,0.0:False}); M.insert(1,'evidence_tier','A*'); first=['region_id','repro_status','repro_z','repro_same_sign','expected_z','repro_one_sided_p']; M=M[first+[c for c in M.columns if c not in first]]; M.to_csv(O+'region_repro.csv.gz',index=False)
# ---- pair file
cols=['pair_id','pair_type','status','analysable','region_id','exposure_id','gene','phenotype','idp','family','is_qc','hit_labels','matched_hit','match_level','z_gap','assay_has_any_hit','variant','chrom','pos37','locus','ea','oa','a1','a2','catalog_row','n_catalog_rows_at_pos','catalog_match','allele_check','sign_a2_to_ea','palindromic_flag','discovery_orientation','eaf','maf','n_disc','n_rep','beta_disc','se_disc','z_disc_ea','beta_native_disc_a2','se_native_disc','z_disc_used','absz_disc','beta_rep_a2','se_rep','mlogp_rep','z_rep_a2','beta_rep_ea','z_rep_ea','z_rep_signed','z_rep_signed_lo','z_rep_signed_hi','same_sign','one_sided_p','z_rep_signed_gt_1p645','expected_z','expected_z_n','expected_z_se','mu_tweedie_family','mu_tweedie_global','expected_z_tweedie','expected_z_tweedie_global','wc_cut_any','digit_status','z_rep_width','n_inferred_chars','n_feasible','read_tries','fetch_cluster_rows','fetch_cluster_bytes','note','log_p_mr','F_x','H4_max_p12_1e5','H4_500k_p12_1e5','H4_1M_p12_1e5','model_linked','read_in_earlier_F2r_or_s07','q_bh_global__drop','q_bh_global__infer','q_by_global__drop','q_by_global__infer','q_bh_family__drop','q_bh_family__infer','q_bh_cis_family__drop','q_bh_cis_family__infer','q_by_cis_family__drop','q_by_cis_family__infer','q_bh_estimated_m__drop','q_bh_estimated_m__infer','best_pair_of_assay__drop','best_pair_of_assay__infer']
T[[c for c in cols if c in T.columns]].to_csv(O+'repro_pairs.csv.gz',index=False)
_c=T[T.pair_type=='region_center'][['exposure_id','phenotype','z_disc_ea']].merge(T[T.pair_type.isin(['hit','null_matched','null_random_unselected'])&~T.palindromic_flag][['exposure_id','phenotype','z_disc_ea']],on=['exposure_id','phenotype'],suffixes=('_native','_atlas')); disc_chk=dict(pairs_in_both_native_window_and_atlas=len(_c),max_abs_z_difference=float((_c.z_disc_ea_native-_c.z_disc_ea_atlas).abs().max()) if len(_c) else None,n_sign_disagree=int((np.sign(_c.z_disc_ea_native)!=np.sign(_c.z_disc_ea_atlas)).sum()))
json.dump(dict(discovery_native_vs_atlas=disc_chk,n_pairs=len(T),status=T.status.value_counts().to_dict(),by_type=T.groupby('pair_type').status.value_counts().unstack(fill_value=0).to_dict('index'),digit_status=T.digit_status.value_counts().to_dict(),discovery_orientation=T[T.palindromic_flag].groupby(['pair_type','discovery_orientation']).size().rename('n').reset_index().to_dict('records'),allele_check=T.allele_check.value_counts().to_dict(),overlap_random_panel=ov2,tweedie_families=sorted(bf),tweedie_tail_points=tails,cmle_cutoffs_global={k:v for k,v in lab_t.items() if not isinstance(v,dict)}),open(O+'work/analysis_meta.json','w'),indent=1,default=str)
cc=['group','scope','n_analysable','n_assays','n_loci','median_absz_discovery','same_sign_frac','same_sign_flip_p_by_locus','one_sided_p05_frac','mean_z_rep_signed','expected_z','ratio','ratio_ci_lo_assay','ratio_ci_hi_assay','ratio_after_wc_tweedie','expected_z_after_wc_tweedie','design_effect_twoway_locus_phenotype']
print(Y[Y.scope=='overall'][[c for c in cc if c in Y.columns]].round(3).to_string()); log('DONE')
