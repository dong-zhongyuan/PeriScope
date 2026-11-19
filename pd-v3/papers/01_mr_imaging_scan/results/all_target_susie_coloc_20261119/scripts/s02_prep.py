"""Per-group inputs for SuSiE: panel-aligned protein and imaging z, panel genotypes, ST16 signals, and the single-signal (ABF) H4 where production has none."""
from common import *
import sys, time, shutil, traceback
from multiprocessing import Pool
from scipy.special import logsumexp
P=pd.read_csv(WORK+'plan_pairs.csv.gz',keep_default_na=False,low_memory=False,dtype={'idp':str})
Q=P[P.pre_status.eq('')].copy()
GD=WORK+'groups/'; os.makedirs(GD,exist_ok=True)
LUT=np.zeros((256,4),dtype=np.uint8)
for b in range(256):
    for k in range(4): LUT[b,k]={0:2,1:255,2:1,3:0}[(b>>(2*k))&3]
PH=pd.read_csv(ROOT+'data/big40_atlas/phenotypes_3935.csv',usecols=['phenotype','n_discovery'],low_memory=False).set_index('phenotype').n_discovery.to_dict()
QSM_N=19720
BIM={}; ST={}
def key_arr(ch,pos,a,b): return [f'{ch}:{p}:'+('_'.join((x,y)) if x<y else '_'.join((y,x))) for p,x,y in zip(pos,a,b)]
def load():
    b=pd.read_csv(BED+'.bim',sep='\t',header=None,usecols=[0,3,4,5],names=['chr','pos','a1','a2'],dtype={'chr':'int8','pos':'int32','a1':str,'a2':str}); b['row']=np.arange(len(b),dtype='int32')
    for c,g in b.groupby('chr',sort=False):
        g=g.sort_values('pos',kind='stable'); BIM[int(c)]=(g.pos.values,g.a1.values,g.a2.values,g.row.values)
    s=pd.read_csv(WORK+'st16_cis_signals.csv.gz',low_memory=False)
    for e,g in s.groupby('exposure_id'): ST[e]=g.to_dict('records')
def panel(ch,lo,hi):
    p,a,b,r=BIM[ch]; i,j=np.searchsorted(p,lo,'left'),np.searchsorted(p,hi,'right'); d=pd.DataFrame(dict(pos=p[i:j],a1=a[i:j],a2=b[i:j],row=r[i:j]))
    d=d[(d.a1.str.len()==1)&(d.a2.str.len()==1)]; d['k']=key_arr(ch,d.pos,d.a1,d.a2); return d.drop_duplicates('k',keep=False).reset_index(drop=True)
def std(d,ch,c,span):
    d=d.rename(columns={'base_pair_location':'pos','effect_allele':'ea','other_allele':'oa','standard_error':'se','p_value':'p','effect_allele_frequency':'eaf'})
    d=d[np.abs(d.pos-c)<=span]; d=d[(d.ea.astype(str).str.fullmatch('[ACGT]'))&(d.oa.astype(str).str.fullmatch('[ACGT]'))&(d.ea!=d.oa)]
    d=d[np.isfinite(d.beta.astype(float))&np.isfinite(d.se.astype(float))&(d.se.astype(float)>0)].copy(); d['k']=key_arr(ch,d.pos,d.ea,d.oa); return d.drop_duplicates('k',keep=False).reset_index(drop=True)
def read_img(r):
    p=r.img_path if r.img_src.startswith('window') else f'{WORK}imaging/{r.family}_{r.idp}/chr{r.chrom}_{r.pos37}.csv.gz'
    if not os.path.exists(p): return None,None
    d=pd.read_csv(p,low_memory=False); sp=FULL if (not r.img_src.startswith('window') or r.img_src=='window_exact') else HALF
    return std(d,r.chrom,r.pos37,FULL),sp
def labf(z,V,W): r=W/(W+V); return .5*(np.log1p(-r)+r*z*z)
def comb(l1,l2,p1,p2,p12):
    s1,s2,s12=logsumexp(l1),logsumexp(l2),logsumexp(l1+l2); dif=s1+s2+np.log(max(-np.expm1(min(s12-s1-s2,-1e-300)),1e-300))
    x=np.array([0,np.log(p1)+s1,np.log(p2)+s2,np.log(p1)+np.log(p2)+dif,np.log(p12)+s12]); return np.exp(x-logsumexp(x))
def abf_single(a,b,c):
    """Production recipe (all_target_coloc.R), one window."""
    m=a[['k','ea','beta','se','eaf','n']].merge(b[['k','ea','beta','se']],on='k',suffixes=('_x','_y'))
    if not len(m): return None
    bx=m.beta_x.astype(float).values; sx=m.se_x.astype(float).values; sy=m.se_y.astype(float).values; by=np.where(m.ea_x==m.ea_y,1,-1)*m.beta_y.astype(float).values
    n=float(m['n'].iloc[0]); f=m.eaf.astype(float).values; maf=np.minimum(f,1-f); ov=1/sx**2; nvx=2*n*maf*(1-maf); cf=(ov*nvx).sum()/(ov*ov).sum()
    if not np.isfinite(cf) or cf<=0: return None
    sd=np.sqrt(cf); ns=len(m); p1=1/(ns+1) if ns*1e-4>=1 else 1e-4
    l1=labf(bx/sx,sx**2,(.15*sd)**2); l2=labf(by/sy,sy**2,.15**2); out={}
    for lab,p in [('1e5',1e-5),('1e6',1e-6)]:
        pp=comb(l1,l2,p1,p1,1/(ns+1) if ns*p>=1 else p); out['H4_'+lab]=pp[4]; out['H3_'+lab]=pp[3]
    out['nsnps']=ns; return out
def ukb_single(sig,w,ch,inst_k):
    """s23 recipe: every ST16 cis signal against the imaging ABF over all window SNVs; 'primary' is the sentinel-carrying signal with the largest sentinel PIP."""
    w=w.copy(); ix={k:i for i,k in enumerate(w.k)}; z=w.beta.astype(float).values/w.se.astype(float).values; l2=labf(z,w.se.astype(float).values**2,.0225); best=None; rows=[]
    for s in sig:
        ids=str(s['cs_ids']).replace(' ','').split(','); pp=[float(x) for x in str(s['cs_pips']).replace(' ','').split(',')]; pl=float(s['pip']); ln=float(s['log10bf'])*np.log(10); l1=np.full(len(w),-50.0); mass=0.0
        for v,p in zip(ids,pp):
            g=v.split(':'); k=snv_key(ch,g[1],g[2],g[3]) if len(g)>=4 else None; j=ix.get(k)
            if j is not None: l1[j]=ln+np.log(max(p,1e-6)/max(pl,1e-6)); mass+=p
        if mass<=0: rows.append((s,mass,None)); continue
        h5=comb(l1,l2,1e-4,1e-4,1e-5); h6=comb(l1,l2,1e-4,1e-4,1e-6); rows.append((s,mass,(h5[4],h5[3],h6[4],h6[3])))
    prim=[r for r in rows if bool(r[0]['has_sentinel'])]
    prim=max(prim,key=lambda r:(0 if pd.isna(r[0]['sentinel_pip']) else float(r[0]['sentinel_pip']))) if prim else None
    anyv=[r[2][0] for r in rows if r[2] is not None]
    return dict(H4_primary=prim[2][0] if prim and prim[2] else np.nan,H3_primary=prim[2][1] if prim and prim[2] else np.nan,H4_primary_1e6=prim[2][2] if prim and prim[2] else np.nan,H3_primary_1e6=prim[2][3] if prim and prim[2] else np.nan,cs_mass_primary=prim[1] if prim else np.nan,H4_any_signal=max(anyv) if anyv else np.nan,n_signals=len(sig),nsnp_single=len(w))
def group(a):
    gid,g=a; d=GD+gid
    if os.path.exists(d+'/DONE') : return gid,'cached'
    tmp=d+f'.tmp{os.getpid()}'; shutil.rmtree(tmp,ignore_errors=True); os.makedirs(tmp)
    try:
        r0=g.iloc[0]; ch=int(r0.chrom); c=int(r0.pos37); plat=r0.platform; ikey=snv_key(ch,c,r0.ea,r0.oa)
        b=panel(ch,c-HALF-1,c+HALF+1); b=b[np.abs(b.pos-c)<=HALF]; bk=b.set_index('k'); prot=None
        if plat=='AGES':
            x=pd.read_csv(r0.protein_path,low_memory=False).rename(columns={'n_metadata':'n'}); prot=std(x,ch,c,FULL); pw=prot[np.abs(prot.pos-c)<=HALF]
        Pairs=[]; imgs={}
        for r in g.itertuples():
            y,sp=read_img(r); m=dict(pair_id=r.pair_id,exposure_id=r.exposure_id,phenotype=r.phenotype,platform=plat,accession=r.accession,region_id=r.region_id,prep_status='ok',prep_reason='')
            if y is None: m.update(prep_status='not_evaluable',prep_reason='imaging window file missing after extraction'); Pairs.append(m); continue
            yw=y[np.abs(y.pos-c)<=HALF]; Ny=float(yw.n_y.median()) if 'n_y' in yw and yw.n_y.notna().any() else float(PH.get(r.phenotype,np.nan)) if r.family=='BIG40' else float(QSM_N)
            m.update(n_imaging=Ny,n_imaging_source='regional_file' if 'n_y' in yw and yw.n_y.notna().any() else ('catalog_n_discovery' if r.family=='BIG40' else 'qsm_cohort_19720_not_idp_specific'),n_imaging_window=len(yw),img_span=sp)
            s=set(yw.k)&set(b.k)
            if plat=='AGES':
                s&=set(pw.k); m.update(n_protein_window=len(pw),n_protein_imaging_shared=len(set(pw.k)&set(yw.k)))
                if r.region_id=='' :
                    o={}
                    for lab,span in [('500k',HALF),('1000k',FULL)]:
                        if span==FULL and sp!=FULL: continue
                        q=abf_single(prot[np.abs(prot.pos-c)<=span],y[np.abs(y.pos-c)<=span],c)
                        if q: o.update({f'{k}_{lab}':v for k,v in q.items()})
                    m.update(o)
            else:
                m.update(ukb_single(ST[r.exposure_id],yw,ch,ikey))
            m['n_panel_shared_preMAF']=len(s); imgs[r.pair_id]=(yw,s,Ny); Pairs.append(m)
        U=set().union(*[v[1] for v in imgs.values()]) if imgs else set(); bb=b[b.k.isin(U)].sort_values(['pos','k']).reset_index(drop=True)
        bed=np.memmap(BED+'.bed',dtype=np.uint8,mode='r'); nb=NB
        if len(bb):
            B=np.stack([bed[3+int(q)*nb:3+(int(q)+1)*nb] for q in bb.row.values]); G=LUT[B].reshape(len(bb),-1)[:,:NIND]; miss=G==255; n_ok=(~miss).sum(1); af=np.where(miss,0,G).sum(1)/(2*n_ok); maf=np.minimum(af,1-af); sdg=np.array([g_[~m_].std() for g_,m_ in zip(G,miss)])
            keep=(maf>=.01)&(sdg>0); bb=bb[keep].reset_index(drop=True); G=G[keep]; af=af[keep]; bb['af_a1']=af
        else: G=np.zeros((0,NIND),dtype=np.uint8); bb['af_a1']=[]
        bb[['k','pos','a1','a2','row','af_a1']].to_csv(tmp+'/snps.csv',index=False); G.tofile(tmp+'/G.bin'); A=bb.set_index('k').a1; kept=set(bb.k)
        def align(fr,cols):
            f=fr[fr.k.isin(kept)].copy(); sg=np.where(f.ea.values==A.reindex(f.k).values,1.0,-1.0); f['b']=sg*f.beta.astype(float); f['s']=f.se.astype(float); f['z']=f.b/f.s; f['p']=f.p.astype(float)
            if 'eaf' in f and f.eaf.notna().any(): f['eaf_a1']=np.where(sg>0,f.eaf.astype(float),1-f.eaf.astype(float))
            else: f['eaf_a1']=np.nan
            return f[['k','pos','b','s','z','p','eaf_a1']]
        if plat=='AGES': align(pw,None).to_csv(tmp+'/protein.csv',index=False); nx=float(prot.n.astype(float).iloc[0])
        else: nx=np.nan
        for m in Pairs:
            if m['prep_status']=='ok': m['n_protein']=nx
        for pid,(yw,s,Ny) in imgs.items(): align(yw,None).to_csv(f'{tmp}/img_{pid}.csv',index=False)
        pd.DataFrame(Pairs).to_csv(tmp+'/pairs.csv',index=False)
        # instrument genotype
        ib=None; p_,a_,b_,r_=BIM[ch]; i0,j0=np.searchsorted(p_,c,'left'),np.searchsorted(p_,c,'right'); inst=dict(key=ikey,ea=r0.ea,oa=r0.oa,pos=c,variant=r0.instrument,in_panel=False,panel_row=-1,panel_a1='',panel_af=np.nan)
        for t in range(i0,j0):
            if {a_[t],b_[t]}=={r0.ea,r0.oa}:
                q=bed[3+int(r_[t])*nb:3+(int(r_[t])+1)*nb]; gv=LUT[q].reshape(-1)[:NIND]; ok=gv!=255; af_=gv[ok].sum()/(2*ok.sum()); inst.update(in_panel=True,panel_row=int(r_[t]),panel_a1=a_[t],panel_af=float(af_)); gv.tofile(tmp+'/inst.bin'); break
        json.dump(inst,open(tmp+'/instrument.json','w'))
        if plat=='UKB':
            sr=[]
            for s in ST[r0.exposure_id]:
                ids=str(s['cs_ids']).replace(' ','').split(','); pp=[float(x) for x in str(s['cs_pips']).replace(' ','').split(',')]
                for v,p in zip(ids,pp):
                    q=v.split(':'); k=snv_key(ch,q[1],q[2],q[3]) if len(q)>=4 else None
                    sr.append(dict(signal=int(s['signal']),lead=s['variant_id'],lead_pip=float(s['pip']),log10bf=float(s['log10bf']),cs_id=v,k=k or '',cs_pip=p,has_sentinel=bool(s['has_sentinel'])))
            pd.DataFrame(sr).to_csv(tmp+'/signals.csv',index=False)
        open(tmp+'/DONE','w').write('ok'); shutil.rmtree(d,ignore_errors=True); os.replace(tmp,d); return gid,'ok',len(g),len(bb)
    except Exception:
        shutil.rmtree(tmp,ignore_errors=True); open(GD+gid+'.error','w').write(traceback.format_exc()[-1500:]); return gid,'error'
if __name__=='__main__':
    import random
    n=int(sys.argv[1]); only=set(open(sys.argv[2]).read().split()) if len(sys.argv)>2 else None
    load(); T=[(g,x) for g,x in Q.groupby('gid') if (only is None or g in only)]; random.Random(1).shuffle(T); print('groups',len(T),flush=True); t=time.time(); n_e=0
    with Pool(n) as X:
        for i,r in enumerate(X.imap_unordered(group,T)):
            n_e+=r[1]=='error'
            if i%100==0 or r[1]=='error': print(i,r,round(time.time()-t),flush=True)
    print('DONE errors',n_e,flush=True)
