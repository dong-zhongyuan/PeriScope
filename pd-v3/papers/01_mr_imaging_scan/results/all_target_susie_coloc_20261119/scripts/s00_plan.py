"""Pair universe, source availability, not_evaluable reasons. One row per platform x exposure x phenotype x accession."""
from common import *
os.makedirs(WORK,exist_ok=True)
t=pd.read_csv(RES+'all_target_delivery_20261119/all_target_region_table.csv.gz',usecols=['region_id','accession','gene','chromosome','center_grch37','phenotype','exposure_id','H4_w1000k_1e5','H4_w500k_1e5','H4_w1000k_1e6','H4_w500k_1e6'],low_memory=False)
h=pd.read_csv(RM+'hits_q05_drop.csv.gz',usecols=['exposure_id','phenotype','q_bh_global'])
u=pd.read_csv(RM+'ukb_cs_coloc_results.csv.gz',usecols=['exposure_id','phenotype','kinds','H4_primary','H3_primary','max_H4_any_signal','n_signals'])
i=pd.read_csv(RES+'full_atlas_20261116/full/selected_cis_instruments.csv',usecols=['exposure_id','gene','chrom','pos37','ea_x','oa_x','variant','MHC'],low_memory=False)
m=pd.read_csv(ROOT+'data/ages_catalog/ukb_all_assays_to_ages.csv',usecols=['exposure_id','match_status','accession'],low_memory=False); m=m[m.match_status.eq('exact_reported_target_match')&m.accession.notna()].drop_duplicates(['exposure_id','accession'])
c=pd.read_pickle('/tmp/mr_imaging_cache/st16_cis.pkl')
keep=['exposure_id','signal','variant_id','rsid','pip','log10bf','cs_size','cs_ids','cs_pips','chrom','pos','has_sentinel','sentinel_pip']
c[keep].to_csv(WORK+'st16_cis_signals.csv.gz',index=False)
st=set(c.exposure_id)
U=pd.concat([t[['exposure_id','phenotype']],h[['exposure_id','phenotype']],u[['exposure_id','phenotype']]]).drop_duplicates().merge(i,on='exposure_id',how='left')
assert U.pos37.notna().all()
H=h.set_index(['exposure_id','phenotype']).q_bh_global.to_dict(); Uk={(a,b) for a,b in zip(u.exposure_id,u.phenotype)}; Tp={(a,b) for a,b in zip(t.exposure_id,t.phenotype)}
U['family']=U.phenotype.str.split(':').str[0]; U['idp']=U.phenotype.str.split(':').str[-1]
def scan(fam):
    r={}
    for d in os.listdir(IM+fam):
        q=IM+fam+'/'+d
        if not os.path.isdir(q) or not os.path.exists(q+'/COMPLETE.json'): continue
        win=[]; whole=None
        for f in os.listdir(q):
            x=re.match(r'chr(\w+)_center(\d+)\.csv\.gz$',f)
            if x: win.append((x.group(1),int(x.group(2)),q+'/'+f))
            elif f==d+'.txt.gz': whole=q+'/'+f
        r[d]=(win,whole)
    return r
B=scan('BIG40'); Q=scan('OxfordQSM2022')
stats={f[:4]:ROOT+'data/big40/download/stats/'+f for f in os.listdir(ROOT+'data/big40/download/stats') if f.endswith('.txt.gz')}
TAR=ROOT+'data/oxford_qsm/download/QSMandT2s_IDPs.tar.gz'
def img(r):
    ch=str(int(r.chrom)); p=int(r.pos37); v=(B if r.family=='BIG40' else Q).get(r.idp)
    if v:
        ws=[w for w in v[0] if w[0]==ch]
        for w in ws:
            if w[1]==p: return 'window_exact',w[2],w[1]
        for w in ws:
            if abs(w[1]-p)<=HALF: return 'window_overlap',w[2],w[1]
        if v[1]: return 'whole',v[1],p
    if r.family=='BIG40' and r.idp in stats: return 'whole',stats[r.idp],p
    if r.family!='BIG40': return 'qsm_tar',TAR,p
    return 'none','',np.nan
x=[img(r) for r in U.itertuples()]; U['img_src']=[a[0] for a in x]; U['img_path']=[a[1] for a in x]; U['img_center']=[a[2] for a in x]
acc=m.groupby('exposure_id').accession.apply(list).to_dict(); rg={}
for r in t.itertuples(): rg.setdefault((r.exposure_id,r.accession),[]).append(r)
rows=[]
for r in U.itertuples():
    base=dict(exposure_id=r.exposure_id,gene=r.gene,phenotype=r.phenotype,family=r.family,idp=r.idp,chrom=int(r.chrom),pos37=int(r.pos37),ea=r.ea_x,oa=r.oa_x,instrument=r.variant,MHC=bool(r.MHC),img_src=r.img_src,img_path=r.img_path,img_center=r.img_center,
              in_region_table=(r.exposure_id,r.phenotype) in Tp,in_ukb_results=(r.exposure_id,r.phenotype) in Uk,mr_hit=(r.exposure_id,r.phenotype) in H,q_bh_global=H.get((r.exposure_id,r.phenotype),np.nan))
    com=None
    if r.chrom>22: com='chrX: 1000G EUR panel on disk is autosomes only'
    elif r.img_src=='none': com='imaging window not on disk (no regional extract covering instrument +/-500 kb, no whole file)'
    # UKB
    s,why=('',''); 
    if com: s,why='not_evaluable',com
    elif r.exposure_id not in st: s,why='no_protein_cs','no ST16 cis signal for this assay'
    rows.append(dict(base,platform='UKB',accession='',region_id='',protein_path='',pre_status=s,reason=why,gid=h12('UKB',r.exposure_id),pair_id=h12('UKB',r.exposure_id,r.phenotype)))
    # AGES
    al=acc.get(r.exposure_id,[])
    if not al:
        rows.append(dict(base,platform='AGES',accession='',region_id='',protein_path='',pre_status='not_evaluable',reason='no exact AGES source for this assay',gid='',pair_id=h12('AGES',r.exposure_id,r.phenotype,'')))
    for a in al:
        rr=rg.get((r.exposure_id,a),[]); rid=next((q.region_id for q in rr if q.phenotype==r.phenotype),'')
        pp=f'{REG}sources/{a}/{rr[0].region_id}/protein_regional_statistics.csv.gz' if rr else ''
        s,why='',''
        if com: s,why='not_evaluable',com
        elif not pp or not os.path.exists(pp): s,why='not_evaluable','no AGES regional protein window on disk'
        rows.append(dict(base,platform='AGES',accession=a,region_id=rid,protein_path=pp,pre_status=s,reason=why,gid=h12('AGES',r.exposure_id,a),pair_id=h12('AGES',r.exposure_id,r.phenotype,a)))
P=pd.DataFrame(rows); assert not P.pair_id.duplicated().any()
assert P.region_id[P.region_id!=''].nunique()==len(t)
P.to_csv(WORK+'plan_pairs.csv.gz',index=False)
print('rows',len(P)); print(P.groupby(['platform','pre_status'],dropna=False).size().to_string()); print(P[P.pre_status==''].groupby(['platform','img_src']).size().to_string()); print('groups',P[P.pre_status==''].groupby('platform').gid.nunique().to_dict())
print(P[P.pre_status=='not_evaluable'].groupby(['platform','reason']).size().to_string())
