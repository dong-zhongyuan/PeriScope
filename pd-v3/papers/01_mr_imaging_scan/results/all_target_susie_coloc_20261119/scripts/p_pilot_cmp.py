from common import *
import sys, glob
s=pd.read_csv(WORK+'pilot/sample.csv'); r=pd.read_csv(RM+'susie_results.csv')
def load(tag):
    R=WORK+'pilot/'+tag+'/'; o=[]
    for p in glob.glob(R+'*/*.csv'):
        if p.endswith('.sig.csv'): continue
        o.append(pd.read_csv(p))
    d=pd.concat(o); L=[pd.read_csv(p) for p in glob.glob(R+'*/*.sig.csv')]; L=pd.concat(L) if L else pd.DataFrame(columns=['pair_id','PP.H3','PP.H4'])
    x=L.groupby('pair_id').agg(mh3=('PP.H3','max'),mh4=('PP.H4','max'),np=('PP.H4','size')).reset_index(); return d.merge(x,how='left',on='pair_id')
for tag in sys.argv[1:]:
    d=load(tag).merge(s,on='pair_id').merge(r,on='region_id',suffixes=('','_p'))
    print('==',tag,len(d),d.status.value_counts().to_dict(), 'seconds med',d.seconds.median())
    m=d[['gene','n_snps','m','n_protein_cs','cs_x','n_imaging_cs','cs_y','lambda_protein','lambda_x','lambda_imaging','lambda_y','H4_susie_any','mh4','max_H4','mh3','max_H3','status','err']].round(4)
    print(m.sort_values('gene').to_string(index=False,max_colwidth=22))
    ok=d[d.err.isna()]
    print('n match m',(ok.n_snps==ok.m).sum(),'/',len(ok),'cs_x',(ok.n_protein_cs==ok.cs_x).sum(),'lambda_x maxdiff',np.nanmax(np.abs(ok.lambda_protein-ok.lambda_x)),'lambda_y maxdiff',np.nanmax(np.abs(ok.lambda_imaging-ok.lambda_y)),'H4 maxdiff',np.nanmax(np.abs(ok.mh4-ok.max_H4)) if ok.max_H4.notna().any() else None,'H3 maxdiff',np.nanmax(np.abs(ok.mh3-ok.max_H3)) if ok.max_H3.notna().any() else None)
