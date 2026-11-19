from common import *
os.makedirs(WORK+'pilot',exist_ok=True)
r=pd.read_csv(RM+'susie_results.csv'); P=pd.read_csv(WORK+'plan_pairs.csv.gz',keep_default_na=False,low_memory=False,dtype={'idp':str}); A=P[(P.platform=='AGES')&(P.region_id!='')]
a=r[(r.cs_x>=1)&(r.cs_y>=1)].region_id.tolist(); b=r[r.err.notna()].region_id.tolist(); o=r[~r.region_id.isin(a+b)].sample(10,random_state=20261119).region_id.tolist(); s=a+b+o
x=A[A.region_id.isin(s)]; print(len(s),len(x), x.pre_status.value_counts().to_dict())
x[['region_id','pair_id','gid']].to_csv(WORK+'pilot/sample.csv',index=False); open(WORK+'pilot/groups.txt','w').write('\n'.join(sorted(set(x.gid)))+'\n'); open(WORK+'pilot/pairs.txt','w').write('\n'.join(x.pair_id)+'\n')
