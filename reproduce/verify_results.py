from pathlib import Path
import json,math,csv,hashlib,statistics,collections
W=Path(__file__).resolve().parent;R=W.parent;D=R/'results';O=D/'optimization_20261003';F=R/'figure_source_data'
checks=[]
def ck(name,actual,expected):
 ok=actual==expected;checks.append(dict(item=name,actual=actual,expected=expected,passed=ok));assert ok,(name,actual,expected)
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
p=read(O/'mr_parallel_gate_20261003/priority_all_programs.csv');g=json.loads((O/'joint_analysis/programs.json').read_text())
unions=collections.defaultdict(set);recipients=collections.Counter();directions=collections.defaultdict(set)
for r in p:
 gs=set(g[r['combo']]);ck('program size '+r['combo'],len(gs),int(r['n_genes']))
 unions[r['protein']]|=gs;recipients[(r['protein'],r['brain_state'])]+=1;directions[(r['protein'],r['brain_state'])].add(r['response_direction'])
ck('target gene union sizes rebuilt from raw programs',{k:len(v) for k,v in sorted(unions.items())},{'CD22':103,'CD40':170,'CD71':196})
ck('union genes',len(set.union(*unions.values())),337)
ck('triple shared genes',len(set.intersection(*unions.values())),23)
ck('exclusive genes',{k:len(v-set.union(*(vv for kk,vv in unions.items() if kk!=k))) for k,v in sorted(unions.items())},{'CD22':48,'CD40':88,'CD71':92})
expected={'CD22/astro':2,'CD22/microglia_mhc2':1,'CD40/astro':3,'CD40/microglia_mhc2':4,'CD71/astro':2,'CD71/microglia_mhc2':4}
ck('Table 2 programs per recipient',{'/'.join(k):v for k,v in recipients.items()},expected)
ck('Table 2 program direction',{'/'.join(k):sorted(v) for k,v in directions.items()},{'CD22/astro':['up'],'CD22/microglia_mhc2':['down'],'CD40/astro':['up'],'CD40/microglia_mhc2':['up'],'CD71/astro':['down'],'CD71/microglia_mhc2':['down']})
fm=read(F/'priority_gene_membership.csv');fs={k:{r['gene'] for r in fm if r[k].lower() in ('true','1')} for k in unions}
ck('figure membership equals raw program unions',fs==unions,True)
m=read(D/'blood_brain_MR_20261003/all_brain_MR_estimates.csv');m=[r for r in m if r.get('selected','True').lower()=='true'];pv=[]
maxdiff=collections.defaultdict(float)
for r in m:
 bx,sx,by,sy,sign=[float(r[k]) for k in ['beta_exposure','se_exposure','beta_outcome','se_outcome','sign']]
 b=by*sign/bx;s=sy/abs(bx);prob=math.erfc(abs(b/s)/math.sqrt(2));pv.append(prob)
 calc={'beta_MR':b,'se_MR':s,'p_MR':prob,'F_exposure':(bx/sx)**2,'se_second_order':math.sqrt(sy**2/bx**2+by**2*sx**2/bx**4),'CI_low':b-1.96*s,'CI_high':b+1.96*s}
 for k,v in calc.items():
  maxdiff[k]=max(maxdiff[k],abs(v-float(r[k])));assert math.isclose(v,float(r[k]),rel_tol=1e-8,abs_tol=1e-12),(r['gene'],k,v,r[k])
order=sorted(range(len(pv)),key=lambda i:pv[i]);q=[0.]*len(pv);running=1.
for rank in range(len(order),0,-1):
 i=order[rank-1];running=min(running,pv[i]*len(order)/rank);q[i]=running
ck('MR family size',len(m),576)
ck('all raw MR estimates and uncertainty verified',True,True)
ck('BH independently recomputed with sorted running minima',all(math.isclose(v,float(r['q_BH']),rel_tol=1e-8,abs_tol=1e-12) for r,v in zip(m,q)),True)
sig=[(r,v) for r,v in zip(m,q) if v<.05]
ck('significant MR target-outcome pairs',len(sig),5)
ck('significant genes',sorted({r['gene'] for r,v in sig}),['CD22','CD40','TFRC'])
tfrc=next(v for r,v in sig if r['gene']=='TFRC' and 'left' in r['outcome_trait']);ck('TFRC left q three significant digits',float(f'{tfrc:.3g}'),8.31e-5)
pa=read(F/'pig_animal_sensitivity.csv');v=[float(r['spearman']) for r in pa if r['method']=='periscope'];ck('pig omission range independently rounded',[round(min(v),3),round(max(v),3)],[.238,.317])
b=read(F/'benchmark_raw.csv');grp=collections.defaultdict(list)
for r in b:grp[(r['axis'],r['method'])].append(float(r['condition_aware_rho']))
means={k:statistics.mean(v) for k,v in grp.items()};axes=sorted({r['axis'] for r in b if r['axis'].endswith('__astro')})
ck('all astrocyte direction routes lead',all(max((v,method) for (a,method),v in means.items() if a==axis)[1]=='periscope' for axis in axes),True)
ck('PeriScope route mean correlations',[round(means[(a,'periscope')],3) for a in axes],[.284,.344,.331,.316])
ck('benchmark rows',len(b),168)

print('PASS',len(checks),'independent checks; all 576 MR rows recalculated; gene membership rebuilt from original programs')
