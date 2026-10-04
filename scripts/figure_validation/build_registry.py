from pathlib import Path
import numpy as np,json,hashlib
A=Path('/public/home/mengxl/dzy/pd_product_assets');C=A/'interim/randko_cache';O=A/'results/figure_program_registry';ga=json.loads((C/'gsea_disc344.json').read_text());gb=json.loads((C/'gsea_conf256.json').read_text());expected=json.loads((O/'programs.json').read_text());rebuilt={};rankings={}
for f in sorted(C.glob('*.npz')):
 z=np.load(f,allow_pickle=True);genes=z['genes'].astype(str).tolist();proteins=z['proteins'].astype(str).tolist();seeds=list(map(int,z['seeds']));q=np.arange(5)-2;D=np.tensordot(q,z['C'][[seeds.index(s) for s in range(42,47)]],axes=(0,2))/sum(q*q);avg=D.mean(0)
 for p in z['top3'].astype(str):
  c=p+'__'+f.stem;gs=sorted({g for terms in [ga.get(c,[]),gb.get(c,[])] for t in terms for g in t['genes']} & set(genes));rebuilt[c]=gs
  if c=='CD22__cDC__x__astro':rankings[c]=sorted(gs,key=lambda g:-abs(avg[proteins.index(p),genes.index(g)]))
assert expected==rebuilt
old=json.loads((O/'discovery_rankings.json').read_text());assert old==rankings
(O/'discovery_rankings.json').write_text(json.dumps(rankings,indent=2))
(O/'server_verification.json').write_text(json.dumps(dict(registry_sha256=hashlib.sha256((O/'programs.json').read_bytes()).hexdigest(),n_programs=len(rebuilt),source_fold_A=str(C/'gsea_disc344.json'),source_fold_B=str(C/'gsea_conf256.json'),rank_method='Absolute mean gene slope over seeds 42–46; ordinal quantile slope',all_members_equal=True,ranking_equal=True),indent=2))
print('24 canonical programs and discovery ranking verified from server caches')
