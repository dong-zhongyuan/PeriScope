import json, csv, numpy as np, sys
O='/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/random_pairing_20261119/training'
A='/public/home/mengxl/dzy/pd_product_assets/results'; R='/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/research_20261118/ablation/retrain'
pairs={'full47':f'{A}/benchmark_20261002/seed_repair/models/seed_47','rp42':f'{R}/random_pairing/seed_42','co44':f'{R}/coupling_off/seed_44'}
g=sys.argv[1] if len(sys.argv)>1 else 'repro'
out={}
for k,p in pairs.items():
    n=f'{O}/checks/{g}/{k}'
    a=json.load(open(n+'/training_inputs.json')); b=json.load(open(p+'/training_inputs.json'))
    ign={'out_dir','run_dir','run_id','ablate','label_shuffle','save_epochs','epochs'}
    aa={x:y for x,y in a['args'].items() if x not in ign}; bb={x:y for x,y in b['args'].items() if x not in ign}
    cfg_same=aa==bb and all(a[x]==b[x] for x in ['strata','n_genes','n_proteins','ot_scale','protein_variance','brain_train_donors','blood_train_donors','cite_train_donors'])
    va={x['epoch']:x for x in json.load(open(n+'/validation.json'))}; vb={x['epoch']:x for x in json.load(open(p+'/validation.json'))}
    rows=[]
    for e in sorted(va):
        rows.append(dict(epoch=e,new=va[e]['selection_score'],stored=vb[e]['selection_score'],abs_diff=abs(va[e]['selection_score']-vb[e]['selection_score']),new_train_loss=va[e]['train_loss'],stored_train_loss=vb[e]['train_loss']))
    la=list(csv.DictReader(open(n+'/train_log.csv'))); lb=list(csv.DictReader(open(p+'/train_log.csv')))
    m=min(len(la),len(lb)); d=[abs(float(la[i]['loss'])-float(lb[i]['loss'])) for i in range(m) if la[i]['step']==lb[i]['step']]
    out[k]=dict(config_identical=cfg_same,validation=rows,train_log_points=len(d),train_log_max_abs_loss_diff=max(d),train_log_exact=sum(x==0 for x in d))
    print(k,cfg_same,[(r['epoch'],round(r['new'],6),round(r['stored'],6),f"{r['abs_diff']:.1e}") for r in rows],'log',len(d),'exact',out[k]['train_log_exact'],'maxdiff',max(d))
json.dump(out,open(f'{O}/checks/{g}/repro_compare.json','w'),indent=1)
